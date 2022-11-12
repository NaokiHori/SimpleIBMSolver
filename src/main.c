#include <stdbool.h>
#include <mpi.h>
#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "temperature.h"
#include "ib.h"
#include "statistics.h"
#include "save.h"
#include "logging.h"
#include "decide_dt.h"
#include "config.h"
#include "fileio.h"
#include "tdm.h"


/**
 * @brief main function
 * @return : error code
 */
int main(void){
  /* ! launch MPI, start timer ! 3 ! */
  MPI_Init(NULL, NULL);
  double wtimes[2] = {0., 0.};
  wtimes[0] = common_get_wtime();
  // load environmental variables
  config.load();
  /* ! initialise time step and current (simulation) time ! 15 ! */
  int step;
  double time;
  {
    const bool restart_sim = config.get_bool("restart_sim");
    if(restart_sim){
      // directory name from which information are loaded
      const char *dirname = config.get_string("restart_dir");
      fileio_r_0d_serial(dirname, "step", NPYIO_INT,    sizeof(   int), &step);
      fileio_r_0d_serial(dirname, "time", NPYIO_DOUBLE, sizeof(double), &time);
    }else{
      // start from 0
      step = 0;
      time = 0.;
    }
  }
  /* ! initialise structures ! 4 ! */
  domain_t      *domain      = domain_init();
  fluid_t       *fluid       = fluid_init(domain);
  temperature_t *temperature = temperature_init(domain);
  ib_t          *ib          = ib_init(domain);
  statistics_t  *statistics  = statistics_init(domain);
  /* main loop */
  for(;;){
    /* ! decide time step size and diffusive term treatment ! 1 ! */
    const double dt = decide_dt(domain, fluid);
    /* ! integrate mass, momentum, and internal energy balances in time ! 17 ! */
    for(int rkstep = 0; rkstep < RKSTEPMAX; rkstep++){
      /* ! update temperature and thermal forcing ! 5 ! */
      if(config.get_bool("solve_temp")){
        temperature_compute_force(domain, temperature);
        temperature_compute_rhs(domain, rkstep, fluid, temperature);
        temperature_update_temp(domain, rkstep, dt, temperature);
      }
      ib_compute_inertia(domain, 0, fluid, ib);
      /* ! update velocity by integrating momentum equation ! 2 ! */
      fluid_compute_rhs(domain, rkstep, fluid, temperature);
      fluid_update_velocity(domain, rkstep, dt, fluid);
      ib_exchange_momentum(domain, fluid, dt, ib);
      ib_update_momentum_fleid(domain, fluid, ib);
      /* ! compute scalar potential ! 1 ! */
      fluid_compute_potential(domain, rkstep, dt, fluid);
      /* ! correct velocity to be solenoidal ! 1 ! */
      fluid_correct_velocity(domain, rkstep, dt, fluid);
      /* ! update pressure ! 1 ! */
      fluid_update_pressure(domain, rkstep, dt, fluid);
      /*** update particles iteratively ***/
      ib_compute_collision_force(domain, 0, ib);
      for(int substep = 0; ; substep++){
        ib_compute_inertia(domain, 1, fluid, ib);
        ib_compute_collision_force(domain, 1, ib);
        const double residual_max = 1.e-8;
        double residual;
        ib_increment_particles(rkstep, dt, ib, &residual);
        if(residual < residual_max){
          break;
        }
        const int substep_max = 100;
        if(substep > substep_max){
          break;
        }
      }
      ib_update_particles(domain, ib);
    }
    /* ! step and time are incremented ! 3 ! */
    step += 1;
    time += dt;
    wtimes[1] = common_get_wtime();
    /* ! output log ! 3 ! */
    if(logging_next < time){
      logging(domain, step, time, dt, wtimes[1]-wtimes[0], fluid, temperature);
    }
    /* ! save flow fields ! 3 ! */
    if(save_next < time){
      save(domain, step, time, fluid, temperature, ib);
    }
    /* ! collect statistics ! 3 ! */
    if(stat_next < time){
      statistics_collect(domain, time, fluid, temperature, statistics);
    }
    /* ! terminate when the simulation is finished ! 3 ! */
    if(time > config.get_double("timemax")){
      break;
    }
    /* ! terminate when wall time limit is reached ! 3 ! */
    if(wtimes[1]-wtimes[0] > config.get_double("wtimemax")){
      break;
    }
  }
  /* ! save restart file and statistics at last ! 2 ! */
  save(domain, step, time, fluid, temperature, ib);
  statistics_output(domain, step, time, statistics);
  /* ! finalise structures ! 4 ! */
  statistics_finalise(statistics);
  ib_finalise(ib);
  temperature_finalise(temperature);
  fluid_finalise(fluid);
  domain_finalise(domain);
  // free internal memories used by tdm
  tdm_cleanup();
  // free memory used to store environmental variables
  config.unload();
  /* ! finalise MPI ! 1 ! */
  MPI_Finalize();
  return 0;
}

