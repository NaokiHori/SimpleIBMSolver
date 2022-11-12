#include <stdio.h>
#include <stdbool.h>
#include "common.h"
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "temperature.h"
#include "fileio.h"
#include "config.h"
#include "logging.h"
#include "arrays/fluid.h"
#include "arrays/temperature.h"
#include "internal.h"


static bool is_scheduled = false;
static double logging_rate;
double logging_next;

/**
 * @brief show current step, time, time step size, diffusive treatments
 * @param[in] fname : file name to which the log is written
 * @param[in] step  : current time step
 * @param[in] time  : current simulation time
 * @param[in] dt    : time step size
 * @param[in] wtime : current wall time
 * @return          : error code
 */
static int show_progress(const char fname[], const int step, const double time, const double dt, const double wtime){
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if(myrank == 0){
    FILE *fp = fileio_fopen(fname, "a");
    if(fp != NULL){
      const int ndigits_timemax  = common_get_ndigits((int)config.get_double("timemax"));
      const int ndigits_wtimemax = common_get_ndigits((int)config.get_double("wtimemax"));
      /* ! show progress to standard output and file ! 12 ! */
      // output to stdout and file
#define MPRINT(...) { \
      fprintf(fp,     __VA_ARGS__); \
      fprintf(stdout, __VA_ARGS__); \
}
      MPRINT("step %*d, time %*.1f, dt %.2e, elapsed %*.1f [sec]\n",
          NDIGITS_STEP, step,
          ndigits_timemax+2, time,
          dt,
          ndigits_wtimemax+2, wtime
      );
#undef MPRINT
      fileio_fclose(fp);
    }
  }
  return 0;
}

/**
 * @brief output log files to be monitored during simulation
 * @param[in] domain      : information related to MPI domain decomposition
 * @param[in] step        : current time step
 * @param[in] time        : current simulation time
 * @param[in] dt          : time step size
 * @param[in] wtime       : current wall time
 * @param[in] fluid       : velocity
 * @param[in] temperature : temperature
 * @return                : error code
 */
int logging(const domain_t *domain, const int step, const double time, const double dt, const double wtime, const fluid_t *fluid, const temperature_t *temperature){
  if(!is_scheduled){
    // timing has not been scheduled yet
    // load from configuration
    logging_rate = config.get_double("log_rate");
    double after = config.get_double("log_after");
    // find time to trigger next event
    logging_next = logging_rate * ceil( fmax(time, after) / logging_rate );
    is_scheduled = true;
  }else{
    show_progress   (DIRNAME_LOG "/progress.dat",  step, time, dt, wtime);
    check_divergence(DIRNAME_LOG "/divergence.dat", domain, time, fluid);
    check_momentum  (DIRNAME_LOG "/momentum.dat",   domain, time, fluid);
    check_energy    (DIRNAME_LOG "/energy.dat",     domain, time, fluid, temperature);
    check_nusselt   (DIRNAME_LOG "/nusselt.dat",    domain, time, fluid, temperature);
    logging_next += logging_rate;
  }
  return 0;
}

