#include <math.h>
#include <time.h>
#include "common.h"
#include "domain.h"
#include "temperature.h"
#include "fileio.h"
#include "config.h"
#include "arrays/temperature.h"



/**
 * @brief allocate temperature_t
 * @param[in] domain : information about domain decomposition and size
 * @return           : structure being allocated
 */
static temperature_t *allocate(const domain_t * restrict domain){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  /* ! structure is allocated ! 1 ! */
  temperature_t * restrict temperature = common_calloc(1, sizeof(temperature_t));
  /* ! temperature and buoyancy force fields are allocated ! 2 ! */
  temperature->temp       = common_calloc(      TEMP_SIZE_0 *       TEMP_SIZE_1 *       TEMP_SIZE_2, sizeof(double));
  temperature->tempforcex = common_calloc(TEMPFORCEX_SIZE_0 * TEMPFORCEX_SIZE_1 * TEMPFORCEX_SIZE_2, sizeof(double));
  /* ! Runge-Kutta source terms are allocated ! 3 ! */
  temperature->srctempa = common_calloc(SRCTEMPA_SIZE_0 * SRCTEMPA_SIZE_1 * SRCTEMPA_SIZE_2, sizeof(double));
  temperature->srctempb = common_calloc(SRCTEMPB_SIZE_0 * SRCTEMPB_SIZE_1 * SRCTEMPB_SIZE_2, sizeof(double));
  temperature->srctempg = common_calloc(SRCTEMPG_SIZE_0 * SRCTEMPG_SIZE_1 * SRCTEMPG_SIZE_2, sizeof(double));
  return temperature;
}

/**
 * @brief load temperature field
 * @param[in ] domain      : information related to MPI domain decomposition
 * @param[out] temperature : temperature
 * @return                 : structure whose members are initialised
 */
static int load(const domain_t * restrict domain, temperature_t * restrict temperature){
  const int glisize = domain->glsizes[0];
  const int gljsize = domain->glsizes[1];
  const int glksize = domain->glsizes[2];
  const int   isize = domain->mysizes[0];
  const int   jsize = domain->mysizes[1];
  const int   ksize = domain->mysizes[2];
  const int ioffset = domain->offsets[0];
  const int joffset = domain->offsets[1];
  const int koffset = domain->offsets[2];
  /* ! temp is loaded ! 21 ! */
  {
    // temp: [0:isize+1] x [1:jsize] x [1:ksize]
    double * restrict temp = temperature->temp;
    const int glsizes[NDIMS] = {glksize, gljsize, glisize+2};
    const int mysizes[NDIMS] = {  ksize,   jsize,   isize+2};
    const int offsets[NDIMS] = {koffset, joffset, ioffset  };
    double *buf = common_calloc(mysizes[0]*mysizes[1]*mysizes[2], sizeof(double));
    const char *dirname = config.get_string("restart_dir");
    fileio_r_nd_parallel(dirname, "temp", NDIMS, glsizes, mysizes, offsets, buf);
    for(int cnt = 0, k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 0; i <= isize+1; i++){
          TEMP(i, j, k) = buf[cnt];
          cnt++;
        }
      }
    }
    common_free(buf);
    // update boundary and halo values
    temperature_update_boundaries_temp(domain, temp);
  }
  return 0;
}

/**
 * @brief initialise temperature field
 * @param[in ] domain      : information related to MPI domain decomposition
 * @param[out] temperature : temperature
 * @return                 : structure whose members are initialised
 */
static int init(const domain_t * restrict domain, temperature_t * restrict temperature){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  /* ! temp is initialised ! 18 ! */
  {
    int mpirank;
    MPI_Comm_rank(domain->sdecomp->comm_cart, &mpirank);
    srand(time(NULL)+mpirank);
    double * restrict temp = temperature->temp;
    // random field
    const double factor = 0.01;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          double r = -0.5+1.*rand()/RAND_MAX;
          TEMP(i, j, k) = factor*r;
        }
      }
    }
    // update boundary and halo values
    temperature_update_boundaries_temp(domain, temp);
  }
  return 0;
}

/**
 * @brief constructor of the structure
 * @param[in] domain : information about domain decomposition and size
 * @return           : structure being allocated and initalised
 */
temperature_t *temperature_init(const domain_t * restrict domain){
  /* ! allocate structure and its members ! 1 ! */
  temperature_t *temperature = allocate(domain);
  /* ! initialise or load temperature ! 1 ! */
  const bool restart_sim = config.get_bool("restart_sim");
  if(restart_sim){
    load(domain, temperature);
  }else{
    init(domain, temperature);
  }
  return temperature;
}

