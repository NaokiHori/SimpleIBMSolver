#include <math.h>
#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "temperature.h"
#include "ib.h"
#include "save.h"
#include "fileio.h"
#include "config.h"
#include "arrays/fluid.h"
#include "arrays/temperature.h"


static bool is_scheduled = false;
static double save_rate;
double save_next;


/**
 * @brief save members in fluid_t
 * @param[in] dirname : name of directory to which *.npy files will be written
 * @param[in] domain  : information related to MPI domain decomposition
 * @param[in] fluid   : velocity, pressure
 * @return            : error code
 */
int save_fluid(const char dirname[], const domain_t *domain, const fluid_t *fluid){
  // pack arrays (ux, uy, p) into buffer without halo and write them to the corresponding file
  const int glisize = domain->glsizes[0];
  const int gljsize = domain->glsizes[1];
  const int   isize = domain->mysizes[0];
  const int   jsize = domain->mysizes[1];
  const int ioffset = domain->offsets[0];
  const int joffset = domain->offsets[1];
  // ux: [1:isize+1] x [1:jsize]
  {
    const double *ux = fluid->ux;
    const int glsizes[NDIMS] = {gljsize, glisize+1};
    const int mysizes[NDIMS] = {  jsize,   isize+1};
    const int offsets[NDIMS] = {joffset, ioffset  };
    double *buf = common_calloc(mysizes[0]*mysizes[1], sizeof(double));
    for(int cnt = 0, j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize+1; i++){
        buf[cnt] = UX(i, j);
        cnt++;
      }
    }
    fileio_w_nd_parallel(dirname, "ux", NDIMS, glsizes, mysizes, offsets, buf);
    common_free(buf);
  }
  // uy: [0:isize+1] x [1:jsize]
  {
    const double *uy = fluid->uy;
    const int glsizes[NDIMS] = {gljsize, glisize+2};
    const int mysizes[NDIMS] = {  jsize,   isize+2};
    const int offsets[NDIMS] = {joffset, ioffset  };
    double *buf = common_calloc(mysizes[0]*mysizes[1], sizeof(double));
    for(int cnt = 0, j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize+1; i++){
        buf[cnt] = UY(i, j);
        cnt++;
      }
    }
    fileio_w_nd_parallel(dirname, "uy", NDIMS, glsizes, mysizes, offsets, buf);
    common_free(buf);
  }
  // p: [0:isize+1] x [1:jsize]
  {
    const double *p = fluid->p;
    const int glsizes[NDIMS] = {gljsize, glisize+2};
    const int mysizes[NDIMS] = {  jsize,   isize+2};
    const int offsets[NDIMS] = {joffset, ioffset  };
    double *buf = common_calloc(mysizes[0]*mysizes[1], sizeof(double));
    for(int cnt = 0, j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize+1; i++){
        buf[cnt] = P(i, j);
        cnt++;
      }
    }
    fileio_w_nd_parallel(dirname, "p", NDIMS, glsizes, mysizes, offsets, buf);
    common_free(buf);
  }
  return 0;
}

/**
 * @brief save members in temperature_t
 * @param[in] dirname     : name of directory to which *.npy files will be written
 * @param[in] domain      : information related to MPI domain decomposition
 * @param[in] temperature : temperature
 * @return                : error code
 */
int save_temperature(const char dirname[], const domain_t *domain, const temperature_t *temperature){
  const int glisize = domain->glsizes[0];
  const int gljsize = domain->glsizes[1];
  const int   isize = domain->mysizes[0];
  const int   jsize = domain->mysizes[1];
  const int ioffset = domain->offsets[0];
  const int joffset = domain->offsets[1];
  // temp: [0:isize+1] x [1:jsize]
  {
    const double *temp = temperature->temp;
    const int glsizes[NDIMS] = {gljsize, glisize+2};
    const int mysizes[NDIMS] = {  jsize,   isize+2};
    const int offsets[NDIMS] = {joffset, ioffset  };
    double *buf = common_calloc(mysizes[0]*mysizes[1], sizeof(double));
    for(int cnt = 0, j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize+1; i++){
        buf[cnt] = TEMP(i, j);
        cnt++;
      }
    }
    fileio_w_nd_parallel(dirname, "temp", NDIMS, glsizes, mysizes, offsets, buf);
    common_free(buf);
  }
  return 0;
}


/**
 * @brief save structures which contains essential information
 *          for restart and post-processing
 * @param[in] domain      : information related to MPI domain decomposition
 * @param[in] step        : time step
 * @param[in] time        : current time units
 * @param[in] fluid       : velocity
 * @param[in] temperature : temperature
 * @return                : error code
 */
int save(const domain_t *domain, const int step, const double time, const fluid_t *fluid, const temperature_t *temperature, const ib_t *ib){
  if(!is_scheduled){
    // timing has not been scheduled yet
    // load from configuration
    save_rate    = config.get_double("save_rate");
    double after = config.get_double("save_after");
    // find time to trigger next event
    save_next = save_rate * ceil( fmax(time, after) / save_rate);
    is_scheduled = true;
  }else{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    /* ! create directory from main process ! 6 ! */
    char *dirname = fileio_concatenate(DIRNAME_SAVE "/step", step);
    if(myrank == 0){
      fileio_mkdir(dirname);
    }
    // wait for the completion of mkdir
    MPI_Barrier(MPI_COMM_WORLD);
    /* ! save parameters ! 5 ! */
    if(myrank == 0){
      fileio_w_0d_serial(dirname, "step", NPYIO_INT,    sizeof(int),    &step);
      fileio_w_0d_serial(dirname, "time", NPYIO_DOUBLE, sizeof(double), &time);
      config.save(dirname);
    }
    /* save domain */
    if(myrank == 0){
      domain_save(dirname, domain);
    }
    /* ! save flow fields ! 4 ! */
    save_fluid(dirname, domain, fluid);
    if(config.get_bool("solve_temp")){
      save_temperature(dirname, domain, temperature);
    }
    if(myrank == 0){
      ib_save(dirname, ib);
    }
    common_free(dirname);
    save_next += save_rate;
  }
  return 0;
}

