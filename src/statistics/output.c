#include "common.h"
#include "config.h"
#include "domain.h"
#include "statistics.h"
#include "fileio.h"
#include "arrays/statistics.h"



/**
 * @brief save statistical data having ux^1, ux^2, uy^1, and uy^2
 * @param[in] dirname    : name of directory to which *.npy files will be written
 * @param[in] domain    : information related to MPI domain decomposition
 * @param[in] statistics : ux^1, ux^2, uy^1, and uy^2
 * @return               : error code
 */
static int save_fluid(const char dirname[], const domain_t *domain, const statistics_t *statistics){
  const int glisize = domain->glsizes[0];
  const int gljsize = domain->glsizes[1];
  const int glksize = domain->glsizes[2];
  const int   isize = domain->mysizes[0];
  const int   jsize = domain->mysizes[1];
  const int   ksize = domain->mysizes[2];
  const int ioffset = domain->offsets[0];
  const int joffset = domain->offsets[1];
  const int koffset = domain->offsets[2];
  // ux1, ux2: [1:isize+1] x [1:jsize] x [1:ksize]
  {
    const double *ux1 = statistics->ux1;
    const double *ux2 = statistics->ux2;
    const int glsizes[NDIMS] = {glksize, gljsize, glisize+1};
    const int mysizes[NDIMS] = {  ksize,   jsize,   isize+1};
    const int offsets[NDIMS] = {koffset, joffset, ioffset  };
    double *buf = common_calloc(mysizes[0]*mysizes[1]*mysizes[2], sizeof(double));
    // ux1
    for(int cnt = 0, k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize+1; i++){
          buf[cnt] = UX1(i, j, k);
          cnt++;
        }
      }
    }
    fileio_w_nd_parallel(dirname, "ux1", NDIMS, glsizes, mysizes, offsets, buf);
    // ux2
    for(int cnt = 0, k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize+1; i++){
          buf[cnt] = UX2(i, j, k);
          cnt++;
        }
      }
    }
    fileio_w_nd_parallel(dirname, "ux2", NDIMS, glsizes, mysizes, offsets, buf);
    common_free(buf);
  }
  // uy1, uy2: [0:isize+1] x [1:jsize] x [1:ksize]
  {
    const double *uy1 = statistics->uy1;
    const double *uy2 = statistics->uy2;
    const int glsizes[NDIMS] = {glksize, gljsize, glisize+2};
    const int mysizes[NDIMS] = {  ksize,   jsize,   isize+2};
    const int offsets[NDIMS] = {koffset, joffset, ioffset  };
    double *buf = common_calloc(mysizes[0]*mysizes[1]*mysizes[2], sizeof(double));
    // uy1
    for(int cnt = 0, k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 0; i <= isize+1; i++){
          buf[cnt] = UY1(i, j, k);
          cnt++;
        }
      }
    }
    fileio_w_nd_parallel(dirname, "uy1", NDIMS, glsizes, mysizes, offsets, buf);
    // uy2
    for(int cnt = 0, k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 0; i <= isize+1; i++){
          buf[cnt] = UY2(i, j, k);
          cnt++;
        }
      }
    }
    fileio_w_nd_parallel(dirname, "uy2", NDIMS, glsizes, mysizes, offsets, buf);
    common_free(buf);
  }
  // uz1, uz2: [0:isize+1] x [1:jsize] x [1:ksize]
  {
    const double *uz1 = statistics->uz1;
    const double *uz2 = statistics->uz2;
    const int glsizes[NDIMS] = {glksize, gljsize, glisize+2};
    const int mysizes[NDIMS] = {  ksize,   jsize,   isize+2};
    const int offsets[NDIMS] = {koffset, joffset, ioffset  };
    double *buf = common_calloc(mysizes[0]*mysizes[1]*mysizes[2], sizeof(double));
    // uz1
    for(int cnt = 0, k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 0; i <= isize+1; i++){
          buf[cnt] = UZ1(i, j, k);
          cnt++;
        }
      }
    }
    fileio_w_nd_parallel(dirname, "uz1", NDIMS, glsizes, mysizes, offsets, buf);
    // uz2
    for(int cnt = 0, k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 0; i <= isize+1; i++){
          buf[cnt] = UZ2(i, j, k);
          cnt++;
        }
      }
    }
    fileio_w_nd_parallel(dirname, "uz2", NDIMS, glsizes, mysizes, offsets, buf);
    common_free(buf);
  }
  return 0;
}

/**
 * @brief save statistical data having T^1 and T^2
 * @param[in] dirname    : name of directory to which *.npy files will be written
 * @param[in] domain    : information related to MPI domain decomposition
 * @param[in] statistics : T^1 and T^2
 * @return               : error code
 */
static int save_temperature(const char dirname[], const domain_t *domain, const statistics_t *statistics){
  const int glisize = domain->glsizes[0];
  const int gljsize = domain->glsizes[1];
  const int glksize = domain->glsizes[2];
  const int   isize = domain->mysizes[0];
  const int   jsize = domain->mysizes[1];
  const int   ksize = domain->mysizes[2];
  const int ioffset = domain->offsets[0];
  const int joffset = domain->offsets[1];
  const int koffset = domain->offsets[2];
  // temp1, temp2: [0:isize+1] x [1:jsize] x [1:ksize]
  {
    const double *temp1 = statistics->temp1;
    const double *temp2 = statistics->temp2;
    const int glsizes[NDIMS] = {glksize, gljsize, glisize+2};
    const int mysizes[NDIMS] = {  ksize,   jsize,   isize+2};
    const int offsets[NDIMS] = {koffset, joffset, ioffset  };
    double *buf = common_calloc(mysizes[0]*mysizes[1]*mysizes[2], sizeof(double));
    // temp1
    for(int cnt = 0, k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 0; i <= isize+1; i++){
          buf[cnt] = TEMP1(i, j, k);
          cnt++;
        }
      }
    }
    fileio_w_nd_parallel(dirname, "temp1", NDIMS, glsizes, mysizes, offsets, buf);
    // temp2
    for(int cnt = 0, k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 0; i <= isize+1; i++){
          buf[cnt] = TEMP2(i, j, k);
          cnt++;
        }
      }
    }
    fileio_w_nd_parallel(dirname, "temp2", NDIMS, glsizes, mysizes, offsets, buf);
    common_free(buf);
  }
  return 0;
}


/**
 * @brief save structures which contains collected statistical data
 * @param[in] domain     : information related to MPI domain decomposition
 * @param[in] step       : current time step
 * @param[in] time       : current time units
 * @param[in] statistics : collected statistical data
 * @return               : error code
 */
int statistics_output(const domain_t *domain, const int step, const double time, const statistics_t *statistics){
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  /* create directory from main process */
  char *dirname = fileio_concatenate(DIRNAME_STAT "/step", step);
  if(myrank == 0){
    fileio_mkdir(dirname);
  }
  // wait for the completion of mkdir
  MPI_Barrier(MPI_COMM_WORLD);
  /* save parameters */
  if(myrank == 0){
    fileio_w_0d_serial(dirname, "num",  NPYIO_INT,    sizeof(int),    &statistics->num);
    fileio_w_0d_serial(dirname, "step", NPYIO_INT,    sizeof(int),    &step);
    fileio_w_0d_serial(dirname, "time", NPYIO_DOUBLE, sizeof(double), &time);
    config.save(dirname);
  }
  /* save domain info (coordinates) */
  if(myrank == 0){
    domain_save(dirname, domain);
  }
  /* save collected statistics */
  save_fluid(dirname, domain, statistics);
  if(config.get_bool("solve_temp")){
    save_temperature(dirname, domain, statistics);
  }
  // don't forget to free memory holding directory name
  common_free(dirname);
  return 0;
}

