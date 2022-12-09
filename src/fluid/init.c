#include "common.h"
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "config.h"
#include "arrays/fluid.h"



/**
 * @brief allocate fluid_t
 * @param[in] domain : information about domain decomposition and size
 * @return           : structure being allocated
 */
static fluid_t *allocate(const domain_t * restrict domain){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  /* ! structure is allocated ! 1 ! */
  fluid_t * restrict fluid = common_calloc(1, sizeof(fluid_t));
  /* ! velocity, pressure, scalar potential are allocated ! 5 ! */
  fluid->ux  = common_calloc( UX_SIZE_0 *  UX_SIZE_1 *  UX_SIZE_2, sizeof(double));
  fluid->uy  = common_calloc( UY_SIZE_0 *  UY_SIZE_1 *  UY_SIZE_2, sizeof(double));
  fluid->uz  = common_calloc( UZ_SIZE_0 *  UZ_SIZE_1 *  UZ_SIZE_2, sizeof(double));
  fluid->p   = common_calloc(  P_SIZE_0 *   P_SIZE_1 *   P_SIZE_2, sizeof(double));
  fluid->psi = common_calloc(PSI_SIZE_0 * PSI_SIZE_1 * PSI_SIZE_2, sizeof(double));
  /* ! Runge-Kutta source terms are allocated ! 9 ! */
  fluid->srcuxa = common_calloc(SRCUXA_SIZE_0 * SRCUXA_SIZE_1 * SRCUXA_SIZE_2, sizeof(double));
  fluid->srcuxb = common_calloc(SRCUXB_SIZE_0 * SRCUXB_SIZE_1 * SRCUXB_SIZE_2, sizeof(double));
  fluid->srcuxg = common_calloc(SRCUXG_SIZE_0 * SRCUXG_SIZE_1 * SRCUXG_SIZE_2, sizeof(double));
  fluid->srcuya = common_calloc(SRCUYA_SIZE_0 * SRCUYA_SIZE_1 * SRCUYA_SIZE_2, sizeof(double));
  fluid->srcuyb = common_calloc(SRCUYB_SIZE_0 * SRCUYB_SIZE_1 * SRCUYB_SIZE_2, sizeof(double));
  fluid->srcuyg = common_calloc(SRCUYG_SIZE_0 * SRCUYG_SIZE_1 * SRCUYG_SIZE_2, sizeof(double));
  fluid->srcuza = common_calloc(SRCUZA_SIZE_0 * SRCUZA_SIZE_1 * SRCUZA_SIZE_2, sizeof(double));
  fluid->srcuzb = common_calloc(SRCUZB_SIZE_0 * SRCUZB_SIZE_1 * SRCUZB_SIZE_2, sizeof(double));
  fluid->srcuzg = common_calloc(SRCUZG_SIZE_0 * SRCUZG_SIZE_1 * SRCUZG_SIZE_2, sizeof(double));
  return fluid;
}

/**
 * @brief load flow field
 * @param[in ] domain : information about domain decomposition and size
 * @param[out] fluid  : velocity and pressure to be loaded
 * @return            : structure whose members are initialised
 */
static int load(const domain_t * restrict domain, fluid_t * restrict fluid){
  const int glisize = domain->glsizes[0];
  const int gljsize = domain->glsizes[1];
  const int glksize = domain->glsizes[2];
  const int   isize = domain->mysizes[0];
  const int   jsize = domain->mysizes[1];
  const int   ksize = domain->mysizes[2];
  const int ioffset = domain->offsets[0];
  const int joffset = domain->offsets[1];
  const int koffset = domain->offsets[2];
  /* ! ux is loaded ! 21 ! */
  {
    // ux: [1:isize+1] x [1:jsize] x [1:ksize]
    double * restrict ux = fluid->ux;
    const int glsizes[NDIMS] = {glksize, gljsize, glisize+1};
    const int mysizes[NDIMS] = {  ksize,   jsize,   isize+1};
    const int offsets[NDIMS] = {koffset, joffset, ioffset  };
    double *buf = common_calloc(mysizes[0]*mysizes[1]*mysizes[2], sizeof(double));
    const char *dirname = config.get_string("restart_dir");
    fileio_r_nd_parallel(dirname, "ux", NDIMS, glsizes, mysizes, offsets, buf);
    for(int cnt = 0, k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize+1; i++){
          UX(i, j, k) = buf[cnt];
          cnt++;
        }
      }
    }
    common_free(buf);
    // update boundary and halo values
    fluid_update_boundaries_ux(domain, ux);
  }
  /* ! uy is loaded ! 21 ! */
  {
    // uy: [0:isize+1] x [1:jsize] x [1:ksize]
    double * restrict uy = fluid->uy;
    const int glsizes[NDIMS] = {glksize, gljsize, glisize+2};
    const int mysizes[NDIMS] = {  ksize,   jsize,   isize+2};
    const int offsets[NDIMS] = {koffset, joffset, ioffset  };
    double *buf = common_calloc(mysizes[0]*mysizes[1]*mysizes[2], sizeof(double));
    const char *dirname = config.get_string("restart_dir");
    fileio_r_nd_parallel(dirname, "uy", NDIMS, glsizes, mysizes, offsets, buf);
    for(int cnt = 0, k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 0; i <= isize+1; i++){
          UY(i, j, k) = buf[cnt];
          cnt++;
        }
      }
    }
    common_free(buf);
    // update boundary and halo values
    fluid_update_boundaries_uy(domain, uy);
  }
  /* ! uz is loaded ! 21 ! */
  {
    // uz: [0:isize+1] x [1:jsize] x [1:ksize]
    double * restrict uz = fluid->uz;
    const int glsizes[NDIMS] = {glksize, gljsize, glisize+2};
    const int mysizes[NDIMS] = {  ksize,   jsize,   isize+2};
    const int offsets[NDIMS] = {koffset, joffset, ioffset  };
    double *buf = common_calloc(mysizes[0]*mysizes[1]*mysizes[2], sizeof(double));
    const char *dirname = config.get_string("restart_dir");
    fileio_r_nd_parallel(dirname, "uz", NDIMS, glsizes, mysizes, offsets, buf);
    for(int cnt = 0, k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 0; i <= isize+1; i++){
          UZ(i, j, k) = buf[cnt];
          cnt++;
        }
      }
    }
    common_free(buf);
    // update boundary and halo values
    fluid_update_boundaries_uz(domain, uz);
  }
  /* ! p is loaded ! 21 ! */
  {
    // p: [0:isize+1] x [1:jsize] x [1:ksize]
    double * restrict p = fluid->p;
    const int glsizes[NDIMS] = {glksize, gljsize, glisize+2};
    const int mysizes[NDIMS] = {  ksize,   jsize,   isize+2};
    const int offsets[NDIMS] = {koffset, joffset, ioffset  };
    double *buf = common_calloc(mysizes[0]*mysizes[1]*mysizes[2], sizeof(double));
    const char *dirname = config.get_string("restart_dir");
    fileio_r_nd_parallel(dirname, "p", NDIMS, glsizes, mysizes, offsets, buf);
    for(int cnt = 0, k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 0; i <= isize+1; i++){
          P(i, j, k) = buf[cnt];
          cnt++;
        }
      }
    }
    common_free(buf);
    // update boundary and halo values
    fluid_update_boundaries_p(domain, p);
  }
  return 0;
}

/**
 * @brief initialise flow field
 * @param[in ] domain : information about domain decomposition and size
 * @param[out] fluid  : velocity and pressure to be initialised
 * @return            : structure whose members are initialised
 */
static int init(const domain_t * restrict domain, fluid_t * restrict fluid){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  /* ! ux is initialised ! 12 ! */
  {
    double * restrict ux = fluid->ux;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 2; i <= isize; i++){
          UX(i, j, k) = 0.;
        }
      }
    }
    // update boundary and halo values
    fluid_update_boundaries_ux(domain, ux);
  }
  /* ! uy is initialised ! 12 ! */
  {
    double * restrict uy = fluid->uy;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          UY(i, j, k) = 0.;
        }
      }
    }
    // update boundary and halo values
    fluid_update_boundaries_uy(domain, uy);
  }
  /* ! uz is initialised ! 12 ! */
  {
    double * restrict uz = fluid->uz;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          UZ(i, j, k) = 0.;
        }
      }
    }
    // update boundary and halo values
    fluid_update_boundaries_uz(domain, uz);
  }
  /* ! p is initialised ! 12 ! */
  {
    double * restrict p = fluid->p;
    for(int k = 1; k <= ksize; k++){
      for(int j = 1; j <= jsize; j++){
        for(int i = 1; i <= isize; i++){
          P(i, j, k) = 0.;
        }
      }
    }
    // update boundary and halo values
    fluid_update_boundaries_p(domain, p);
  }
  return 0;
}


/**
 * @brief constructor of the structure
 * @param[in] domain : information about domain decomposition and size
 * @return           : structure being allocated and initalised
 */
fluid_t *fluid_init(const domain_t * restrict domain){
  /* ! allocate structure and its members ! 1 ! */
  fluid_t * restrict fluid = allocate(domain);
  /* ! initialise or load velocity and pressure ! 6 ! */
  const bool restart_sim = config.get_bool("restart_sim");
  if(restart_sim){
    load(domain, fluid);
  }else{
    init(domain, fluid);
  }
  return fluid;
}

