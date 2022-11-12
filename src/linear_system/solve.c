#include "sdecomp.h"
#include "tdm.h"
#include "linear_system.h"


#if NDIMS == 2

/**
 * @brief solve linear system in x direction
 * @param[inout] linear_system : structure storing all information about A x = b
 *                                in particular x1pncl contains right-hand-side (in) and answer (out)
 * @return                     : error code
 */
int linear_system_solve_in_x(linear_system_t * restrict linear_system){
  const int * restrict glsizes = linear_system->glsizes;
  const sdecomp_t * restrict sdecomp = linear_system->sdecomp;
  const int x1pncl_isize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_XDIR, glsizes[0]);
  const int x1pncl_jsize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_YDIR, glsizes[1]);
  tdm_solve_double(
      /* size of system  */ x1pncl_isize,
      /* number of rhs   */ x1pncl_jsize,
      /* periodic        */ false,
      /* lower- diagonal */ linear_system->tdm_l,
      /* center-diagonal */ linear_system->tdm_c,
      /* upper- diagonal */ linear_system->tdm_u,
      /* input / output  */ linear_system->x1pncl
  );
  return 0;
}

/**
 * @brief solve linear system in y direction
 * @param[inout] linear_system : structure storing all information about A x = b
 *                                in particular y1pncl contains right-hand-side (in) and answer (out)
 * @return                     : error code
 */
int linear_system_solve_in_y(linear_system_t * restrict linear_system){
  const int * restrict glsizes = linear_system->glsizes;
  const sdecomp_t * restrict sdecomp = linear_system->sdecomp;
  const int y1pncl_isize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_XDIR, glsizes[0]);
  const int y1pncl_jsize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_YDIR, glsizes[1]);
  tdm_solve_double(
      /* size of system  */ y1pncl_jsize,
      /* number of rhs   */ y1pncl_isize,
      /* periodic        */ true,
      /* lower- diagonal */ linear_system->tdm_l,
      /* center-diagonal */ linear_system->tdm_c,
      /* upper- diagonal */ linear_system->tdm_u,
      /* input / output  */ linear_system->y1pncl
  );
  return 0;
}

#else // NDIMS == 3

/**
 * @brief solve linear system in x direction
 * @param[inout] linear_system : structure storing all information about A x = b
 *                                in particular x1pncl contains right-hand-side (in) and answer (out)
 * @return                     : error code
 */
int linear_system_solve_in_x(linear_system_t * restrict linear_system){
  const int * restrict glsizes = linear_system->glsizes;
  const sdecomp_t * restrict sdecomp = linear_system->sdecomp;
  const int x1pncl_isize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_XDIR, glsizes[0]);
  const int x1pncl_jsize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_YDIR, glsizes[1]);
  const int x1pncl_ksize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_ZDIR, glsizes[2]);
  tdm_solve_double(
      /* size of system  */ x1pncl_isize,
      /* number of rhs   */ x1pncl_jsize * x1pncl_ksize,
      /* periodic        */ false,
      /* lower- diagonal */ linear_system->tdm_l,
      /* center-diagonal */ linear_system->tdm_c,
      /* upper- diagonal */ linear_system->tdm_u,
      /* input / output  */ linear_system->x1pncl
  );
  return 0;
}

/**
 * @brief solve linear system in y direction
 * @param[inout] linear_system : structure storing all information about A x = b
 *                                in particular y1pncl contains right-hand-side (in) and answer (out)
 * @return                     : error code
 */
int linear_system_solve_in_y(linear_system_t * restrict linear_system){
  const int * restrict glsizes = linear_system->glsizes;
  const sdecomp_t * restrict sdecomp = linear_system->sdecomp;
  const int y1pncl_isize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_XDIR, glsizes[0]);
  const int y1pncl_jsize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_YDIR, glsizes[1]);
  const int y1pncl_ksize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_ZDIR, glsizes[2]);
  tdm_solve_double(
      /* size of system  */ y1pncl_jsize,
      /* number of rhs   */ y1pncl_ksize * y1pncl_isize,
      /* periodic        */ true,
      /* lower- diagonal */ linear_system->tdm_l,
      /* center-diagonal */ linear_system->tdm_c,
      /* upper- diagonal */ linear_system->tdm_u,
      /* input / output  */ linear_system->y1pncl
  );
  return 0;
}

/**
 * @brief solve linear system in z direction
 * @param[inout] linear_system : structure storing all information about A x = b
 *                                in particular z1pncl contains right-hand-side (in) and answer (out)
 * @return                     : error code
 */
int linear_system_solve_in_z(linear_system_t * restrict linear_system){
  const int * restrict glsizes = linear_system->glsizes;
  const sdecomp_t * restrict sdecomp = linear_system->sdecomp;
  const int z1pncl_isize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Z1PENCIL, SDECOMP_XDIR, glsizes[0]);
  const int z1pncl_jsize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Z1PENCIL, SDECOMP_YDIR, glsizes[1]);
  const int z1pncl_ksize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Z1PENCIL, SDECOMP_ZDIR, glsizes[2]);
  tdm_solve_double(
      /* size of system  */ z1pncl_ksize,
      /* number of rhs   */ z1pncl_isize * z1pncl_jsize,
      /* periodic        */ true,
      /* lower- diagonal */ linear_system->tdm_l,
      /* center-diagonal */ linear_system->tdm_c,
      /* upper- diagonal */ linear_system->tdm_u,
      /* input / output  */ linear_system->z1pncl
  );
  return 0;
}

#endif // NDIMS
