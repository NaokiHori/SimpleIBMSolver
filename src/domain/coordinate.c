#include <math.h>
#include "common.h"
#include "config.h"
#include "domain.h"
#include "internal.h"


/**
 * @brief define cell-face positions in x direction
 * @param[in] isize : number of cell-centers in x direction (boundary excluded)
 * @param[in] lx    : physical domain size in x direction
 *                      normally 1 since this is the reference length scale
 *                      to non-dimensionalise equations
 * @return : cell-face positions in x direction
 */
double *allocate_and_init_xf(const int isize, const double lx){
  /* ! xf: cell face coordinates ! 11 ! */
  double *xf = common_calloc(XF_SIZE, sizeof(double));
  // uniform grid, which enables us to use
  //   efficient DCT-based Poisson solver
  //   see src/fluid/compute_potential.c
  const double dx = lx/isize;
  for(int i = 1; i <= isize+1; i++){
    XF(i) = 1.*(i-1)*dx;
  }
  // force boundary values just in case
  XF(      1) = 0.;
  XF(isize+1) = lx;
  return xf;
}

/**
 * @brief define cell-center positions in x direction
 * @param[in] isize : number of cell-centers in x direction (boundary excluded)
 * @param[in] xf    : cell-face positions in x direction
 * @return : cell-center positions in x direction
 */
double *allocate_and_init_xc(const int isize, const double *xf){
  /* ! xc: cell center coordinates ! 8 ! */
  double *xc = common_calloc(XC_SIZE, sizeof(double));
  // center between two XFs
  for(int i = 1; i <= isize; i++){
    XC(i) = 0.5*(XF(i  )+XF(i+1));
  }
  // at boundaries, face positions are assigned
  XC(      0) = XF(      1);
  XC(isize+1) = XF(isize+1);
  return xc;
}

/**
 * @brief define face-to-face distances in x direction
 * @param[in] isize : number of cell-centers in x direction (boundary excluded)
 * @param[in] xf    : cell-face positions in x direction
 * @return : face-to-face distances in x direction
 */
double *allocate_and_init_dxf(const int isize, const double *xf){
  /* ! dxf: distance from cell face to cell face ! 4 ! */
  double *dxf = common_calloc(DXF_SIZE, sizeof(double));
  for(int i = 1; i <= isize; i++){
    DXF(i) = XF(i+1)-XF(i  );
  }
  return dxf;
}

/**
 * @brief define center-to-center distances in x direction
 * @param[in] isize : number of cell-centers in x direction (boundary excluded)
 * @param[in] xf    : cell-face positions in x direction
 * @return : center-to-center distances in x direction
 */
double *allocate_and_init_dxc(const int isize, const double *xc){
  /* ! dxc: distance from cell center to cell center (generally) ! 4 ! */
  double *dxc = common_calloc(DXC_SIZE, sizeof(double));
  for(int i = 1; i <= isize+1; i++){
    DXC(i) = XC(i  )-XC(i-1);
  }
  return dxc;
}

