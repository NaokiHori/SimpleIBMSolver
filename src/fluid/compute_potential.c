#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "common.h"
#include "config.h"
#include "domain.h"
#include "sdecomp.h"
#include "fluid.h"
#include "tdm.h"
#include "arrays/fluid.h"


static int tell_fftw_plan_creation_failure(const int line){
  // function to just dump error message and abort
  printf("Line %d: FFTW plan creation failed\n", line);
  printf("A possible reason is you link Intel-MKL lib\n");
  printf("Make sure you use FFTW3, NOT its wrapper by MKL\n");
  MPI_Abort(MPI_COMM_WORLD, 0);
  return 0;
}


// structure only used in this source,
//   which contains variables to solve Poisson equation
// NOTE: this is shared among normal and efficient solvers
//   thus some variables may not be used
typedef struct {
  double *x1pncl_r, *y1pncl_r;
  fftw_complex *x1pncl_c, *y1pncl_c;
  fftw_plan fftw_plan_fwrd_x, fftw_plan_bwrd_x;
  fftw_plan fftw_plan_fwrd_y, fftw_plan_bwrd_y;
  double *tdm_l, *tdm_c, *tdm_u;
  double *eigenvalues_x, *eigenvalues_y;
  sdecomp_transpose_t *transposer_x1_to_y1_r, *transposer_y1_to_x1_r;
  sdecomp_transpose_t *transposer_x1_to_y1_c, *transposer_y1_to_x1_c;
} poisson_solver_t;

static poisson_solver_t *poisson_solver = NULL;

/**
 * @brief same as init_poisson_solver_normal, but more efficient (limited to uniform x grids)
 * @param[in] domain : information about domain decomposition and size
 * @return           : error code
 */
static int init_poisson_solver_efficient(const domain_t * restrict domain){
  poisson_solver = common_calloc(1, sizeof(poisson_solver_t));
  const sdecomp_t *sdecomp = domain->sdecomp;
  /* compute size of each pencil in each direction */
  const int glisize = domain->glsizes[0];
  const int gljsize = domain->glsizes[1];
  // x1 pencil, real
  const int x1pncl_isize_r = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_XDIR, glisize);
  const int x1pncl_jsize_r = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_X1PENCIL, SDECOMP_YDIR, gljsize);
  // y1 pencil, real
  const int y1pncl_isize_r = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_XDIR, glisize);
  const int y1pncl_jsize_r = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_YDIR, gljsize);
  // tri-diagonal matrix solver in y
  {
    poisson_solver->tdm_l = common_calloc(y1pncl_jsize_r, sizeof(double));
    poisson_solver->tdm_c = common_calloc(y1pncl_jsize_r, sizeof(double));
    poisson_solver->tdm_u = common_calloc(y1pncl_jsize_r, sizeof(double));
    /* set lower and upper diagonal components */
    // NOTE: since lower- and upper-diagonal components are
    //   independent to x direction and time,
    //   we compute here and re-use them
    const double dy = domain->dy;
    for(int j = 0; j < y1pncl_jsize_r; j++){
      poisson_solver->tdm_l[j] = 1./dy/dy;
      poisson_solver->tdm_u[j] = 1./dy/dy;
    }
  }
  /* parallel matrix transpose */
  {
    int sizes[NDIMS];
    // between x1 and y1 (real)
    sizes[0] = glisize;
    sizes[1] = gljsize;
    poisson_solver->transposer_x1_to_y1_r = sdecomp_transpose_fwrd_init(sdecomp, SDECOMP_X1PENCIL, sizes, sizeof(double), MPI_DOUBLE);
    poisson_solver->transposer_y1_to_x1_r = sdecomp_transpose_bwrd_init(sdecomp, SDECOMP_Y1PENCIL, sizes, sizeof(double), MPI_DOUBLE);
  }
  /* allocate pencils */
  {
    int sizes[NDIMS];
    // real x1 pencil
    sizes[0] = x1pncl_isize_r;
    sizes[1] = x1pncl_jsize_r;
    poisson_solver->x1pncl_r = fftw_malloc(sizes[0]*sizes[1]*sizeof(double));
    // real y1 pencil
    sizes[0] = y1pncl_isize_r;
    sizes[1] = y1pncl_jsize_r;
    poisson_solver->y1pncl_r = fftw_malloc(sizes[0]*sizes[1]*sizeof(double));
  }
  /* create x forward fftw plan */
  {
    // forward transform in x direction, real to real
    // information of each FFT
    const fftw_iodim dims[1] = {{
      .n  = glisize, // size of FFT: number of elements in x
      .is = 1, // distance between two input  elements: 1 (contiguous)
      .os = 1 // distance between two output elements: 1 (contiguous)
    }};
    // information of consecutive FFTs
    const fftw_iodim hdims[1] = {{
      .n  = x1pncl_jsize_r, // how many times do we want to repeat the same FFT?
      .is = x1pncl_isize_r, // distance between each FFT input
      .os = x1pncl_isize_r // distance between each FFT output
    }};
    // "the" DCT
    const fftw_r2r_kind kind[1] = {FFTW_REDFT10};
    // create plan
    poisson_solver->fftw_plan_fwrd_x = fftw_plan_guru_r2r(1, dims, 1, hdims, poisson_solver->x1pncl_r, poisson_solver->x1pncl_r, kind, FFTW_PATIENT);
    // error check
    if(poisson_solver->fftw_plan_fwrd_x == NULL){
      tell_fftw_plan_creation_failure(__LINE__);
    }
  }
  /* create x backward fftw plan */
  {
    // backward transform in x direction, real to real
    // information of each FFT
    const fftw_iodim dims[1] = {{
      .n  = glisize, // size of FFT: number of elements in x
      .is = 1, // distance between two input  elements: 1 (contiguous)
      .os = 1 // distance between two output elements: 1 (contiguous)
    }};
    // information of consecutive FFTs
    const fftw_iodim hdims[1] = {{
      .n  = x1pncl_jsize_r, // how many times do we want to repeat the same FFT?
      .is = x1pncl_isize_r, // distance between each FFT input
      .os = x1pncl_isize_r // distance between each FFT output
    }};
    // "the" iDCT
    const fftw_r2r_kind kind[1] = {FFTW_REDFT01};
    // create plan
    poisson_solver->fftw_plan_bwrd_x = fftw_plan_guru_r2r(1, dims, 1, hdims, poisson_solver->x1pncl_r, poisson_solver->x1pncl_r, kind, FFTW_PATIENT);
    // error check
    if(poisson_solver->fftw_plan_bwrd_x == NULL){
      tell_fftw_plan_creation_failure(__LINE__);
    }
  }
  /* compute eigenvalues in x */
  {
    const double dx = domain->dx;
    const int glisize = domain->glsizes[0];
    const int myisize = sdecomp_get_pencil_mysize(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_XDIR, glisize);
    const int ioffset = sdecomp_get_pencil_offset(sdecomp, SDECOMP_Y1PENCIL, SDECOMP_XDIR, glisize);
    poisson_solver->eigenvalues_x = common_calloc(myisize, sizeof(double));
    for(int local_i = 0; local_i < myisize; local_i++){
      const int global_i = local_i + ioffset;
      double val = -4./pow(dx, 2.)*pow(
          sin( M_PI * global_i / (2.*glisize) ),
          2.
      );
      poisson_solver->eigenvalues_x[local_i] = val;
    }
  }
  return 0;
}

#define X1PNCL_R(I, J) (x1pncl_r[((J)-1)*(isize)+((I)-1)])

/**
 * @brief compute scalar potential \psi to correct velocity
 * @param[in   ] domain : information about domain decomposition and size
 * @param[in   ] rkstep : Runge-Kutta step
 * @param[in   ] dt     : time step size
 * @param[inout] fluid  : velocity (in), scalar potential \psi (out)
 * @return              : error code
 */
int fluid_compute_potential(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid){
  // initalise solver, whose buffer "x1pncl_r" is necessary below
  if(poisson_solver == NULL){
    init_poisson_solver_efficient(domain);
  }
  /* ! compute right-hand-side ! 22 ! */
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const double gamma = RKCOEFS[rkstep].gamma;
    const double * restrict dxf = domain->dxf;
    const double            dy  = domain->dy;
    const double * restrict ux = fluid->ux;
    const double * restrict uy = fluid->uy;
    double * restrict x1pncl_r = poisson_solver->x1pncl_r;
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        double ux_xm = UX(i  , j  );
        double ux_xp = UX(i+1, j  );
        double uy_ym = UY(i  , j  );
        double uy_yp = UY(i  , j+1);
        X1PNCL_R(i, j) = 1./(gamma*dt)*(
           +(ux_xp-ux_xm)/DXF(i)
           +(uy_yp-uy_ym)/dy
        );
      }
    }
  }
  // efficient solver with DCT
  /* project to wave space (x) */
  fftw_execute(poisson_solver->fftw_plan_fwrd_x);
  /* transpose real x1pencil to y1pencil */
  sdecomp_transpose_execute(
      poisson_solver->transposer_x1_to_y1_r,
      poisson_solver->x1pncl_r,
      poisson_solver->y1pncl_r
  );
  /* solve linear systems in y */
  {
    const int glisize = domain->glsizes[0];
    const int gljsize = domain->glsizes[1];
    const int isize = sdecomp_get_pencil_mysize(domain->sdecomp, SDECOMP_Y1PENCIL, SDECOMP_XDIR, glisize);
    const int jsize = sdecomp_get_pencil_mysize(domain->sdecomp, SDECOMP_Y1PENCIL, SDECOMP_YDIR, gljsize);
    const double * restrict tdm_l = poisson_solver->tdm_l;
    const double * restrict tdm_u = poisson_solver->tdm_u;
    double * restrict tdm_c = poisson_solver->tdm_c;
    for(int i = 0; i < isize; i++){
      const double eval_x = poisson_solver->eigenvalues_x[i];
      /* set center diagonal components */
      for(int j = 0; j < jsize; j++){
        tdm_c[j] = -tdm_l[j]-tdm_u[j]+eval_x;
      }
      double * restrict rhs = poisson_solver->y1pncl_r + i * jsize;
      /* solve linear system */
      tdm_solve_double(
          /* size of system  */ jsize,
          /* number of rhs   */ 1,
          /* periodic        */ true,
          /* lower- diagonal */ tdm_l,
          /* center-diagonal */ tdm_c,
          /* upper- diagonal */ tdm_u,
          /* input / output  */ rhs
      );
      /* normalise FFT */
      for(int j = 0; j < jsize; j++){
        // for DCT, factor is doubled
        rhs[j] /= 2. * glisize;
      }
    }
  }
  /* transpose real y1pencil to x1pencil */
  sdecomp_transpose_execute(
      poisson_solver->transposer_y1_to_x1_r,
      poisson_solver->y1pncl_r,
      poisson_solver->x1pncl_r
  );
  /* project to physical space (x) */
  fftw_execute(poisson_solver->fftw_plan_bwrd_x);
  /* ! store result ! 13 ! */
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const double * restrict x1pncl_r = poisson_solver->x1pncl_r;
    double * restrict psi = fluid->psi;
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        PSI(i, j) = X1PNCL_R(i, j);
      }
    }
    // NOTE: psi and p are assumed to have the same array shape
    fluid_update_boundaries_p(domain, psi);
  }
  return 0;
}

#undef X1PNCL_R

/**
 * @brief destruct a structure poisson_solver_t
 * @return : error code
 */
int fluid_compute_potential_finalise(void){
  // use fftw_free since they are allocated using fftw_malloc
  fftw_free(poisson_solver->x1pncl_r);
  fftw_free(poisson_solver->y1pncl_r);
  fftw_free(poisson_solver->x1pncl_c);
  fftw_free(poisson_solver->y1pncl_c);
  fftw_destroy_plan(poisson_solver->fftw_plan_fwrd_y);
  fftw_destroy_plan(poisson_solver->fftw_plan_bwrd_y);
  fftw_cleanup();
  common_free(poisson_solver->eigenvalues_y);
  sdecomp_transpose_finalise(poisson_solver->transposer_x1_to_y1_r);
  sdecomp_transpose_finalise(poisson_solver->transposer_y1_to_x1_r);
  sdecomp_transpose_finalise(poisson_solver->transposer_x1_to_y1_c);
  sdecomp_transpose_finalise(poisson_solver->transposer_y1_to_x1_c);
  common_free(poisson_solver);
  return 0;
}

