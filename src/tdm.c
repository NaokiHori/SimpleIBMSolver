#include <math.h>
#include <float.h>
#include <complex.h>
#include <fftw3.h>
#include "common.h"
#include "tdm.h"


// buffer to store updated upper-diagonal components,
//   which is needed to keep the original matrix given by the user
// this will be allocated when needed, and re-allocated
//   when nitems_buf is smaller than requested size
// N.B. NOT thread-safe
static int nitems_buf = 0;
static double *buf = NULL;

// buffer to store additional right-hand-side
//   arising in periodic systems (Sherman-Morrison formula)
// this will be allocated when needed, and re-allocated
//   when nitems_q1 is smaller than requested size
// N.B. NOT thread-safe
static int nitems_q1 = 0;
static double *q1 = NULL;

/**
 * @brief kernel function to solve a linear system whose right-hand-side is "double"
 * @param[in   ] n : matrix size
 * @param[in   ] l : lower  diagonal part
 * @param[in   ] c : center diagonal part
 * @param[in   ] u : upper  diagonal part
 * @param[inout] q : right-hand-side & answers
 * @return         : error code
 */
static int dgtsv(const int n, const double * restrict l, const double * restrict c, const double * restrict u, double * restrict q){
  /* ! divide first row by center-diagonal term ! 4 ! */
  // buf: static variable to store updated u
  //   (whose scope is inside this source)
  buf[0] = u[0]/c[0];
  q[0]   = q[0]/c[0];
  /* ! forward sweep ! 12 ! */
  for(int i = 1; i < n; i++){
    double val = c[i]-l[i]*buf[i-1];
    if(fabs(val) > DBL_EPSILON){
      // normal
      val = 1./val;
      buf[i] = val*(u[i]            );
      q[i]   = val*(q[i]-l[i]*q[i-1]);
    }else{
      // singular
      q[i] = 0.;
    }
  }
  /* ! backward sweep ! 3 ! */
  for(int i = n-2; i >= 0; i--){
    q[i] -= buf[i]*q[i+1];
  }
  return 0;
}

/**
 * @brief kernel function to solve a linear system whose right-hand-side is "fftw_complex"
 * @param[in   ] n : matrix size
 * @param[in   ] l : lower  diagonal part
 * @param[in   ] c : center diagonal part
 * @param[in   ] u : upper  diagonal part
 * @param[inout] q : right-hand-side & answers
 * @return         : error code
 */
static int zgtsv(const int n, const double * restrict l, const double * restrict c, const double * restrict u, fftw_complex * restrict q){
  /* ! divide first row by center-diagonal term ! 4 ! */
  // buf: static variable to store updated u
  //   (whose scope is inside this source)
  buf[0] = u[0]/c[0];
  q[0]   = q[0]/c[0];
  /* ! forward sweep ! 12 ! */
  for(int i = 1; i < n; i++){
    double val = c[i]-l[i]*buf[i-1];
    if(fabs(val) > DBL_EPSILON){
      // normal
      val = 1./val;
      buf[i] = val*(u[i]            );
      q[i]   = val*(q[i]-l[i]*q[i-1]);
    }else{
      // singular
      q[i] = 0.;
    }
  }
  /* ! backward sweep ! 3 ! */
  for(int i = n-2; i >= 0; i--){
    q[i] -= buf[i]*q[i+1];
  }
  return 0;
}

/**
 * @brief solve linear system whose right-hand-side is "double"
 * @param[in   ] n           : size of tri-diagonal matrix
 * @param[in   ] m           : how many right-hand-sides do you want to solve?
 * @param[in   ] is_periodic : periodic boundary condition is imposed
 *                               (Sherman-Morrison formula is used)
 *                             or not (normal Thomas algorithm is used)
 * @param[in   ] l           : pointer to lower- diagonal components
 * @param[in   ] c           : pointer to center-diagonal components
 * @param[in   ] u           : pointer to upper- diagonal components
 * @param[inout] q           : right-hand-sides (size: "n", repeat for "m" times) & answers
 *                               N.B. memory is contiguous in "n" direction, sparse in "m" direction
 * @return                   : error code
 */
int tdm_solve_double(const int n, const int m, const bool is_periodic, const double * restrict l, const double * restrict c, const double * restrict u, double * restrict q){
  // reallocate buffer to store updated "u"
  //   when the size is insufficient
  if(nitems_buf < n){
    if(buf != NULL){
      common_free(buf);
    }
    buf = common_calloc(n, sizeof(double));
    nitems_buf = n;
  }
  if(is_periodic){
    /* ! solve additional system arising from periodicity ! 17 ! */
    // q1: static variable to store the right-hand-side of Sherman-Morrison formula
    //   (whose scope is inside this source)
    // reallocate if size is insufficient
    if(nitems_q1 < n){
      if(q1 != NULL){
        common_free(q1);
      }
      q1 = common_calloc(n, sizeof(double));
      nitems_q1 = n;
    }
    for(int i = 0; i < n-1; i++){
      q1[i]
        = i ==   0 ? -l[i]
        : i == n-2 ? -u[i]
        : 0.;
    }
    dgtsv(n-1, l, c, u, q1);
    // for each "m" independent input
    for(int j = 0; j < m; j++){
      // pointer to the head of this system
      double *q0 = q + j * n;
      /* ! solve small system ! 1 ! */
      dgtsv(n-1, l, c, u, q0);
      /* ! find n-1-th component ! 3 ! */
      double num = q0[n-1]-u[n-1]*q0[0]-l[n-1]*q0[n-2];
      double den = c [n-1]+u[n-1]*q1[0]+l[n-1]*q1[n-2];
      q0[n-1] = fabs(den) < DBL_EPSILON ? 0. : num / den;
      /* ! reflect periodicity ! 3 ! */
      for(int i = 0; i < n-1; i++){
        q0[i] = q0[i]+q0[n-1]*q1[i];
      }
    }
  }else{
    // solve "m" independent systems
    for(int j = 0; j < m; j++){
      dgtsv(n, l, c, u, q + j * n);
    }
  }
  return 0;
}

/**
 * @brief solve linear system whose right-hand-side is "fftw_complex"
 * @param[in   ] n           : size of tri-diagonal matrix
 * @param[in   ] m           : how many right-hand-sides do you want to solve?
 * @param[in   ] is_periodic : periodic boundary condition is imposed
 *                               (Sherman-Morrison formula is used)
 *                             or not (normal Thomas algorithm is used)
 * @param[in   ] l           : pointer to lower- diagonal components
 * @param[in   ] c           : pointer to center-diagonal components
 * @param[in   ] u           : pointer to upper- diagonal components
 * @param[inout] q           : right-hand-sides (size: "n", repeat for "m" times) & answers
 *                               N.B. memory is contiguous in "n" direction, sparse in "m" direction
 * @return                   : error code
 */
int tdm_solve_fftw_complex(const int n, const int m, const bool is_periodic, const double * restrict l, const double * restrict c, const double * restrict u, fftw_complex * restrict q){
  // reallocate buffer to store updated "u"
  //   when the size is insufficient
  if(nitems_buf < n){
    if(buf != NULL){
      common_free(buf);
    }
    buf = common_calloc(n, sizeof(double));
    nitems_buf = n;
  }
  if(is_periodic){
    /* ! solve additional system arising from periodicity ! 17 ! */
    // q1: static variable to store the right-hand-side of Sherman-Morrison formula
    //   (whose scope is inside this source)
    // reallocate if size is insufficient
    if(nitems_q1 < n){
      if(q1 != NULL){
        common_free(q1);
      }
      q1 = common_calloc(n, sizeof(double));
      nitems_q1 = n;
    }
    for(int i = 0; i < n-1; i++){
      q1[i]
        = i ==   0 ? -l[i]
        : i == n-2 ? -u[i]
        : 0.;
    }
    dgtsv(n-1, l, c, u, q1);
    // for each "m" independent input
    for(int j = 0; j < m; j++){
      // pointer to the head of this system
      fftw_complex *q0 = q + j * n;
      /* ! solve small system ! 1 ! */
      zgtsv(n-1, l, c, u, q0);
      /* ! find n-1-th component ! 3 ! */
      fftw_complex num = q0[n-1]-u[n-1]*q0[0]-l[n-1]*q0[n-2];
      double       den = c [n-1]+u[n-1]*q1[0]+l[n-1]*q1[n-2];
      q0[n-1] = fabs(den) < DBL_EPSILON ? 0. : num / den;
      /* ! reflect periodicity ! 3 ! */
      for(int i = 0; i < n-1; i++){
        q0[i] = q0[i]+q0[n-1]*q1[i];
      }
    }
  }else{
    // solve "m" independent systems
    for(int j = 0; j < m; j++){
      zgtsv(n, l, c, u, q + j * n);
    }
  }
  return 0;
}

/**
 * @brief clean-up internal buffers
 * @return : error code
 */
int tdm_cleanup(void){
  common_free(buf);
  common_free(q1);
  buf = NULL;
  q1 = NULL;
  nitems_buf = 0;
  nitems_q1 = 0;
  return 0;
}

/*
 * Followed by test functions to confirm APIs
 *   1. tdm_solve_double
 *   2. tdm_solve_fftw_complex
 * are working fine
 */

#if defined(DEBUG_TEST_TDM)

#include <assert.h>
#include <limits.h>

#define M 6

int test1(void){
  // solve M tri-diagonal matrices
  //   d2p/dx2 = q, where p = j * sin(2 pi x), q = -4 pi^2 p
  // only n < USHRT_MAX (usually 65535) is considered,
  //   which is sufficiently large for fluid simulations
  // Dirichlet boundary conditions: p(0) = p(1) = 0
  FILE *fp = fopen("tdm1.log", "w");
  for(int iter = 0; ; iter++){
    const int n = (4 << iter) + 1;
    if(n > USHRT_MAX){
      break;
    }
    double *x   = common_calloc(n  , sizeof(double));
    double *res = common_calloc(n*M, sizeof(double));
    double *ans = common_calloc(n*M, sizeof(double));
    const double lx = 1.;
    const double dx = lx/(n-1);
    // init systems
    for(int i = 0; i < n; i++){
      x[i] = 1.*i*dx;
    }
    for(int j = 0; j < M; j++){
      for(int i = 0; i < n; i++){
        ans[j*n+i] =  sin(2.*M_PI*x[i])*j;
        res[j*n+i] = -pow(2.*M_PI, 2.)*ans[j*n+i];
      }
    }
    // init matrix
    double *l = common_calloc(n, sizeof(double));
    double *c = common_calloc(n, sizeof(double));
    double *u = common_calloc(n, sizeof(double));
    for(int i = 0; i < n; i++){
      if(i == 0 || i == n-1){
        l[i] = 0.;
        u[i] = 0.;
        c[i] = 1.;
      }else{
        l[i] = +1./dx/dx;
        u[i] = +1./dx/dx;
        c[i] = -2./dx/dx;
      }
    }
    // solve
    tdm_solve_double(n, M, false, l, c, u, res);
    // clean-up solver
    tdm_cleanup();
    double maxdif = 0.;
    for(int j = 0; j < M; j++){
      for(int i = 0; i < n; i++){
        maxdif = fmax(maxdif, fabs(res[j*n+i]-ans[j*n+i]));
      }
    }
    // resolution residual
    fprintf(fp, "%10d % .7e\n", n-1, maxdif);
    common_free(x);
    common_free(res);
    common_free(ans);
    common_free(l);
    common_free(c);
    common_free(u);
  }
  fclose(fp);
  return 0;
}

int test2(void){
  // solve M tri-diagonal matrices
  //   d2p/dx2 = q, where p = sin(2 pi x + \phi), q = -4 pi^2 p
  //   where \phi is argument 2 pi / M * j
  // only n < USHRT_MAX (usually 65535) is considered,
  //   which is sufficiently large for fluid simulations
  // Periodic boundary conditions: p(0) = p(N)
  FILE *fp = fopen("tdm2.log", "w");
  for(int iter = 0; ; iter++){
    const int n = 4 << iter;
    if(n > USHRT_MAX){
      break;
    }
    double *x   = common_calloc(n  , sizeof(double));
    double *res = common_calloc(n*M, sizeof(double));
    double *ans = common_calloc(n*M, sizeof(double));
    const double lx = 1.;
    const double dx = lx/n;
    // init systems
    for(int i = 0; i < n; i++){
      x[i] = 1.*i*dx;
    }
    for(int j = 0; j < M; j++){
      const double phi = 2.*M_PI/M*j;
      for(int i = 0; i < n; i++){
        ans[j*n+i] =  sin(2.*M_PI*x[i]+phi);
        res[j*n+i] = -pow(2.*M_PI, 2.)*ans[j*n+i];
      }
    }
    // init solver
    double *l = common_calloc(n, sizeof(double));
    double *c = common_calloc(n, sizeof(double));
    double *u = common_calloc(n, sizeof(double));
    for(int i = 0; i < n; i++){
      l[i] = +1./dx/dx;
      u[i] = +1./dx/dx;
      c[i] = -2./dx/dx;
    }
    // solve
    tdm_solve_double(n, M, true, l, c, u, res);
    // clean-up solver
    tdm_cleanup();
    // correct gauge (we know average leads 0)
    for(int j = 0; j < M; j++){
      double ave = 0.;
      for(int i = 0; i < n; i++){
        ave += res[j*n+i];
      }
      ave /= n;
      for(int i = 0; i < n; i++){
        res[j*n+i] -= ave;
      }
    }
    double maxdif = 0.;
    for(int j = 0; j < M; j++){
      for(int i = 0; i < n; i++){
        maxdif = fmax(maxdif, fabs(res[j*n+i]-ans[j*n+i]));
      }
    }
    // resolution residual
    fprintf(fp, "%10d % .7e\n", n, maxdif);
    common_free(x);
    common_free(res);
    common_free(ans);
    common_free(l);
    common_free(c);
    common_free(u);
  }
  fclose(fp);
  return 0;
}

int main(void){
  test1();
  test2();
  return 0;
}

#endif // DEBUG_TEST_TDM

