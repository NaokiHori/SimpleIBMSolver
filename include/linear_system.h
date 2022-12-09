#if !defined(LINEAR_SYSTEM_H)
#define LINEAR_SYSTEM_H

#include "sdecomp.h"
#include "domain.h"


/** @struct linear_system_t
 *  @brief structure storing buffers and plans to solve tri-diagonal linear systems in each dimension A x = b
 *  @var glsizes                   : GLOBAL sizes of the system to be solved
 *  @var sdecomp                   : structure storing MPI domain decomposition
 *  @var x1pncl, y1pncl            : buffers to store x1- and y1-pencils (in/out of system)
 *  @var tdm_[lcu]                 : coefficients of tri-diagonal matrices (shared among all directions)
 *  @var transposer_[xy]1_to_[yx]1 : plans to transpose between x1-pencil and y1-pencil
 */
typedef struct {
  int *glsizes;
  sdecomp_t *sdecomp;
  double *x1pncl, *y1pncl;
  double *tdm_l, *tdm_c, *tdm_u;
  sdecomp_transpose_t *transposer_x1_to_y1, *transposer_y1_to_x1;
} linear_system_t;

extern int linear_system_solve_in_x(linear_system_t * restrict linear_system);
extern int linear_system_solve_in_y(linear_system_t * restrict linear_system);


extern linear_system_t *init_linear_system(sdecomp_t * restrict sdecomp, const int glsizes[NDIMS]);
extern int linear_system_finalise(linear_system_t * restrict linear_system);

#endif // LINEAR_SYSTEM_H
