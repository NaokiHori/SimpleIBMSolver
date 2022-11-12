#include "common.h"
#include "ib.h"
#include "fileio.h"


#if NDIMS == 2

int ib_save(const char dirname[], const ib_t *ib){
  const int n_p = ib->n_particles;
  double *rs = common_calloc(n_p, sizeof(double));
  double *xs = common_calloc(n_p, sizeof(double));
  double *ys = common_calloc(n_p, sizeof(double));
  for(int n = 0; n < n_p; n++){
    const particle_t *p = ib->particles[n];
    rs[n] = p->r;
    xs[n] = p->x;
    ys[n] = p->y;
  }
  fileio_w_1d_serial(dirname, "p_rs", NPYIO_DOUBLE, sizeof(double), n_p, rs);
  fileio_w_1d_serial(dirname, "p_xs", NPYIO_DOUBLE, sizeof(double), n_p, xs);
  fileio_w_1d_serial(dirname, "p_ys", NPYIO_DOUBLE, sizeof(double), n_p, ys);
  common_free(rs);
  common_free(xs);
  common_free(ys);
  return 0;
}

#else // NDIMS == 3

int ib_save(const char dirname[], const ib_t *ib){
  const int n_p = ib->n_particles;
  double *rs = common_calloc(n_p, sizeof(double));
  double *xs = common_calloc(n_p, sizeof(double));
  double *ys = common_calloc(n_p, sizeof(double));
  double *zs = common_calloc(n_p, sizeof(double));
  for(int n = 0; n < n_p; n++){
    const particle_t *p = ib->particles[n];
    rs[n] = p->r;
    xs[n] = p->x;
    ys[n] = p->y;
    zs[n] = p->z;
  }
  fileio_w_1d_serial(dirname, "p_rs", NPYIO_DOUBLE, sizeof(double), n_p, rs);
  fileio_w_1d_serial(dirname, "p_xs", NPYIO_DOUBLE, sizeof(double), n_p, xs);
  fileio_w_1d_serial(dirname, "p_ys", NPYIO_DOUBLE, sizeof(double), n_p, ys);
  fileio_w_1d_serial(dirname, "p_zs", NPYIO_DOUBLE, sizeof(double), n_p, zs);
  common_free(rs);
  common_free(xs);
  common_free(ys);
  common_free(zs);
  return 0;
}

#endif // NDIMS
