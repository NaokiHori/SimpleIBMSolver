#include "common.h"
#include "ib.h"
#include "fileio.h"
#include "arrays/ib.h"



static ib_t *allocate(const domain_t *domain, const int n_particles){
  ib_t *ib = common_calloc(1, sizeof(ib_t));
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  ib->dux = common_calloc(DUX_SIZE_0 * DUX_SIZE_1, sizeof(double));
  ib->duy = common_calloc(DUY_SIZE_0 * DUY_SIZE_1, sizeof(double));
  ib->n_particles = n_particles;
  ib->particles = common_calloc(n_particles, sizeof(particle_t*));
  for(int n = 0; n < n_particles; n++){
    ib->particles[n] = common_calloc(1, sizeof(particle_t));
  }
  return ib;
}

static int init(ib_t *ib){
  const int n_particles = ib->n_particles;
  const char dirname[] = {"init_p"};
  // prepare buffers
  double *rs   = common_calloc(n_particles, sizeof(double));
  double *xs   = common_calloc(n_particles, sizeof(double));
  double *ys   = common_calloc(n_particles, sizeof(double));
  double *uxs  = common_calloc(n_particles, sizeof(double));
  double *uys  = common_calloc(n_particles, sizeof(double));
  double *vzs  = common_calloc(n_particles, sizeof(double));
  fileio_r_1d_serial(dirname, "p_rs",  NPYIO_DOUBLE, sizeof(double), n_particles,  rs);
  fileio_r_1d_serial(dirname, "p_xs",  NPYIO_DOUBLE, sizeof(double), n_particles,  xs);
  fileio_r_1d_serial(dirname, "p_ys",  NPYIO_DOUBLE, sizeof(double), n_particles,  ys);
  fileio_r_1d_serial(dirname, "p_uxs", NPYIO_DOUBLE, sizeof(double), n_particles, uxs);
  fileio_r_1d_serial(dirname, "p_uys", NPYIO_DOUBLE, sizeof(double), n_particles, uys);
  fileio_r_1d_serial(dirname, "p_vzs", NPYIO_DOUBLE, sizeof(double), n_particles, vzs);
  for(int n = 0; n < n_particles; n++){
    particle_t *p = ib->particles[n];
    p->r  = rs[n];
    p->x  = xs[n];
    p->y  = ys[n];
    p->ux = uxs[n];
    p->uy = uys[n];
    p->vz = vzs[n];
  }
  common_free(rs);
  common_free(xs);
  common_free(ys);
  common_free(uxs);
  common_free(uys);
  common_free(vzs);
  return 0;
}


ib_t *ib_init(const domain_t *domain){
  const char dirname[] = {"init_p"};
  int n_particles = 0;
  fileio_r_0d_serial(dirname, "n_particles", NPYIO_INT, sizeof(int), &n_particles);
  ib_t *ib = allocate(domain, n_particles);
  init(ib);
  return ib;
}

