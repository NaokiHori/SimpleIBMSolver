#if !defined(IB_H)
#define IB_H

#include "domain.h"
#include "structure.h"


typedef struct particle_t_ {
  double r;
  double x, y, z;
  double dx, dy, dz;
  double ux, uy, uz;
  double vx, vy, vz;
  double dux, duy, duz;
  double dvx, dvy, dvz;
  double fux, fuy, fuz;
  double tvx, tvy, tvz;
  double iux[2], iuy[2], iuz[2];
  double ivx[2], ivy[2], ivz[2];
  double cfx[2], cfy[2], cfz[2];
} particle_t;

struct ib_t_ {
  int n_particles;
  particle_t **particles;
  double *dux, *duy, *duz;
};

extern double ib_s_weight(const double grid_size, const double r, const double px, const double py, const double pz, const double x, const double y, const double z);
extern double ib_v_weight(const double grid_size, const double r, const double px, const double py, const double pz, const double x, const double y, const double z);


/* constructor and destructor */
extern ib_t *ib_init(const domain_t *domain);
extern int ib_finalise(ib_t *ib);

/* called by main update routine */
extern int ib_reset_particle_increments(ib_t *ib);
extern int ib_compute_inertia(const domain_t *domain, const int cnstep, const fluid_t *fluid, ib_t *ib);
extern int ib_compute_collision_force(const domain_t *domain, const int cnstep, ib_t *ib);
extern int ib_exchange_momentum(const domain_t *domain, const fluid_t *fluid, const double dt, ib_t *ib);
extern int ib_increment_particles(const int rkstep, const double dt, ib_t *ib, double *residual);
extern int ib_update_momentum_fleid(const domain_t *domain, fluid_t *fluid, const ib_t *ib);
extern int ib_update_particles(const domain_t *domain, ib_t *ib);

/* other utility functions */
extern int ib_decide_loop_size(const int lbound, const int ubound, const double grid_size, const double radius, const double grav_center, int *min, int *max);
extern double ib_compute_volume(const double r);
extern double ib_compute_mass(const double den, const double r);
extern double ib_compute_moment_of_inertia(const double den, const double r);

extern int ib_save(const char dirname[], const ib_t *ib);

#endif // IB_H
