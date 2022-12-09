#if !defined(IB_H)
#define IB_H

#include "domain.h"
#include "structure.h"


typedef struct particle_t_ {
  // fixed parameters
  // radius
  double r;
  // gravity center,
  // translational (x, y)
  double x, y;
  double dx, dy;
  // translational and rotational velocities
  double ux, uy, vz;
  double dux, duy, dvz;
  // surface forces and torque at k+1/2 step
  // 1/m * \int_{S} a_i dS or 1/I * \int_{S} \epsilon_{ijk} \omega_j a_k dS
  // a_i^k = \alpha^k ( U_i^k - u_i^* ) / ( \gamma \Delta t )
  double fux, fuy, tvz;
  // internal inertia, e.g.,
  //   0: 1/C * \int_{V^{k  }} u_i^{k  } dV^{k  }
  //   1: 1/C * \int_{V^{k+1}} u_i^{k+1} dV^{k+1}
  // similar to the rotational component,
  //   where C is a pre-factor, mass or moment of inertia
  double iux[2], iuy[2], ivz[2];
  // collision force based on k (0) and k+1 (1) step particle positions and velocities
  double cfx[2], cfy[2];
} particle_t;

struct ib_t_ {
  int n_particles;
  particle_t **particles;
  // responses of surface forces and torque on the momentum fields
  double *dux, *duy;
};

extern double ib_s_weight(const double grid_size, const double r, const double px, const double py, const double x, const double y);
extern double ib_v_weight(const double grid_size, const double r, const double px, const double py, const double x, const double y);


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
