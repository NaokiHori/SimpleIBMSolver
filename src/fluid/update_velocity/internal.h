#if !defined(FLUID_UPDATE_VELOCITY_INTERNAL_H)
#define FLUID_UPDATE_VELOCITY_INTERNAL_H


extern int fluid_update_velocity_ux(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid);
extern int fluid_update_velocity_uy(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid);
extern int fluid_update_velocity_uz(const domain_t * restrict domain, const int rkstep, const double dt, fluid_t * restrict fluid);

extern int fluid_update_velocity_finalise_ux(void);
extern int fluid_update_velocity_finalise_uy(void);
extern int fluid_update_velocity_finalise_uz(void);


#endif // FLUID_UPDATE_VELOCITY_INTERNAL_H
