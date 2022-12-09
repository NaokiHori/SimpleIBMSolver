#if !defined(DOMAIN_H)
#define DOMAIN_H

#include "arrays/domain.h"
#include "sdecomp.h"


/* ! definition of a structure domain_t ! 21 !*/
/** @struct domain_t
 *  @brief struct storing parameters relevant to spatial domain
 *  @var sdecomp  : MPI domain decomposition
 *  @var glsizes  : global     number of grid points in each direction
 *  @var mysizes  : local (my) number of grid points in each direction
 *  @var offsets  : offsets to my starting index in each direction
 *  @var lengths  : domain size in each direction
 *  @var xf, xc   : cell-face and cell-center locations in x direction
 *  @var dxf, dxc : face-to-face and center-to-center distances in x direction
 *  @var d[xy]    : representative grid sizes (dx is dummy unless uniform in x)
 */
typedef struct {
  sdecomp_t *sdecomp;
  int *glsizes;
  int *mysizes;
  int *offsets;
  double *lengths;
  double *xf, *xc;
  double *dxf, *dxc;
  double dx, dy;
} domain_t;


// constructor and destructor
extern domain_t *domain_init(void);
extern int domain_finalise(domain_t *domain);

// save members which will be used for post-processing
extern int domain_save(const char dirname[], const domain_t *domain);

// halo communication for 1. ux-like, 2. uy-like, 3. uz-like, and 4: p-like arrays
int domain_communicate_halo_ux_like(const domain_t *domain, double *ux_like);
int domain_communicate_halo_uy_like(const domain_t *domain, double *uy_like);
int domain_communicate_halo_p_like (const domain_t *domain, double  *p_like);

#endif // DOMAIN_H
