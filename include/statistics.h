#if !defined(STATISTICS_H)
#define STATISTICS_H

#include "domain.h"
#include "fluid.h"
#include "temperature.h"


/* ! definition of a structure statistics_t_ ! 13 ! */
/** @struct statistics_t
 *  @brief struct storing statistics-related variables
 *  @var num          : number of samples which have been summed
 *  @var ux1, ux2     : mean and squared ux
 *  @var uy1, uy2     : mean and squared uy
 *  @var temp1, temp2 : mean and squared temp
 */
typedef struct {
  int num;
  double *ux1, *ux2;
  double *uy1, *uy2;
  double *temp1, *temp2;
} statistics_t;


// next time to trigger collect
extern double stat_next;

extern statistics_t *statistics_init(const domain_t *domain);
extern int statistics_finalise(statistics_t *statistics);

extern int statistics_collect(const domain_t *domain, const double time, const fluid_t *fluid, const temperature_t *temperature, statistics_t *statistics);
extern int statistics_output(const domain_t *domain, const int step, const double time, const statistics_t *statistics);

#endif // STATISTICS_H
