#include "common.h"
#include "fluid.h"


/**
 * @brief destruct a structure fluid_t
 * @param[inout] fluid : structure to be cleaned-up
 * @return             : error code
 */
int fluid_finalise(fluid_t * restrict fluid){
  common_free(fluid->ux);
  common_free(fluid->uy);
#if NDIMS == 3
  common_free(fluid->uz);
#endif
  common_free(fluid->p);
  common_free(fluid->psi);
  common_free(fluid->srcuxa);
  common_free(fluid->srcuxb);
  common_free(fluid->srcuxg);
  common_free(fluid->srcuya);
  common_free(fluid->srcuyb);
  common_free(fluid->srcuyg);
#if NDIMS == 3
  common_free(fluid->srcuza);
  common_free(fluid->srcuzb);
  common_free(fluid->srcuzg);
#endif
  common_free(fluid);
  return 0;
}

