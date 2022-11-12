#include "common.h"
#include "sdecomp.h"
#include "temperature.h"


/**
 * @brief destruct a structure temperature_t
 * @param[inout] temperature : structure to be cleaned-up
 * @return                   : error code
 */
int temperature_finalise(temperature_t * restrict temperature){
  temperature_update_temp_finalise();
  common_free(temperature->temp);
  common_free(temperature->tempforcex);
  common_free(temperature->srctempa);
  common_free(temperature->srctempb);
  common_free(temperature->srctempg);
  common_free(temperature);
  return 0;
}

