#include <math.h>
#include "fluid.h"
#include "temperature.h"
#include "statistics.h"
#include "config.h"
#include "arrays/fluid.h"
#include "arrays/temperature.h"
#include "arrays/statistics.h"


static bool is_scheduled = false;
static double stat_rate;
double stat_next;


/**
 * @brief compute ux^1 and ux^2 and add results to the arrays
 * @param[in   ] domain    : information related to MPI domain decomposition
 * @param[in   ] fluid      : x velocity
 * @param[inout] statistics : arrays of the collected statistics
 * @return                  : error code
 */
static int collect_mean_ux(const domain_t *domain, const fluid_t *fluid, statistics_t *statistics){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *ux = fluid->ux;
  double *ux1 = statistics->ux1;
  double *ux2 = statistics->ux2;
  /* ! ux and its square are added ! 8 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize+1; i++){
        UX1(i, j, k) += pow(UX(i, j, k), 1.);
        UX2(i, j, k) += pow(UX(i, j, k), 2.);
      }
    }
  }
  return 0;
}

/**
 * @brief compute uy^1 and uy^2 and add results to the arrays
 * @param[in   ] domain    : information related to MPI domain decomposition
 * @param[in   ] fluid      : y velocity
 * @param[inout] statistics : arrays of the collected statistics
 * @return                  : error code
 */
static int collect_mean_uy(const domain_t *domain, const fluid_t *fluid, statistics_t *statistics){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *uy = fluid->uy;
  double *uy1 = statistics->uy1;
  double *uy2 = statistics->uy2;
  /* ! uy and its square are added ! 8 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize+1; i++){
        UY1(i, j, k) += pow(UY(i, j, k), 1.);
        UY2(i, j, k) += pow(UY(i, j, k), 2.);
      }
    }
  }
  return 0;
}

/**
 * @brief compute uz^1 and uz^2 and add results to the arrays
 * @param[in   ] domain    : information related to MPI domain decomposition
 * @param[in   ] fluid      : z velocity
 * @param[inout] statistics : arrays of the collected statistics
 * @return                  : error code
 */
static int collect_mean_uz(const domain_t *domain, const fluid_t *fluid, statistics_t *statistics){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *uz = fluid->uz;
  double *uz1 = statistics->uz1;
  double *uz2 = statistics->uz2;
  /* ! uz and its square are added ! 8 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize+1; i++){
        UZ1(i, j, k) += pow(UZ(i, j, k), 1.);
        UZ2(i, j, k) += pow(UZ(i, j, k), 2.);
      }
    }
  }
  return 0;
}

/**
 * @brief compute T^1 and T^2 and add results to the arrays
 * @param[in   ] domain     : information related to MPI domain decomposition
 * @param[in   ] temperature : temperature
 * @param[inout] statistics  : arrays of the collected statistics
 * @return                   : error code
 */
static int collect_mean_temp(const domain_t *domain, const temperature_t *temperature, statistics_t *statistics){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double *temp = temperature->temp;
  double *temp1 = statistics->temp1;
  double *temp2 = statistics->temp2;
  /* ! temp and its square are added ! 8 ! */
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize+1; i++){
        TEMP1(i, j, k) += pow(TEMP(i, j, k), 1.);
        TEMP2(i, j, k) += pow(TEMP(i, j, k), 2.);
      }
    }
  }
  return 0;
}

/**
 * @brief accumulate statistical data
 * @param[in   ] domain      : information related to MPI domain decomposition
 * @param[in   ] time        : current simulation time units
 * @param[in   ] fluid       : velocity
 * @param[in   ] temperature : temperature
 * @param[inout] statistics  : arrays of the collected statistics
 * @return                   : error code
 */
int statistics_collect(const domain_t *domain, const double time, const fluid_t *fluid, const temperature_t *temperature, statistics_t *statistics){
  if(!is_scheduled){
    // timing has not been scheduled yet
    // load from configuration
    stat_rate    = config.get_double("stat_rate");
    double after = config.get_double("stat_after");
    // find time to trigger next event
    stat_next = stat_rate * ceil( fmax(time, after) / stat_rate );
    is_scheduled = true;
  }else{
    /* ! collect temporally-averaged quantities ! 4 ! */
    collect_mean_ux(domain, fluid, statistics);
    collect_mean_uy(domain, fluid, statistics);
    collect_mean_uz(domain, fluid, statistics);
    collect_mean_temp(domain, temperature, statistics);
    /* ! number of samples is incremented ! 1 ! */
    statistics->num += 1;
    stat_next += stat_rate;
  }
  return 0;
}

