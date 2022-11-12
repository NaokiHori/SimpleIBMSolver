#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <mpi.h>
#include "common.h"


const rkcoef_t RKCOEFS[RKSTEPMAX] = {
  {.alpha =  32./60., .beta  =   0./60., .gamma = 32./60.},
  {.alpha =  25./60., .beta  = -17./60., .gamma =  8./60.},
  {.alpha =  45./60., .beta  = -25./60., .gamma = 20./60.}
};

/**
 * @brief memory allocation with error handler
 * @param[in] count : number of elements
 * @param[in] size  : size of each element
 * @return          : pointer to the allocated buffer
 */
void *common_calloc(const size_t count, const size_t size){
  void *ptr = calloc(count, size);
  if(ptr == NULL){
    fprintf(stderr, "memory allocation error\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  return ptr;
}

/**
 * @brief memory deallocation with error handler
 * @param[in] ptr : pointer to the allocated buffer
 */
void common_free(void *ptr){
  free(ptr);
}

/**
 * @brief get current wall time (NOT simulation time units)
 * @return : current time (synchronised among all processes)
 */
double common_get_wtime(void){
  double now = MPI_Wtime();
  // share value of main process
  MPI_Bcast(&now, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  return now;
}

/**
 * @brief compute number of digits of the given integer
 * @param[in] num : a target integer whose number of digits are computed
 * @return        : number of digits
 */
int common_get_ndigits(int num){
  /*
   * E.g., num =    3 -> return 1
   * E.g., num =   13 -> return 2
   * E.g., num = 1234 -> return 4
   * N.B. negative values are not considered
   */
  assert(num >= 0);
  int retval = 1;
  while(num /= 10){
    retval++;
  }
  return retval;
}

