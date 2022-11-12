#include <mpi.h>
#include "sdecomp.h"


/**
 * @brief destruct a structure sdecomp_t
 * @param[in,out] sdecomp : structure to be cleaned up
 * @return                : error code
 */
int sdecomp_finalise(sdecomp_t *sdecomp){
  MPI_Comm_free( &(sdecomp->comm_cart) );
  sdecomp_free(sdecomp);
  return 0;
}

