#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <mpi.h>
#include "sdecomp.h"


void *sdecomp_calloc(const size_t count, const size_t size){
  void *ptr = calloc(count, size);
  errno = 0;
  if(ptr == NULL){
    perror(__func__);
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  return ptr;
}

void sdecomp_free(void *ptr){
  free(ptr);
}

