#if !defined(FILEIO_H)
#define FILEIO_H

#include <stdio.h>  // FILE, size_t

/* general file opener / closer */
extern FILE *fileio_fopen(const char * restrict path, const char * restrict mode);
extern int fileio_fclose(FILE *stream);

/* prepare directory to be stored */
extern char *fileio_concatenate(const char prefix[], const int step);
extern int fileio_mkdir(const char dirname[]);

/* for simple_npyio lib */
// datatypes, which are immersed in NPY files (argument: dtype)
// only two types are considered for now
// N.B. include single quotations
#define NPYIO_INT    "'int32'"
#define NPYIO_DOUBLE "'float64'"
// wrapper functions of simple_npyio header writer / reader
// scalar io (only by single process)
extern int fileio_r_0d_serial(const char dirname[], const char dsetname[], const char dtype[], const size_t size,       void *data);
extern int fileio_w_0d_serial(const char dirname[], const char dsetname[], const char dtype[], const size_t size, const void *data);
// vector io (only by single process)
extern int fileio_r_1d_serial(const char dirname[], const char dsetname[], const char dtype[], const size_t size, const size_t nitems,       void *data);
extern int fileio_w_1d_serial(const char dirname[], const char dsetname[], const char dtype[], const size_t size, const size_t nitems, const void *data);
// multi-dimensional array (by all processes)
extern int fileio_r_nd_parallel(const char dirname[], const char dsetname[], const int ndims, const int *array_of_sizes, const int *array_of_subsizes, const int *array_of_starts,       double *data);
extern int fileio_w_nd_parallel(const char dirname[], const char dsetname[], const int ndims, const int *array_of_sizes, const int *array_of_subsizes, const int *array_of_starts, const double *data);

#endif // FILEIO_H
