#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <mpi.h>
#include "common.h"
#include "fileio.h"
#include "simple_npyio.h"


/**
 * @brief allocate and initialise string having an error message
 * @param[in] fname   : name of file now being handled andan error is detected
 * @param[in] line    : line number in this script
 * @param[in] message : additional message to be attached
 * @return            : message
 */
static char *generate_error_message(const char fname[], const int line, const char message[]){
  /* format is sprintf("%s:%d, %s", filename, line, message) */
  size_t nchars[3] = {0};
  nchars[0] = strlen(fname);
  nchars[1] = common_get_ndigits(line);
  nchars[2] = strlen(message);
  /* allocate and copy contents */
  // number of all characters + separators (':' and ", ") and NUL, 4 additional characters in total
  char *retval = common_calloc(nchars[0]+nchars[1]+nchars[2]+4, sizeof(char));
  sprintf(retval, "%s:%d, %s", fname, line, message);
  return retval;
}

/**
 * @brief allocate and initialise string having a npy file name
 * @param[in] dirname  : name of directory, e.g., "output/save/step0000000000"
 * @param[in] dsetname : name of the dataset, e.g., "ux"
 * @return             : file name
 */
static char *generate_npy_filename(const char dirname[], const char dsetname[]){
  size_t nchars[2] = {0};
  nchars[0] = strlen(dirname);
  nchars[1] = strlen(dsetname);
  // allocate with "/" and ".npy" and NUL, 6 additional characters in total
  char *fname = common_calloc(nchars[0]+nchars[1]+6, sizeof(char));
  sprintf(fname, "%s/%s.npy", dirname, dsetname);
  return fname;
}

/**
 * @brief file opener
 * @param[in] path : pointer to the file name to be opened
 * @param[in] mode : file open mode
 * @return         : file pointer
 */
FILE *fileio_fopen(const char * restrict path, const char * restrict mode){
  FILE *stream = fopen(path, mode);
  if(stream == NULL){
    char *error_message = generate_error_message(__FILE__, __LINE__, path);
    perror(error_message);
    common_free(error_message);
  }
  return stream;
}

/**
 * @brief file closer
 * @param[in] stream : file pointer to be closed
 * @return           : error code
 */
int fileio_fclose(FILE *stream){
  int retval = fclose(stream);
  if(retval != 0){
    char *error_message = generate_error_message(__FILE__, __LINE__, "");
    perror(error_message);
    common_free(error_message);
    return 1;
  }
  stream = NULL;
  return 0;
}

/**
 * @brief create directory name from prefix and current time step
 * @param[in] prefix : e.g., "output/save/step"
 * @param[in] step   : current time step, e.g., 0
 * @return           : resulting directory name, e.g., "output/save/step0000000000"
 */
char *fileio_concatenate(const char prefix[], const int step){
  char *dirname = common_calloc(
      strlen(prefix)+NDIGITS_STEP+1, // + NUL
      sizeof(char)
  );
  sprintf(dirname, "%s%0*d", prefix, NDIGITS_STEP, step);
  return dirname;
}

/**
 * @brief create directory
 * @param[in] dirname : name of the directory to be created
 * @return            : error code
 */
int fileio_mkdir(const char dirname[]){
  // NOTE: call this function ONLY from the main process
  // NOTE: continue even if failed,
  //   since we want to override previous data (errorcode: EEXIST)
  if(mkdir(dirname, 0777) != 0){
    // failed to create directory
    char *error_message = generate_error_message(__FILE__, __LINE__, dirname);
    perror(error_message);
    common_free(error_message);
    return 1;
  }
  return 0;
}

/**
 * @brief wrapper function of simple_npyio_r_header with error handling
 * @param[in] fname            : name of the file from which data is loaded
 * @param[in] ndims            : number of dimensions of the data to be loaded
 * @param[in] shape            : sizes of the data in each dimension to be loaded
 * @param[in] dtype            : datatype of the data to be loaded
 * @param[in] is_fortran_order : memory contiguous direction, normally false
 * @return                     : size of npy header in byte
 */
static size_t fileio_r_npy_header(const char fname[], const size_t ndims, const size_t *shape, const char *dtype, const bool is_fortran_order){
  // success: return npy file header size
  // failure: return 0
  FILE *fp = fileio_fopen(fname, "r");
  if(fp == NULL){
    // failed to open file
    char *error_message = generate_error_message(__FILE__, __LINE__, fname);
    perror(error_message);
    common_free(error_message);
    return 0;
  }
  // load header, return header size when succeeded, return 0 otherwise
  size_t ndims_ = 0;
  size_t *shape_ = NULL;
  char *dtype_ = NULL;
  bool is_fortran_order_ = false;
  size_t header_size = simple_npyio_r_header(&ndims_, &shape_, &dtype_, &is_fortran_order_, fp);
  fclose(fp);
  // check arguments, return header size when all OK, return 0 otherwise
  bool success = true;
  // ndims
  if(ndims != ndims_){
    fprintf(stderr, "%s:%d npyio header read failed\n", __FILE__, __LINE__);
    fprintf(stderr, "ndims: %zu expected, %zu obtained\n", ndims, ndims_);
    success = false;
  }
  // shape (for each dimension)
  for(size_t n = 0; n < (ndims < ndims_ ? ndims : ndims_); n++){
    if(shape[n] != shape_[n]){
      fprintf(stderr, "%s:%d npyio header read failed\n", __FILE__, __LINE__);
      fprintf(stderr, "shape[%zu]: %zu expected, %zu obtained\n", n, shape[n], shape_[n]);
      success = false;
    }
  }
  // dtype
  if(strcmp(dtype, dtype_) != 0){
    fprintf(stderr, "%s:%d npyio header read failed\n", __FILE__, __LINE__);
    fprintf(stderr, "dtype: %s expected, %s obtained\n", dtype, dtype_);
    success = false;
  }
  // is_fortran_order
  if(is_fortran_order != is_fortran_order_){
    fprintf(stderr, "%s:%d npyio header read failed\n", __FILE__, __LINE__);
    fprintf(stderr, "is_fortran_order: %s expected, %s obtained\n", is_fortran_order ? "true" : "false", is_fortran_order_ ? "true" : "false");
    success = false;
  }
  common_free(shape_);
  common_free(dtype_);
  if(!success){
    header_size = 0;
  }
  return header_size;
}

/**
 * @brief wrapper function of simple_npyio_w_header with error handling
 * @param[in] fname            : name of the file to which data is dumped
 * @param[in] ndims            : number of dimensions of the data to be written
 * @param[in] shape            : sizes of the data in each dimension to be written
 * @param[in] dtype            : datatype of the data to be written
 * @param[in] is_fortran_order : memory contiguous direction, normally false
 * @return                     : size of npy header in byte
 */
static size_t fileio_w_npy_header(const char fname[], const size_t ndims, const size_t *shape, const char dtype[], const bool is_fortran_order){
  // success: return npy file header size
  // failure: return 0
  size_t header_size = 0;
  FILE *fp = fileio_fopen(fname, "w");
  if(fp != NULL){
    header_size = simple_npyio_w_header(ndims, shape, dtype, is_fortran_order, fp);
    if(header_size == 0){
      fprintf(stderr, "%s:%d npyio header write failed\n", __FILE__, __LINE__);
    }
    fclose(fp);
  }
  return header_size;
}

/**
 * @brief read / write data from / to a npy file
 * @param[in   ] is_read          : mode, read (true) or write (false)
 * @param[in   ] dirname          : name of directory in which a target npy file is contained
 * @param[in   ] dsetname         : name of dataset
 * @param[in   ] dtype            : datatype written in npy file, e.g., '<f8', 'float64'
 * @param[in   ] ndims            : number of dimensions of the dataset
 * @param[in   ] shape            : shape of the dataset
 * @param[in   ] is_fortran_order : memory order, normally false
 * @param[in   ] totalsize        : total size of the data
 * @param[inout] data             : pointer to the data, (in) when write, (out) when read
 */
static int kernel_serial(const bool is_read, const char dirname[], const char dsetname[], const char dtype[], const size_t ndims, const size_t *shape, const bool is_fortran_order, const size_t totalsize, void *data){
  // whether the whole io procedure is successful or not
  bool success = true;
  // create a string "xxx/yyy.npy"
  char *fname = generate_npy_filename(dirname, dsetname);
  // load (and sanitise) / write header
  size_t header_size = is_read
    ? fileio_r_npy_header(fname, ndims, shape, dtype, is_fortran_order)
    : fileio_w_npy_header(fname, ndims, shape, dtype, is_fortran_order);
  // load data
  if(header_size != 0){
    // open file with read-only mode when read,
    // while with append mode (since header is already there) when write
    FILE *fp = is_read
      ? fileio_fopen(fname, "r")
      : fileio_fopen(fname, "a");
    if(fp != NULL){
      // skip header, whose size is obtained above
      fseek(fp, header_size, SEEK_SET);
      // load / write data
      size_t retval = is_read
        ? fread (data, totalsize, 1, fp)
        : fwrite(data, totalsize, 1, fp);
      // fread / fwrite should return the number of loaded / dumped elements
      // since we load / dump 1 element whose size is totalsize,
      //   the number of elements is 1
      if(retval != 1){
        char *error_message = generate_error_message(__FILE__, __LINE__, fname);
        if(is_read){
          fprintf(stderr, "%s fread  failed\n", error_message);
        }else{
          fprintf(stderr, "%s fwrite failed\n", error_message);
        }
        common_free(error_message);
        success = false;
      }
      fclose(fp);
    }else{
      // failed to open file
      char *error_message = generate_error_message(__FILE__, __LINE__, fname);
      perror(error_message);
      common_free(error_message);
      success = false;
    }
  }else{
    // header size is 0, simple_npyio returns an error
    success = false;
  }
  // free heap, which was allocated by generate_npy_filename
  // abort program for read process, since the loaded data must be used
  common_free(fname);
  if(is_read && !success){
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  return 0;
}

/**
 * @brief read scalar (0d) data from a npy file, by one process
 * @param[in ] dirname          : name of directory in which a target npy file is contained
 * @param[in ] dsetname         : name of dataset
 * @param[in ] dtype            : datatype written in npy file, e.g., '<f8', 'float64'
 * @param[out] data             : pointer to the data to be loaded
 */
int fileio_r_0d_serial(const char dirname[], const char dsetname[], const char dtype[], const size_t size, void *data){
  const size_t ndims = 0;
  const size_t *shape = NULL;
  const size_t totalsize = size;
  kernel_serial( true, dirname, dsetname, dtype, ndims, shape, false, totalsize, (void *)data);
  return 0;
}

/**
 * @brief write scalar (0d) data to a npy file, by one process
 * @param[in] dirname          : name of directory in which a target npy file is contained
 * @param[in] dsetname         : name of dataset
 * @param[in] dtype            : datatype written in npy file
 * @param[in] data             : pointer to the data to be written
 */
int fileio_w_0d_serial(const char dirname[], const char dsetname[], const char dtype[], const size_t size, const void *data){
  const size_t ndims = 0;
  const size_t *shape = NULL;
  const size_t totalsize = size;
  kernel_serial(false, dirname, dsetname, dtype, ndims, shape, false, totalsize, (void *)data);
  return 0;
}

/**
 * @brief read vector (1d) data from a npy file, by one process
 * @param[in ] dirname          : name of directory in which a target npy file is contained
 * @param[in ] dsetname         : name of dataset
 * @param[in ] dtype            : datatype written in npy file
 * @param[in ] nitems           : number of items (i.e., vector length)
 * @param[out] data             : pointer to the data to be loaded
 */
int fileio_r_1d_serial(const char dirname[], const char dsetname[], const char dtype[], const size_t size, const size_t nitems, void *data){
  /* read 1d (vector) data from a npy file */
  const size_t ndims = 1;
  const size_t shape[1] = {nitems};
  const size_t totalsize = size * nitems;
  kernel_serial( true, dirname, dsetname, dtype, ndims, shape, false, totalsize, (void *)data);
  return 0;
}

/**
 * @brief write vector (1d) data to a npy file, by one process
 * @param[in] dirname          : name of directory in which a target npy file is contained
 * @param[in] dsetname         : name of dataset
 * @param[in] dtype            : datatype written in npy file
 * @param[in] nitems           : number of items (i.e., vector length)
 * @param[in] data             : pointer to the data to be written
 */
int fileio_w_1d_serial(const char dirname[], const char dsetname[], const char dtype[], const size_t size, const size_t nitems, const void *data){
  /* write 1d (vector) data to a npy file */
  const size_t ndims = 1;
  const size_t shape[1] = {nitems};
  const size_t totalsize = size * nitems;
  kernel_serial(false, dirname, dsetname, dtype, ndims, shape, false, totalsize, (void *)data);
  return 0;
}

/**
 * @brief read / write N-dimensional data from / to a npy file, by all processes
 * @param[in   ] is_read           : mode, read (true) or write (false)
 * @param[in   ] dirname           : name of directory in which a target npy file is contained
 * @param[in   ] dsetname          : name of dataset
 * @param[in   ] ndims             : number of dimensions of the array
 * @param[in   ] array_of_sizes    : global sizes   of the dataset
 * @param[in   ] array_of_subsizes : local  sizes   of the dataset
 * @param[in   ] array_of_starts   : local  offsets of the dataset
 * @param[inout] data              : pointer to the data to be loaded or written
 */
static int kernel_parallel(const bool is_read, const char dirname[], const char dsetname[], const int ndims, const int *array_of_sizes, const int *array_of_subsizes, const int *array_of_starts, double *data){
  // NOTE: to define subarray, data is limited to double
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  char *fname = generate_npy_filename(dirname, dsetname);
  // check header
  size_t header_size;
  if(myrank == 0){
    // values expected
    size_t *shape = common_calloc(ndims, sizeof(size_t));
    for(int dim = 0; dim < ndims; dim++){
      shape[dim] = array_of_sizes[dim];
    }
    const char dtype[] = {NPYIO_DOUBLE};
    const bool is_fortran_order = false;
    header_size = is_read
      ? fileio_r_npy_header(fname, ndims, shape, dtype, is_fortran_order)
      : fileio_w_npy_header(fname, ndims, shape, dtype, is_fortran_order);
    common_free(shape);
  }
  // abort when an error detected when the header is loaded
  MPI_Bcast(&header_size, sizeof(size_t)/sizeof(uint8_t), MPI_BYTE, 0, MPI_COMM_WORLD);
  if(header_size == 0){
    common_free(fname);
    return 1;
  }
  // open file
  MPI_File fh = NULL;
  {
    const int amode = is_read
      ? MPI_MODE_RDONLY
      : MPI_MODE_CREATE | MPI_MODE_RDWR;
    int mpi_error_code = MPI_File_open(MPI_COMM_WORLD, fname, amode, MPI_INFO_NULL, &fh);
    if(MPI_SUCCESS != mpi_error_code){
      char string[MPI_MAX_ERROR_STRING];
      int resultlen;
      MPI_Error_string(mpi_error_code, string, &resultlen);
      fprintf(stderr, "%s:%d %s\n", __FILE__, __LINE__, string);
      common_free(fname);
      return 1;
    }
  }
  // load / dump
  {
    // create data type and set file view
    const MPI_Offset disp = header_size;
    const MPI_Datatype etype = MPI_DOUBLE;
    MPI_Datatype filetype;
    MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_C, MPI_DOUBLE, &filetype);
    MPI_Type_commit(&filetype);
    MPI_File_set_view(fh, disp, etype, filetype, "native", MPI_INFO_NULL);
    // number of elements
    int count = 1;
    for(int n = 0; n < ndims; n++){
      count *= array_of_subsizes[n];
    }
    // load / dump
    if(is_read){
      MPI_File_read_all(fh, data, count, etype, MPI_STATUS_IGNORE);
    }else{
      MPI_File_write_all(fh, data, count, etype, MPI_STATUS_IGNORE);
    }
    // clean-up
    MPI_Type_free(&filetype);
  }
  // clean-up
  MPI_File_close(&fh);
  common_free(fname);
  return 0;
}

/**
 * @brief read N-dimensional data from a npy file, by all processes
 * @param[in ] dirname           : name of directory in which a target npy file is contained
 * @param[in ] dsetname          : name of dataset
 * @param[in ] ndims             : number of dimensions of the array
 * @param[in ] array_of_sizes    : global sizes   of the dataset
 * @param[in ] array_of_subsizes : local  sizes   of the dataset
 * @param[in ] array_of_starts   : local  offsets of the dataset
 * @param[out] data              : pointer to the data to be loaded
 */
int fileio_r_nd_parallel(const char dirname[], const char dsetname[], const int ndims, const int *array_of_sizes, const int *array_of_subsizes, const int *array_of_starts,       double *data){
  kernel_parallel( true, dirname, dsetname, ndims, array_of_sizes, array_of_subsizes, array_of_starts, (double *)data);
  return 0;
}

/**
 * @brief write N-dimensional data to a npy file, by all processes
 * @param[in] dirname           : name of directory in which a target npy file is contained
 * @param[in] dsetname          : name of dataset
 * @param[in] ndims             : number of dimensions of the array
 * @param[in] array_of_sizes    : global sizes   of the dataset
 * @param[in] array_of_subsizes : local  sizes   of the dataset
 * @param[in] array_of_starts   : local  offsets of the dataset
 * @param[in] data              : pointer to the data to be written
 */
int fileio_w_nd_parallel(const char dirname[], const char dsetname[], const int ndims, const int *array_of_sizes, const int *array_of_subsizes, const int *array_of_starts, const double *data){
  kernel_parallel(false, dirname, dsetname, ndims, array_of_sizes, array_of_subsizes, array_of_starts, (double *)data);
  return 0;
}

