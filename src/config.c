#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <errno.h>
#include <mpi.h>
#include "common.h"
#include "fileio.h"
#include "config.h"


typedef struct dict_t_ {
  const char *key;
  char *value;
} dict_t;

static dict_t **dict = NULL;

#if NDIMS == 2

/* ! available environmental variables ! 26 ! */
static const char *keys[] = {
  // name of environmental variables
  "restart_sim",
  "restart_dir",
  "solve_temp",
  "add_buoyancy",
  "timemax",
  "wtimemax",
  "log_rate",
  "log_after",
  "save_rate",
  "save_after",
  "stat_rate",
  "stat_after",
  "lx",
  "ly",
  "glisize",
  "gljsize",
  "implicitx",
  "implicity",
  "coef_dt_adv",
  "coef_dt_dif",
  "Ra",
  "Pr"
};

#else // NDIMS == 3

/* ! available environmental variables ! 29 ! */
static const char *keys[] = {
  // name of environmental variables
  "restart_sim",
  "restart_dir",
  "solve_temp",
  "add_buoyancy",
  "timemax",
  "wtimemax",
  "log_rate",
  "log_after",
  "save_rate",
  "save_after",
  "stat_rate",
  "stat_after",
  "lx",
  "ly",
  "lz",
  "glisize",
  "gljsize",
  "glksize",
  "implicitx",
  "implicity",
  "implicitz",
  "coef_dt_adv",
  "coef_dt_dif",
  "Ra",
  "Pr"
};

#endif

static int get_nitems(void){
  // get number of keys,
  //   i.e., size of dictionary defined above
  return sizeof(keys)/sizeof(char *);
}

static int find_key_index(const char *key){
  // return index of the given key in the dictionary
  //   which is used to access the corresponding value
  const int nitems = get_nitems();
  for(int n = 0; n < nitems; n++){
    if(0 == strcmp(key, dict[n]->key)){
      return n;
    }
  }
  printf("unknown envname: %s\n", key);
  MPI_Abort(MPI_COMM_WORLD, 0);
  return -1;
}

/**
 * @brief load environmental variables and create dictionary to store them
 * @return : error code
 */
static int load(void){
  // allocate pointers to nitems key-value pairs
  const int nitems = get_nitems();
  dict = common_calloc(nitems, sizeof(dict_t *));
  for(int n = 0; n < nitems; n++){
    // allocate space for one key-value
    dict[n] = common_calloc(1, sizeof(dict_t));
    // assign ptr to key
    const char *key = keys[n];
    dict[n]->key = key;
    // copy value instead of assign ptr,
    //   since ptr to env var might be changed by system
    char *value = getenv(key);
    if(value != NULL){
      size_t nchar = strlen(value);
      dict[n]->value = common_calloc(nchar+1, sizeof(char));
      memcpy(dict[n]->value, value, sizeof(char)*nchar);
    }else{
      dict[n]->value = NULL;
    }
  }
  // print the final dictionary
  {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    const int maxchar = 18;
    for(int n = 0; n < nitems; n++){
      if(myrank == 0){
        const char *key = dict[n]->key;
        char *value = dict[n]->value;
        if(value == NULL){
          printf("#. ENV %*s is NOT found\n",     maxchar, key);
        }else{
          printf("#. ENV %*s is     found: %s\n", maxchar, key, value);
        }
      }
    }
  }
  return 0;
}

/**
 * @brief clean-up memories used by storing dictionary
 * @return : error code
 */
static int unload(void){
  const int nitems = get_nitems();
  for(int n = 0; n < nitems; n++){
    common_free(dict[n]->value);
    common_free(dict[n]);
  }
  common_free(dict);
  return 0;
}

/**
 * @brief load ENV and return it
 * @param[in] envname : name of the ENV
 * @return            : value
 */
static char *get_string(const char envname[]){
  /*
   * return "value" of key-value pair,
   *   where "key" is the given name of ENV
   */
  const int index = find_key_index(envname);
  char *value = dict[index]->value;
  if(value == NULL){
    printf("ERROR: %s is not given\n", envname);
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  return value;
}

/**
 * @brief load ENV and return it as a boolean value
 * @param[in] envname : name of the ENV
 * @return            : value
 */
static bool get_bool(const char envname[]){
  /*
   * return true  if true  is given
   * return false if false is given, or "envname" is not given
   */
  const int index = find_key_index(envname);
  const char *value = dict[index]->value;
  if(value == NULL){
    // not given: regard it as false
    return false;
  }
  if(0 == strcmp(value, "true")){
    return true;
  }
  if(0 == strcmp(value, "false")){
    return false;
  }
  printf("ERROR: %s cannot be interpreted as bool\n", value);
  MPI_Abort(MPI_COMM_WORLD, 0);
  // return something to avoid compiler warning
  return false;
}

/**
 * @brief load ENV and return it as an integer value
 * @param[in] envname : name of the ENV
 * @return            : value
 */
static int get_int(const char envname[]){
  const int index = find_key_index(envname);
  const char *value = dict[index]->value;
  if(value == NULL){
    printf("ERROR: %s is not given\n", envname);
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  // try to convert value to an integer
  errno = 0;
  long retval_ = strtol(value, NULL, 10);
  if(errno != 0 || INT_MIN > retval_ || INT_MAX < retval_){
    printf("ERROR: over/underflow is detected: %s\n", value);
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  // conert long to int
  return (int)retval_;
}

/**
 * @brief load ENV and return it as a double value
 * @param[in] envname : name of the ENV
 * @return            : value
 */
static double get_double(const char envname[]){
  const int index = find_key_index(envname);
  const char *value = dict[index]->value;
  if(value == NULL){
    printf("ERROR: %s is not given\n", envname);
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  errno = 0;
  double retval = strtod(value, NULL);
  if(errno != 0){
    printf("ERROR: over/underflow is detected: %s\n", value);
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  return retval;
}

/**
 * @brief save parameters
 * @param[in] envname : name of directory to which data is saved
 * @return            : error code
 */
static int save(const char dirname[]){
  const double Ra = get_double("Ra");
  const double Pr = get_double("Pr");
  fileio_w_0d_serial(dirname, "Ra", NPYIO_DOUBLE, sizeof(double), &Ra);
  fileio_w_0d_serial(dirname, "Pr", NPYIO_DOUBLE, sizeof(double), &Pr);
  return 0;
}

const config_t config = {
  .load   = load,
  .unload = unload,
  .get_string = get_string,
  .get_bool   = get_bool,
  .get_int    = get_int,
  .get_double = get_double,
  .save       = save
};

