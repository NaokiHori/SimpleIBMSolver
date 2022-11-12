#if !defined(CONFIG_H)
#define CONFIG_H

#include <stdbool.h>

typedef struct {
  // constructor
  int (*  load)(void);
  // destructor
  int (*unload)(void);
  // getter
  char * (*get_string)(const char envname[]);
  bool   (*get_bool  )(const char envname[]);
  int    (*get_int   )(const char envname[]);
  double (*get_double)(const char envname[]);
  // save information
  int (*save)(const char dirname[]);
} config_t;

extern const config_t config;

#endif // CONFIG_H
