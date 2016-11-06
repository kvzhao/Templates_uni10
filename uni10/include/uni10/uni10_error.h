#ifndef __UNI10_ERROR_H__
#define __UNI10_ERROR_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define uni10_error_msg(is_true, msg){error_msg ( (is_true), (msg), __PRETTY_FUNCTION__, __FILE__, __LINE__ );}           
static inline void error_msg(bool is_true, char const* msg, char const *const func, const char *const file, int const line)
{
  if (is_true)
  {
    fprintf(stderr, "\n# Uni10 error occur at %s\n# error: %s\n# file : %s (%d)\n\n", func, msg, file, line) ;
    exit(EXIT_FAILURE);
  }

}

#endif
