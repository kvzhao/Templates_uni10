#ifndef __UNI10_ERROR_H__
#define __UNI10_ERROR_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#define uni10_error_msg(is_true, format, ...) {error_msg ( __PRETTY_FUNCTION__,  __FILE__, __LINE__, (is_true), (format), __VA_ARGS__);}           
static inline void error_msg( char const *const func, const char *const file, int const line, bool is_true, char const* format, ...){
  if (is_true)
  {
    va_list args;
    char msg[512];
    va_start(args, format);
    vsprintf(msg, format, args);
    fprintf(stderr, "\n# Uni10 error occur at %s\n# error: %s\n# file : %s (%d)\n\n", func, msg,file, line) ;
    va_end(args);
    exit(EXIT_FAILURE);
  }

}

#endif
