#ifndef __UNI10_MEM_INFO_H__
#define __UNI10_MEM_INFO_H__

#include <stdio.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <sys/vmmeter.h>

#include "uni10/uni10_type.h"

#define REAL_MEMSIZE  get_real_memsize; 

inline void get_real_memsize(uni10_unit64& real_memsize){

  

}

#endif
