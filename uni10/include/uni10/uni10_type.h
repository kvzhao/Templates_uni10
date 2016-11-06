#ifndef __UNI10_TYPE_H__
#define __UNI10_TYPE_H__

#include <stdio.h>
#include <string.h>

#include <complex>

typedef double uni10_double;
typedef size_t uni10_int;  // MKL uses long long int, not int64_t
typedef std::complex<double> uni10_complex;

typedef double* uni10_double_ptr;
typedef std::complex<double>* uni10_complex_ptr;


#endif
