#ifndef __UNI10_ELEM_RNG_LAPACK_CPU_H__
#define __UNI10_ELEM_RNG_LAPACK_CPU_H__

#include "uni10/uni10_type.h"

#if defined(BOOST)

#include <boost/random.hpp>

// engines
typedef boost::random::mt19937                 uni10_mt19937;
typedef boost::random::mt19937_64              uni10_mt19937_64;
typedef boost::random::ranlux24                uni10_ranlux24;
typedef boost::random::ranlux48                uni10_ranlux48;
typedef boost::random::minstd_rand             uni10_minstd;
typedef boost::random::minstd_rand0            uni10_minstd0;

// distributions
typedef boost::random::uniform_real_distribution<double>    uni10_uniform_real;
typedef boost::random::uniform_int_distribution<int>        uni10_uniform_int;
typedef boost::random::normal_distribution<double>          uni10_normal;
typedef boost::random::student_t_distribution<double>       uni10_student_t;
typedef boost::random::lognormal_distribution<double>       uni10_lognormal;

#else

#include <random>

// engines
typedef std::mt19937                 uni10_mt19937;
typedef std::mt19937_64              uni10_mt19937_64;
typedef std::ranlux24                uni10_ranlux24;
typedef std::ranlux48                uni10_ranlux48;
typedef std::minstd_rand             uni10_minstd;
typedef std::minstd_rand             uni10_minstd0;

// distributions
typedef std::uniform_real_distribution<double>    uni10_uniform_real;
typedef std::uniform_int_distribution<int>        uni10_uniform_int;
typedef std::normal_distribution<double>          uni10_normal;
typedef std::student_t_distribution<double>       uni10_student_t;
typedef std::lognormal_distribution<double>       uni10_lognormal;

#endif

#define uni10_elem_rng(M, elemNum, eng, dis, up, dn, seed)\
do{ eng generator(seed);\
  dis distribution(up, dn);\
  for(uni10_uint64 i = 0; i < elemNum; i++)\
    M[i] = distribution(generator);}\
while(0);

#endif
