/*
   $Id$

   ADMB adaptations Copyright (c) 2009-2012 ADMB Foundation

   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)

   Modified for AD Model Builder by Kasper Kristensen <kkr@aqua.dtu.dk> and
   Anders Nielsen <an@aqua.dtu.dk> 2009
*/
/**
 * \file
 * Description not yet available.
 */
#include "fvar.h"
#ifdef DEBUG
  #include <cassert>
#endif

#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

/**
  \ingroup RNG
  Constructor for random_number_generator class.
  Based on the C-program for MT19937, originally coded by
  Takuji Nishimura and Makoto Matsumoto.
  \param seed Integer used to initialize the random number generator.
  Using different values of seed will generat different series of random numbers.
*/
random_number_generator::random_number_generator(const int seed)
{
#ifdef DEBUG
  assert(seed >= 0);
#endif
  unsigned long s = static_cast<unsigned long>(seed);
  mt=new unsigned long [N]; /* the array for the state vector  */
  mti=N+1; /* mti==N+1 means mt[N] is not initialized */

  mt[0]= s & 0xffffffffUL;
  for (mti=1; mti<N; mti++) {
    mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + static_cast<unsigned long>(mti));
    /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
    /* In the previous versions, MSBs of the seed affect   */
    /* only MSBs of the array mt[].                        */
    /* 2002/01/09 modified by Makoto Matsumoto             */
    mt[mti] &= 0xffffffffUL;
    /* for >32 bit machines */
  }
  better_rand();
}

/**
Destructor
*/
random_number_generator::~random_number_generator()
{
   delete [] mt;
   mt = 0;
}

/**
  \ingroup RNG
  Reinitialize random number seed.
  Based on the C-program for MT19937, originally coded by
  Takuji Nshimura and Makoto Matsumoto.
  \param seed Integer used to initialize the random number generator.
  Using different values of seed will generat different series of random numbers.
*/
void random_number_generator::reinitialize(int seed)
{
#ifdef DEBUG
  assert(seed >= 0);
#endif
  unsigned long s = static_cast<unsigned long>(seed);
  mt[0]= s & 0xffffffffUL;
  for (mti=1; mti<N; mti++) {
      mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + static_cast<unsigned long>(mti));
      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
      /* In the previous versions, MSBs of the seed affect   */
      /* only MSBs of the array mt[].                        */
      /* 2002/01/09 modified by Makoto Matsumoto             */
      mt[mti] &= 0xffffffffUL;
      /* for >32 bit machines */
  }
  better_rand();
}

/**
  \ingroup RNG
  Random number generator.
  Based on the Mersenne twister alorithm, MT19937, originally coded by
  Takuji Nishimura and Makoto Matsumoto.\n\n
  See Nishimura, T. and M. Matusomoto (1998) Mersenne twister:
  a 623-dimensionally equidistributed uniform pseudo-random number generator.
  ACM Transactions on Modeling and Computer Simulation (TOMACS) 8(1):3-30.

  \returns double containing uniformly distributed pseudorandom number
   between zero and one.
*/
double random_number_generator::better_rand()
{
  unsigned long y;
  static unsigned long mag01[2]={0x0UL, MATRIX_A};
  /* mag01[x] = x * MATRIX_A  for x=0,1 */

  if (mti >= N) { /* generate N words at one time */
      int kk = 0;

      //if (mti == N+1)   /* if init_genrand() has not been called, */
      //    init_genrand(5489UL); /* a default initial seed is used */

      for (;kk<N-M;kk++) {
          y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
          mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      for (;kk<N-1;kk++) {
          y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
          mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
      mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

      mti = 0;
  }

  y = mt[mti++];

  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  return (((double)y) + 0.5)*(1.0/4294967296.0);
}

#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK

/**
  \ingroup RNG
  Normal number generator.
  \returns N(0,1) double containing Normally distributed pseudorandom number
   with mean zero and standard deviation one.
*/
double randn(const random_number_generator& rng)
{
  double x=((random_number_generator&) rng).better_rand();
  double y=((random_number_generator&) rng).better_rand();
  double u=sqrt(-2*log(x))*cos(2*PI*y);
  return u;
}

/**
  \ingroup RNG
  Uniform random number generator.

  \returns double containing uniformly distributed pseudorandom number
   between zero and one.
*/
double randu(const random_number_generator& rng)
{
  return ((random_number_generator&)rng).better_rand();
}

/**
 * Description not yet available.
 * \param
 */
void dvector::fill_randbi(double p, const random_number_generator& rng)
{
  if ( p<0 || p>1)
  {
    cerr << "Error in dvar_vector::fill_randbi proportions of"
     " successes must lie between 0 and 1\n";
    ad_exit(1);
  }
  for (int i=indexmin(); i<=indexmax(); i++)
  {
    if ( ((random_number_generator&) rng).better_rand()<=p)
    {
      elem(i)=1;
    }
    else
    {
      elem(i)=0;
    }
  }
}

/**
 * Description not yet available.
 * \param
 */
void dvector::fill_randu(const random_number_generator& rng)
{
  for (int i=indexmin(); i<=indexmax(); i++)
  {
    elem(i)=((random_number_generator&) rng).better_rand();
  }
}

/**
 * Description not yet available.
 * \param
 */
void dmatrix::colfill_randu(const int&j, const random_number_generator& rng)
{
  for (int i=rowmin(); i<=rowmax(); i++)
  {
    elem(i,j)=((random_number_generator&) rng).better_rand();
  }
}

/**
 * Description not yet available.
 * \param
 */
void dmatrix::rowfill_randu(const int& i, const random_number_generator& rng)
{
  for (int j=colmin(); j<=colmax(); j++)
  {
    elem(i,j)=((random_number_generator&) rng).better_rand();
  }
}

/**
 * Description not yet available.
 * \param
 */
void dvector::fill_randn(const random_number_generator& rng)
{
  for (int i=indexmin(); i<=indexmax(); i++)
  {
    (*this)(i)=randn(rng);
  }
}

/**
 * Description not yet available.
 * \param
 */
void dmatrix::fill_randn(const random_number_generator& rng)
{
  for (int i=rowmin(); i<=rowmax(); i++)
  {
    elem(i).fill_randn(rng);
  }
}

/**
 * Description not yet available.
 * \param
 */
void d3_array::fill_randn(const random_number_generator& rng)
{
  for (int i=slicemin(); i<=slicemax(); i++)
  {
    elem(i).fill_randn(rng);
  }
}

/**
 * Description not yet available.
 * \param
 */
void d3_array::fill_randu(const random_number_generator& rng)
{
  for (int i=slicemin(); i<=slicemax(); i++)
  {
    elem(i).fill_randu(rng);
  }
}

/**
 * Description not yet available.
 * \param
 */
void dmatrix::fill_randu(const random_number_generator& rng)
{
  for (int i=rowmin(); i<=rowmax(); i++)
  {
    elem(i).fill_randu(rng);
  }
}

/**
 * Description not yet available.
 * \param
 */
void dmatrix::colfill_randn(const int&j, const random_number_generator& rng)
{
  for (int i=rowmin(); i<=rowmax(); i++)
  {
    elem(i,j)=randn(rng);
  }
}

/**
 * Description not yet available.
 * \param
 */
void dmatrix::rowfill_randn(const int& i, const random_number_generator& rng)
{
  for (int j=colmin(); j<=colmax(); j++)
  {
    elem(i,j)=randn(rng);
  }
}
