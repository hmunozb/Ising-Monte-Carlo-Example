/* 
 * File:   RandGenerators.h
 * Author: hmunozbauza
 *
 * Created on April 22, 2016, 1:11 PM
 */

#ifndef RANDGENERATORS_H
#define	RANDGENERATORS_H

#include "randutils.hpp"
using RNG_ENGINE = mt19937_64;
using RAND_UL = uniform_int_distribution<unsigned long>;
using RAND_REAL = uniform_real_distribution<double>;

class RandomUnitRealGenerator{
public:
	RandomUnitRealGenerator( RNG_ENGINE& RNG) :
			 rng(RNG), rand_real(0.0, 1.0){ }
	double operator()(){
		return rand_real(rng);
	}
private:
	RNG_ENGINE& rng;
	RAND_REAL rand_real;
};

class RandomUlongGenerator{
public:
	RandomUlongGenerator( 
			unsigned long min,
			unsigned long max, 
			RNG_ENGINE& RNG) 
			: rng(RNG), Rand_UL(min, max){ }
	unsigned long operator()(){
		return Rand_UL(rng);
	}
	
private:
	RNG_ENGINE& rng;
	RAND_UL Rand_UL;
};
#endif	/* RANDGENERATORS_H */

