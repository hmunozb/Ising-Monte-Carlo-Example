/* 
 * File:   StatMechFunctionals.h
 * Author: hmunozbauza
 *
 * Created on April 22, 2016, 1:10 PM
 */

#ifndef STATMECHFUNCTIONALS_H
#define	STATMECHFUNCTIONALS_H

#include <map>
#include "MCMC.h"

using namespace std;

template< typename NumericType >
class BoltzmannTPF : public TransitionProbabilityFn< NumericType > {
	using TransitionProbabilityFn< NumericType >::it;
	using TransitionProbabilityFn< NumericType >::stored_results;
public:
	
	BoltzmannTPF(double beta_temp) : beta(beta_temp){ }

	inline double operator()(NumericType& delta_h) {
		it = stored_results.find(delta_h);
		if (it != stored_results.end()) {
			return it->second;
		} else {
			double d = exp(-beta * delta_h);
			stored_results[delta_h] = d;
			return d;
		}
	}
	inline double operator()(NumericType&& delta_h) {
		it = stored_results.find(delta_h);
		if (it != stored_results.end()) {
			return it->second;
		} else {
			double d = exp(-beta * delta_h);
			stored_results[delta_h] = d;
			return d;
		}
	}
	
	double get_beta(){ return beta;}
private:
	double beta;
};

/*
 * Generates the probability that an effective field in a temperature beta
 * will result in a positive ising spin
 */
template< typename NumericType >
class HeatbathTPF : public TransitionProbabilityFn< NumericType > {
	using TransitionProbabilityFn< NumericType >::it;
	using TransitionProbabilityFn< NumericType >::stored_results;
public:
	
	HeatbathTPF(double beta_temp) : beta(beta_temp) { }

	inline double operator()(NumericType fld) {
		it = stored_results.find(fld);
		if (it != stored_results.end()) {
			return it->second;
		} else {
			double d = 1.0 / (1 + exp(-2 * beta * fld));
			stored_results[fld] = d;
			return d;
		}
	}
	
private:
	double beta;
};

#endif	/* STATMECHFUNCTIONALS_H */

