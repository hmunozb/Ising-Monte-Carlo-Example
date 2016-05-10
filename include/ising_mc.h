/* 
 * File:   ising_mc.h
 * Author: hmunozbauza
 *
 * Created on December 20, 2015, 8:06 AM
 */

#ifndef ISING_MC_H
#define	ISING_MC_H
#include "ising.h"
#include "MCMC.h"
#include "monte_carlo.h"
#include "topologies.h"

void ising_basic_mc(int ising_l, unsigned long sweeps, double beta);

void ising_pt(int ising_l, unsigned long sweeps, vector<double> beta_arr, string output_file);

vector<double> ising_pt_temps(unsigned int l,
		double beta_0, 
		double max_beta, 
		unsigned num_temps,
		double acceptance_goal = 0.5,
		unsigned long time_out = 20);

#endif	/* ISING_MC_H */

