/* 
 * File:   monte_carlo.h
 * Author: hmunozbauza
 *
 * Created on June 21, 2015, 3:18 AM
 */

#ifndef MONTE_CARLO_H
#define	MONTE_CARLO_H


#include "concurrentqueue.h"
#include "blockingconcurrentqueue.h"
#include "EMDataQueue.h"

#include <chrono>
#include <random>
#include <deque>
#include <tuple>
#include <list>
#include <thread>
#include <future>
#include <vector>
#include <functional>
#include <unordered_map>
#include <map>
#include <sstream>
#include <algorithm>

#include "readerwriterqueue.h"

#include "ising.h"
#include "Log.h"
#include "MCMC.h"
#include "RandGenerators.h"
#include "StatMechFunctionals.h"



using namespace std;


/*
 * Runs Hastings-Metropolis on a spin glass ISM for the indicated number of sweeps
 */

//Overload for advancing an rng engine
template< typename NumericType, bool use_queue = true >
void Ising_HM(
		double beta,
		unsigned long sweeps,
		ising_model< NumericType >& ISM,
		RNG_ENGINE& eng,
		EMDataQueue& Q,
		long replica_id = 0){ 
	//Data packet with energy, magnetization, and acceptance count
	EMData em;
	em.replica_id = replica_id;
	em.error_code = EMDataCode::OK;
	//random number generators
	RandomUnitRealGenerator rand_real( eng );
	RandomUlongGenerator rand_spin(0, ISM.size() - 1, eng);
	unsigned long size = ISM.size();
	
	//Transition probability function
	BoltzmannTPF<NumericType> transition_p(beta);
	//Move generator lambda: draw random spin and calculate its acceptance probability
	auto MoveFN = [&](size_t& spin, double& prob){
		spin = rand_spin();
		prob = transition_p( ISM.ProposeFlip_H(spin) );
	};
	//Move acceptor lambda: flip spin and increment acceptance count
	auto AcceptFN = [&](size_t spin) {
		ISM.flip_delta_H(spin);
		++(em.acceptances);
	};
	//Instantiate the HM advancer
	auto HAST_METR = instantiate_HM< size_t >(
			MoveFN,
			AcceptFN,
			rand_real);	
	
	//Run the required sweeps
	for(unsigned long i = 0; i < sweeps; ++i){
		em.acceptances = 0;
		
		HAST_METR(size);
		em.E = ISM.get_H();
		em.M = ISM.get_M();
		if(use_queue){
			if(!(Q.enqueue(em))){
				LOG(ERROR) << "Data queue full. Skipping one sweep of data.";
			}
		}
	}
	em.error_code = EMDataCode::SIM_DONE; //End of simulation
	if(use_queue){
		if(!(Q.enqueue(em))){
			LOG(ERROR) << "Seriously what?";
		}
	}
}



//Overload for using running without transfering data through a queue
template< typename NumericType >
void Ising_HM_Run_Only(double beta,
		unsigned long sweeps,
		ising_model< NumericType >& ISM,
		RNG_ENGINE& eng) {
	EMDataQueue dummy;
	return Ising_HM<NumericType, false>(beta, sweeps, ISM, eng, dummy);
}



//todo: add ability to seed
template<typename NumericType, bool use_queue = true>
vector<unsigned long> Ising_parallel_tempering_Threaded(const vector<double>& beta_arr,
		unsigned long MCS,
		unsigned long eql_sweeps,
		ising_model<NumericType>& ISM,
		//RNG_ENGINE& mt,
		EMDataQueue& Q) {
	EMDataQueue dummy;
	unsigned int ensemble_size = beta_arr.size();
		
	LOG(INFO) << "Ensemble size = " << ensemble_size;
	unsigned long size = ISM.size();
	LOG(INFO) << " Ising size = " << size;
	
	unsigned no_thrds = max(std::thread::hardware_concurrency(),(unsigned)1);
//PREPARATORY OBJECTS:

	//Array of ising model objects
	vector< ising_model<NumericType> > ism_arr; 
	ism_arr.reserve(ensemble_size);
	for (unsigned int i = 0; i < ensemble_size; ++i)
		ism_arr.push_back(ising_model<NumericType>(ISM));
	
	// Temperature n mapping to ising model ism_arr[ beta_to_ism[n] ]
	vector< size_t > beta_to_ism; 
	for (unsigned int i = 0; i < ensemble_size; ++i)
		//Initial mapping is straightforward
		beta_to_ism.push_back(i);
	
	//Transition probability calculators for each temperature
	vector< BoltzmannTPF< NumericType > > BTPV;
	BTPV.reserve( ensemble_size );
	for (unsigned int i = 0; i < ensemble_size; ++i) 
		BTPV.push_back( BoltzmannTPF<NumericType>( beta_arr[i] ) );
	
	//Assign acceptance counts vector passed in args
	vector<unsigned long> acc_counts(ensemble_size - 1, 0);
	
	//Vector of beta differences
	vector<double> beta_diffs;
	for (unsigned long i = 0; i < ensemble_size - 1; ++i) {
		beta_diffs.push_back(beta_arr[i + 1] - beta_arr[i]);
	}
	
	//Random Number Generators
	vector< std::mt19937_64 > rngs;
	for(unsigned i  = 0; i < ensemble_size + 1; ++i){
		rngs.push_back(mt19937_64(randutils::auto_seed_256{}.base()));
	}
	
	RAND_REAL RandUnitReal(0.0, 1.0);
	RAND_UL RandSpin(0, size - 1);
	
	//Initial Equilibration
	for( unsigned long i = 0; i < ensemble_size; ++i){
		Ising_HM<NumericType, false>(
				beta_arr[i], eql_sweeps, 
				ism_arr[ beta_to_ism[i] ], 
				rngs[ i ],
				dummy);
	}
	
	moodycamel::BlockingConcurrentQueue<unsigned long> in_q;
	moodycamel::BlockingConcurrentQueue<unsigned long> out_q;
	//Thread procedure
	auto sim_thread = [&](unsigned thrd_id) {
		unsigned long j;
		double Gm, MvProb;
		long spin;
		while(true){
			//Wait for assignment
			in_q.wait_dequeue(j);
			if(j == ensemble_size) return;	//PT is done

			//Monte carlo sweeps
			for(unsigned long i = 0; i < 5*size; ++i){
				spin = RandSpin(rngs[j]);
				MvProb = BTPV[j](ism_arr[beta_to_ism[j]].ProposeFlip_H(spin) );
				
				if(MvProb >= 1){
					ism_arr[ beta_to_ism[j] ].flip_delta_H( spin);
				} else {
					Gm = RandUnitReal(rngs[j]);
					if(Gm < MvProb)
						ism_arr[ beta_to_ism[j] ].flip_delta_H( spin);
				}
			}
			out_q.enqueue(j);
		}
	};
	
	//loop variables
	double dlt;
	size_t n;
	EMData em;
	unsigned long sims_done = 0;
	em.error_code = EMDataCode::OK;
	 
	//Start simulation worker threads
	vector<thread> thrds;
	for(unsigned i = 0; i < no_thrds; ++i) {
		thrds.push_back(thread(sim_thread, i));
	}
	
	LOG(INFO) << "Starting MCS loop";
	
	for( unsigned long i = 0; i < MCS; ++i){
		//Run sweep tasks for all temperatures
		//LOG(INFO) << "Iteration " << i;
		sims_done = 0;
		for( unsigned long j = 0; j < ensemble_size; ++j){
			//thrds[j] = thread(hm0, 10*size, j);
			in_q.enqueue(j);
		}
		unsigned long item;
		while(sims_done < ensemble_size){
			out_q.wait_dequeue(item);
			++sims_done;
		}
		if(use_queue){
			for(unsigned long j = 0; j < ensemble_size; ++j){
				em.replica_id = j;
				em.acceptances = 0;
				em.E = ism_arr[ beta_to_ism[j] ].get_H();
				em.M = ism_arr[ beta_to_ism[j] ].get_M();
				Q.enqueue(em);
			}
		}
		//LOG(INFO) << "Performing exchange move.";
		//For all temperature boundaries
		for (unsigned long j = 0; j < ensemble_size - 1; ++j) {
			dlt = beta_diffs[j]
					*(ism_arr[ beta_to_ism[j + 1] ].get_H() 
						- ism_arr[ beta_to_ism[j] ].get_H());
			
			if (dlt >= 0 || (RandUnitReal(rngs[ensemble_size]) < exp(dlt))) {
				//swap isms in temps j and j+1
				n = beta_to_ism[j];
				beta_to_ism[j] = beta_to_ism[j + 1];
				beta_to_ism[j + 1] = n;
				//Increase acceptance count
				acc_counts[j]++;
			}
		}
		
		if(i % (MCS/10) == 0){
			LOG(INFO) << double(i)/MCS*100.0 << "%";
		}
	}
	
	LOG(DEBUG) << "MCS Loop End";

	for(int i = 0; i < no_thrds; ++i) in_q.enqueue(ensemble_size);
	
	em.replica_id = ensemble_size;
	em.error_code = EMDataCode::SIM_DONE;
	Q.enqueue(em);
	
	for (unsigned long i = 0; i < ensemble_size - 1; ++i) {
		cout << beta_arr[i]<< ":\t\t" << ((double)(acc_counts[i])) / (MCS) << endl;
	}
	LOG(INFO) << "PT done.";
	
	for(unsigned i = 0; i < no_thrds; ++i) {
		thrds[i].join();
	}
	
	return acc_counts;
}

#endif	/* MONTE_CARLO_H */

