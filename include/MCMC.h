/* 
 * File:   MCMC.h
 * Author: hmunozbauza
 *
 * Created on June 25, 2015, 3:27 PM
 */

#ifndef MCMC_H
#define	MCMC_H

#include <random>
#include <map>

using namespace std;
/*
 * Hast_Metr_Advance()
 * Generic advancer for the Hastings-Metropolis algorithm 
 * 
 * Explicit Template Parameters:
 * * MvType - The data type that represents a "move" in the HM algorithm. e.g. 
 *	A simple spin flip in an Ising spin glass entails a Move type of unsigned long,
 *	the vertex in the graph that is to be flipped. The Potts model would require
 *	a pair/struct with the vertex number and the new spin vector.
 * * CntRejects - Whether to count the number of rejected moves, by incrementing
 *	the contents of the rejections argument (see rejections).
 * 
 * Implicit Template Parameters:
 * * MvFn - see getMove, acceptFn - see acceptMove, 
 * * URNG - type of random number generator from, or similar to, <random> generator types
 * 
 * Arguments:
 * * steps - Number of steps to advance
 * * getMove - A void function object that takes, as arguments, a MvType& Move and
 *	a double& MvProb and shall assign respectively a new randomly generated move
 *	and its probability of acceptance.
 * * acceptFn - A void function object that will take a const MvType& Move and shall
 *	change any necessary object state to "accept" Move. This should also, if
 *	necessary, change any state so that the moves generated by getMove will
 *	remain legal to accept.
 * * g - A pseudorandom number generator. Shall be passable to a uniform_real_distribution
 *	to generate reals between 0 and 1.
 * * rejections - Optional. A pointer to a ulong that will be incremented every time
 *	a move is rejected, if the CntRejects parameter is true.
 */

template<typename MvType,
	bool CntRejects,
	typename MvFn,
	typename acceptFn,
	typename URNG>
inline void Hast_Metr_Advance(
		unsigned long steps,
		MvFn getMove,
		acceptFn acceptMove,
		URNG& g) {

	uniform_real_distribution<double> RandReal(0.0, 1.0);
	double Gm, MvProb;
	MvType Move;
	for (unsigned long int i = 0; i < steps; ++i) {
		getMove(Move, MvProb);
		if (MvProb - 1 >= 0) {
			acceptMove(Move);
		} else {
			Gm = RandReal(g);
			if (Gm < MvProb) acceptMove(Move);
		}
	}   
}


template<typename MvType,
	typename MvFn,
	typename acceptFn,
	typename randRealFn>
class HM_Advancer{
public:
    HM_Advancer(MvFn get_move_fn, acceptFn accept_move_fn, randRealFn get_rand_real_fn)
            : getMove(get_move_fn), acceptMove(accept_move_fn), 
			random_real( get_rand_real_fn ){
        
    }
    inline void operator()(const unsigned long& steps){
        for( unsigned long int i = 0; i < steps; ++i){
            getMove( Move, MvProb );
            if( MvProb  >= 1){
                acceptMove(Move);
            } else {
                if( random_real() < MvProb)
					acceptMove(Move);
            }
        }  
    }
    template< typename... MvArgTs >
    inline void operator()(const unsigned long& steps, MvArgTs... Mv_args ){
        for( unsigned long int i = 0; i < steps; ++i){
            getMove( Move, MvProb, Mv_args... );
            if( MvProb  >= 1){
                acceptMove(Move, Mv_args...);
            } else {
                if( random_real() < MvProb) 
					acceptMove(Move, Mv_args... );
            }
        }  
    }

private:
    MvFn getMove;
    acceptFn acceptMove;
    randRealFn random_real;
//Private temporaries
    double MvProb;
    MvType Move;
};

template
<	typename MvType,
	typename MvFn,
	typename acceptFn,
	typename randRealFn >
HM_Advancer< MvType, MvFn, acceptFn, randRealFn >
instantiate_HM(	
		MvFn get_move_fn,
		acceptFn accept_move_fn,
		randRealFn get_rand_real_fn){
    return HM_Advancer< MvType,  MvFn, acceptFn, randRealFn>( 
			get_move_fn, 
			accept_move_fn, 
			get_rand_real_fn );
}

template<typename NumericType>
class TransitionProbabilityFn {
protected:
	TransitionProbabilityFn() = default;
	map<NumericType, double> stored_results;
	typename map<NumericType, double>::iterator it;
};

#endif	/* MCMC_H */

