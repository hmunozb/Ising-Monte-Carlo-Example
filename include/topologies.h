/* 
 * File:   topologies.h
 * Author: hmunozbauza
 *
 * Created on December 20, 2015, 7:56 AM
 */

#ifndef TOPOLOGIES_H
#define	TOPOLOGIES_H

#include <boost/graph/adjacency_list.hpp>
#include "_ising/ising_instance.h"
#include "randutils.hpp"

ising_graph<int> make_2D_ising(int l, int J);

ising_graph<int> make_2D_ising_periodic(int l, int J) ;

ising_graph<int> make_2D_ising_spin_glass(int l, int m, unsigned int sd);
#endif	/* TOPOLOGIES_H */

