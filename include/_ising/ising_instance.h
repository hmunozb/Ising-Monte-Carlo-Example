/* 
 * File:   ising_instance.h
 * Author: hmunozbauza
 *
 * Created on December 18, 2015, 7:07 PM
 */

#ifndef ISING_INSTANCE_H
#define	ISING_INSTANCE_H

#include <vector>
#include <random>
#include <boost/graph/adjacency_list.hpp>

using namespace std;
using namespace boost;


		
template<class BondNumeric>
using ising_graph = 
		adjacency_list<	setS,	//neighbors in a set
						vecS,	//vertices in a vector. *note: this means vertex_descriptor == integer type*
						undirectedS,	//Undirected graph
						no_property,
						property<edge_weight_t, BondNumeric, no_property> >; //edges have "weight" property for bond strength



//generate a random array of spins of length size
template<typename URNG>
vector<int> RandomSpinVector(size_t size, int magnitude, URNG& gen){
	vector<int> v;
	uniform_int_distribution<unsigned char> rand_int(0, 1);
	for(size_t i = 0; i < size; ++i){
		if(rand_int(gen)){
			v.push_back(magnitude);
		} else {
			v.push_back(-magnitude);
		}
	}
	return v;
}


using namespace std;

/*Sharable, read-only accessor class to a vector of ising fields */
template< typename FieldNumeric >
class ising_fields{
public:
	ising_fields( size_t size, FieldNumeric str) : fields(size, str){ }
	ising_fields( const vector<FieldNumeric>& field_strs) : fields(field_strs){ }
	
	inline FieldNumeric operator[](size_t n) const{
		return fields[n];
	}
	inline size_t size() const{
		return fields.size();
	}
private:
	vector<FieldNumeric> fields;
};

/*Sharable, read-only, and condensed accessor to the graph of an ising instance*/
template< typename NumericType >
class ising_instance{
public:
	explicit ising_instance(const ising_graph<NumericType>& graph) : 
			G(graph), 
			J(get(edge_weight, G))
			{ }
	
	inline size_t size() const{
		return num_vertices(G);
	}
	
	inline pair<
		typename ising_graph<NumericType>::adjacency_iterator ,
		typename ising_graph<NumericType>::adjacency_iterator  >
	neighbors_of( size_t i) const{
		return adjacent_vertices(vertex(i,G), G);
	}
	inline const NumericType& get_edge(size_t i, size_t j) const{
		std::pair<typename ising_graph<NumericType>::edge_descriptor, bool> e = edge(i, j, G);
		return J[e.first];
	}
	
	inline const ising_graph<NumericType>& graph() const{
		return G;
	}
private:
	
	ising_graph<NumericType> G;
	typename property_map<ising_graph<NumericType>, edge_weight_t>::type J;

	
};
#endif	/* ISING_INSTANCE_H */

