/* 
 * File:   ising_model.h
 * Author: hmunozbauza
 *
 * Created on December 18, 2015, 11:49 PM
 */

#ifndef ISING_MODEL_H
#define	ISING_MODEL_H
#include "ising_instance.h"

using namespace std;
template< typename NumericType>
class ising_model {
public:
	/*constructs an ising model object from the given instance. Spin nodes take on
	 *+-spin_str values.
	 * 
	 * mem_eff_fields: 
	 *	Option to calculate Hamiltonian by maintaining a local field array.
	 *	That is, H = sum_i F_i s_i, where F depends on the external fields, bond strengths,
	 *	and neighboring spins.
	 *
	 */
	

	ising_model(	const ising_instance<NumericType>& ising_instance,
					const vector<NumericType>& initial_spin_arr,
					NumericType spin_str,
					const ising_fields<NumericType>& field_strs): 
			instance(ising_instance), 
			G(instance.graph()),
			spin_strength(spin_str), 
			external_fields(field_strs),
			spin_array(initial_spin_arr){
		
		//Normalize array to spin strength
		for (NumericType& n : spin_array) {
			n = (n > 0 ? spin_strength : -spin_strength);
		}
		//Clearly these must be the same size
		if (external_fields.size() != instance.size() ) {
			throw runtime_error("Ising size mismatch");
		}
		//These will clean up the state of the class
		clear_proposal();
		recalculate_Hamiltonian();
	}

	//Copy constructor. Simply refers to the same ising instance as "that"
	//References save space in tempering
	ising_model(const ising_model& that) 
			: ising_model(	that.instance, 
							that.spin_array,
							that.spin_strength, 
							that.external_fields){
		if(that.mem_eff_fields){
			memorize_ef_mode(true);
		}
	}

	//Default move constructor is fine.
	ising_model( ising_model&& that) = default;
	
    template<typename URNG>
    void RandomizeSpins(URNG& r){
        spin_array = RandomSpinVector(spin_array.size(), spin_strength, r);
        Magnetization = 0;
		clear_proposal();
		recalculate_Hamiltonian();
    }
	
	//Only calculates the change in hamiltonian from a spin flip, and caches the proposal
	inline NumericType ProposeFlip_H(const size_t spin) {
		NumericType local_change = calculateFlip(spin);
		set_proposal(spin, local_change);
		return local_change;
	}

	//Flips the spin and return the Hamiltonian change
	NumericType flip_delta_H(const size_t spin) {
		delta_h = 0;
		if (proposalExists) {
			//If this spin was previously proposed, simply update the system
			if (proposalFlip == spin) {
				NumericType old_spin = spin_array[spin];
				spin_array[spin] = -old_spin;
				delta_h = proposedChange;
				Hamiltonian += delta_h;
				Magnetization += -2 * old_spin;
			}
			clear_proposal();
		} else {
			//Not previously proposed: calculate the Hamiltonian change
			delta_h = calculateFlip(spin);

			NumericType old_spin = spin_array[spin];
			spin_array[spin] = -old_spin;

			Hamiltonian += delta_h;
			Magnetization += -2 * old_spin;
		}
		if (mem_eff_fields) {
			//MEF mode: recalculate neighboring fields due to flipping this spin
			adj_it it_begin, it_end;
			tie(it_begin,it_end) = instance.neighbors_of(spin);
			for(auto it = it_begin; it != it_end; ++it ) {
				effective_fields[*it] = calculateField(*it);
			}
		}
		return delta_h;
	}

	void suspend() {
		if (proposalExists) clear_proposal();
		hamiltonian_valid = false;
	}

	void suspended_flip(const size_t spin) {
		spin_array[spin] *= -1;
	}

	void unsuspend() {
		recalculate_Hamiltonian();
		if (mem_eff_fields) memorize_ef_mode(true); //recalculate fields	
	}

	inline NumericType get_H() const {
		return Hamiltonian;
	}

	inline NumericType get_M() const {
		return Magnetization;
	}

	inline size_t size() const {
		return instance.size();
	}

	inline NumericType calculateField(const size_t spin) {
		local_fld = 0;
		adj_it it_begin, it_end;
		tie(it_begin,it_end) = instance.neighbors_of(spin);
		for (auto it = it_begin; it != it_end; ++it) {
			local_fld += spin_array[*it] 
					* instance.get_edge(spin, *it);
		}
		local_fld += external_fields[spin];

		return local_fld;
	}

	NumericType memorizedField(const size_t spin) {
		return effective_fields[spin];
	}

	NumericType calculateFlip(const size_t spin) {
		if (mem_eff_fields) {
			local_fld = effective_fields[spin];
		} else {
			local_fld = calculateField(spin);
		}
		local_fld *= 2;
		local_fld *= spin_array[spin];
		return local_fld;
	}

	void memorize_ef_mode(bool mode) {
		mem_eff_fields = mode;

		if (mem_eff_fields) {
			size_t sz = size();
			effective_fields.resize(sz);
			for (size_t i = 0; i < sz; ++i) {
				effective_fields[i] = calculateField(i);
			}
		}
	}

	const ising_graph<NumericType>& access_graph() const {
		return instance;
	}

	NumericType get_spin(size_t n) const{
		return spin_array[n];
	}
	void recalculate_Hamiltonian() {
		Hamiltonian = 0;
		Magnetization = 0;
		size_t size = instance.size();
		
		for (size_t i = 0; i < size; ++i) {
			//sum bond strengths over i < j to avoid double counting
			adj_it it_begin, it_end;
			tie(it_begin,it_end) = instance.neighbors_of(i);
			
			auto it = upper_bound(it_begin, it_end, i);
			while (it != it_end) {
				NumericType jj;
				jj = instance.get_edge(i, *it);
				Hamiltonian +=
						spin_array[i] * spin_array[*it] * jj;
				++it;
			}
			Hamiltonian += external_fields[i] * spin_array[i];
			Magnetization += spin_array[i];
		}
		Hamiltonian = -Hamiltonian;
		hamiltonian_valid = true;
	}

	ising_fields< NumericType > get_ext_fields() const{
		return external_fields;
	}
	NumericType spin_magnitude() const{
		return spin_strength;
	}
	virtual ~ising_model() {

	}

private:
	typedef typename ising_graph<NumericType>::adjacency_iterator adj_it;
	
	const ising_instance<NumericType>& instance;
	const ising_graph<NumericType>& G;
	NumericType spin_strength;
	const ising_fields< NumericType >& external_fields;
	
	vector<NumericType> spin_array;
	
	vector< NumericType > effective_fields;
	bool mem_eff_fields = false;
	NumericType Hamiltonian;
	NumericType Magnetization;

	bool proposalExists;
	size_t proposalFlip;
	NumericType proposedChange;

	bool hamiltonian_valid = false;

	//used as temporary variables
	
	NumericType local_fld;
	NumericType delta_h;
	
	inline void set_proposal(const size_t spin, const NumericType delta_H) {
		proposalExists = true;
		proposalFlip = spin;
		proposedChange = delta_H;
	}

	inline void clear_proposal() {
		proposalExists = false;
		proposedChange = 0;
		proposalFlip = size();
	}


};


#endif	/* ISING_MODEL_H */

