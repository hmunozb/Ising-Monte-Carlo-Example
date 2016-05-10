#include "ising_mc.h"
#include <thread>
#include <fstream>
#include <gsl/gsl_rstat.h>

void ising_basic_mc(int ising_l, unsigned long sweeps, double beta){

    randutils::random_generator<std::mt19937_64> rng;
    uniform_int_distribution<int> dist(0, 1);
    
    
    cout << ising_l << endl;
    cout << beta << endl;
    ising_graph<int> G = make_2D_ising_periodic(ising_l, 1);
    ising_instance<int> instance(G);
    //Random initial configuration
	vector<int> spins = RandomSpinVector(num_vertices(G), 1, rng.engine());
    //No fields
	ising_fields<int> fields(num_vertices(G), 0);
	//ising model object
    ising_model<int> ISM(instance, spins, 1, fields);

    LOG(INFO) << " Ising size: " << ising_l << '\n' 
	    << "Sweeps = " << sweeps << '\n'
	    << " Beta = " << beta << "\n";
    ofstream ofs("out.txt", fstream::binary);
    auto get_data = [&](EMDataQueue& Q){
		EMData em;
		EMDataCode C;
		unsigned long i = 0;
		while(true){
			if(!(Q.try_dequeue(em))){
				continue;
			}
			C = em.error_code;
			if(C == EMDataCode::SIM_DONE){
				break;
			}
			
			ofs << i <<',' 
					<< em.E << ',' 
					<< em.M << ',' 
					<< em.acceptances << "\n" ;
			++i;
		}
	};
    ISM.memorize_ef_mode(true);
	EMDataQueue Q;
	thread t1(get_data, std::ref(Q));
    Ising_HM(beta, sweeps, ISM, rng.engine(), Q); 
	t1.join();
    
}

vector<double> ising_pt_temps(unsigned int l,
		double beta_0, 
		double max_beta, 
		unsigned num_temps,
		double acceptance_goal,
		unsigned long time_out){
	unsigned N = l*l;
	vector<double> beta_arr;
	vector<double> new_beta_arr;
	vector<unsigned long> acceptances;
	
	double beta_step = (max_beta-beta_0)/(num_temps-1);
	for(int i = 0; i < num_temps; ++i){
		beta_arr.push_back(beta_0+beta_step*i);
		new_beta_arr.push_back(beta_arr[i]);
	}
	randutils::random_generator<std::mt19937_64> rng;
	ising_graph<int> G = make_2D_ising_periodic(l, 1);
	
	ising_instance<int> instance(G);
	vector<int> spins = RandomSpinVector(num_vertices(G), 1, rng.engine());
	ising_fields<int> fields(num_vertices(G), 0);
	ising_model<int> ISM(instance, spins, 1, fields);
	
	EMDataQueue dummy;
	vector<double> acc_ratios(num_temps-1);

	unsigned long sweeps = 1000;
	int n;
	for(unsigned long i = 0; i < time_out; ++i){
		acceptances = Ising_parallel_tempering_Threaded<int, false>(
			beta_arr,
			sweeps,
			50,
			ISM,
			dummy);
		
		for(unsigned long j = 1; j < num_temps; ++j){
			acc_ratios[j-1] = max<unsigned long>(1, acceptances[j-1])/(sweeps*acceptance_goal);
			new_beta_arr[j] = min(100.0*(1.0-j/(1000.0*num_temps)),
					new_beta_arr[j-1]
					+(beta_arr[j] - beta_arr[j-1])*acc_ratios[j-1]);
		}
		beta_arr = new_beta_arr;
		bool done = true;
		for(unsigned long j = 0; j < num_temps - 1;++j){
			if(acc_ratios[j] > 1.1 || acc_ratios[j] < .91)
				done = false;
		}
		if(done){
			n = i+1;
			break;
		}
	}
	cout << "Temperatures found after " << n << "steps\n";
	for( double b : beta_arr){
		cout << "\t " << b << ",\n";
	}
	
	return beta_arr;
	
	
	
}

void ising_pt(int ising_l, unsigned long sweeps, vector<double> beta_arr,string output_file){
    int n = ising_l*ising_l;
    
    vector<unsigned long> acc_counts;
	
     randutils::random_generator<std::mt19937_64> rng;
    ising_graph<int> G = make_2D_ising_periodic(ising_l, 1);
    
    ising_instance<int> instance(G);
    //Random initial configuration
	vector<int> spins = RandomSpinVector(num_vertices(G), 1, rng.engine());
    //No fields
	ising_fields<int> fields(num_vertices(G), 0);
	//ising model object
    ising_model<int> ISM(instance, spins, 1, fields);
    ofstream ofs(output_file, ofstream::binary);
	
	auto get_data = [&](EMDataQueue& Q){
		EMData em;
		EMDataCode C;
		unsigned long i = 0;
		std::streamsize lsize = sizeof(long);
		while(true){
			Q.wait_dequeue(em);
			C = em.error_code;
			if(C == EMDataCode::SIM_DONE){
				break;
			}
			
			ofs.write(reinterpret_cast<char*>(&i), lsize);
			ofs.write(reinterpret_cast<char*>(&(em.replica_id)),lsize);
			ofs.write(reinterpret_cast<char*>(&(em.E)),lsize);
			ofs.write(reinterpret_cast<char*>(&(em.M)),lsize);
			
			++i;
		}
	};
	EMDataQueue Q;
	thread t1(get_data, std::ref(Q));
	Ising_parallel_tempering_Threaded(
			beta_arr,
			sweeps,
			1000,
			ISM,
			//rng.engine(),
			Q);
    
	t1.join();

	ofs.close();	

}