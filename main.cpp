/* 
 * File:   main.cpp
 * Author: hmunozbauza
 *
 * Created on December 18, 2015, 6:38 PM
 */

#include <iostream>
#include "ising.h"
#include "Log.h"
#include "ising_mc.h"

using namespace std;

/*
 * 
 */

INITIALIZE_EASYLOGGINGPP

unsigned long get_current_time(){
    return chrono::system_clock::now().time_since_epoch().count();
}

void init_EasyLogger(int argc, char** argv, string log_suffix = "") {
/*********Set up the Logging system********/

    //Filename for logfile
    //Get current time string
    char time_str[64];
    std::time_t t = std::time(NULL);
    if (std::strftime(time_str, sizeof (time_str), "%F_%H-%M-%S", std::localtime(&t))) {
        //You are spared from harsh judgment
    } else {
        strcpy(time_str, "T_ERR\0");
    }
    //Log file and log configurations
    string log_file = "./logs/ISING" + string(time_str) + log_suffix+".log";
    const string elpp_config
            = R"(* GLOBAL:
   FORMAT               =  "%datetime{%b %d, %Y  %H:%m:%s} [%level] %msg"
   FILENAME             =  ")" + log_file + R"("
   ENABLED              =  true
   TO_FILE              =  true
   TO_STANDARD_OUTPUT   =  true
   MILLISECONDS_WIDTH   =  3
   PERFORMANCE_TRACKING =  true
   MAX_LOG_FILE_SIZE    =  2097152 ## 2MB
   LOG_FLUSH_THRESHOLD  =  100 ## Flush after every 100 logs
* DEBUG:
   FORMAT               = "%datetime{%b %d, %Y  %H:%m:%s.%g} [%level]   %func %msg")";
    //EL startup
    START_EASYLOGGINGPP(argc, argv);
    el::Configurations lcnf;
    lcnf.setToDefault();
    lcnf.parseFromText(elpp_config);
    el::Loggers::reconfigureAllLoggers(lcnf);
    LOG(INFO) << "Info Logger Ready!";
    LOG(DEBUG) << "Debug Logger Ready!";
    /****************************************/
}


int main(int argc, char** argv) {
    init_EasyLogger(argc, argv);
    int l1 = 30, l2 = 60, l3 = 90,
	l4 = 20, l5 = 15,l6 = 10, l7 = 5;
    unsigned long swps = (2 << 16) - 1;

	//The temp array was calculated with the following:
	//auto beta_arr1 = ising_pt_temps(30, 0.3, 0.6, 40, .7);
	vector<double> beta_arr1 = 
	{   0.3,
         0.309791,
         0.319931,
         0.329522,
         0.338636,
         0.348023,
         0.357037,
         0.366331,
         0.374694,
         0.383238,
         0.391816,
         0.39997,
         0.407163,
         0.414324,
         0.420851,
         0.42717,
         0.432786,
         0.438486,
         0.444416,
         0.450974,
         0.458129,
         0.465803,
         0.473854,
         0.483133,
         0.492918,
         0.503714,
         0.515219,
         0.527376,
         0.540349,
         0.554813,
         0.570533,
         0.587131,
         0.60605,
         0.626146,
         0.64906,
         0.672977,
         0.700566,
         0.730375,
         0.764164,
         0.800962
	};
	
	ising_pt(l1, swps, beta_arr1, "out30_1.bin");
	ising_pt(l2, swps, beta_arr1, "out60_1.bin");
	ising_pt(l3, swps, beta_arr1, "out90_1.bin");
	
	vector<double> beta_arr2;
	double beta1 = 1.0/2.35, beta2 = 1.0/2.2; unsigned steps = 35;
	double dbeta = (beta2-beta1)/(steps-1);
	
	for(unsigned long i = 0; i < steps; ++i){
		beta_arr2.push_back(beta1 + i*dbeta);
	}
	ising_pt(l1, swps, beta_arr2, "out30_2.bin");
	ising_pt(l2, swps, beta_arr2, "out60_2.bin");
	ising_pt(l3, swps, beta_arr2, "out90_2.bin");
	
	ising_pt(l4, swps, beta_arr1, "out20_1.bin");
	ising_pt(l5, swps, beta_arr1, "out15_1.bin");
	ising_pt(l6, swps, beta_arr1, "out10_1.bin");
	ising_pt(l7, swps, beta_arr1, "out05_1.bin");
	
	ising_pt(l4, swps, beta_arr2, "out20_2.bin");
	ising_pt(l5, swps, beta_arr2, "out15_2.bin");
	ising_pt(l6, swps, beta_arr2, "out10_2.bin");
	ising_pt(l7, swps, beta_arr2, "out05_2.bin");
	
    return 0;
}

