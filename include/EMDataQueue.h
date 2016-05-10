/* 
 * File:   EMDataQueue.h
 * Author: hmunozbauza
 *
 * Created on April 23, 2016, 1:23 PM
 */

#ifndef EMDATAQUEUE_H
#define	EMDATAQUEUE_H

#include "readerwriterqueue.h"

using namespace std;
using namespace moodycamel;

enum EMDataCode{
	OK,
	PHASE_DONE,
	SIM_DONE
};

struct EMData {
	long replica_id;
	EMDataCode error_code;
	long E;
	long M;
	unsigned long acceptances;
	
};


typedef BlockingReaderWriterQueue<EMData> EMDataQueue;

class EMStatistics{
	double avg_e;
	double avg_e_sq;
	double avg_m;
	double avg_m_sq;
};

#endif	/* EMDATAQUEUE_H */
