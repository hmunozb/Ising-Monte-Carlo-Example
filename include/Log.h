/* 
 * File:   Log.h
 * Author: hmunozbauza
 *
 * Created on June 7, 2015, 11:21 PM
 */

#ifndef LOG_H
#define	LOG_H

#define ELPP_DEBUG_ASSERT_FAILURE
#define ELPP_NO_DEFAULT_LOG_FILE
#define ELPP_THREAD_SAFE
#if defined(__GNUC__)
    #define ELPP_STACKTRACE_ON_CRASH
#endif
#include "easylogging++.h"


#endif	/* LOG_H */

