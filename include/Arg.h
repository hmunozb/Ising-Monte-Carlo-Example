/* 
 * File:   Arg.h
 * Author: hmunozbauza
 *
 * Created on June 14, 2015, 3:06 PM
 */

#ifndef ARG_H
#define	ARG_H

/* 
 * Arg.h
 * Some structures to build on optionparser to parse
 * the options passed to the program.
 */
#ifndef LIBCOMPILE

#include "optionparser.h"

#include <cstdlib>
#include <iostream>
using namespace std;


/*
 * Argument struct, based on example file for optionparser.h 
 */
struct Arg: public option::Arg
{
  static void printError(const char* msg1, const option::Option& opt, const char* msg2);

  static option::ArgStatus Unknown(const option::Option& option, bool msg);

  static option::ArgStatus Required(const option::Option& option, bool msg);

  static option::ArgStatus NonEmpty(const option::Option& option, bool msg);
  
  static option::ArgStatus Numeric(const option::Option& option, bool msg);
  
  static option::ArgStatus Directory(const option::Option& option, bool msg);
};

enum OptionIndex{
    HELP,
    UNKNOWN,
    MISC,
    N1,
    BETA,
    STEPS,
    SWEEPS,
    VERB
};
enum  OptionType{
    N_A
};

extern const option::Descriptor opts[];

#endif


#endif	/* ARG_H */

