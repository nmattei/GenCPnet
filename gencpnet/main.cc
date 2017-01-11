/* 
 *  main.cc: Command line interaction for GenCPnet

 *  Copyright (c) 2015 Thomas E. Allen
 *
 *  This file is part of GenCPnet, the Uniform CP-net Generator.
 *
 *  GenCPnet is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published
 *  by the Free Software Foundation, either version 3 of the License,
 *  or (at your option) any later version.
 *
 *  GenCPnet is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with GenCPnet.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "help.h"
#include "netcount.h" 
#include "outcomes.h"
#include "tables.h" 

// GLOBAL VARIABLES
unsigned int DOM_SIZE = 2; // domain size (homogeneous)
double DEG_INC = 0.00;  // degree of incompleteness (also homogeneous)
bool VERBOSE = false; // output extra details to clog when verbose mode is set
bool QUIET = false; // output few if any details to clog when quiet mode is set
bool COUNT_CPNETS = false; // only counts CP-nets; no generation
bool COUNT_DAGS = false; // only counte the dependency graphs; no generation

// Is there sufficient memory for this problem to be reasonable?
// This is just a quick check to avoid infeasible problems given system resources.
// In boundary cases, enough_memory will return true, but program will
// exhaust memory several minutes later--hopefully with error message, not crash.
bool enough_memory(int c, int d)
{
  // First check for truly ridiculous attempts that we can't even represent as integers
  double d_to_c_power = pow(d, c);
  if (pow(d,c) > std::numeric_limits<std::size_t>::max())
  {
    std::cerr << "Are you joking?\n"
	      << "No, there isn't enough memory for CPTs with "
	      << d << '^' << c << " = " << d_to_c_power << " entries!\n";
    return false;
  }
  
  bool result = true;
  size_t sizeCPT = size_t(d_to_c_power);
  volatile int** dummy = nullptr; // "volatile" tells compiler not to optimize this away

  try
  {
    // Allocate and promptly deallocate dummy table
    dummy = new volatile int*[sizeCPT];
    for (size_t r = 0; r < sizeCPT; ++r)
      dummy[r] = new volatile int[c+1];
    for (size_t r = 0; r < sizeCPT; ++r)
      delete[] dummy[r];
    delete[] dummy;
  }
  catch (std::bad_alloc& ba)
  {
    std::cerr << "Sorry, there is not enough memory for CPTs with "
	      << sizeCPT << " (" << d_to_c_power << ") entries\n";
    result = false;
  } // end try-catch

  return result;
}

int main(int argc, char** argv)
{
  // parse input parameters
  if (argc < 2) {
    std::cerr << "usage: " << argv[0] << " [OPTIONS] [DIRECTORY]: "
	      << "type " << argv[0] << " --help for details\n";
    return EXIT_FAILURE;
  } // end if

  if (argc >= 2) {
    std::string argv1 = argv[1];
    if (argv1 == "--help") {
      std::cout << HELPMESSAGE << std::endl;
      return EXIT_SUCCESS;
    } 
    else if (argv1 == "--version") {
      std::cout << VERSIONMESSAGE << std::endl;
      return EXIT_SUCCESS;
    }// end if argv1
    
  } // end if argc

  // default parameters
  int n = -1; // n and c are initially "undefined"
  int c = -1; // note that d defaults to 2 with DOM_SIZE=2 above
  int num_instances = 1; // defaults to a single instance
  int test_pairs = 0; // defaults to only a CP-net, not DT problems
  int HD = 0; // Hamming distance defaults to unconstrained ("0")
  std::string directory = "."; // defaults to current directory

  // This should work but is not very robust.  For example, if the
  // same parameter is entered twice, only the last one will be
  // applied, but no error message will be generated.
  {
    int i = 1; // indexes command line arguments
    for (;;) {
      if (i >= argc) break;

      std::string args = argv[i];
      if (args == "-n") { // number of nodes
	if (i <= argc) {
	  n = atoi(argv[i+1]);
	  i += 2;
	}
      } else if (args == "-c") { // bound on indegree
	if (i <= argc) {
	  c = atoi(argv[i+1]);
	  i += 2;
	}
      } else if (args == "--count") { // number of CP-nets
	COUNT_CPNETS = true;
	++i;
      } else if (args == "--countdags") { // number of DAGs
	COUNT_DAGS = true;
	++i;
      } else if (args == "-d") { // homogeneous domain size
	if (i <= argc) {
	  DOM_SIZE = atoi(argv[i+1]);
	  i += 2;
	}
      } else if (args == "-g") { // number of CP-nets to generate
	if (i <= argc) {
	  num_instances = atoi(argv[i+1]);
	  i += 2;
	}
      } else if (args == "-h") { // Hamming distance of DT instances (if -t used)
	if (i <= argc) {
	  HD = atoi(argv[i+1]);
	  i += 2;
	}
      } else if (args == "-i") { // degree of incompleteness
	if (i <= argc) {
	  DEG_INC = atof(argv[i+1]);
	  i += 2;
	}
      } else if (args == "--quiet" || args == "-q") { // for running in batch mode
	QUIET = true;
	++i;
      } else if (args == "-t") { // number of DT problems to generate
	if (i <= argc) {
	  test_pairs = atoi(argv[i+1]);
	  i += 2;
	}
      } else if (args == "--verbose" || args == "-V") { // debugging output
	VERBOSE = true;
	++i;
      } else if (directory == ".") { // where to output CP-nets and DT instances
	directory = argv[i++];
      } else {
	std::cerr << argv[0] << ": invalid option " << argv[i] << '\n'
		  << "Type " << argv[0] << " --help for details.\n";
	return EXIT_FAILURE;
      } // end if else block
    } // end infinite for loop
  }

  // Apply defaults for bound on indegree c
  if (c < 0) {
    c = (n > 6) ? 5 : n - 1;
  } else if (c >= n) {
    c = n - 1; // Surely this is what the user intended, right?
    // std::cerr << "Error: Bound on indegree c must be less than n.\n";
    // return EXIT_FAILURE;
  }

  // Check for error conditions in parameters
  if (n < 1) {
    std::cerr << "Error: Number of nodes n > 0 must be specified: e.g., -n 10.\n";
    return EXIT_FAILURE;
  } else if (n > 63) {
    std::cerr << "Error: Number of nodes n must be less than 64.\n";
    return EXIT_FAILURE;
  } else if (HD > n || HD < 0) {
    std::cerr << "Error: Hamming distance must be the range 0..n\n";
    return EXIT_FAILURE;
  }

  if (DEG_INC < 0.0 || DEG_INC >= 1.0) {
    std::cerr << "Error: degree of incompleteness must be in range [0.0, 1.0).\n";
    return EXIT_FAILURE;
  }

  // Output model specifications (parameters)
  if (!QUIET)
    std::clog << "Building distribution tables for CP-nets with the following specs:\n"
	      << "Number of nodes: n = " << n << std::endl
	      << "Bound on in-degree c = " << c << std::endl
	      << "Homogeneous domains of size d = " << DOM_SIZE << std::endl
	      << "Probability of incompleteness i = " << DEG_INC << std::endl;

  // Check for generation problem feasibility
  if (!enough_memory(c, DOM_SIZE))
  {
    std::cerr << "Aborting generation.\n";
    return EXIT_FAILURE;
  }

  // Initialize the CP-net counting object
  Netcount nc(n + 1, c + 1); 
  nc.init();
  for (int nn = 1; nn <= n; ++nn)
    for (int cc = 0; cc <= c; ++cc)
      nc.prob_cpnet(nn, cc);

  // If only counting the instances, not actually generating anything
  if (COUNT_DAGS)
  {
    if (QUIET)
      std::cout << nc.count_ldag(n, c);
    else      
      std::cout << "Number of DAGs: " << nc.count_ldag(n, c) << '\n';
  }

  if (COUNT_CPNETS)
  {
    if (QUIET)
      std::cout << nc.count_cpnet(n, c);
    else
      std::cout << "Number of CP-nets: " << nc.count_cpnet(n, c) << '\n';
  }

  if (COUNT_CPNETS || COUNT_DAGS)
    return EXIT_SUCCESS; // immediately terminate

  // Output model specifications (parameters)
  if (!QUIET)
    std::clog << "Generating " << num_instances << " random CP-nets with these specs.\n";

  // Init dagcode and CPTs
  uint64_t dc[n];
  unsigned long int* cpt[n];
  for (int i = 0; i < n; ++i)
    cpt[i] = nullptr;

  // File namer
  for (unsigned int counter = 0; counter < num_instances; ++counter)
  {
    // Generate dagcode
    nc.cdist->generate_random_cpnet(n, c, dc, cpt);

    // Output file count as a string 0000, 0001, etc.
    std::stringstream counts {};
    counts << std::setfill('0') << std::setw(4) << counter;

    // File name NOT including path
    std::string fname = "cpnet_n" + std::to_string(n) 
      + "c" + std::to_string(c)
      + "d" + std::to_string(DOM_SIZE)
      + (DEG_INC > 0.0 ? "i" + std::to_string(int(100*DEG_INC)) : "")
      + "_" + counts.str()
      + ".xml";
    if (VERBOSE) {
      std::clog << "Directory = \"" << directory << '"' << std::endl;
      std::clog << "Filename = \"" << fname << '"' << std::endl;
    } // end if

    // Details when in verbose mode
    if (VERBOSE) {
      std::clog << "Generating CP-net " << counter
		<< " (" << fname << ")\n";
      for (int i=0; i<n; ++i) {
	std::clog << i << ": ";
	if (i == 0)
	  print_cpt(cpt[i], 1);
	else {
	  int q = __builtin_popcount(dc[i]);
	  int r = int(pow(DOM_SIZE, q));
	  // std::clog << "d=" << DOM_SIZE << " q=" << q << " r=" << r << endl;
	  print_cpt(cpt[i], pow(DOM_SIZE, __builtin_popcount(dc[i])));
	} // end if
      } // end for 
    } // end if verbose

    // Output CP-net as its XML specification
    nc.cdist->dc_and_cpts_to_xml(n, dc, cpt, directory, fname);

    // Generate corresponding DT problems if specified
    if (test_pairs > 0)
    {
      if (!QUIET)
	std::clog << "Generating " << test_pairs
		  << " corresponding DT problems for each CP-net.\n";

      for (int k = 0; k < test_pairs; ++k)
      {
	// Output file count as a string 0000, 0001, etc.
	std::stringstream pair_counts {};
	pair_counts << std::setfill('0') << std::setw(4) << k;
	std::string ofname = "dt_n" + std::to_string(n)
	  + "c" + std::to_string(c)
	  + "d" + std::to_string(DOM_SIZE)
	  + (DEG_INC > 0.0 ? "i" + ((DEG_INC < 0.10 ? "0" : "") + std::to_string(int(100*DEG_INC))) : "")
	  + "_" + counts.str()
	  + '_' + pair_counts.str()
	  + ".xml";

	// Generate pair of outcomes
	Outcomes pair (n, DOM_SIZE, HD);
	
	if (VERBOSE)
	{
	  std::clog << ofname << "\t";
	  pair.dump();
	} // if verbose
	
	// Write to the XML file for DT problems
	pair.xmlout(directory, fname, ofname);
      } // end for
    } // end if test_pairs
    
    // Deallocate resources
    nc.cdist->free_cpt(n, cpt);
  } // for counter

  // Cleanly terminate program
  if (!QUIET)
    std::clog << "Generation complete." << std::endl;
  return EXIT_SUCCESS;
}

