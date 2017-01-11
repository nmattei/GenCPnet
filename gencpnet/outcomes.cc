/* 
 *  outcomes.cc: Generates XML files consisting of DT problems

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

// Note: DT = dominance testing.  We include this code to generate DT
// problems uniformly at random, not just CP-nets.

#include <iostream>
#include <fstream>
#include <random>
#include "outcomes.h"
#include "tables.h" // random_k_subset

extern std::mt19937_64 mt_rand; // random number generator
extern bool VERBOSE; // output details to clog when verbose mode is set

Outcomes::Outcomes( const unsigned int n,  // number of features
		    const unsigned int d,  // domain size (homogeneous)
		    const unsigned int h ) // Hamming distance (0 indicates unconstrained)
  : VARS(n), VALS(d), HD(h)
{
  init_pair();
  random_pair();
} // end constructor (HD specified)

Outcomes::Outcomes( const unsigned int n,  // number of features
		    const unsigned int d ) // domain size (homogeneous)
  : VARS(n), VALS(d), HD(0) // implicit HD=0 i.e. unconstrained
{
  init_pair();
  random_pair();
} // end constructor (HD unconstrained)

void Outcomes::init_pair()
{
  o1 = new unsigned int[VARS];
  o2 = new unsigned int[VARS];
} // end init_pair

void Outcomes::random_pair()
{
  std::uniform_int_distribution<int> randi(1, VALS);
  std::uniform_int_distribution<int> randflip(1, VALS - 1);
  
  if (HD == 0) // unconstrained Hamming distance
  {
    bool are_different = false;
    do {
      for (int i = 0; i < VARS; ++i)
      {
	o1[i] = randi(mt_rand);
	o2[i] = randi(mt_rand);
	if (o1[i] != o2[i])
	  are_different = true;
      } // end for 
    } while (!are_different);
  } // end if HD == 0
  else
  {
    // generate first outcome
    for (int i = 0; i < VARS; ++i)
      o1[i] = randi(mt_rand);

    // select a random subset of the variables to flip
    uint64_t T = (01ULL << VARS) - 1;
    uint64_t flip = random_k_subset(T, VARS, HD);

    // std::clog << std::hex << "0x" << T << std::dec << std::endl;
    // std::clog << std::hex << "0x" << flip << std::dec << std::endl;
    
    // flip the selected subset to different values
    for (int i = 0; i < VARS; ++i)
    {
      if (flip & 01ULL) // if rightmost bit is set
      {
	unsigned int newvalue = randflip(mt_rand);
	if (newvalue >= o1[i])
	  ++newvalue; // this ensures we get a different value
	o2[i] = newvalue;
      }
      else // copy from o1
      {
	o2[i] = o1[i];
	// std::clog << "o1[" << i << "]=" << o1[i]
	// 	  << ", o2[" << i << "]=" << o2[i] << '\n';
      } // end if rightmost bit is set
      flip >>= 1; // move on to the next bit and its variable
    } // end while flip set not empty
    
  } // end else HD > 0
} // end method random_pair

void Outcomes::dump()
{
  std::clog << "o1 = [";
  for (int i = 0; i < VARS; ++i)
    std::clog << o1[i] << ',';
  std::clog << "]; ";

  std::clog << "o2 = [";
  for (int i = 0; i < VARS; ++i)
    std::clog << o2[i] << ',';
  std::clog << "]\n";
} // end method dump

void Outcomes::xmlout( const std::string& working_directory, // folder
		       const std::string& cpnet_fname, // preference specification file
		       const std::string& dt_fname )   // preference query file
{
  // Open file--but first make sure a file with this name does not already exist!
  std::fstream xout (working_directory + '/' + dt_fname, std::ios_base::in);
  if (xout) { // file exists, report error and die!
    std::cerr << "Error: filename " << dt_fname << " already exists.\n"
	      << "Delete file(s) first or output to another directory.\n";
    exit(EXIT_FAILURE);
  } else { // a file with this name does not already exist, so attempt to create it
    xout.open(working_directory + '/' + dt_fname, std::ios_base::out);
    if (!xout.is_open()) // Verify file is open.  If not, perhaps the wrong directory?
    {
      std::cerr << "Error: Cannot open file "
		<< dt_fname << " for output.\n"
		<< "Check that directory " + working_directory
		<< " exists and that you have the necessary permissions.\n";
      exit(EXIT_FAILURE);
    } // end if unable to create file
  } // end if file already exists

  // Header
  // NOTE: Is the file naming convention correct?! Do we need the directory here? A full path?
  xout << "<PREFERENCE-QUERY>\n"
       << "  <PREFERENCE-SPECIFICATION-FILENAME>" << cpnet_fname << "</PREFERENCE-SPECIFICATION-FILENAME>\n"
       << "  <QUERY-TYPE>DOMINANCE</QUERY-TYPE>\n";

  // First outcome
  xout << "  <OUTCOME>\n"
       << "    <LABEL>BETTER</LABEL>\n";
  for (int i = 0; i < VARS; ++i)
  {
    xout << "    <ASSIGNMENT>\n"
	 << "      <PREFERENCE-VARIABLE>x" << i + 1 << "</PREFERENCE-VARIABLE>\n"
	 << "      <VALUATION>" << o1[i] << "</VALUATION>\n"
	 << "    </ASSIGNMENT>\n";
  } // end for
  xout << "  </OUTCOME>\n";

  // Second outcome
  xout << "  <OUTCOME>\n"
       << "    <LABEL>WORSE</LABEL>\n";
  for (int i = 0; i < VARS; ++i)
  {
    xout << "    <ASSIGNMENT>\n"
	 << "      <PREFERENCE-VARIABLE>x" << i + 1 << "</PREFERENCE-VARIABLE>\n"
	 << "      <VALUATION>" << o2[i] << "</VALUATION>\n"
	 << "    </ASSIGNMENT>\n";
  } // end for
  xout << "  </OUTCOME>\n";

  // Footer 
  xout << "</PREFERENCE-QUERY>\n";
  
  // Close the file
  xout.close(); 
} // end method xmlout

// End of file: outcomes.cc
