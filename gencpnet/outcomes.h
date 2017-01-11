/* 
 *  outcomes.cc: Interface to generate XML files with DT problems
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

#ifndef OUTCOMES_H
#define OUTCOMES_H

class Outcomes
{
 public:
  Outcomes( const unsigned int n,   // number of features
	    const unsigned int d,   // domain size (homogeneous)
	    const unsigned int h ); // Hamming distance (0 indicates unconstrained)

  Outcomes( const unsigned int n,   // number of features
	    const unsigned int d ); // domain size (homogeneous)

  ~Outcomes() { delete[] o1; delete[] o2; }; // destructor

  void xmlout( const std::string& working_directory, // folder
	       const std::string& cpnet_fname, // preference specification file
	       const std::string& dt_fname );  // preference query file
  
  unsigned int hd() { return HD; }; // returns HD

  void dump(); // output pair of outcomes for debugging purposes

private:
  const unsigned int VARS;
  const unsigned int VALS;
  const unsigned int HD;

  unsigned int* o1;
  unsigned int* o2;

  void init_pair();
  void random_pair();
}; // end class Outcomes;

#endif // OUTCOMES_H

// End of file: outcomes.h
