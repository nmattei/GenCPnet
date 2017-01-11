/* 
 *  degen_multi.h: Interface for degeneracy testing in multivalued
 *  CPTs

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

#ifndef DEGEN_MULTI_H
#define DEGEN_MULTI_H

bool degen_multi( unsigned int m, // number of inputs (indegree)
		  unsigned int nAsst, // number of rows
		  unsigned int* C ); // array of outputs defining function

void rand_cpt( unsigned long int cpt[], // Return value
	       const int k,              // Number of parents
	       const long unsigned int size_cpt ); // Number of rows = d^k

void print_cpt( const unsigned long int cpt[], // CPT to output
	        const unsigned long int size_cpt ); // Number of rows = d^k

#endif // DEGEN_MULTI_H

// End of file: degen_multi.h
