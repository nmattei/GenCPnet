/* 
 *  findperm.h: Interface for mapping integers to permutations

 *  Copyright (c) 2015 Cory Siler and Thomas E. Allen
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

#ifndef FINDPERM_H
#define FINDPERM_H

#include <string>
#include <vector>

// Generate a permutation from a number using the Lehmer code mapping.

std::string perm_num_to_xml(unsigned long int num, int factorial);
std::string perm_to_ranking_xml(std::vector<int> perm);

// Convert an integer in 1..d! to a ranking on {1,2,,...,d}
std::string num_to_ranking_1(unsigned long int num, int factorial);
// Convert an integer in 0..d!-1 to a ranking on {0,1,...,d-1}
std::string num_to_ranking_0(unsigned long int num, int factorial);

// Helper functions
std::string perm_to_ranking(std::vector<int> perm);
std::vector<int> num_to_perm(unsigned long int num, int factorial);
std::vector<unsigned long int> num_to_lehmer(unsigned long int num, int factorial);
std::vector<int> lehmer_to_perm(std::vector<unsigned long int> lehmer, int factorial);
unsigned long int fac(int num);
std::string int_to_string(int num);

#endif // FINDPERM_H

// End of file: findperm.h
