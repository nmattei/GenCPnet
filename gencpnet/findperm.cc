/* 
 *  findperm.cc: Maps integers to rankings

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

// If you are unfamiliar with factoradic number representations, we
// suggest starting with: 
// https://en.wikipedia.org/wiki/Factorial_number_system

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>
#include "findperm.h"
// #include <iostream>

using std::string;
using std::ostringstream;
using std::vector;

// Converts an integer in 1..d! to an XML formatted ranking
// Inputs:
//		num - The integer to be converted.
//		factorial - A value d such that num < d!
// Outputs:
//		Returns an XML string specifying the ranking
// Note:
//              Adapted from num_to_ranking_1 
string perm_num_to_xml(unsigned long int num, int factorial)
{
  if ( num < 1 || num > fac(factorial) )
    return string("");

  // Convert to form expected by other functions
  num -= 1;
  vector<int> perm = num_to_perm(num, factorial);
  for (int i = 0; i < perm.size(); ++i)
    perm[i] += 1;

  return perm_to_ranking_xml(perm);
}

// Convert a ranking from vector form to XML string form
// Inputs:
//		perm - A vector representing a ranking.
// Outputs:
//		Returns a string in the specified XML format
// Note:
//              Adapted from perm_to_ranking
string perm_to_ranking_xml(vector<int> perm)
{
  string ranking{};
  for (int i = 0; i < perm.size() - 1; ++i)
  {
    // std::clog << "i=" << i << ", perm=" << int_to_string(perm[i])
    // 	      << ", next-perm=" << int_to_string(perm[i+1]) << std::endl;
    ranking += "  <PREFERENCE>" + int_to_string(perm[i]) + ':'
      + int_to_string(perm[i + 1]) + "</PREFERENCE>\n";
  }
  return ranking;
}

// Convert an integer in 1..d! to a ranking on {1,2,,...,d}
// Inputs:
//		num - The integer to be converted.
//		factorial - A value d such that num < d!
// Outputs:
//		Returns a string with the ranking separated by colons.
string num_to_ranking_1(unsigned long int num, int factorial)
{
  if ( num < 1 || num > fac(factorial) )
    return string("");

  // Convert to form expected by other functions
  num -= 1;
  vector<int> perm = num_to_perm(num,factorial);
  for (int i = 0; i < perm.size(); i+=1)
  {
    perm[i] += 1;
  }

  return perm_to_ranking(perm);
}

// Convert an integer in 0..d!-1 to a ranking on {0,1,...,d-1}
// Inputs:
//		num - The integer to be converted.
//		factorial - A value d such that num < d!
// Outputs:
//		Returns a string with the ranking separated by colons.
string num_to_ranking_0(unsigned long int num, int factorial)
{
  if ( num < 0 || num >= fac(factorial) )
    return string("");

  return perm_to_ranking( num_to_perm(num,factorial) );
}

// Convert a ranking from vector form to string form, e.g. "2:4:1:3"
// Inputs:
//		perm - A vector representing a ranking.
// Outputs:
//		Returns a string with the ranking separated by colons.
string perm_to_ranking(vector<int> perm)
{
  if ( perm.size() == 0 )
    return string("");

  string ranking = string(int_to_string(perm[0]));
  for (int i = 1; i < perm.size(); i+=1)
  {
    ranking = ranking + ":" + string(int_to_string(perm[i]));
  }
  return ranking;
}

// Find a permutation encoded by a number.
// Returns a vector permutation of the digits 0..factorial-1
vector<int> num_to_perm(unsigned long int num, int factorial)
{
  return lehmer_to_perm(num_to_lehmer(num,factorial),factorial);
}

// Find the Lehmer code of a number.
// Returns a vector of Lehmer code digits, most significant first.
vector<unsigned long int> num_to_lehmer(unsigned long int num, int factorial)
{
  vector<unsigned long int> lehmer(factorial);
  unsigned long int quotient = num;
  for (unsigned long int divisor = 1; divisor <= factorial; divisor+=1)
  {
    lehmer[factorial-divisor] = (quotient % divisor);
    quotient = quotient / divisor;
  }

  return lehmer;
}

// Convert a Lehmer code to a permutation.
// Algorithm from Joel G. Silva Filho's "A Tutorial on Generating Permuations
//		using Factoradic Number Representation"
// Returns a vector permutation of the digits 0..factorial-1

vector<int> lehmer_to_perm(vector<unsigned long int> lehmer, int factorial)
{
  reverse(lehmer.begin(),lehmer.end()); // they do it backwards
  vector<int> S;
  for (int i = 0; i < factorial; i+=1)
  {
    S.push_back(i);
  }
  vector<int> perm(factorial);
	

  for (int i = 0; i < factorial; i+=1)
  {
    unsigned long int j = factorial-1-i;
    unsigned long int r = lehmer[j];
    perm[i] = S[r];
    for (j = r; j <=factorial-2; j+=1)
    {
      S[j] = S[j+1];
    }
  }

  return perm;
}

// Compute a factorial.
// Seriously, C++? You don't have this already?
unsigned long int fac(int num)
{
  if (num <= 1)
    return 1;

  return num * fac(num-1);
}

// Integer to string conversion.
// Seriously, C++? Again?
string int_to_string(int num)
{
  ostringstream strm;
  strm << num;
  return strm.str();
}

// End of file: findperm.cc
