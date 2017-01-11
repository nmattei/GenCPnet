/* 
 *  degen_multi.cc: Tests for degeneracy in multivalued CPTs

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

#include <iostream> // print_cpt
#include <cmath> // for pow
#include <random> // Mersenne Twister
#include <functional> // We need bind for MT
#include <new> // Possibly needed for bad_alloc
#include "degen_multi.h"
#include "tables.h"
#include "findperm.h"

extern unsigned int DOM_SIZE; // domain size (homogeneous)
extern double DEG_INC;  // degree of incompleteness (also homogeneous)
extern std::mt19937_64 mt_rand; // random number generator

bool degen_multi( unsigned int m, // number of inputs (indegree)
  		 // unsigned int d, // size of domains (homogeneous)  
		 unsigned int nAsst, // number of rows
		 long unsigned int C[] ) // array of outputs defining function
{
  // Returns true iff degenerate
  // Note: This is not the most efficient approach, but should be adequate.

  // Build CPT containing the assignments to parents as well as preference
  // The resulting matrix is a truth table if d=2 and tables are complete

  // Allocate cpt array dynamically to avoid stack overflow from:
  // unsigned int cpt[nAsst][m+1]; //!stack overflow 
  unsigned int** cpt = nullptr;
  try
  {
    cpt = new unsigned int*[nAsst];
    for (int r = 0; r < nAsst; ++r) {
      cpt[r] = new unsigned int[m+1];
      cpt[r][m] = C[r];
    } // end for
  }
  catch (std::bad_alloc& ba)
  {
    std::cerr << "Out of memory [degen_multi:cpt]: Insufficient memory for CPTs with "
	      << nAsst << " entries!\n";
    exit(EXIT_FAILURE);
  } // end try-catch
  
  for (int c = 0; c < m; ++c) {
    int r = 0; // init row
    while (r < nAsst) {
      for (int x = 1; x <= DOM_SIZE; ++x) {
	for (int copies = 0; copies < int(pow(DOM_SIZE, m - c - 1)); ++copies)
	  cpt[r++][c] = x;
      } // end for copies
    } // end while r
  } // end for c

  // Determine degeneracy
  for (int y = 0; y < m; ++y) { // Iterate over Y in Pa(X)
    bool dependsOnY = false; // set flag: does function depend on Y?

    for (int a = 1; a <= DOM_SIZE; ++a) { // Iterate over pairs of values (y, y'), 
      for (int b = 1; b <= DOM_SIZE; ++b) { // such that y, y' in dom(Y)
	if (a == b) continue; // and y != y' of course

	// Allocate u1 and u2 arrays dynamically to avoid stack overflow from:
	// unsigned int u1[nAsst / DOM_SIZE]; // cpt(cpt(:,Y) == a, end)
	// unsigned int u2[nAsst / DOM_SIZE]; // cpt(cpt(:,Y) == b, end)
	unsigned int* u1 = nullptr;
	unsigned int* u2 = nullptr;
	try
	{
	  // Translating MATLAB code very literally 
	  u1 = new unsigned int[nAsst / DOM_SIZE];
	  u2 = new unsigned int[nAsst / DOM_SIZE];
	  int u1x = 0; // indexes u1
	  int u2x = 0; // indexes u2
	  int ax = 0; // indexes rows with cpt(:,Y) == a
	  int bx = 0; // indexes rows with cpt(:,Y) == b

	  // Iterate over the whole CPT but of course not beyond it
	  while (ax < nAsst && bx < nAsst) {

	    // Get the next pair of rows s.t. *cpt(:,Y) == a and b resp.
	    if (cpt[ax][y] != a) {
	      ++ax;
	      continue;
	    } // end if
	    if (cpt[bx][y] != b) {
	      ++bx;
	      continue;
	    } // end if

	    // copy values and proceed to the next pair of rows
	    u1[u1x++] = cpt[ax++][m];
	    u2[u2x++] = cpt[bx++][m];
	  
	  } // end while

	  // MATLAB: if ~all(cpt(cpt(:,Y) == a, end) == cpt(cpt(:,Y) == b, end))
	  for (int rr = 0; rr < nAsst / DOM_SIZE; ++rr) {
	    if (u1[rr] != u2[rr]) { 
	      dependsOnY = true; // that is, if outputs for some y, y' are different
	      delete[] u1; // darned goto statement complicates everything :(
	      delete[] u2; // darned goto statement complicates everything :(
	      goto nextY; // then continue with next Y
	    } // end if
	  } // end for rr

	  delete[] u1;
	  delete[] u2;
	}
	catch (std::bad_alloc& ba)
	{
	  std::cerr << "Out of memory [degen_multi:u1]: Insufficient memory for CPTs with "
		    << nAsst << " entries!\n";
	  exit(EXIT_FAILURE);
	} // end try-catch

      } // end for b
    } // end for a
	
  nextY:
    if (!dependsOnY)
    {
      // deallocate possibly huge array
      for (int r = 0; r < nAsst; ++r)
	delete[] cpt[r];
      delete[] cpt;
      return true; // if vacuous in Y, no need to look at other inputs
    } // end if
  } // end for Y in Pa(X)

  // deallocate possibly huge array
  for (int r = 0; r < nAsst; ++r)
    delete[] cpt[r];
  delete[] cpt;
  return false; // function is non-degenerate iff we reach this line
  
} // function degen_multi



// Generate a random, nondegenerate CPT (multi, incompleteness)
void rand_cpt( unsigned long int cpt[], // Return value
	       const int k,              // Number of parents
	       const long unsigned int size_cpt ) // Number of rows = d^k
{
  // Use rejection sampling to obtain a nondegenerate CPT.  Note: with
  // high probability this will occur on the first attempt--unless the
  // degree of incompleteness is set fairly high (at 1.0, it would
  // always be degenerate except for root nodes!)
  do {
    std::uniform_int_distribution<int> randomranking(1, factorial[DOM_SIZE]); 
    std::uniform_real_distribution<double> random01(0.0, 1.0);
    for (int i=0; i < size_cpt; ++i) // iterate over the CPRs of the CPT
      if (random01(mt_rand) >= DEG_INC) // roll the incompleteness dice
	cpt[i] = randomranking(mt_rand); // assign random permuation number
      else
	cpt[i] = 0; // no rule for this assignment to parents
  } while (degen_multi(k, size_cpt, cpt)); // reject if degenerate
}

// Print (to standard error) the permutations but without parent labels
void print_cpt(const unsigned long int cpt[], // CPT to output
	       const unsigned long int size_cpt ) // Number of rows = d^k
{
  // using namespace std; //!D
  // cout << "sizeof(cpt) = " << sizeof(cpt)
  //      << " size_cpt = " << size_cpt 
  //      << " cpt = " << cpt << endl;

  std::clog << "[ ";
  for (int i=0; i < size_cpt; ++i) {
    if (cpt[i])
      std::clog << num_to_ranking_1(cpt[i], DOM_SIZE) << ' ';
    else
      std::clog << "* "; // missing rule (values incomparable for this assignment)
  } //end for i
  std::clog << "]\n";
}

// End of file degen_multi.cc
