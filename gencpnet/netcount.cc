/* 
 *  netcount.cc: Counts labeled DAGs and CP-nets 

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

// Uses GMP library to implement counting graphs and CP-nets with
// multi-valued variables and incomplete tables

#include <iostream>
#include <algorithm> // max()
#include <utility> // ordered pair
#include <map> // multimap
#include <gmpxx.h> // Project makes heavy use of the Gnu MP library
#include <cmath> // pow()
#include "netcount.h" 
#include "tables.h" 

extern unsigned int DOM_SIZE;

Netcount::~Netcount()
{
  delete_pascal();
  delete_ldag();
  delete_bldag();
  delete_cpnet();
  delete_gamma();
  delete cdist;
  cdist = nullptr;
}

void Netcount::init()
{
  if (initialized_)
  {
    delete_bldag();
    delete_ldag();
    delete_cpnet();
    delete_gamma();
    delete_pascal();
    delete cdist;
    cdist = nullptr;
  }
  initialized_ = true;
  init_pascal();
  init_gamma();
  cdist = new cpnet_dist(MAX_N_, MAX_GAMMA_);
}

mpz_class Netcount::binom(const int n, const int k)
{
  mpz_class result = 0;
  if (n < 0 || k < 0)
    return result; // Knuth's convention?
  if (n >= MAX_N_ || k >= MAX_K_) { // range check
    std::cerr << "Binomial out of range\n";
    throw "Binomial out of range";
  }
  return PASCAL[n][k];
}
  
// Return the total number of multi-valued functions phi assuming that
// the output (CPR) is a ranking over the d values, or no entry if
// incompleteness is allowed
mpz_class Netcount::phi(const int n, const bool winc)
{
  // Compute the exponent part: d^n
  long unsigned int p = (long unsigned int) pow(DOM_SIZE, n); // p = d^n
 
  // Compute the factorial (base) part: d! (or d!+1 if incomplete)
  long unsigned int fac;
  if (DOM_SIZE <= 16) {
    fac = factorial[DOM_SIZE]; // table look up in netcount.h
  } else { // not sure the correct range on this--depends on long unsigned maxint
    fac = 1;
    for (int i = 1; i <= DOM_SIZE; ++i)
      fac *= i;
  }
  if (winc) 
    ++fac;

  // Compute and return the total number of functions phi
  mpz_class result;
  mpz_ui_pow_ui(result.get_mpz_t(), fac, p); // phi = (d!)^(d^n) or (d!+1)^(d^n)
  return result;
}

// Generalization allowing multivalued domains and possibly incomplete tables
mpz_class Netcount::gamma(const int n)
{
  mpz_class result = GAMMA[n];
  if (!result) {
    for (int k = 0; k <= n; ++k)
    {
      int sign = (n - k) % 2 ? -1 : 1;
      result += sign * binom(n, k) * 
	( DEG_INC * phi(k, true) 
	  + (1.0 - DEG_INC) * phi(k, false) );
      // WAS BROKEN:
      // result += sign * binom(n, k) * phi(k, DOM_SIZE, true); //!D
    } // for 
  } // if
  return result;
}

void Netcount::print_pascal()
{
  for (int n = 0; n < MAX_N_; ++n)
  {
    for (int k = 0; k <= n && k < MAX_K_; ++k)
      std::cout << PASCAL[n][k] << ' ';
    std::cout << std::endl;
  }
}

// Initialization function that builds Pascal's triangle up to the
// size specified by the parameter
void Netcount::init_pascal()
{
  PASCAL = new mpz_class*[MAX_N_];
  mpz_class **p = PASCAL; // alias
  p[0] = new mpz_class[MAX_K_];
  for (int k = 1; k < MAX_K_; ++k)
    p[0][1] = 0;
  p[0][0] = 1;
  for (int n = 1; n < MAX_N_; ++n)
  {
    p[n] = new mpz_class[MAX_K_];
    p[n][0] = 1;
    for (int k = 1; k < MAX_K_; ++k)
      p[n][k] = p[n-1][k-1] + p[n-1][k];
  }
}

void Netcount::init_gamma()
{
  if (GAMMA != nullptr)
    delete[] GAMMA;
  GAMMA = new mpz_class[MAX_GAMMA_] {};
  mpz_class tmp;
  for (int k = 0; k < MAX_GAMMA_; ++k) {
    tmp = gamma(k);
  }
}

void Netcount::init_ldag()
{
  if (LDAG != nullptr) return;
  LDAG = new mpz_class**[MAX_N_];
  mpz_class ***L = LDAG; // alias for convenience
  for (int i = 0; i < MAX_N_; ++i)
  {
    L[i] = new mpz_class*[i + 1];
    for (int j = 0; j <= i; ++j) 
    {
      L[i][j] = new mpz_class[j + 1];
      for (int k = 0; k <= j; ++k)
	L[i][j][k] = 0;
    }
  }  
}

void Netcount::init_bldag()
{
  if (BLDAG != nullptr) return;
  BLDAG = new mpz_class***[MAX_N_];
  mpz_class ****B = BLDAG; // alias for convenience
  for (int n = 0; n < MAX_N_; ++n)
  {
    B[n] = new mpz_class**[MAX_K_ + 1];
    for (int c = 0; c <= MAX_K_; ++c)
    {
      B[n][c] = new mpz_class*[n + 1];
      for (int j = 0; j <= n; ++j) 
      {
	B[n][c][j] = new mpz_class[j + 1];
	for (int q = 0; q <= j; ++q)
	{
	  B[n][c][j][q] = 0;
	} // for q
      } // for j
    } // for c
  } // for i
}

void Netcount::init_cpnet()
{
  if (CPNET != nullptr) return;
  CPNET = new mpz_class***[MAX_N_];
  mpz_class ****C = CPNET; // alias for convenience
  for (int n = 0; n < MAX_N_; ++n)
  {
    C[n] = new mpz_class**[MAX_K_ + 1];
    for (int c = 0; c <= MAX_K_; ++c)
    {
      C[n][c] = new mpz_class*[n + 1];
      for (int j = 0; j <= n; ++j) 
      {
	C[n][c][j] = new mpz_class[j + 1];
	for (int q = 0; q <= j; ++q)
	{
	  C[n][c][j][q] = 0;
	} // for q
      } // for j
    } // for c
  } // for i
}

void Netcount::delete_gamma()
{
  if (GAMMA == nullptr) return;
  delete[] GAMMA;
  GAMMA = nullptr;
}

void Netcount::delete_pascal()
{
  if (PASCAL == nullptr) return;
  for (int n = 0; n < MAX_N_; ++n)
    delete[] PASCAL[n];
  delete[] PASCAL;
  PASCAL = nullptr;
}

void Netcount::delete_ldag()
{
  if (LDAG == nullptr) return;
  for (int i = 0; i < MAX_N_; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      delete[] LDAG[i][j];
      LDAG[i][j] = nullptr;
    }
    delete[] LDAG[i];
  }
  delete[] LDAG;
  LDAG = nullptr;
}

void Netcount::delete_bldag()
{
  if (BLDAG == nullptr) return;
  for (int n = 0; n < MAX_N_; ++n)
  {
    for (int c = 0; c <= MAX_K_; ++c)
    {
      for (int j = 0; j <= n; ++j)
      {
	delete[] BLDAG[n][c][j];
      }
      delete[] BLDAG[n][c];
    }
    delete[] BLDAG[n];
  }
  delete[] BLDAG;
  BLDAG = nullptr;
}

void Netcount::delete_cpnet()
{
  if (CPNET == nullptr) return;
  mpz_class ****C = CPNET; // alias for convenience
  for (int n = 0; n < MAX_N_; ++n)
  {
    for (int c = 0; c <= MAX_K_; ++c)
    {
      for (int j = 0; j <= n; ++j)
      {
	delete[] C[n][c][j];
      }
      delete[] C[n][c];
    }
    delete[] C[n];
  }
  delete[] C;
  C = nullptr;
}

// Tom's LDAG recurrence based on DAG codes
mpz_class Netcount::count_ldag(int n, int j, int q)
{
  if (j >= n) return 1; // base case
  mpz_class count = 0;
  for (int s = 0; s <= q; ++s)
    for (int t = 0; q + t <= j; ++t)
    {
      mpz_class rr;
      if (!(rr = LDAG[n][j+1][q+t]))
	rr = LDAG[n][j+1][q+t] = count_ldag(n, j + 1, q + t);
      count += binom(q, s) * binom(n - q, t) * rr;
    } // end for t
  return count;
}

mpz_class Netcount::count_bounded_ldag(int n, int c, int j, int q)
{
  if (j >= n) return 1; // base case
  mpz_class count = 0;
  for (int s = 0; s <= c && s <= q; ++s)
    for (int t = 0; s + t <= c && q + t <= j; ++t)
    {
      mpz_class rr;
      if (!(rr = BLDAG[n][c][j+1][q+t]))
	rr = BLDAG[n][c][j+1][q+t] = count_bounded_ldag(n, c, j + 1, q + t);
      count += binom(q, s) * binom(n - q, t) * rr;
    } // end for t
  return count;
}

mpz_class Netcount::count_cpnet(int n, int c, int j, int q)
{
  if (j >= n) return 1; // base case
  if (j <= 0) return gamma(0) * count_cpnet(n, c, 1, 0);
  mpz_class count = CPNET[n][c][j][q];
  if (!count)
  {
    for (int s = 0; s <= c && s <= q; ++s)
      for (int t = 0; s + t <= c && q + t <= j; ++t)
	count += gamma(s + t) * binom(q, s) * binom(n - q, t) * count_cpnet(n, c, j + 1, q + t);
    CPNET[n][c][j][q] = count;
  } // end if
  return count;
}

mpz_class Netcount::get_cpnet_cdf(int n, int c, int j, int q)
{
  if (j >= n) return 1; // base case
  // FOR IID-DAG REPLACE 2 WITH 1:
  if (j <= 0) return gamma(0) * get_cpnet_cdf(n, c, 1, 0); // node 0 has no parents but d! or d!+1 choices of CPT
  mpz_class count = CPNET[n][c][j][q]; // see if we have already evaluated "node" (n,c,j,q)
  if (!count) // if not, then let's evaluate it
  {
    // We are actually implementing a (temporary) ordered triple st_counts = (GMP-float, int, int) 
    std::multimap <mpz_class, std::pair <int, int>> st_counts;

    // Implements Notes 1:48.1 
    for (int s = 0; s <= c && s <= q; ++s)
      for (int t = 0; s + t <= c && q + t <= j; ++t)
      {
	// TO CREATE TABLES FOR IID-DAG, PUT COMMENT BEFORE GAMMA:
	mpz_class count_st = gamma(s + t) *
	  binom(q, s) * binom(n - q, t) *
	  get_cpnet_cdf(n, c, j + 1, q + t);
	count += count_st;
	st_counts.insert(std::make_pair(count_st, std::make_pair(s, t)));
      } // end for t

    // Create a distribution table of the proper size in the node (n,c,j,q)
    cdist->init(n, c, j, q, int(st_counts.size())); 

    // Convert counts to cumulative probabilities such that the last one has P = 1.0
    mpf_class prob = 0.0; // ranges from [0.0, 1.0]
    for (auto rit = st_counts.rbegin(); rit != st_counts.rend(); rit++) { // sorted least to greatest
      mpf_class f_count_st = rit->first; // convert from GMP int to GMP float
      mpf_class f_count = count;
      prob += f_count_st / f_count; // scale to [0,1]

      // store (prob, s, t) tuple in table
      cdist->dist(n, c, j, q)->store(mpf_get_d(prob.get_mpf_t()), rit->second.first, rit->second.second);
    } // end if
    if (abs(1.0 - prob) < 1e-6) // we allow for numerical issues up to 1 in a million
      cdist->dist(n, c, j, q)->adjust_last(); // adjust last item to probability 1
    else
      throw 0; // probabilities should total 1.0 to machine precision
    CPNET[n][c][j][q] = count; // store count (GMP int) in table for efficient computation
  }
  return count; 
}

