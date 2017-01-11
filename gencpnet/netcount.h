/* 
 *  netcount.h: Interface for counting labeled DAGs and CP-nets 

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

// Uses GMP library to implement counting graphs and CP-nets
#ifndef COUNT_H
#define COUNT_H

#include <gmpxx.h>
#include "tables.h"

extern unsigned int DOM_SIZE; // domain size (homogeneous)
extern double DEG_INC;  // degree of incompleteness (also homogeneous)

class Netcount
{
 public:
  // Limiting to one constructor due to time constraints
  Netcount( const int n, // number of nodes in model
	    const int c ) // bound c on indegree (assume c <= n-1)
    : MAX_N_(n), 
      MAX_K_(std::min(n-1, c) + 1), 
      MAX_GAMMA_(std::min(n-1, c) + 1) {};

  ~Netcount();

  int MAX_N() { return MAX_N_; }
  int MAX_K() { return MAX_K_; }
  int MAX_GAMMA() { return MAX_GAMMA_; }

  void init();

  mpz_class binom(const int n, const int k);
  mpz_class phi(const int n, const bool winc);
  mpz_class gamma(const int k);
  
  mpz_class count_ldag(int n) { init_ldag(); return count_ldag(n, 1, 0); };
  mpz_class count_ldag(int n, int c) { init_bldag(); return count_bounded_ldag(n, c, 1, 0); };
  mpz_class count_cpnet(int n) { init_cpnet(); return count_cpnet(n, n - 1, 0, 0); };
  mpz_class count_cpnet(int n, int c) { init_cpnet(); return count_cpnet(n, c, 0, 0); };
  mpz_class prob_cpnet(int n) { init_cpnet(); return get_cpnet_cdf(n, n - 1, 0, 0); };
  mpz_class prob_cpnet(int n, int c) { init_cpnet(); return get_cpnet_cdf(n, c, 0, 0); };

  cpnet_dist* cdist = nullptr;

  void print_pascal();

 private:
  bool initialized_ = false;
  const int MAX_N_;
  const int MAX_K_;
  const int MAX_GAMMA_;

  mpz_class **PASCAL = nullptr;
  mpz_class *GAMMA = nullptr;
  mpz_class ***LDAG = nullptr;
  mpz_class ****BLDAG = nullptr;
  mpz_class ****CPNET = nullptr;

  void init_pascal();
  void init_gamma();
  void init_ldag();
  void init_bldag();
  void init_cpnet();

  void delete_gamma();
  void delete_pascal();
  void delete_ldag();
  void delete_bldag();
  void delete_cpnet();

  mpz_class count_ldag(int n, int j, int q);
  mpz_class count_bounded_ldag(int n, int c, int j, int q);
  mpz_class count_cpnet(int n, int c, int j, int q);
  mpz_class get_cpnet_cdf(int n, int c, int j, int q);
};

#endif // COUNT_H

// End of file: netcount.h
