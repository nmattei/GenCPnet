// Implements the joint conditional cummulative density functions for CP-nets
// For node j in the lexicographic topological sort as defined by the dag code:
// the table gives P( s, t | n, c, j ).
// The number of rows in each table is the triangular number of min(c, j).

#ifndef TABLES_H
#define TABLES_H

#include <random> // Mersenne Twister
#include <chrono> // Used to initialize MT

const unsigned long int factorial[17] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000};

void random_outcome_pair(uint64_t& o, uint64_t& oo, const unsigned int n, const unsigned int hd = 0);

uint64_t random_k_subset(uint64_t T, const int N, const int n); 

class cpnet_ccdf
{
 public:
  cpnet_ccdf(const int length);
  cpnet_ccdf(cpnet_ccdf &obj);
  ~cpnet_ccdf();
  int length() { return length_; };
  void store(double p, const int s, const int t) { P[row_] = p; S[row_] = s; T[row_] = t; ++row_; };
  void adjust_last() { P[length_ - 1] = 1.0; };
  void print();
  friend std::ostream& operator<<(std::ostream& os, const cpnet_ccdf& obj);
  friend std::istream& operator>>(std::istream& is, cpnet_ccdf& obj); // test
  void random_st(int& s, int& t);
  // static uint64_t random_cpt(const int k); // use rand_cpt in degen_multi.cc
  void random_node(const int n, int& q, uint64_t& U, uint64_t& A, unsigned long int** cpt, int j);
  // void random_node(const int n, int& q, uint64_t& U, uint64_t& A, uint64_t& cpt);
  // void random_toporder_node(const int n, const int j, int& q, uint64_t& U, uint64_t& A, uint64_t& cpt);
  
 private:
  int length_;
  int row_;
  double *P;
  int *S;
  int *T;
};

// all tables for all available values of (n, c, j, q) go here
class cpnet_dist
{
 public:
  cpnet_dist(int maxn, int maxk);
  ~cpnet_dist();
  cpnet_ccdf* dist(int n, int c, int j, int q) { return DIST[n][c][j][q]; }; 
  void init(int n, int c, int j, int q, int len);
  void print();
  void dump();
  void dump(std::string fname);
  static void dagcode_to_dag(const unsigned int n, const uint64_t dc[], unsigned int ch_ix[], char ch[], unsigned long imp[]); 
  static void dc_and_cpts_to_xml(const unsigned int n, const uint64_t dc[], unsigned long int** cpts, std::string working_directory, std::string fname); 
  void free_cpt(const int n, long unsigned int** cpt);
  void generate_random_cpnet(const int n, const int c, uint64_t dc[], long unsigned int* cpt[]);
  // void random_toporder_cpnet(const int n, const int c, uint64_t dc[], uint64_t cpt[]);
  
 private:
  const int MAX_N;
  const int MAX_K;
  cpnet_ccdf***** DIST = nullptr;
};

bool degen_multi( unsigned int m, // number of inputs (indegree)
		  unsigned int nAsst, // number of rows
		  unsigned int* C ); // array of outputs defining function

void rand_cpt( unsigned long int cpt[], // Return value
	       const int k,              // Number of parents
	       const long unsigned int size_cpt ); // Number of rows = d^k

void print_cpt( const unsigned long int cpt[], // CPT to output
	        const unsigned long int size_cpt ); // Number of rows = d^k

#endif // TABLES_H
