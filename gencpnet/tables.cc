/* 
 *  tables.cc: Constructs distribution tables for CP-net generation

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

// Implements the tables representing the distribution over CP-nets
// For node j in the lexicographic topological sort as defined by the dag code:
// the table gives P( s, t | n, c, j ).
// The number of rows in each table is the triangular number of min(c, j).

#include <cstdint>
#include <cstring>
#include <fstream> // import-export
#include <functional> // We need bind for MT
#include <iostream>
#include <new> // Possibly needed for bad_alloc
#include <random> // Mersenne Twister
#include <string>
#include "bitwise.h"
#include "degen_multi.h" 
#include "findperm.h"
#include "tables.h"

extern unsigned int DOM_SIZE;
extern bool VERBOSE;

// x86 assembly language linkage
// extern "C" uint64_t vacuous(uint64_t x, uint64_t n);

// Initialize random number generator
std::mt19937_64 mt_rand(std::chrono::high_resolution_clock::now().time_since_epoch().count());

// Given the number of features (bits) n and an optional Hamming
// distance constraint (the number of features in which the two
// outcomes differ), return a randomly generated pair of outcomes.
void random_outcome_pair(uint64_t& o, uint64_t& oo, const unsigned int n, const unsigned int hd)
{
  // Define distribution for the random number generator.
  // U64(mt_rand) will then call the Mersenne Twister generator to
  // obtain an n-bit uniformly generated random integer (n <= 64).
  std::uniform_int_distribution<uint64_t> U64(0, (01ULL << n) - 1); 

  o = U64(mt_rand); // generate first outcome
  if (hd)
    oo = o ^ random_k_subset((01ULL << n) - 1, n, hd);
  else
    do oo = U64(mt_rand); while (!(o ^ oo)); // generate a distinct second outcome 
}

// Knuth's algorithm 3.4.2S: Select a subset of size n from a set of size N.
uint64_t random_k_subset(uint64_t T, const int N, const int n)
{
  // S1. [Initialize.]
  std::uniform_real_distribution<double> random01(0.0, 1.0);
  uint64_t S = 0; // result subset
  int t = 0, m = 0; // t things dealt with, m selected so far

  while (m < n) 
  {
    // Select and remove the next element from set T
    uint64_t I = T & -T; // extract rightmost bit
    T ^= I; // remove bit from "set"

    // S2. [Generate U.]
    double U = random01(mt_rand);

    // S3. [Test.]
    if ((N - t) * U < n - m)
    {
      // S4. [Select.]
      S |= I;
      ++m;
    } // else S5. [Skip.]
    ++t; // NOW (not before S3!) count it among the things dealt with 
  } // end while

  return S;
}

cpnet_ccdf::cpnet_ccdf(int length)
{
  length_ = length;
  row_ = 0;
    
  P = new double [length_];
  S = new int [length_];
  T = new int [length_];
}

cpnet_ccdf::cpnet_ccdf(cpnet_ccdf &obj)
{
  // std::cerr << "copying object with length " << obj.length_ << "\n";
  length_ = obj.length_;
  row_ = obj.row_;

  P = new double [length_];
  S = new int [length_];
  T = new int [length_];
  for (int i = 0; i < length_; ++i)
  {
    P[i] = obj.P[i];
    S[i] = obj.S[i];
    T[i] = obj.T[i];
  }
}

cpnet_ccdf::~cpnet_ccdf()
{
  delete[] P;
  delete[] S;
  delete[] T;
}

void cpnet_ccdf::print()
{
  std::cout << std::dec;
  for (int i = 0; i < length_; ++i)
    std::cout << P[i] << '\t' << S[i] << '\t' << T[i] << std::endl;
}

std::ostream& operator<<(std::ostream& os, const cpnet_ccdf& obj)
{
  for (int i=0; i < obj.length_; ++i)
  {
    std::uint64_t u;
    std::memcpy(&u, &obj.P[i], sizeof(obj.P[i]));
    os << std::hex << u << ' ' << obj.S[i] << ' ' << obj.T[i] << std::endl;
  }
  return os;
}

// Reads a complete line of text from the input stream.
// We expect three integers in hex, such as: 3fe60dca404a8b61 0 1
// The first is the probability (a 64-bit double masquerading as
// a 64-bit int), the next two are the integer values of s and t.
std::istream& operator>>(std::istream& is, cpnet_ccdf& obj)
{
  std::uint64_t u;
  is >> std::hex >> u;
  std::memcpy(&obj.P[obj.row_], &u, sizeof(obj.P[obj.row_]));
  is >> obj.S[obj.row_];
  is >> obj.T[obj.row_];

  // std::cout << obj.P[obj.row_] << ' ' << obj.S[obj.row_] << obj.T[obj.row_];
  return is;
}

// Selects a pair (s, t) iid from the distribution defined in the table
void cpnet_ccdf::random_st(int &s, int& t)
{
  std::uniform_real_distribution<double> U(0.0, 1.0);
  double p = U(mt_rand);
  // std::cout << p << '\t';
  int i = 0;
  while (p > P[i]) ++i; // Linear search, but we expect p to be close to the top
  s = S[i];
  t = T[i]; 
}

// Returns a random node consisting of the values of (s, t) as well as a random CPT
// void cpnet_ccdf::random_node(int& s, int& t, uint64_t& cpt)
void cpnet_ccdf::random_node( const int n, 
			      int& q, 
			      uint64_t& U, 
			      uint64_t& A, 
			      unsigned long int** cpt, 
			      int j )
{
  // using namespace std;
  // cout << "Entering random_node\n";

  int s, t; // subset sizes s elements selected from U, t from V \ U
  random_st(s, t); // choose (s, t) from the distribution defined by the table

  unsigned long int size_cpt = int(pow(DOM_SIZE, s + t)); // get size of this node's CPT
  
  cpt[j] = new long unsigned int[size_cpt]; // allocate space for the array of permutation numbers

  rand_cpt(cpt[j], s + t, size_cpt); // choose a cpt of the appropriate size (indegree)
  
  // print_cpt(cpt[j], size_cpt);

  // cpt = random_cpt(s + t); // choose a cpt of the appropriate size (indegree)
  uint64_t S = random_k_subset(U, q, s); // random subset of U of size s
  uint64_t T = random_k_subset(~U, n - q, t); // random subset of V\U of size t

  if (popcount(S) > n || popcount(T) > n)
  {
    std::cerr << "Stop!!\n";
    exit(EXIT_FAILURE);
  }
  
  U |= T; // U' = U union T
  q += t; // q' = |U'| = q + t
  A = S | T; // Dagcode element 

  // cout << "Exit random_node\n";
}

cpnet_dist::cpnet_dist(int maxn, int maxk) : MAX_N(maxn), MAX_K(maxk)
{
  DIST = new cpnet_ccdf****[MAX_N];
  cpnet_ccdf *****C = DIST; // alias for convenience
  for (int n = 0; n < MAX_N; ++n)
  {
    C[n] = new cpnet_ccdf***[MAX_K + 1];
    for (int c = 0; c <= MAX_K; ++c)
    {
      C[n][c] = new cpnet_ccdf**[n + 1];
      for (int j = 0; j <= n; ++j) 
      {
	C[n][c][j] = new cpnet_ccdf*[j + 1];
	for (int q = 0; q <= j; ++q)
	{
	  C[n][c][j][q] = nullptr;
	} // for q
      } // for j
    } // for c
  } // for i
}

cpnet_dist::~cpnet_dist()
{
  if (DIST == nullptr) return;
  cpnet_ccdf *****C = DIST; // alias for convenience
  for (int n = 0; n < MAX_N; ++n)
  {
    for (int c = 0; c <= MAX_K; ++c)
    {
      for (int j = 0; j <= n; ++j)
      {
	for (int q = 0; q <= j; ++q)
	{
	  delete C[n][c][j][q]; // created with new not new[]
	}
	delete[] C[n][c][j];
      }
      delete[] C[n][c];
      // std::cout << ".\n";
    }
    delete[] C[n];
  }
  delete[] C;
  C = nullptr;
}

void cpnet_dist::init(int n, int c, int j, int q, int len)
{
  // Initialize data structure
  DIST[n][c][j][q] = new cpnet_ccdf(len);
}

// Convert dagcode to dag (adapted from Steinsky p 2)
void cpnet_dist::dagcode_to_dag(
  const unsigned int n, // number of nodes
  const uint64_t dc[],  // the DagCode 
  unsigned int ch_ix[], // the index of the child set for each node in cch[]
  char cch[],           // the index and parent magnitude of each child, indexed by ch_ix
  unsigned long imp[] ) // importance vector -- needed for DT* and variants (see Li et al.)
{
  // std::clog << dc << '\t';
  // for (int i=0; i < n; ++i) std::clog << dc[i] << ' ';
  // std::clog << '\n';

  // Stack overflow? Seems most unlikely for n <= 64 which is already stipulated
  char ch[2*n*n]; // temp structure that we will compact and return as cch[]
  for (int i = 0; i < 2*n*n; ++i) // initialize temporary child storage structure
    ch[i] = 0;
  
  uint64_t S = (01ULL << n) - 1; // bit mask with 1...111 (n ones) = {1, ..., n}
  for (int k = n - 1; k > 0; --k)
  {
    uint64_t U = 0; // set union
    for (int j = 1; j <= k; ++j)
      U |= dc[j];
    uint64_t unseen = S ^ U;
    uint64_t s = bsr(unseen);
    // __asm__ ( "bsrq\t%1, %0" // gcc doesn't give us an intrinsic for _bitscan_reverse 
    // 	      : "=r" (s)
    // 	      : "rm" (unseen) );
    S ^= 01ULL << s; // remove s from set

    uint64_t Pa = dc[k]; // copy parent set, since we need to remove from it
    int mag = 1;
    while (Pa)
    {
      // int u = __builtin_ffs(Pa) - 1; // get first parent // 32-bit limit
      int u = bsf(Pa) - 1; // get first parent
      if (u >= n)
      {
	std::cerr << "Range error in dagcode_to_dag (u=" << u << " exceeds n=" << n << ")\n";
	throw "range error"; // range error somewhere
      } // end if u out of range
      Pa ^= 01ULL << u; // remove from set
      int nCh = ch[2*n*u];
      ch[2*n*u] = nCh + 1;
      ch[2*n*u + 2*nCh + 1] = char(s + 1);
      ch[2*n*u + 2*nCh + 2] = char(mag);
      mag <<= 1;
    } // while Pa nonempty

    imp[s] = 1; // give one to the node for itself
    int child = 2*n*s + 1; // find its children info in ch
    for (int i = 0; i < ch[2*n*s]; ++i)
    {
      imp[s] += imp[ch[child] - 1]; // and add the weight/importance of each child
      child += 2; // skip over magnitude data since not needed fof this
    } // for i
  } // for k
  
  // We still need to compute the importance of the "root" node corresponding to dc[0]
  // int u = __builtin_ffs(S) - 1; // recover the label of dagcode "dc[0]" (the "root" node)
  int u = bsf(S) - 1; // recover the label of dagcode "dc[0]" (the "root" node)
  imp[u] = 1; // counting process same as described above for each node s
  int child = 2*n*u + 1;
  for (int i = 0; i < ch[2*n*u]; ++i)
  {
    imp[u] += imp[ch[child] - 1];
    child += 2;
  } // for i

  // Pack ch to return as cch indexed by ch_ix
  int i = 0; // index into cch
  for (int j = 0; j < n; ++j) // iterate over the n nodes
  {
    int nCh = ch[2*n*j]; // j has nCh children
    ch_ix[j] = i;
    for (int k = 0; k < nCh; ++k) // iterate over children of node j
    {
      cch[i++] = ch[2*n*j + 2*k + 1]; // index of child k
      cch[i++] = ch[2*n*j + 2*k + 2]; // magnitude of node j as the parent of child k
    }
    cch[i++] = char(-1);
    cch[i++] = char(-1);
    // if (i > 2*n + 2*c*n - c*c - 2*c) throw "Index out of range"; // good idea but c hidden from us
  } // for j
}

// Format described (in part) at http://www.ece.iastate.edu/~gsanthan/crisner.html
void cpnet_dist::dc_and_cpts_to_xml(
  const unsigned int n, // number of nodes
  const uint64_t dc[],  // the DagCode 
  unsigned long int** cpts, // the CPTs
  std::string working_directory, // file folder
  std::string fname ) // file name NOT including path
{
  if (VERBOSE) {
    std::clog << "Writing to " << working_directory << '/' << fname << std::endl;
  } // end if

  // Open file--but first make sure a file with this name does not already exist!
  std::fstream xout (working_directory + '/' + fname, std::ios_base::in);
  if (xout) { // file exists, report error and die!
    std::cerr << "Error: filename " << fname << " already exists.\n"
	      << "Delete file(s) first or output to another directory.\n";
    exit(EXIT_FAILURE);
  } else { // a file with this name does not already exist, so attempt to create it
    xout.open(working_directory + '/' + fname, std::ios_base::out);
    if (!xout) { // we were unable to create the file 
      std::cerr << "Error: cannot open output file " << fname << "\n"
		<< "Make sure specified directory " << working_directory << " is accessible.\n";
      exit(EXIT_FAILURE);
    } // end if unable to create file
  } // end if file already exists

  // Outer tag
  xout << "<PREFERENCE-SPECIFICATION>\n\n";

  // Output variables and domains
  for (int i = 1; i <= n; ++i) {

    xout << "<PREFERENCE-VARIABLE>" << std::endl;
    xout << " <VARIABLE-NAME>x" << i << "</VARIABLE-NAME>\n";
    for (int j = 1; j <= DOM_SIZE; ++j)
      xout << " <DOMAIN-VALUE>" << j << "</DOMAIN-VALUE>\n";
    xout << "</PREFERENCE-VARIABLE>\n\n";
  } // end for i

  // Initialize data structures for mappings from labels to CPTs and parents
  unsigned int map[n]; // the mapping from node labels to dagcode elements
  unsigned int parents[n][n]; // the parents stored by their child's node label
  unsigned int nPa[n]; // the number of parents 
  unsigned long int sizeCPT[n]; // the size of each CPT
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j)
      parents[i][j] = 0;
    nPa[i] = 0;
    sizeCPT[i] = 1;
  }

  // Adapted from dagcode_to_dag 
  uint64_t S = (01ULL << n) - 1; // bit mask with 1...111 (n ones) = {1, ..., n}
  for (int k = n - 1; k > 0; --k) // iterate through the dagcode from right to left
  {
    uint64_t U = 0; // set union
    for (int j = 1; j <= k; ++j)
      U |= dc[j];
    uint64_t unseen = S ^ U;
    uint64_t s = bsr(unseen); // the label that corresponds to dagcode element k
    // __asm__ ( "bsrq\t%1, %0" // gcc doesn't give us intrinsic for _bitscan_reverse 
    // 	      : "=r" (s)
    // 	      : "rm" (unseen) );

    if (s >= n) // check for range error
    {
      std::cerr << "Range error in dc_and_cpts_to_xml "
		<< "(s should be in range 0.." << n-1 << ", but s=" << s << ")\n";
      //!D
      std::cerr << "unseen = 0x" << std::hex << unseen << std::endl;
      throw "range error"; // range error somewhere
    } // end if u out of range

    if (VERBOSE)
    {
      std::cerr << "Decoding node s=" << (s+1) << '\n';
    }
    S ^= 01ULL << s; // remove s from set S

    map[s] = k; // the node with label s corresponds to dagcode k and CPT k

    // std::clog << "k=" << k << " s=" << s << '\n'; //!D

    uint64_t Pa = dc[k]; // copy parent set, since we need to remove from it
    // nPa[s] = __builtin_popcount(Pa); // get number of parents of node s
    nPa[s] = popcount(Pa); // get number of parents of node s
    sizeCPT[s] = int(pow(DOM_SIZE, nPa[s])); // and size of CPT
    int j = 0; // indexes the parent labels of s encoded in dc[k]
    while (Pa)
    {
      // uint64_t u = __builtin_ffs(Pa); // get parent with smallest (?) label // 32 bit limit?!
      uint64_t u = bsf(Pa); // get parent with smallest (?) label
      // uint64_t u; // the label that corresponds to parent with smallest (?) label
      if (!u || u > n)
      {
	std::cerr << "Range error in dc_and_cpts_to_xml "
		  << "(u should be in range 1.." << n << ", but u=" << u << ")\n";
	//!D
	std::cerr << "Pa = 0x" << std::hex << Pa << std::endl;
	throw "range error"; // range error somewhere
      } // end if u out of range
      Pa ^= 01ULL << (u - 1); // remove u from temporary parent set
      if (VERBOSE)
	std::clog << "Node x" << u << " is a parent of x" << (s + 1) << ".\n"; //!D
      parents[s][j++] = u; // index labels from 0, but not arrays
    } // while Pa nonempty
  } // for dagcode element k
  // End of code adapted from Steinsky
  
  // Process root node corresponding to dc[0]
  // int u = __builtin_ffs(S) - 1; // recover the label of root  // 32 bit limit
  int u = bsf(S) - 1; // recover the label of root 
  if (u >= n)
  {
    std::cerr << "Range error in dc_and_cpts_to_xml "
	      << "(u should be in range 1.." << n << ", but u=" << u << ")\n";
    // std::cerr << "S = 0x" << std::hex << S << std::endl;
    throw "range error"; // range error somewhere
  } // end if u out of range
  map[u] = 0;
  // Note that parents[0][j] = 0 for j = 0..n-1
  // std::clog << "k=" << 0 << " s=" << u << '\n'; //!D

  // Output mapping and parents by node labels
  if (VERBOSE) {
    for (int i = 0; i < n; ++i) {
      std::clog << "Node x" << (i + 1) << ": "
		<< "map[" << i << "]=" << map[i] 
		<< std::hex << " <" << (map[i] ? dc[map[i]] : 0) << ">: " << std::dec;

      int j = 0;
      while (parents[i][j] && j < n) {
	std::clog << parents[i][j++] << ' ';
      }
      std::clog << "nPa=" << nPa[i] << ' '
		<< "sizeCPT=" << sizeCPT[i] << " CPT=";
      // std::clog << std::endl;
      print_cpt(cpts[map[i]], sizeCPT[i]);
    } // end for i
  }

  // Output conditions in XML format
  for (int i = 0; i < n; ++i) // Iterate over nodes in order of their labels
  {
    // Build table containing the assignments to parents for CPT entries
    // [Adapted from degen_multi.cc]

    // Dynamically allocate fullCPT to avoid stack overflow:
    // unsigned int fullCPT[sizeCPT[i]][nPa[i]+1]; // m? i.e., nPa[i]??
    unsigned int** fullCPT = nullptr;
    try
    {
      fullCPT = new unsigned int*[sizeCPT[i]];
      for (int r = 0; r < sizeCPT[i]; ++r)
	fullCPT[r] = new unsigned int[nPa[i]+1]; 
    }
    catch (std::bad_alloc& ba)
    {
      std::cerr << "Out of memory: Insufficient memory for CPTs with "
		<< sizeCPT[i] << " entries!\n";
      exit(EXIT_FAILURE);
    } // end try-catch
    
    //!DEBUG
    for (int r = 0; r < sizeCPT[i]; ++r) 
      fullCPT[r][nPa[i]] = cpts[map[i]][r]; // may not need this

    for (int c = 0; c < nPa[i]; ++c) {
      int r = 0; // init row
      while (r < sizeCPT[i]) {
	for (int x = 1; x <= DOM_SIZE; ++x) {
	  for (int copies = 0; copies < int(pow(DOM_SIZE, nPa[i] - c - 1)); ++copies)
	    fullCPT[r++][c] = x;
	} // end for copies
      } // end while r
    } // end for c

    for (int r = 0; r < sizeCPT[i]; ++r)
    {
      if (fullCPT[r][nPa[i]]) // if CPR r specifies an ordering over x_i
      {
	xout << "<PREFERENCE-STATEMENT>\n"
		  << "  <STATEMENT-ID>p" << (i + 1) << "_" << (r + 1) 
		  << "</STATEMENT-ID>\n"
		  << "  <PREFERENCE-VARIABLE>x" << (i + 1) 
		  << "</PREFERENCE-VARIABLE>\n";

	// Whether the if is needed as well as the for statement
	// depends on whether multiple conditions (parents) are
	// expressed as separate lines or together with commas.
	// Assuming the former for now.
	if (nPa[i] > 0) // if there are parents, hence conditions to express
	{
	  for (int c = 0; c < nPa[i]; ++c) {
	    xout << "  <CONDITION>"
		      << "x" << parents[i][c] // parent label
		      << "=" << fullCPT[r][c] // parent value
		      << "</CONDITION>\n";
	  } // end for column c of table (== parent of i)
	} // end if there are parents

	xout << perm_num_to_xml(fullCPT[r][nPa[i]], DOM_SIZE);

	xout << "</PREFERENCE-STATEMENT>\n\n";
      } // end if CPR r specifies an ordering over variable x_i      
    } // for row r in CPT of node with label i

    // Deallocate big array
    for (int r = 0; r < sizeCPT[i]; ++r)
      delete[] fullCPT[r];
    delete[] fullCPT;
    
  } // end for node label i

  // Outer tag close
  xout << "</PREFERENCE-SPECIFICATION>\n";

}

void cpnet_dist::free_cpt(const int n, long unsigned int** cpt)
{
  for (int i=0; i<n; ++i) {
    delete[] cpt[i];
    cpt[i] = nullptr;
  } // end for
}

void cpnet_dist::generate_random_cpnet(const int n, const int c, uint64_t dc[], long unsigned int** cpt)
{
  // using namespace std;
  // cout << "Entering generate_random_cpnet\n";

  uint64_t U = 00; // initially the union is the empty set
  int q = 0; // hence q = 0 = |U| = POPCNT(U)
  for (int j = 1; j < n; ++j) {
    dist(n, c, j, q)->random_node(n, q, U, dc[j], cpt, j);
    // long unsigned int size_cpt = pow(DOM_SIZE, q);
    // cout << "Node generated: q=" << q 
    // 	 << " size_cpt=" << size_cpt 
    // 	 << " &cpt=" << cpt[j] << endl;
    // print_cpt(cpt[j], size_cpt);
  }
  
  // cout << "Generating root node cpt\n";
  cpt[0] = new long unsigned int[1]; // allocate space for CPT of root node
  rand_cpt(cpt[0], 0, 1); // obtain the CPT entries
  // print_cpt(cpt[0], 1);
}

void cpnet_dist::print()
{
  std::cout << MAX_N << ',' << MAX_K << std::endl;
  for (int n = 1; n < MAX_N; ++n)
    for (int c = 0; c < MAX_K; ++c)
      for (int j = 0; j <= n; ++j)
	for (int q = 0; q <= j; ++q)
	  if (DIST[n][c][j][q] != nullptr)
	  {
	    std::cout << n << ',' << c << ',' << j << ',' << q << ','
		      << DIST[n][c][j][q]->length() << '\n';
	    DIST[n][c][j][q]->print();
	  }
}

void cpnet_dist::dump()
{
  std::cout << MAX_N << ',' << MAX_K << std::endl;
  for (int n = 1; n < MAX_N; ++n)
    for (int c = 0; c < MAX_K; ++c)
      for (int j = 0; j <= n; ++j)
	for (int q = 0; q <= j; ++q)
	  if (DIST[n][c][j][q] != nullptr)
	    std::cout << n << ',' << c << ',' << j << ',' << q << ','
		      << DIST[n][c][j][q]->length() << '\n'
		      << *DIST[n][c][j][q];
}

void cpnet_dist::dump(std::string fname)
{
  std::ofstream outf;
  outf.open(fname);
  outf << MAX_N << ',' << MAX_K << std::endl;
  for (int n = 1; n < MAX_N; ++n)
    for (int c = 0; c < MAX_K; ++c)
      for (int j = 0; j <= n; ++j)
	for (int q = 0; q <= j; ++q)
	  if (DIST[n][c][j][q] != nullptr)
	    outf << n << ',' << c << ',' << j << ',' << q << ','
	       << DIST[n][c][j][q]->length() << '\n'
	       << *DIST[n][c][j][q];
  outf.close();
}

// End of file: tables.cc
