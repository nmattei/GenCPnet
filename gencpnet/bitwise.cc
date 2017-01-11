/* 
 *  bitwise.cc: Bitwise operations for scanning and counting

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

#include <cstdint>
#include "bitwise.h"

// Bit scan reverse 
unsigned int bsr(uint64_t word)
{
  if (!word) throw 0;
  uint64_t i; 
  __asm__ ( "bsrq\t%1, %0" // gcc doesn't give us intrinsic for _bitscan_reverse 
	    : "=r" (i)
	    : "rm" (word) );
  return i; // implicit type conversion 
} // end function bsr

// Bit scan forward
unsigned int bsf(uint64_t word)
{
  if (!word) throw 0;
  unsigned i = __builtin_ffs(word);
  if (i == 0)
  {
    word >>= 32;
    i = __builtin_ffs(word) + 32;
  }
  return i;
} // end function bsf

// Computing the Hamming weight of a 64-bit word
unsigned int popcount(uint64_t word)
{
  return __builtin_popcount(word) + __builtin_popcount(word >> 32);
} // end function popcount

// End of file: bitwise.cc
