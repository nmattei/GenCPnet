/* 
 *  bitwise.h: Interface for bitwise helper functions.

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



// Better bit_scan_forward and bit_scan_reverse

#ifndef BITWISE_H
#define BITWISE_H

unsigned int bsf(uint64_t word);
unsigned int bsr(uint64_t word);
unsigned int popcount(uint64_t word);

#endif // BITWISE_H
 
// End of file: bitwise.h
