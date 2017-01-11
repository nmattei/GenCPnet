#######################################################################
#
#   Makefile for GenCPnet 
#   Copyright (c) 2015 Thomas E. Allen
# 
#   This file is part of GenCPnet, the Uniform CP-net Generator.
# 
#   GenCPnet is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published
#   by the Free Software Foundation, either version 3 of the License,
#   or (at your option) any later version.
# 
#   GenCPnet is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   General Public License for more details.
# 
#   You should have received a copy of the GNU General Public License
#   along with GenCPnet.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################


# QUICK OVERVIEW:
#
# make 		compiles only the code
# make all	compiles
# make clean	removes binaries to facilitate a clean rebuild
# make install	copies executable file 'gencpnet' to directory 'BIN'
#
# If speed is really important to you, see below about CXXFLAGS


# Current software release number
VERSION = 0.70

# Directory in which to install executable; used with 'make install'
BIN = ~/bin

# Directory in which documentation sources are installed
# Used with 'make all' and 'make doc'
DOC = doc

# C++ compiler and make tool to use
CXX = g++
MAKE = make

# C++ compilation flags
# NOTE: Choose one of the following and comment out the other
CXXFLAGS = -Wfatal-errors --std=c++11 -g # Choose this for debugging
#CXXFLAGS = --std=c++11 -O3 # Choose this if speed is really important

# By default, 'make' checks that GMP is installed and then compiles the code
code:	check gencpn 
	@echo "Success!"

# Use 'make all' to compile the code and also the documentation 
all:	code doc

# Compiles the code for GenCPnet
gencpn:	main.o bitwise.o netcount.o tables.o degen_multi.o findperm.o outcomes.o 
	$(CXX) $(CXXFLAGS) $^ -o gencpnet -lgmpxx -lgmp 

doc:
	$(MAKE) -C $(DOC) all

# Code dependencies 
main.o:	main.cc help.h netcount.h outcomes.h tables.h
bitwise.o:	bitwise.cc bitwise.h
netcount.o:	netcount.cc netcount.h tables.h
tables.o:	tables.cc tables.h degen_multi.h findperm.h bitwise.h
degen_multi.o:	degen_multi.cc degen_multi.h tables.h findperm.h
findperm.o:	findperm.cc findperm.h
outcomes.o:	outcomes.cc outcomes.h tables.h

# Checks that the GMP library is installed
check:
	@echo "Checking if GMP library is installed..."
	@if ldconfig -p | grep libgmp; then \
		echo "GMP library found."; \
	else \
		echo "GMP LIBRARY NOT FOUND"; \
		echo "Gnu Multiprecision Library must be installed first"; \
		echo "See: https://gmplib.org/manual/Installing-GMP.html"; \
		exit 1; \
	fi 

# Phantom target to clean up object and executable files (only in this directory!)
clean:
	@echo "Removing binaries."
	rm -f *.o gencpnet

# Creates a zipped archive for distribution the source code
dist:
	rm -f gencpnet.zip 
	rm -rf gencpnet_$(VERSION)
	mkdir -p gencpnet_$(VERSION)
	cp *.cc gencpnet_$(VERSION)
	cp *.h gencpnet_$(VERSION)
	cp Makefile gencpnet_$(VERSION)
	cp gpl-3.0.txt gencpnet_$(VERSION)
	cp README.txt gencpnet_$(VERSION)
	cp $(DOC)/*.pdf gencpnet_$(VERSION)
	cp -r $(DOC) gencpnet_$(VERSION)
	zip -r gencpnet_$(VERSION).zip gencpnet_$(VERSION)

# Installs the executable file to a suitable destination (e.g., ~/bin)
install:
	cp gencpnet $(BIN)


# End of Makefile for GenCPnet
