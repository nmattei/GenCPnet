# CP-Net Generator

The CP-net Generator (GenCPnet) generates acyclic conditional preference networks (CP-nets) [2] uniformly at random with respect to a specified set. It is possible to specify parameters such as the number of nodes, bound on indegree, and the size of domains.

GenCPnet implements the method described in our paper, “Generating CP-nets Uniformly at Random”. If you use or adapt this software, we kindly ask that you would cite our paper. GenCPnet is free software, released under the GNU Public License version 3. GenCPNet is written in C++ and designed to run on a GNU Linux system. 

This code is described in detail in our AAAI 2016 Paper, [Generating CP-nets Uniformly at Random](http://www.nickmattei.net/docs/gencp.pdf).  Please cite our paper if you use this code for academic purposes.

For comprehensive instructions please [see the guide](/doc/Gencpnet_guide.pdf).

# Basic Instructions

CP-net generator (gencpnet) generates CP-nets uniformly at random parameterized by number of nodes, bound on indegree, and size of domains.

Usage: `gencpnet <OPTIONS> <DIRECTORY>`

OPTIONS (can be applied in any order):
```
-c              bound on indegree for all nodes (default -c 5)
--count         outputs number of CP-nets (given n, c, d) [no generation]
--countdags     outputs number of graphs (given n, c) [no generation]
-d              domain size (homogeneous for all features) (default -d 2)
-g              number of CP-nets to generate (default -g 1)
-h              Hamming distance of outcome pairs (optional and only 
                    used in conjunction with the -t option)
-i              probability that a given rule is missing (default -i 0.0)
-n              number of features/nodes
-t              also generates XML files each with a pair of outcomes for
                    dominance testing experiments (default -o 0)
--help          display this README.txt file
-V, --verbose   output generation details to standard error for debugging
--version       display version information
```

EXAMPLES:

```gencpnet -n 10 -c 5 -d 3 -g 100 testdirec```

Generates 100 CP-nets specifications in the testdirec subfolder, each with 10 nodes and indegree not greater than 5 parents for any node. The size of the domains is 3 for all 10 variables.  The files will be named cpnet_n10c5d3_0000.xml to cpnet_n10c5d3_0099.xml.

```gencpnet -n 5```

Generates just 1 CP-net with 5 nodes, unbounded indegree, to current directory.

```gencpnet -g 5 -n 7 -i 0.40 testdirec```

Generates 5 CP-nets with 7 nodes, maximum indegree 5 (default), binary domains, and on average 40% of the conditional preference rules left unspecified.

```gencpnet -n 10 -c 2 -d 3 -i 0.20 -t 20 experiments```

Generates one CP-net with 10 nodes, bound 2 on indegree, expected 20% incompleteness, stored as "experiments/cpnet_n10c2d3i20_0000.xml". Also generates 10 pairs of outcomes for DT experiments, each with 10 features with domains of size 3, stored as files dt_n10d3_0000_0000.xml dt_n10d3_0000_0019.xml in subdirectory "experiments/".  The outcome pairs are not constrained by Hamming distance since -h is unspecified.

The XML formats for CP-nets and DT problems are described at the site:
<http://www.ece.iastate.edu/~gsanthan/crisner.html>.
