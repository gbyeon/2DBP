# Dedicated Benders Decomposition for Bilevel Problems with Convex Follower
  
This repository provides an implementation of a dedicated Benders decomposition method for bilevel optimization problems with convex follower, proposed in the following paper:

- Geunyeong Byeon and Pascal Van Hentenryck. [Benders Subproblem Decomposition for Bilevel Problems with Convex Follower.](https://arxiv.org/pdf/1801.03520) *Arxiv*, 2019

The method features a Benders subproblem decomposition scheme and a variety of acceleration schemes (e.g., numerically stable Benders cuts and a heuristic for finding an incumbent solution). The implementation is written in C++ language and utilizes an off-the-shelf solver, IBM CPLEX Optimizer. 

## Requirements
We recommend using Mac or Linux machines. To build the source code, users need to install IBM CPLEX Optimizer and provide paths to the libraries and header files of CPLEX using `UserConfig.cmake` file.

In addition, the following libraries are required and they can be installed via `apt-get` on Linux machines or `brew` on Mac machines, if they are not on your machine already:
- CMake for building the source code
- BLAS/LAPACK/bzip2/zlib: These packages are requirements for an external software package used in this implementation

## How to Install
1. Clone the repository on your machine (e.g., type `git clone --recursive https://github.com/gbyeon/bilevel.git` in a terminal window)
2. Provide paths to CPLEX libraries and header files using `UserConfig.cmake`
3. Build the source code by typing the following commands from the root directory in a terminal window:
```
mkdir build
cd build
cmake ..
make
```

## How to Use
The installation steps will create an executable file `./src/runbilevel` (assuming that you are in the `build` directory).
The executable file can be run as a command-line tool: 

Usage: `./runbilevel -f <instance path> [-t <time limit> -s <whether to use the numerically stable cut> -u <dual variable upper bound> -h <whether to use the heuristic for finding an incumbent solution> -ht <time limit on the heuristic>]
`

```
A required parameter
      -f        a path to a MPS file without file extension
```
```
Optional parameters
      -t        time limit in seconds (Default: 3600 sec)
      -s        1: use the numerically stable cut (Default); 0: otherwise
      -u        an upper bound on dual variables when you are using the standard Benders cut, not the numerically stable version (i.e., when you set -s=0). The standard Benders cut is sensitive to the dual variable bound, so you may adjust this parameter when the solution seems incorrect. It is set as 100000 by default.
      -h        1: use the heuristic for finding the incumbent solutions; 0: otherwise (Default)
      -ht       limit on the time to spend for the heuristic in seconds (Default: 150 seconds)
```
## Instances
Instances are available in the `instances` folder. A problem instance should provide two files with the same file name and different extensions---`mps` and `aux`; This is in accordance with a Bilevel Optimization Problem Library and details on the file format can be found [here.](https://coral.ise.lehigh.edu/data-sets/bilevel-instances/)

## Key Publications
Geunyeong Byeon and Pascal Van Hentenryck. [Benders Subproblem Decomposition for Bilevel Problems with Convex Follower.](https://arxiv.org/pdf/1801.03520) *Arxiv*, 2019
Geunyeong Byeon and Pascal Van Hentenryck. [Unit Commitment with Gas Network Awareness" *IEEE Transactions on Power Systems.](https://ieeexplore.ieee.org/abstract/document/8844828) 35.2 (2019): 1327-1339.

## Acknowledgements
This repository is initially developed as part of a work that was partly supported by an NSF CRISP Award (NSF-1638331).
