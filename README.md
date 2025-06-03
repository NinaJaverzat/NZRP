# Monte Carlo algorithm for central-force rigidity percolation in dimension 2 

This project contains a code that implements the algorithm for rigidity percolation presented in the paper [title of the paper, link to it]. The algorithm is implemented on the triangular lattice. We tried to organize and comment the code in a pedagogical way.

## Files organization

The code that implements the algorithm is in the folder `code`. Output files are written in the folder `res`. The folder `conv` contains a second code, which is used to convolute the output with the binomial distribution. Results of the convolution are stored in the folder `conv_res`.

## Usage

#### Compilation

The code does not require any dependencies but a C++ compiler, the standard C++ library and GNU's `make`. The default compiler, specified in the `makefile` file, is g++, but the code can be compiled with icpc as well. The compilation has been tested with g++ >= 10 and icpc >= 19.

To compile the code enter the folder `code` and launch `make` 

```bash
cd code
make
```

Depending on the platform, it might be necessary to explicitly pass the makefile to `make`, which can be done as 

```bash
make -f makefile
```

To compile with icpc, change the second line of the makefile from

```bash
CXX = g++
```

to

```bash
CXX = icpc
```

By default, the makefile compiles the code without any flag enabled. To enable aggressive optimization, compile as 

```bash
make o=1
```

Note that in this case the `-march=native` flag is used.  If you want to compile and run the code on computers with different architectures, remove that flag from the makefile.

To compile with debug flags on, use
```bash
make g=1
```
and to compile with warnings on compile with
```bash
make w=1
```

Note that options can be combined. For example, 

```bash
make o=1 w=1
```
enables both warnings and optimization flags.

At the end of the compilation, the executable `rp.exe` is created.

#### Input

The program takes three mandatory inputs: `L`, `T` and `num`, all integers.

- `L` is the linear size of the triangular lattice 
- `T` is the number of trials, i.e., the number of realizations of the algorithm
- `num` is simply meant to enumerate simulations, so to perform several parrallel runs with the same settings


#### Output 

A simulation produces exactly two files: `LOG_Lval1_Tval2_numval3.txt` and `Smax_Lval1_Tval2_numval3.txt`, where `val1`, `val2` and `val3` are the values of `L`, `T` and `num` given in input. 

The `LOG` file contains basic information about the simulation: the seed of the random number generator, the linear size of the lattice (`L`), the number of nodes (`N=L*L`), the number of bonds (`M=3*N)`, the number of trials (`T`). It further contains `T` rows where the time (in seconds) to execute each trial is written. The last line contains the maximum memory used during the simulation (in KB).   

The `Smax` file contains two columns, which are the relative size of the largest CP and RP clusters. The file is made of `M` rows and each raw represents a value of the number of occupied bonds $m$, i.e., the first row contains the values for $m=1$, the second for $m=2$ and so on. **Importantly,** note that the values are not averaged: the two colums have to be dived by `T` to obtain the average over realizations. The reason for this choice is that if several parallel simulations are performed, all the data produced must be aggregated and only later normalized. 

#### Convolution

The algorithm produces results that are functions of the number of active bonds. To obtain results as a function of the bond concentration, the data must be convoluted (see [the paper by Newmann and Ziff](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.64.016706)). Results over several parallel simulations are aggregated and then convoluted by a second program, convolutor.exe. Once the simulations are ended, enter the folder `conv` and compile the convolutor:

```bash
cd ../conv
make
```

Note that the default behaviour is to compile with aggressive optimization enabled. Once the program is generated, run it by giving in input `L`, `T` and `TOT`, the latter being the total number of realizations performed. **Importantly**, the convolutor is written under the assumption that the first value of `num` is 1. So, if (say) 100 simulations are performed, the convolutor assumes they are numbered from `num=1` to `num=100` and `TOT=100` must be used. The convolutor outputs one file, named `Smax_Lval1_Tval2_TOTval3_convoluted.txt`, where `val1`, `val2` and `val3` are the three input respectively. The output file is stored in the folder `conv_res`. It contains four columns, which are: 

1) The average size of the largest CP cluster as a function of the number of active bonds. This is basically the aggregated average over all the simulations performed.
2) The average size of the largest CP cluster as a function of the bond concentration. This is the convolution of the first line with the binomial distribution.
3) The average size of the largest RP cluster as a function of the number of active bonds. This is basically the aggregated average over all the simulations performed.
4) The average size of the largest RP cluster as a function of the bond concentration. This is the convolution of the third line with the binomial distribution.

Values of the number of active bonds $m$ correspond to the line number, i.e., the first line is for $m=1$, the second for $m=2$ and so on. Values of the bond concentration are $p=m/M$.

#### runner

For simplicity, we provide with a simple runner that can be executed to perform `TOT` parallel simulations using `num_cores` cores. The runner takes in input `L`, `T`, `TOT` and `num_cores`. The first three inputs have the same meaning as before. The last input, `num_cores`, is used to have exactly `num_cores` simulations running in paraller. For example, executing

```bash
./runner 64 1000 100 10
```

will perform 100 simulations (`num=1` to `num=100`) of `T=1000` trials each, on a lattice with linear size `L=64` while keeping occupied no more than 10 cores with 10 parallel simulations. The runner creates the folders `res` and `conv_res`, enters the folder `code`, compiles the code with optimization enabled, makes the simulations and, when they are terminated, enters the folder `conv`, compiles the convolutor and executes it.

## Code organization

We have tried to logically split the code up in modules. We have also tried to use naming conventions that make it easier to understand the purpose of each module, function and (to a lesser extent) variable. The code is also heavily commented. 

- The `main.cpp` file initializes the simulation, creates output files and performs the main loop: one iteration of the main loop represents one trial, i.e., the lattice goes from empty to filled. It also writes the output of the simulation.

- The general data structures used throughout the code are stored in the `defs.h` file.

- The file `basic_functions.cpp` mainly contains initialization functions. 

- The file `pebble_game.cpp` contains functions used to play the pebble game. The file is divided in two major blocks. In the first block, the pebble game of type I is played: a path to a pebble is identified (we recall that in our algorithm the pebble game of type I always finds a pebble) and is reversed. By default, the pebble is searched using a Breadth First Search over the pebble graph. If Depth First Search is preferred, change line 79 of the main from `scalars.search_type = BFS;` to `scalars.search_type = DFS;`. In the second block of functions, the type II game is played. Also in this case, a path to a floppy node is searched with BFS by default. If a floppy nodes is found, the nodes that belong to the path to the floppy node gets marked floppy, otherwise all the visited nodes are marked rigid and triangulated over the base made by the new bond (where three pebbles are frozen).

- The file `percolation.cpp` contains the main logic of the algorithm. The function `single_trial` executes one trial of the algorithm. In particular, it contains all the functions to implement the Newmann-Ziff algorithm for connectivity percolation, the functions to implent the pivoting, rigidification and overconstraining steps of the algorithm and the functions to detect the wrapping of a cluster.

## Datastructures

- `scalars` contains scalar variables like the number of nodes and the number of bonds. Scalars also contains three sets that are used to record the roots of the bonds that wrap around the lattice.
- `Mixed_FIFO_LIFO` is a struct that can be used both as a stack and as a queue. It contains all the basic operations of stack and queues, but it assumes that the memory is pre-allocated. If the number of elements exceeds the allocated memory, a reallocation is needed.
- `CPSmax` and `RPSmax` are the arrays used to monitor the relative size of the largest CP and RP clusters respectively.
- `NZ` and `RNZ` are the vectors used to store the trees that monitor the cluster stuctures.
- `network` is the neighbor table of the lattice. Node that if the bond uv is not active, the corresponding entries in the neighbor table have negative value.
- `pebble_graph` is the vector where the pebble graph is stored.
- `searched_bonds` is used to reconstruct the path that leads to a pebble (type I game) or to a floppy node (type II game). Specifically, if the path is made of the pebble edge `i->j`, then the entry with index `j` has value `i`. In this way patsh can be backtracked (to reverse them, in type I game, or to mark the nodes as floppy, in type II game).
- `np` contains the number of pebbles of each node.
- `marks` and `marks_indices`: the first contains the marks (FLOPPY, RIGID, UNMARKED) of nodes. Everytime a series of searches of type II needs to be performed, the nodes are initially all UNMARKED, except the two nodes where three pebbles are frozen, i.e., `u` and `v`, which are marked RIGID. As searches are performed, new marks are put on nodes. Indiced of the nodes that get marked are stored in `marked_indices`, starting from its second element. The total number of marked nodes is stored in the first element of `marked_indices`. At the end of the rigidification step, once all the type II searches are done, the `marks` array is reset to all nodes being UNMARKED by iterating over `marked_indices`.
- `visited` and `visited_indices`: this pair of arrays works just like `marks` and `marked_indices`, but `visited` contains the nodes that are visited in a single pebble search (of both type). At the end of each search, visited gets re-initialized opportunely. 
- `enqueued` and `enqueued_indices`: also this pair of arrays works just like `marks` and `marked_indices`, but `enqueued` contains the pivots that enter in the stack of pivots needed for steps 5-8 of the algorithm.
- `P` contains the pivotal class of each node.
- `Prc` (Pivots Rigid Clusters) is the pivot network: each entry is a set and, if the bond `i` is a root, then `PRC[i]` contains the pivots of the rigid cluster with root `i`.
- `dx` and `dy` are used to monitor wrapping of rigid clusters. 
- `ROOTS` is used to check if a rigid clusters containing both the nodes `u` and `v` that make up the new bond exists.

## License

The code has been written by D. Notarmuzi and N. Javerzat. Feel free to contact us if you have questions/comments. If you use the code in your project please cite the publication corresponding to this work: _A fast algorithm for 2D Rigidity Percolation_. 
