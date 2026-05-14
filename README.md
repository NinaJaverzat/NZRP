# Monte Carlo algorithm for central-force rigidity percolation in dimension 2

This project contains the code implementing the Newman–Ziff-based algorithm for rigidity percolation presented in the paper *A fast algorithm for 2D Rigidity percolation*, https://arxiv.org/abs/2507.00741. The algorithm is implemented on the triangular lattice. The code is organised and commented in a pedagogical way.

## Repository layout

```
code/               Main simulator source files → rp.exe
res/                Raw simulation output (one file set per run)
aggregate/          Aggregator source → aggregate.exe
agg_res/            Aggregated output
conv/               Legacy Smax-only convolutor → convolutor.exe
fullrange_conv/     Generic full-range convolutor → fullrange_conv.exe
partialrange_conv/  Partial-range convolutor for Pinf/chi → partialrange_conv.exe
conv_res/           All convolution output
figures/            Plots produced by the notebook
NZRP_pipeline.ipynb Reproducibility notebook (compile → run → aggregate → convolve → plot)
```

---

## Dependencies

| Tool | Required by |
|---|---|
| C++20-capable compiler (g++ ≥ 10 or icpc ≥ 19) | all C++ programs |
| GNU `make` | all C++ programs |
| Boost.Math (`boost/math/distributions/binomial.hpp`) | `conv`, `fullrange_conv`, `partialrange_conv` |
| Python ≥ 3.9 with `scipy`, `pandas`, `matplotlib` | `NZRP_pipeline.ipynb` |

---

## 1 · Simulator (`code/`)

### Compilation

Enter the `code` folder and run `make`:

```bash
cd code
make
```

If `make` does not find the makefile automatically, pass it explicitly:

```bash
make -f makefile
```

Build options are controlled by single-letter variables:

| Variable | Effect |
|---|---|
| `o=1` | Enable aggressive optimisation (`-O3 -march=native -flto -ffast-math`). **Removes portability across CPU architectures** — remove `-march=native` from the makefile when cross-compiling. |
| `g=1` | Enable debug symbols (`-g -g2`) |
| `w=1` | Enable warnings (`-Wall -Wextra -pedantic`) |

Options may be combined, e.g. `make o=1 w=1`. The default compiler is `g++`; override it with `CXX=icpc` or any other C++20 compiler.

Compilation produces the executable `rp.exe`.

### Input

```
rp.exe <L> <T> <num>
```

| Argument | Meaning |
|---|---|
| `L` | Linear size of the triangular lattice. The lattice has `N = L²` nodes and `M = 3N` bonds. |
| `T` | Number of independent Monte Carlo trials (realizations). |
| `num` | Run index used only to label output files, allowing many parallel runs with the same `L` and `T`. |

### Output files

All output files are written to `../res/` relative to the `code` directory, i.e. to the `res/` folder in the repository root. Every file name follows the pattern `<STEM>_L<L>_T<T>_num<num>.txt`. The stems and their content are:

| File stem | Columns / content |
|---|---|
| `LOG` | Header with RNG seed, lattice parameters, and trial count; then `T` lines with per-trial wall time (s); final line with peak memory (KB). |
| `Smax` | `M` rows, one per bond number $m$. Two columns: cumulative (un-normalised) largest CP cluster size and largest RP cluster size. Divide by `T` to obtain averages. |
| `COP` | Header line `# m1 m2` giving the measured bond range. Then one row per $m \in [m_1, m_2]$ with three columns: cumulative CP order parameter (`Pinf` numerator), normalisation count, and cumulative CP susceptibility numerator. |
| `ROP` | Same layout as `COP` but for the rigidity percolation order parameter and susceptibility. |
| `Searches` | `M` rows. One column: cumulative number of pebble-game searches performed at each $m$. |
| `TimePerBond` | `M` rows. One column: cumulative wall time (s) spent processing each bond number $m$. |
| `Events` | `M` rows. Three columns: cumulative counts of pivoting events, rigidification events, and overconstraining events. |
| `Operations` | `M` rows. Three columns: cumulative counts of Type-I visited nodes, Type-II visited nodes, and pivot pushes. |
| `EventTimes` | `M` rows. Three columns: cumulative wall time (s) spent in pivoting, rigidification, and overconstraining routines. |
| `PERFLOG` | `T` rows, one per trial. Eleven columns: trial wall time; pivoting, rigidification, and overconstraining event counts; searches; Type-I visited; Type-II visited; pivot pushes; pivoting, rigidification, and overconstraining times. |
| `CPSxy` | One row per trial that produced a spanning CP cluster. Single column: size of the spanning cluster. |
| `RPSxy` | Same as `CPSxy` but for the rigidity percolation spanning cluster. |
| `WRAPS` | `T` rows. Six columns: bond number at which the CP cluster first wraps in X (`cp_mx`), bond number at which it wraps in XY (`cp_mxy`), corresponding RP values (`rp_mx`, `rp_mxy`), maximum CP susceptibility, and maximum RP susceptibility for that trial. |

**Note on normalisation:** the `Smax`, `COP`, `ROP`, `Searches`, `TimePerBond`, `Events`, `Operations`, and `EventTimes` files store *cumulative* (not averaged) values so that data from several parallel runs can be summed before dividing. The aggregator and convolutors handle this automatically.

---

## 2 · Aggregator (`aggregate/`)

The aggregator combines the raw per-run files from `res/` into a single summarised file written to `agg_res/`. It does not perform convolution.

### Compilation

```bash
cd aggregate
make          # add o=1 / g=1 / w=1 as needed
```

This produces `aggregate.exe`.

### Usage

```
aggregate.exe <L> <T> <FIRST_NUM> <LAST_NUM> <STEM> <MODE> [COLUMN]
```

| Argument | Meaning |
|---|---|
| `L`, `T` | Lattice size and trial count, matching the simulation parameters. |
| `FIRST_NUM`, `LAST_NUM` | Inclusive range of `num` values to aggregate. Missing files in the range are silently skipped. |
| `STEM` | File stem to read, e.g. `Searches`, `WRAPS`, `PERFLOG`. |
| `MODE` | Aggregation mode (see table below). |
| `COLUMN` | Zero-based column index, required for `samples` and `hist` modes. |

| Mode | Output suffix | Description |
|---|---|---|
| `curve` | `_avg.txt` | Averages every column of every row across all runs. Preserves the `# m1 m2` range header if present. |
| `samples` | `_col<C>_samples.txt` | Collects all values from column `C` across all rows and runs. Header line contains the mean and second moment. |
| `hist` | `_col<C>_hist.txt` | Same as `samples` but produces a histogram (value, count) instead of the raw list. |

Output files are named `<STEM>_L<L>_T<T>_num<FIRST>_to_<LAST>[_col<C>]_<suffix>.txt` and are written to `agg_res/`.

---

## 3 · Convolutors

The simulator measures observables as functions of the number of active bonds $m$. The Newman–Ziff convolution (see [Newman & Ziff, Phys. Rev. E 64, 016706 (2001)](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.64.016706)) converts them to functions of the bond concentration $p = m/M$ by convolving with the binomial distribution $B(m; M, p)$.

Three convolutors are provided for different use cases. All three require Boost.Math at compile time and write their output to `conv_res/`.

### 3a · Legacy Smax convolutor (`conv/`)

This is the original, narrowly-scoped convolutor for the `Smax` observable only.

```bash
cd conv
make
```

Produces `convolutor.exe`.

```
convolutor.exe <L> <T> <TOT>
```

`TOT` is the total number of runs. **The convolutor expects runs numbered from 1 to `TOT`**. It reads `res/Smax_L<L>_T<T>_num<i>.txt` for `i = 1 … TOT`, aggregates them, and produces:

```
conv_res/Smax_L<L>_T<T>_TOT<TOT>_convoluted.txt
```

Four columns: average CP $S_\text{max}(m)$, convolved CP $S_\text{max}(p)$, average RP $S_\text{max}(m)$, convolved RP $S_\text{max}(p)$.

### 3b · Full-range convolutor (`fullrange_conv/`)

A general-purpose convolutor that can handle any observable stored over the full bond range $m = 1 \ldots M$.

```bash
cd fullrange_conv
make          # add o=1 / g=1 / w=1 as needed
```

Produces `fullrange_conv.exe`.

```
fullrange_conv.exe <L> <T> <FIRST_NUM> <LAST_NUM> <STEM> <MODE> <NCOLS> <COL0> [COL1 …] [MODEL]
```

| Argument | Meaning |
|---|---|
| `FIRST_NUM`, `LAST_NUM` | Inclusive range of run indices. |
| `STEM` | Input file stem (read from `res/`). |
| `MODE` | Processing mode (see below). |
| `NCOLS` | Number of columns to process. |
| `COL0 …` | Zero-based column indices (exactly `NCOLS` values). |
| `MODEL` | Optional alternative stem used when naming CDF/wrap output files. Defaults to `STEM`. |

| Mode | Description |
|---|---|
| `pair` | Reads each column as a raw cumulative observable, normalises by $T \times \text{found\_files}$, and convolves the full range. Output has columns: $p$, avg($m$), conv($p$), repeated for each requested column. |
| `curve` | Identical to `pair`. |
| `cdf` | Interprets the two requested columns as the bond numbers at which the CP and RP wrapping clusters first appear (per trial). Builds empirical CDFs and also writes an `_averaged.txt` file pairing the two CDFs. |
| `wrap` | Alias for `cdf`. |

Output file: `conv_res/<OUTSTEM>_L<L>_T<T>_num<FIRST>_to_<LAST>_convoluted.txt` (and `_averaged.txt` for `cdf`/`wrap`).

### 3c · Partial-range convolutor (`partialrange_conv/`)

Specialised convolutor for the `COP` and `ROP` files, which only cover the bond range $[m_1, m_2]$ stored in their header. The convolution is performed only over this restricted range, which avoids extrapolation beyond the measured region.

```bash
cd partialrange_conv
make          # add o=1 / g=1 / w=1 as needed
```

Produces `partialrange_conv.exe`.

```
partialrange_conv.exe <L> <T> <FIRST_NUM> <LAST_NUM> <STEM> <PINF_COL> [NORM_COL CHI_COL]
```

| Argument | Meaning |
|---|---|
| `STEM` | File stem, e.g. `COP` or `ROP`. |
| `PINF_COL` | Zero-based column index of the order-parameter numerator. |
| `NORM_COL` | Zero-based column index of the normalisation count (default: 1). |
| `CHI_COL` | Zero-based column index of the susceptibility numerator (default: 2). |

The range $[m_1, m_2]$ is read automatically from the `# m1 m2` header of the first existing input file. The convolution stops before $p$ values where the binomial would require data beyond $m_2$.

Output files (in `conv_res/`):

- `<STEM>_L<L>_T<T>_num<FIRST>_to_<LAST>_averaged.txt` — averaged $P_\infty(m)$ and $\chi(m)$ over the measured range.
- `<STEM>_L<L>_T<T>_num<FIRST>_to_<LAST>_convoluted.txt` — convolved $P_\infty(p)$ and $\chi(p)$.

---

## 4 · Reproducibility notebook (`NZRP_pipeline.ipynb`)

`NZRP_pipeline.ipynb` is a Jupyter notebook that drives the entire workflow from a single place. It must be run from the repository root (the directory containing `code/`).

The notebook is divided into three parts:

**Part 1 – Simulation.** A configuration cell at the top defines `Ls` (list of system sizes), `Ts` (trials per size), `run_numbers` (how many parallel runs per size), and compiler options. Helper functions compile `rp.exe` via `make` and invoke the simulator for each `(L, T, num)` combination. A status table shows which output files have been produced.

**Part 2 – Aggregation and convolution.** Helper functions build `aggregate.exe`, `fullrange_conv.exe`, and `partialrange_conv.exe` on demand (rebuilding only when source files are newer than the executable). Aggregation and convolution specs are assembled automatically from the `Ls` / `Ts` / `run_numbers` lists and dispatched to the appropriate tool. The resulting file paths are collected in dictionaries for the plotting step.

**Part 3 – Display and plots.** Aggregated and convolved data are loaded with `pandas` and displayed with `matplotlib`. Publication-ready figures can be saved to `figures/`.

Python dependencies: `scipy`, `pandas`, `matplotlib`, `numpy`.

---

## 5 · Data structures

The following C++ structs are defined in `code/defs.h` and used throughout the simulator.

| Struct | Purpose |
|---|---|
| `Scalars` | Holds all scalar simulation parameters (`L`, `N`, `M`, `T`, current bond count `m`, seed, RNG state) and running accumulators for the CP and RP wrapping state and susceptibility extrema. |
| `OrderParam` | A bond-indexed observable: a `std::vector<long double> y` of length proportional to `M`, a companion `std::vector<unsigned long long> norm` for per-entry normalisations, and `start`/`stop` indices marking the measured range. |
| `LightOrderParam` | A simplified version of `OrderParam` without the per-entry normalisation vector, used for susceptibility accumulation. |
| `TrialData` | Stores per-trial performance counters (wall time, event counts, operation counts, sub-routine times) for logging to `PERFLOG`. |
| `WrapState` | Tracks wrapping cluster information during a trial: three `std::set<size_t>` for the roots of X-, Y-, and XY-wrapping clusters, and three `std::vector<WrappingTriplet>` holding the detailed records. |
| `WrappingTriplet` | A single record of a wrapping cluster: its root bond index (`root`), the bond number at which wrapping was detected (`bond_num`), and the cluster size at that moment (`size`). |
| `Mixed_FIFO_LIFO` | A flat-array container that can be used as either a FIFO queue (BFS) or a LIFO stack (DFS) via `push`/`pop`/`pop_back`/`front`/`top` operations. Used in the pebble-game search routines. |

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

The code has been written by D. Notarmuzi and N. Javerzat. Feel free to contact us if you have questions/comments. If you use the code in your project please cite 
D. Notarmuzi & N. Javerzat. (2025) https://doi.org/10.5281/zenodo.15584520
as well as the publication corresponding to this work: D. Notarmuzi & N. Javerzat. (2025) _A fast algorithm for 2D Rigidity Percolation_, https://arxiv.org/abs/2507.00741
