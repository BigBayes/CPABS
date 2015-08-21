# Cancer Phylogenetic Reconstruction with Aldous' Beta Splitting (CPABS)


## Installation:

To install CPABS, you will need the Julia technical language.  After installing this, open the Julia shell and type:

```
Pkg.init()
```

This will create a julia package directory, e.g. `~/.julia/v0.3` if you installed julia 0.3.  Next, clone the CPABS git repository into your julia package directory, for example

```
cd ~/.julia/v0.3
git clone https://github.com/leviboyles/cpabs.git
```

Next, you will need to install all required Julia packages, as well as any optional ones you need.  This is done from the Julia shell, e.g.:

```
Pkg.add("ArgParse")
```

### Required Packages

ArgParse
OptionsMod
HDF5
Distributions
NumericExtensions

### Optional Packages

For real time diagnostics plotting: Winston
For phylogeny plotting: Compose, Color

After all packages are installed you should be ready to use CPABS.

## Usage

You can use CPABS from within Julia or from the command line.

### Within Julia
First, you will need to import the CPABS module:
```
using CPABS
```

The main tool for running an MCMC chain is `run_phylo_experiment` 
```
CPABS.run_phylo_experiment(inputfile, alpha, multilocus_filename,
                           wl_K_boundaries=wl_K_boundaries,
                           num_iterations=num_iterations,
                           outputfile=outputfile,
                           init_K=init_K)
```
The first 3 arguments are required, although `nothing` may be entered for `multilocus_filename`.  `alpha` is as described in the CPABS paper.

`num_iterations` is the number of MCMC iterations.

`wl_K_boundaries` is a vector of floats indicating the boundaries of the Wang-Landau partition on the number of clones `K` used for rjmcmc.  For example, `[3,5,Inf]` specifies the partition `(1:2, 3:4, 5:Inf)`.

`outputfile` is the name and path to write the resulting chain.

`init_K` is the number of clones in the initial MCMC state.

#### Evaluation
There are several tools for evaluating the chains produced by CPABS.  CPABS uses the Wang-Landau MCMC scheme, where each sample generated is weighted.  These weights *must* be used in computing expectations or other statistics from the chains in order to provide meaningful results.

The data are required to compute the weights for each sample.  The data can be loaded into the appropriate format using `constructDataState()`:

```
data = CPABS.constructDataState("CLL077.csv")
```

We also require the chain itself:
```
models = CPABS.loadModels("your_chain.models");
```

We can then compute expectations with respect to arbitrary functions of the model states (this automatically computes and uses the WL weights):
```
CPABS.compute_expectation(data, models, foo)
```

We can also simply apply a function to all states of the chain (we will need to explicitly compute the weights in this case):
```
wl_weights = get_chain_WL_weights(data, models)
bars = apply_func_to_chain(models, bar) 
```

The functions `foo` and `bar` take in an object of type `ModelState` can return an object of arbitrary type.  However, to use `compute_expectation`, the returned type needs to have scalar multiplication and addition operators defined in Julia.

`ModelState` has two entries that are particularly relevant for evaluation.  If `model` has type `ModelState`, then
```
model.tree
model.Z
```
are the relevant fields.  `model.tree` is an array of `TreeNode`s, giving a fairly standard binary tree structure whose definition is in `model/tree.jl`, which also includes several utilities for manipulating/analyzing trees.  `model.Z` is the assignment vector, if there are `M` mutations, then `model.Z` has `M` integer entries indexing the node of `model.tree` to which that mutation belongs.

There are also a few shortcuts for computing the cocluster and ancestor/descendant/sibling matrices:
```
ccmat = CPABS.compute_cocluster_matrix("your_chain.models", data)
cadsmat = CPABS.compute_cluster_ancestor_descendent_sibling_matrix("your_chain.models", data)
```
 
### From the Command Line
A command line tool `cpabs_cmd.jl` is also provided for convenience for easily running large batches of MCMC runs.  It is essentially a wrapper around `run_phylo_experiment`.  Usage example, outputting a 30,000 iteration run to ``my_chain.models``:
```
julia cpabs_cmd.jl -o my_chain.models -n 30000 CLL077.csv
```
