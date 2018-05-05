## README

This is code that was used for Experiment 2 in the WWW 2018 paper. We use both the triangle-fixing algorithm as well as the lazy constraints method to get LP lower bounds on a 1000-node BTER graph, a 1000 node LFR graph (not shown in paper), and the largest component of ca-GrQc.


## Method 1: lazy constraints
For using Lazy constraints, run:

Cagrqc_lazyLP.jl 
bter1000_lazyLP.jl

and the results will be stored in

cagrqc_conn_lazyLP_bounds.txt and
bter

btr1000_lazyLP_bounds.txt

respectively.

Note that these start at the highest values of Lambda and work their way down (easiest to hardest), but you should not expect them to converge for the tiniest values of Lambda. We will need to run the triangle fixing algorithm for the smallest values.

## Method 2: triangle fixing algorithm

For the triangle fixing algorithm, we run this for very tiny lambda, and do it one experiment (one lambda) at a time. 

The file to use is cagrqc_weighted_tfix.jl. Do the following:

1. Open julia in the terminal
2. include("cagrqc_TrifixLP_best.jl")
3. run_cagrqc("cagrqc_lam77",0.000077426368268,.1, .1, 20.0)

The above gives the filename for storing output, the lambda value, the change tolerance (we stop the algorithm aftet the solution vector changes by a small enough amount), and triangle tolerance (we also have to satisfy triangle inequality constraints to at least .1), and then set gamma = 20, which should be small enough for any value of lambda we choose.

## Lambda values

We use these lambda values for cagrqc:

   0.000010000000000
   0.000027825594022
   0.000077426368268
   0.000215443469003
   0.000599484250319
   0.001668100537200
   0.004641588833613
   0.012915496650149
   0.035938136638046
   0.100000000000000
   0.150000000000000
   0.250000000000000
   0.350000000000000
   0.450000000000000
   0.550000000000000
   0.650000000000000
   0.750000000000000
   0.850000000000000
   0.950000000000000
   
This is logspace(-5,-1,10) and then .15:.1:.95.

Note that we don't have to run 1e-5, because I confirmed using Gurobi that this is smaller than the minimum scaled sparsest cut of cagrqc (I'll have to get together that code later).

