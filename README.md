## README

This folder contains code for algorithms and experiments associated with the following paper:

A Projection Method for Metric-Constrained Optimization
(Nate Veldt, David Gleich, Anthony Wirth, and James Saunderson)

https://arxiv.org/abs/1806.01678

In particular this repository contains Julia implementations of DykstraCC and DykstraSC: projection-based methods for solving convex relaxations of the correlation clustering and sparsest cut graph clustering objectives respectively.

The original set of experiments was run in Julia 0.6. The Julia-0.6 subdirectory contains code for algorithm implementations and all experiments in Julia 0.6.

The Julia-1.0 folder includes algorithm implementations for DykstraCC and DykstraSC that work with version 1.0 of Julia.

### Graphs

Graphs can be obtained online at the Suitesparse Matrix Collection (https://sparse.tamu.edu/), and stored in the "graphs" folder. You should make graphs unweighted and undirected, and take the largest connected component if you wish to reproduce experiments from the paper. The folder already containes three real-world networks for you to test the algorithms on.

### Important Julia Packages

For the experiments, we use JuMP and Gurobi packages in Julia, which you will need to have installed. Note that you also need to obtain a license at Gurobi.com in order to run Gurobi. 

For some convergence plots we use Plots.jl, but this isn't necessary for running the algorithm itself.

### Matlab Plots

Many of our plots are also produced using Matlab. Examples are given. We include data for running our algorithms on the graphs described in the paper. You can produce convergence plots for any of the experiments in the paper with the code given here.

All of these plots could also be generated using Plots.jl.


### LambdaCC Code

The correlation clustering folder includes code for solving the LP relaxation of the LambdaCC objective function, introduced in the paper "A Correlation Clustering Framework for Community Detection" (Nate Veldt, David Gleich, and Anthony Wirth).