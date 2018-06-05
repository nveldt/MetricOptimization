## README

This folder contains code for algorithms and experiments associated with the following paper:

A Projection Method for Metric-Constrained Optimization
(Nate Veldt, David Gleich, Anthony Wirth, and James Saunderson)

In particular this contains Julia implementations of DykstraCC and DykstraSC: projection-based methods for solving convex relaxations of correlation clustering and sparsest cut objective respectively.

There are two directories, one for experiments involving the sparsest cut relaxation, and one for the correlation clustering relaxations.

These are compared against running Gurobi optimization software. A free academic license for Gurobi can be obtained online at Gurobi.com.

Graphs can be obtained online at the Suitesparse Matrix Collection, and stored in the "Graphs" folder. You should make graphs unweighted and undirected, and take the largest connected component if you wish to reproduce experiments from the paper.

