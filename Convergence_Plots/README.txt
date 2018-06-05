This folder contains all the data from running DykstraSC and DykstraCC on several real world graphs. Matlab code is given for showing convergence plots for all the experiments in our paper on “Efficient Solvers for Metric-Constrained Optimization Problems”.

To see plots, open and run the file Plots_for_CC.m or Plots_for_SCut.m

Notes

* As described in the paper, DykstraSC was run with lambda = 1/n and gamma = 5 for all graphs except Vassar85, for which we used lambda = 1/100 and gamma = 2.
* The text_data folder contains the text output from the experiments
* The mat_data folder contains .mat files storing the progress of the algorithm 
* The CC_ouput contains text files from running DykstraCC

The constraint satisfaction plots all look like step functions. This is because in the algorithm we only perform a full constraint satisfaction check every 10-20 iterations. Performing a full check every iteration is expensive.