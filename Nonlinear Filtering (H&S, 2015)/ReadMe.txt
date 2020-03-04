<Readme file for Nonlinear Filtering>
- Last modified: 2/24/2016


main.m : This file replicates the result in Table 8.1. of Herbst and Schorfheide (2016). 
It approximates the log likelihood of small NK model described in Chapter 1-2 of the textbook. 
 
These files use a number of procedures collected in the folder "Mfiles". The most relevant files are


	a) model_solution.m  : Takes as inputs the vector of structural parameters. Returns the coefficient 
                               matrices of the log-linear approximate solution of the DSGE model. Please, 
                               refer to Sims () for details on the procedure.

	b) sysmat.m          : Takes as inputs the solution of the DSGE model and the vector of structural parameters.
                               Returns the matrices for the state space representation. 

	c) PF_lik.m          : Takes as inputs the matrices for the state space representation, number of particles, data, 
                               initial prior for the filter, indicator for bootstrap particle filter method. Returns the 
                               period log likelihood, updated particles, and effective sample size.   
	


The folder contains a subfolder "LRE" including files to solve the linear rational expectation model. 

us.txt is the data file whose first column is output growth, second column is inflation, third column is federal fund rates. 
The observations used in the estimation range from 1983:I to 2002:IV, giving a total of 80 observations. 
See Appendix for detailed definition of the variables.


	