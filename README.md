# Computational codes for simulation studies in the manuscript "Semiparametric Mixed-effects Ordinary Differential Equation Models with Heavy-Tailed Distributions" by Baisen Liu, Yunlong Nie, Liangliang Wang, and Jiguo Cao. Accepted by Journal of Agricultural, Biological and Environmental Statistics in 2021.

SemiMEODE2021Jan31.pdf is the main manuscript. 

supplement.pdf is the supplementary document.

This foler also contains 5 matlab files.

main.m is the main Matlab program to simulate the Scenario I and II in Section 5 of the article "Bayesian Robust Estimation of Semiparametric Mixed-Effects Ordinary Differential Equation Models Using Heavy-Tailed Distributions" by Liu, Wang, Nie and Cao.

This Matlab file will use the other four files:

Sim1Ode.m is the true semiparametric ODEs function to generate the simulated data.

Sim1WorkingOde.m is the working semiparametric ODEs function where the time-varying function eta(t) is approximated by a linear combination of B-splines basis functions.

Sim1MCMC_funNN.m is the MCMC algorithm assuming that the random-effects and errors follow Gaussian distributions.

Sim1MCMC_funTT.m is the MCMC algorithm assuming that the random-effects and errors follow heavy-tailed distributions.

fdaM.zip is the zipped folder of the FDA package. Our program needs to use this package. Please unzip this file before running our program.
