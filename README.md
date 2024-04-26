A collection of R functions for conditional quasi-likelihood inference for mean residual life regression with clustered survival data.

A set of clustered survival data with dependent censoring is included in "data.csv".

The main function is "CQL_main in R", which compares the estimate of the regression coefficient from the proposed conditional quasi-likelihood inference method 
with Buckley-James estimates to those from three alternative methods, conditional quasi-likelihood inference with Kaplan-Meier estimates, hierarchical 
quasi-likelihood method with lognormal frailty and the conventional method ignoring both within-cluster correlation and dependent censoring.

The mentioned methods above are implemented via Newton-Raphson iteration, corresponding to the functions "NRcqlbj in R", "NRcqlkm in R", "NRhqlln in R" and "NRnon in R".
