# mvot
code for "A Likelihood-based Approach for Multivariate One-Sided Tests With Missing Data"


################
=== Steps to replicate results of real data analysis
################

1. Run the r code in "real_data_analysis.r" to prepare the input data files required by the C++ program.

2. Compile/link "header.h" and "real_data_analysis.cpp" to obtain a C++ program. Make sure the "scythe" library (downloadable at http://scythe.berkeley.edu/) and the "dlib" library (downloadable at http://dlib.net/) are both included in the inclusion path.  

3. Run the C++ program to obtain a text file that contains the results of the real data analysis. 

Note: steps 2 and 3 are repeated for each of the four intervention groups.

################
=== Steps to replicate results of simulation studies
################

1. Set up the simulation configuration file "simulation_configurations.txt" ( or "simulation_configurations_dim5.txt" for p=5)

2. Compile/link "header.h" and "simulation.cpp" ( or "simulation_dim5.cpp" for p=5) to obtain a C++ program. Make sure the "scythe" library (downloadable at http://scythe.berkeley.edu/) and the "dlib" library (downloadable at http://dlib.net/) are both included in the inclusion path.  

3. Run the C++ program to detailed simulation results.


################
=== file for Gaussian quadrature points and weights
################

The file "GH100.txt" contains the quadrature points and weights from the R package fastGHQuad. It is obtained in R by running the following code:

library(fastGHQuad); 
rule100 = gaussHermiteData(100); GH100 = data.frame(x=rule100$x,w=rule100$w)
write.table(GH100,"GH100.txt",row.names=FALSE,col.names=FALSE)
