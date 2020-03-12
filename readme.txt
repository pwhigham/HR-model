This repository contains the basic code required to perform the Hill-Robertson experiments.

The files are:

fixedspace.cpp - the C code to produce the HR model.  Set to N=64.  Needs to be compiled for use on a 
computer, however, fixed64.exe will run on windows machines.

The other C programs implement the supplementary material used for analyzing the steady-state behaviour of the HR model. 

networks.R - R code to build the network matrices.  Assumes that you will place each network in a folder, and these are named 1.txt, 2.txt, etc.  fixedspace.cpp assumes this naming convention for the files.

Fig5*N64.bat - example windows batch files to run a panmictic/ring/lattice/sf1/sf2/star models (N=64) varying Nbeta and Nc.  Note this could have been written using a for loop, but at the time being explicit for each run made sense.  Note that the last argument to fixedspace.cpp indicates whether a header should be written to the file -- this is only set to 1 for the first run.  

For example, if you execute the batch file Fig5panN64.bat it will create the panmictic data for Figure 5 of the HR paper. 

The R file:  paperFigures.R has a large number of example functions that produce the figures (and additional material) for the paper.  Note that the default paths used in these functions assume that the folder containing the data is Spacefixed/output/*.  

** Note ** I will supply all the batch files if this is required for publication, although I do feel that if a researcher wants to validate the model they should build some tests, etc. for themselves.  The main aspect to validate is the c

  