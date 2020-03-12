##################################################
#
# Create networks for HR model.
# Save as files for use with C version of HR paper
#
##################################################
library(MASS)
library(igraph)

################### WEIGHTED GRAPHS ############################
#
# Generate fixed graphs with random outward weights that sum to 1.
# This is equivalent to Liebermann's defn. - the rows should 
# sum to 1 for each location so that it is a stochastic matrix.
#
# For regular graphs (panmictic, ring, lattice) the matrix is
# doubly-stochastic.  For other networks the matrix is right-stochastic.
#
# The networks are created with names 1.txt, 2.txt, ....
# and placed in a folder as specified by the input parameter <folder>.
# NOTE that the code assumes this folder exists.
# This naming is used by the C program that does the HR modelling.
# For regular networks 1 definition is only required.
################################################################

###############
## PANMICTIC
###############
ds.weighted.panmictic <- function(N=64,num.files=1,folder="pan64/")
{
  res <- matrix(nrow=N,ncol=N)
  
  for (fn in 1:num.files)
  {
    for (i in 1:N)  # for each row...
    {
      val <- 1/N
      res[i,] <- rep(val,N)
    }
    filename <- paste(folder,fn,".txt",sep="")
    write.matrix(res,file=filename)
  }
}
ds.weighted.ring <- function(N=64,leftnodes = 2,
                             num.files=1,
                             folder="ring64/")
{
  
  for (fn in 1:num.files)
  {
    g <- watts.strogatz.game(1,N,leftnodes,0.0) # Create ring
    network <- as.matrix(as_adjacency_matrix(g)) # Fixed for now
    diag(network) <- 1  # Note includes self
    
    for (i in 1:N)
    {
      cols <- which(network[i,]==1)
      val <- 1/length(cols)
      network[i,cols] <- rep(val,length(cols))
    }
    filename <- paste(folder,fn,".txt",sep="")
    write.matrix(network,file=filename)
  }
}
ds.weighted.star <- function(N=64,num.files=1,
                             folder="star64/")
{
  
  for (fn in 1:num.files)
  {
    network <- matrix(nrow=N,ncol=N)
    network[,] <- 0
    network[1,] <- 1
    network[,1] <- 1
    diag(network) <- 1
    for (i in 1:N)
    {
      cols <- which(network[i,]==1)
      vals <- 1/length(cols)
      network[i,cols] <- rep(vals,length(cols))
    }
    filename <- paste(folder,fn,".txt",sep="")
    write.matrix(network,file=filename)
  }
  
}

ds.weighted.lattice <- function(dimvector=c(8,8),N=64,
                                num.files=100,
                                folder="lattice64/")
{
  
  for (fn in 1:num.files)
  {
    g <- make_lattice(dimvector=dimvector,
                      directed=FALSE,
                      circular=TRUE)
    network <- as.matrix(as_adjacency_matrix(g)) # Fixed for now
    diag(network) <- 1
    
    for (i in 1:N)
    {
      cols <- which(network[i,]==1)
      vals <- 1/length(cols)
      network[i,cols] <- rep(vals,length(cols))
    }
    filename <- paste(folder,fn,".txt",sep="")
    write.matrix(network,file=filename)
  }
}
######################################################
# This is the scale-free network generator.
#####################################################
# Note we now require many examples since the process is stochastic.
# The linear model is created with power = 1.0
# Defaults to the power=2 model which is the more standard barabasi model.
############################################################################
ds.weighted.sf <- function(power=2.0,N=64,num.files=100,
							folder="sf64p2/")
{
  
  for (fn in 1:num.files)
  {
    g <- sample_pa(n=N,power=power,directed=FALSE)
    while(!is_connected(g)) g <- sample_pa(N,power=power,directed=FALSE)
    
    #plot(degree_distribution(g))
    network <- as.matrix(as_adjacency_matrix(g)) # Fixed for now
    diag(network) <- 1
    
    for (i in 1:N)
    {
      cols <- which(network[i,]==1)
      vals <- 1/length(cols)
      network[i,cols] <- rep(vals,length(cols))
    }
    filename <- paste(folder,fn,".txt",sep="")
    write.matrix(network,file=filename)
  }
}