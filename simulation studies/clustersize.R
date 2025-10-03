generate_increasing_N_size <- function(Ncluster, Ntime, basesize = 10) {
  # Initialize an empty list to store the vectors
  vectors_mtrx <- matrix(NA, nrow = Ncluster, ncol = Ntime)
  
  add.i <- sample(5:10,3)
  # Generate N vectors with a bit of randomness
  for(i in 1:Ncluster) {
    vectors_mtrx[i,] <- basesize  + ifelse(i >= Ncluster/2, 1, 0)*1:Ntime + 
      ifelse(i < 3*Ncluster/4,
             ifelse(i < Ncluster/2,
                    ifelse(i < Ncluster/4, add.i[1],add.i[2]),add.i[3]), 
             0) #+ sample(0:random_size, 1) # Adds a random integer to each element
  }
  
  return(vectors_mtrx)
}
