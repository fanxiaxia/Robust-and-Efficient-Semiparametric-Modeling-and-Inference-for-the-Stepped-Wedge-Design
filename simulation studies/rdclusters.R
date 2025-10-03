randomize_clusters <- function(Ncluster, Time.vec) {
  
  # Shuffle the clusters randomly
  shuffled_clusters <- sample(1:Ncluster)
  
  # Assign clusters to seq_groups
  seq_assignments <- cut(seq_along(shuffled_clusters), length(Time.vec), labels = FALSE)
  
  # Create a dataframe with clusters and their assigned seqs
  df <- data.frame(Cluster = shuffled_clusters, cross =  Time.vec[seq_assignments])
  df <- df[order(df$Cluster), ]
  
  return(df)
}