isomap <- function(data_csv) {
  dist_vec <- dist(data_csv, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
  size = length(data_csv[,1])
  
  matrix <- matrix(Inf, nrow = size, ncol = size)
  for (i in 1:(size - 1)) {
    for (j in (i + 1):size) {
      matrix[i, j] = dist_vec[size*(i - 1) - i*(i - 1)/2 + j - i]
      matrix[j, i] = matrix[i, j]
    }  
  }
  
  knn <- 5
  nn = matrix(0, size, knn)
  mtx <- matrix
  
  for (i in 1:size) {
    matrix[i, i] <- Inf
  }
  
  for(i in 1:size){
    nn[i,] <- order(mtx[i,])[1:knn]
  }
  
  for (i in 1:length(nn[,1])) {
    matrix[i, -nn[i,]] <- Inf
  }
  
  for (i in 1:size) {
    matrix[i, i] <- 0
  }
  
  for (k in 1:size) {
    for (i in 1:size) {
      for (j in 1:size) {
        oldvalue <- matrix[i, j]
        candidate <- matrix[i, k] + matrix[k, j]
        if(oldvalue > candidate) {
          matrix[i, j] <- candidate
          if (matrix[j, i] > candidate) {
            matrix[j, i] <- candidate
          }
        }
      }
    }
  }
  
  maxind <- which(matrix == max(matrix), arr.ind = T)
  matrixmax <- matrix[maxind[1], maxind[2]]
  print(matrixmax)
  
  for (i in 1:size) {
    for (j in 1:size) {
      if (matrix[i,j] == Inf) {
        matrix[i,j] <- matrixmax * 2
      }
    }
  }
  
  fit <- cmdscale(matrix, eig = FALSE, k = 2)
}