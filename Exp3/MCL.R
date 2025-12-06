#' Perform Markov Clustering (MCL) on a Graph
#'
#' @description
#' This function implements the Markov Clustering (MCL) algorithm for detecting
#' communities (clusters) in a graph. MCL simulates random walks within the graph
#' by alternating between two operations: expansion and inflation. It is particularly
#' efficient for biological networks.
#'
#' @param g An \code{igraph} object. The graph to be clustered. It can be directed or undirected.
#' @param inflation Numeric. The inflation parameter. Default is 2.5.
#' @param max_iter Integer. The maximum number of iterations to perform if convergence is not reached. Default is 100.
#' @param pruning Numeric. A threshold for pruning small values in the matrix to zero.
#'   This preserves the sparsity of the matrix and significantly speeds up computation
#'   while saving memory. Default is 1e-5.
#' @importFrom igraph as_adjacency_matrix is_igraph V
#' @importFrom Matrix Diagonal colSums drop0
#' @importFrom methods as
#' @return An igraph object containing MCL clustering labels.
#' @examples
#' \dontrun{
#' library(igraph)
#' g <- make_graph("Zachary")
#' g <- run_MCL(g, inflation = 2)
#' print(head(V(g)$MCL_cluster))
#'
#' # Visualize
#' plot(g, 
#' vertex.color = V(g)$MCL_cluster,
#' vertex.size = 15,
#' vertex.label = V(g)$name)
#' }
#' @export
run_MCL <- function(g,
                    inflation = 2.5,
                    max_iter = 100,
                    pruning = 1e-5){
  
  stopifnot(igraph::is_igraph(g))
  
  ## create adjacency matrix
  adj <- igraph::as_adjacency_matrix(g, sparse = TRUE)
  
  ## add self-loops
  M <- adj + Matrix::Diagonal(nrow(adj))
  
  ## scale initially (Column Normalization)
  col_sum <- Matrix::colSums(M)
  col_sum[col_sum == 0] <- 1 # Avoid division by zero for isolated nodes
  
  inv_diag <- Matrix::Diagonal(x = 1/col_sum)
  M <- M %*% inv_diag # faster
  
  ## MCL iteration
  for (i in 1:max_iter){
    M_prev <- M
    
    # 1. expansion (Matrix Multiplication)
    M <- M %*% M
    
    # 2. inflation (Element-wise power)
    M <- M ^ inflation
    
    # 3. pruning (Keep matrix sparse)
    M[M < pruning] <- 0
    M <- Matrix::drop0(M) # Physically remove zeros from storage
    
    # 4. re-scale (Re-normalize columns)
    col_sum <- Matrix::colSums(M)
    col_sum[col_sum == 0] <- 1
    inv_diag <- Matrix::Diagonal(x = 1/col_sum)
    M <- M %*% inv_diag
   
    # 5. check convergence
    diff <- sum((M - M_prev)^2)
    
    if (diff < 1e-5) {
      cat(sprintf("Converged at iteration %d (Diff: %.2e)\n", i, diff))
      break
    }
  }
  
  ## get results
  clusters <- apply(M, 2, which.max)
  
  # Assign node names
  igraph::V(g)$MCL_cluster <- as.integer(clusters)
  
  return(g)
}
