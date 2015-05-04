#' Find the liquid association scouting gene
#' 
#' @param network.graph An igraph object representing the gene network.
#' @param express.matrix A matrix represeting the express matrix for the genes in gene network.Row names are the gene id in gene network.
#' @param k Integer giving the order of the network.
#' @param n.cores core number used for parallel computing
#' @return A logical matrix representing the LA-scouting genes for each gene. Rows represent the center gene id and columns represents the LA-scouting genes.
#' @export
#' 
lascouting <- function(network.graph, express.matrix, k=2, n.cores=4){
  
  network.node <- V(network.graph)$name
  matrix.node <- row.names(matrix)
  if(!identical(intersect(network.node,matrix.node),union(network.node,matrix.node))){
    common.node <- getCommonNode(network.graph, express.matrix)
    network.graph <- cleanGraph(network.graph, common.node)
    express.matrix <- cleanMatrix(express.matrix, common.node)
  }
  size <- length(common.node)
  normalizeInputMatrix(express.matrix)
  
  if(k!=1)
  {
    graph.connected <- connect.neighborhood(network.graph,k)
    connected.list <- as.matrix(get.edgelist(graph.connected))
  }
  else
  {
    connected.list <- as.matrix(get.edgelist(network.graph))
  }
  row.size <- nrow(connected.list)
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  

  express.matrix.t <- t(express.matrix)/75
  ptm <- proc.time()

  cat("loop begin\n")

  result <- foreach(i=1:row.size) %dopar%
  {

    xy <- express.matrix[connected.list[i,1],]*express.matrix[connected.list[i,2],]
    
    la.vector <- c(xy%*%express.matrix.t)

    lfdr <- fdrtool(la.vector, verbose=FALSE, plot = FALSE)$lfdr

    return(rownames(express.matrix)[which(lfdr<0.2)])

  }
  stopCluster(cl)
  print(proc.time() - ptm)
  
  node.z <- Matrix(0, nrow = size, ncol = size,dimnames=list(rownames(express.matrix),rownames(express.matrix)))
  for(i in 1:row.size)
  {
    if(length(result[[i]])!=0 )
    {
      x = connected.list[i,1]
      y = connected.list[i,2]
      node.z[x,c(result[[i]])] = 1
      node.z[y,c(result[[i]])] = 1
    }
  }
  
  result = node.z#z.kernel.density(node.z, network.graph)
  return(result)
  
}



#' Evaluate the result using kernel density estimation.
#' 
#' @param relate.matrix The matrix returned by lascouting.
#' @param network.graph The igraph object representing the gene network.
#' @param smoothing.normalize Different ways to normalize the result.
#' @return A  matrix representing the kernel density of each gene. Each row is a gene, columns are the weights of scouting genes for the gene.
#' @export
#' 
z.kernel.density <- function(relate.matrix, network.graph, smoothing.normalize=c("one","squareM","none") ) {
  smoothing.normalize <- match.arg(smoothing.normalize)
  
  weight0 = dnorm(0)
  weight1 = dnorm(1)
  weight2 = dnorm(2)
  
  size = nrow(relate.matrix)
  

  relate.matrix = relate.matrix[order(rownames(relate.matrix)), ] 
  relate.matrix = relate.matrix[,order(colnames(relate.matrix)) ] 
  
  



  adjacency1 <- get.adjacency(network.graph, type="both")
  adjacency2 <- get.adjacency(connect.neighborhood(network.graph,2), type="both")-adjacency1

  adjacency1 = adjacency1[order(rownames(adjacency1)), ] 
  adjacency1 = adjacency1[,order(colnames(adjacency1)) ] 
  adjacency2 = adjacency2[order(rownames(adjacency2)), ] 
  adjacency2 = adjacency2[,order(colnames(adjacency2)) ] 
  
  weight.matrix = diag(size)*weight0 +adjacency1*weight1+adjacency2*weight2

  if(smoothing.normalize=="one")
  {
    rsum = rowSums(as.matrix(weight.matrix))
    nmatrix = diag(1/rsum)
    colnames(nmatrix) <- rownames(weight.matrix)
    weight.matrix = nmatrix %*% weight.matrix

  }
  else if(smoothing.normalize=="squareM")
  {
    temp.weight.matrix = as.matrix(weight.matrix)
    rsum = rowSums(temp.weight.matrix)
    m = rowSums(temp.weight.matrix != 0)
    nmatrix = diag(sqrt(m)/rsum)
    colnames(nmatrix) <- rownames(weight.matrix)
    weight.matrix = nmatrix%*%weight.matrix
    
  }

  result <- relate.matrix%*%weight.matrix

  
  return(result)
  
}





getCommonNode <- function(network.graph, matrix){
  network.node <- V(network.graph)$name
  matrix.node <- row.names(matrix) 
  common.node <- intersect(network.node, matrix.node)
  return(common.node)
}

cleanGraph <- function(network.graph, remain){
  network.node <- V(network.graph)$name
  delete <- setdiff(network.node, remain)
  network.graph <- delete.vertices(network.graph, delete)
  return(network.graph)
}

cleanMatrix <- function(express.matrix, remain)
{
  express.matrix <- express.matrix[rownames(b)%in%remain,]
  return(express.matrix)
}

lascore = function(x,y,z)
{
  sum(x*y*z)/length(x)
}








