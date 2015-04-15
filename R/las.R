
getCommonNode <- function(graph, matrix){
  graph.node = V(graph)$name
  matrix.node = row.names(matrix) 
  common.node = intersect(graph.node, matrix.node)
  return(common.node)
}

cleanGraph <- function(graph, remain){
  graph.node = V(graph)$name
  delete = setdiff(graph.node, remain)
  graph = delete.vertices(graph, delete)
  return(graph)
}

cleanMatrix <- function(matrix, remain)
{
  matrix = matrix[rownames(b)%in%remain,]
  return(matrix)
}

las <- function(graph, express.matrix, k=2, width=2, n.cores=1, need.normalize=TRUE){
  graph.node = V(graph)$name
  matrix.node = row.names(matrix)
  if(!identical(intersect(graph.node,matrix.node),union(graph.node,matrix.node))){
    common.node = getCommonNode(graph, express.matrix)
    graph = cleanGraph(graph, common.node)
    express.matrix = cleanMatrix(express.matrix, common.node)
  }
  
  if(need.normalize)
  {
    normalizeMatrix(express.matrix)
    cat("normalized")
  }
  
  
  
}

