#' visualize
#' @export
#'

visualize <- function(g,result, x, k, cutoff=0.8)
{

  X = as.character(x)
  Y = V(g)$name[unlist(igraph::neighborhood(g, 2, nodes=X))]
  
  z = result[X,]
  W = names(z[z>cutoff])
  W1 = V(g)$name[unlist(igraph::neighborhood(g, 1, nodes=W))]
  subg = induced.subgraph(g, unique(c(X,Y,W,W1)))
  
  type <- vector(mode="character", length=length(V(subg)))
  type <- setNames(type, V(subg)$name)
  
  
  for(v in W1)
  {
    type[v] = "W neighbor"
  }
  for(v in W)
  {
    type[v] = "W"
  }
  
  for(v in Y)
  {
    type[v] = "Y"
  }
  type[X] = "X"
  print(type)
  network = asNetwork(subg)
  ggnet(network,node.group=type,segment.size=1,label.nodes=T,col="white")
}

#' visualize
#' @export
#'
visualizewitoutw1 <- function(g,result, x, k, cutoff=0.8)
{
  
  X = as.character(x)
  Y = V(g)$name[unlist(igraph::neighborhood(g, 2, nodes=X))]
  
  z = result[X,]
  W = names(z[z>cutoff])
  #W1 = V(g)$name[unlist(igraph::neighborhood(g, 1, nodes=W))]
  
  #print(c(X,Y,W,W1))
  subg = induced.subgraph(g, unique(c(X,Y,W)))
  
  type <- vector(mode="character", length=length(V(subg)))
  type <- setNames(type, V(subg)$name)
  
  
  
  for(v in W)
  {
    type[v] = "W"
  }
  for(v in Y)
  {
    type[v] = "Y"
  }
  type[X] = "X"
  print(type)
  network = asNetwork(subg)

  ggnet(network,mode="target",node.group=type,size=20,weight.method = 'degree',segment.size=1,label.nodes=T,col="black")
}


