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
  
  #print(c(X,Y,W,W1))
  subg = induced.subgraph(g, unique(c(X,Y,W,W1)))
  
  type <- vector(mode="character", length=length(V(subg)))
  type <- setNames(type, V(subg)$name)
  
  
  for(v in Y)
  {
    type[v] = "Y"
  }
 
  for(v in W1)
  {
    type[v] = "W neighbor"
  }

  for(v in W)
  {
    type[v] = "W"
  }
  
  type[X] = "X"
  network = asNetwork(subg)
  ggnet(network,node.group=type,size=4,segment.size=1)
}