#' visualize
#' @export
#'
visualize <- function(g,result, x, k, cutoff=0.8, path=NULL)
{
  
  X = as.character(x)
  Y = V(g)$name[unlist(igraph::neighborhood(g, 2, nodes=X))]
  
  z = result[X,]
  W = names(z[z>cutoff])

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

  network = asNetwork(subg)
  size = length(V(subg))
  scale = (17/961)*size + 1999/961
  
  output =  ggnet(network,node.group=type,segment.size=1,label.nodes=T,col="black",subset.threshold = 1
                  )
  ggsave(output, file=paste(as.character(x),".jpg",sep = ""), path=path,w=4, h=3, scale=scale,limitsize=FALSE)
}


