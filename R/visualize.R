#' visualize
#' @export
#'
visualize <- function(g, x, k, w)
{
  

  subg <- induced.subgraph(graph=g,vids=x)
  wc = walktrap.community(subg)
  network = asNetwork(subg)
  ggnet(network,size=4,segment.size=1)
  
}