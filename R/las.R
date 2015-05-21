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

getSubNet = function(graph, x, order)
{
  return(neighborhood(graph, order,x))
  
}


getCommunity <- function(z, g,cutoff, community.min)
{
  zcutoff = names(z[z>cutoff]) 
  if(length(zcutoff)==0)
  {
    return(NULL)
  }
  subg <- induced.subgraph(graph=g,vids=zcutoff)
  wc = walktrap.community(subg)

  member = membership(wc)
  w = names(member[member==which.max(sizes(wc))])
  if(length(w)>=community.min)
  {
    return(wc)
  }
  else
  {
    return(NULL)
  }
}

#' get GO
#' 
#' @export
#' 
getGO <- function(sel.entrez,all.entrez)
{
  params <- new("GOHyperGParams", geneIds=sel.entrez, universeGeneIds=all.entrez, ontology="BP", pvalueCutoff=0.01,conditional=F, testDirection="over", annotation="hgu133plus2.db")
  Over.pres = tryCatch({
    Over.pres<-hyperGTest(params)
  }, error = function(e) {
    return(NULL)
  })
  if(is.null(Over.pres))
  {
    return(NULL)
  }
  
  summary <- getGeneric("summary")
  ov<-summary(Over.pres)
  return(ov[ov$Size<1000&&ov$Size>5,])
}



gen.data <- function(ci, member, all.entrez)
{

    
  w = names(member[member==ci])

  sel.entrez<-w
  wgo = getGO(sel.entrez, all.entrez)
  if(is.null(wgo)||is.na(wgo$Pvalue)||length(wgo$Term)==0)
  {
    return(NULL)
  }
  else
  {
    wgo <- paste(wgo$Term, signif(wgo$Pvalue,digits = 5), sep=": ")
    wgo <- paste(wgo, collapse = '\n')
  }
  
  w<-paste(w, collapse = ' ')
  
  return(data.frame(w, wgo))
}






#' get GO
#' 
#' @export
#' 
getgobp <- function(graph, z.matrix, k=2, n.cores=4, cutoff=0.8, community.min=5)
{

  community = apply(z.matrix, 1, getCommunity, graph,cutoff,community.min)
  community = community[!sapply(community, is.null)]
  all.entrez<-colnames(z.matrix)


  cl <- makeCluster(n.cores, outfile="")
  registerDoParallel(cl)
  cat('loop begin\n')
  print(community)
  resulttable <- foreach(i=1:length(names(community)), .combine='rbind') %dopar%
  {
    x = names(community)[i]
    print(x)
    wc = community[[x]]
    member = membership(wc)
    
    community_index = names(sizes(wc)[sizes(wc)>5])
    

    
    
    
    sel.entrez<-x
    xgo = getGO(sel.entrez, all.entrez)
    
    print(length(xgo$Term))
    if(is.null(xgo)||is.na(xgo$Pvalue)||length(xgo$Term)==0)
    {
      return(NULL)
    }  
    else
    {
      xgo <- paste(xgo$Term, signif(xgo$Pvalue,digits = 5), sep=": ")
      xgo <- paste(xgo, collapse = '\n')
    }
    
    xk = neighborhood(graph,k,nodes=x)
    sel.entrez = as.character(unlist(xk))
    xkgo = getGO(sel.entrez, all.entrez)
    
    
    if(is.null(xkgo)||is.na(xkgo$Pvalue)||length(xkgo$Term)==0)
    {
      return(NULL)
    }
    else
    {
      xkgo <- paste(xkgo$Term, signif(xkgo$Pvalue,digits = 5), sep=": ")
      xkgo <- paste(xkgo, collapse = '\n')
    }
    print(community_index)
    w.result = do.call("rbind",lapply(community_index, gen.data, member, all.entrez))
    if(is.null(w.result))
    {
      return(NULL)
    }
    else
    {
      return(cbind(x, xgo, xkgo,w.result))
    }
    
  }


  stopCluster(cl)
  return(resulttable)

}






