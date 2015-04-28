library(igraph)
library(microbenchmark)
library(fdrtool)
network_list = as.matrix(read.table("/Users/danielqiu/Workspace/Bio/HumanBinaryHQ_HINT.txt"))

load("/Users/danielqiu/Workspace/Bio/GSE18864_entrez_norm.bin")

g = graph.data.frame(network_list, directed=FALSE)

las(g,b)

e = E(g)


lasr = function(x,y,z)
{
  sum(x*y*z)
}

x = sample(1:100,100,replace=TRUE)
y = sample(1:100,100,replace=TRUE)
z = sample(1:100,100,replace=TRUE)

microbenchmark(
     lascore(x,y,z),
     lasr(x,y,z)
   )


x = b['1',]
y = b['310']
z = b['780']




for(i in 1:500)
{
  if(length(h[[i]])!=0 )
  {
    print(h[[i]][[1]])
  }
  
}
