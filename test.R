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



library(GGally)
library(intergraph)
library(network)

testdata =relate.matrix['5700',]
target = names(testdata[testdata==1])

terget.graph = induced.subgraph(graph=g,vids=target)
network = asNetwork(terget.graph)
ggnet(network,size=4,segment.size=1)

network = asNetwork(g)
istarget = rownames(result) %in% target
ggnet(network, size=3,node.group=istarget,color="blue" ,)



testdata =newresult['5700',]
target = names(testdata[testdata>0.6])

target = names(testdata)[index.top.N(testdata,N=40)]

target

index.top.N = function(xs, N=10){
  if(length(xs) > 0) {
    o = order(xs, na.last=FALSE)
    o.length = length(o)
    if (N > o.length) N = o.length
    o[((o.length-N+1):o.length)]
  }
  else {
    0
  }
}

