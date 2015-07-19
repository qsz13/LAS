library(igraph)
library(microbenchmark)
library(fdrtool)
network_list = as.matrix(read.table("~/LAS/HumanBinaryHQ_HINT.txt"))

load("/Users/danielqiu/Workspace/Bio/GSE18864_entrez_norm.bin")
load("~/LAS/GSE10255_entrez.bin")
g = graph.data.frame(as.matrix(read.table("~/LAS/HumanBinaryHQ_HINT.txt")), directed=FALSE)

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



testdata =result['5700',]
target = names(testdata[testdata>0.6])
terget.graph = induced.subgraph(graph=g,vids=target)


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


testz = result['1000',]
cutoff=0.8
w = names(testz[testz>cutoff])
subg <- induced.subgraph(graph=g,vids=w)
wc = walktrap.community(subg)
network = asNetwork(subg)
ggnet(network,size=4,segment.size=1)


plot(subg)

library(GOstats)


library(GOstats)

sel.entrez<-rownames(b)[1:5]
all.entrez<-rownames(b)
params <- new("GOHyperGParams", geneIds=sel.entrez, universeGeneIds=all.entrez, ontology="BP", pvalueCutoff=0.01,conditional=F, testDirection="over", annotation="hgu133a.db")
Over.pres<-hyperGTest(params)
ov<-summary(Over.pres)
ov$Term





for(c in com)
{
  print(names(com))
}

  sel.entrez = as.character(unlist(neighborhood(g,2,c('10148'))))
  all.entrez<-rownames(result)
  params <- new("GOHyperGParams", geneIds=sel.entrez, universeGeneIds=all.entrez, ontology="BP", pvalueCutoff=0.01,conditional=F, testDirection="over", annotation="hgu133a.db")
  Over.pres<-hyperGTest(params)
  ov<-summary(Over.pres)
ov$Term

typeof(sel.entrez)
unlist(sel.entrez)
sel.entrez
all.entrez



for(ci in community_index)
{
  print(ci)
}
cluster1 <- c("835", "5261","241", "994")
cluster2 <- c("307", "308", "317", "321", "506", "540", "378", "388", "396")
clusterSim(cluster1, cluster2, ont="MF", organism="human", measure="Wang")
temp <- geneSim("5921", "9046", ont = "BP", organism = "human", measure = "Wang", combine = "max")




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
  
 
  
  for(v in W)
  {
    type[v] = "W"
  }
  for(v in W1)
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
  ggnet(network,node.group=type,size=4,segment.size=1)
}

ifelse(V(g)%in%y, "Y","X")
type <- vector(mode="character", length=length(V(g)))
type <- setNames(type, V(g)$name)
for(v in V(g)$name)
{
  #print(v)
  if(as.integer(v)%in%y[[1]])
  {
    type[v] = 'y'
  }
}

print(type)

network = asNetwork(g)
ggnet(network,node.group=type,size=4,segment.size=1)
