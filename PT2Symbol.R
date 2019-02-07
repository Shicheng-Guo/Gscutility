
PT2Symbol<-function(ReactomePathWayName){
  library("ReactomePA")
  library("graphite")
  library("GOSemSim")
  library("igraph")
org2org <- list(arabidopsis = "athaliana", bovine = "btaurus", 
                canine = "cfamiliaris", chicken = "ggallus", ecolik12 = "ecoli", 
                fly = "dmelanogaster", human = "hsapiens", mouse = "mmusculus", 
                pig = "sscrofa", rat = "rnorvegicus", celegans = "celegans", 
                xenopus = "xlaevis", yeast = "scerevisiae", zebrafish = "drerio")
p <- pathways(org2org[["human"]], "reactome")[[ReactomePathWayName]]
p <- convertIdentifiers(p, "symbol")
g <- pathwayGraph(p)
gg <- igraph.from.graphNEL(g)
gg <- as.undirected(gg)
gg <- set.graph.attribute(gg,"name","RING")
GeneSymbol <- sub("[^:]+:", "", V(gg)$name)
return(GeneSymbol)
}

x1<-PT2Symbol("Interleukin-1 signaling")
x2<-PT2Symbol("Interleukin-2 signaling")
x3<-PT2Symbol("Interleukin-10 signaling")
x4<-PT2Symbol("Interleukin-12 signaling")
x5<-PT2Symbol("Interleukin-6 signaling")
x6<-PT2Symbol("Interleukin-17 signaling")
x7<-PT2Symbol("Regulation of beta-cell development")
x8<-PT2Symbol("Adaptive Immune System")
x9<-PT2Symbol("Immune System")
genelist<-c(x1,x2,x3,x4,x5,x6,x7,x8,x9)
