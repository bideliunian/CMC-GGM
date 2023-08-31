## ploting the connected graph from control group and alcohol group together
network.plotting <- function(AdjMat1, AdjMat2, node.names){
  ## get the intersection of two adjacant matrices/ two edge adjmat
  p <- nrow(AdjMat1)
  AdjMat.and <- matrix(as.numeric(AdjMat1 & AdjMat2), nrow=p, ncol=p)
  AdjMat1.edge <- AdjMat1 - AdjMat.and
  AdjMat2.edge <- AdjMat2 - AdjMat.and
  
  row.names(AdjMat.and) <- node.names
  colnames(AdjMat.and) <- node.names
  row.names(AdjMat1.edge) <- node.names
  colnames(AdjMat1.edge) <- node.names
  row.names(AdjMat2.edge) <- node.names
  colnames(AdjMat2.edge) <- node.names
  
  net.and <- graph_from_adjacency_matrix(AdjMat.and, mode="undirected")
  net1 <- graph_from_adjacency_matrix(AdjMat1.edge, mode="undirected")
  net2 <- graph_from_adjacency_matrix(AdjMat2.edge, mode="undirected")

  #V(net.and)$label.cex <- 1
  plot(net.and, edge.color="#619CFF", edge.width=2, vertex.size=3, margin=0, layout=layMat,  vertex.color='grey90', 
       edge.curved=0, vertex.label.dist=1, vertex.label.cex=1)
  plot(net1, edge.color="#F8766D", edge.width=2, vertex.size=3, margin=0, layout=layMat, vertex.color='grey90',
       edge.curved=0.3, vertex.label.dist=1, vertex.label.cex=1, add=T)
  plot(net2, edge.color="#00BA38", edge.width=2, vertex.size=3, margin=0, layout=layMat, vertex.color='grey90',
       edge.curved=0.3, vertex.label.dist=1, vertex.label.cex=1, add=T)
}