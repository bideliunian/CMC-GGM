network.plotting <- function(edge.matrix, cols){
  ## get the intersection of two adjacant matrices/ two edge adjmat
  p <- nrow(edge.matrix)
  ################### PART 3: Draw network #####################
  layMat <- as.matrix(expand.grid(seq(0, 1, len=sqrt(p)), seq(0, 1, len=sqrt(p))))
  
  rownames(edge.matrix) <- as.character(c(1:p))
  colnames(edge.matrix) <- as.character(c(1:p)) 
  edge.matrix[lower.tri(edge.matrix, diag = TRUE)] <- 0
  
  df <- as.data.frame(edge.matrix)
  df$rowId <- rownames(df)
  df.long <- melt(df, id='rowId', variable.name = 'colId')
  names(df.long)[3] <- 'weight'
  df.long$'edgeId' <- rownames(df.long)
  df.long <- filter(df.long, weight > 0)
  
  df.longlong <- melt(df.long, id.var = c('weight', 'edgeId'), variable.name = 'nodeId')
  df.longlong$value <- as.numeric(df.longlong$value)
  df.longlong$x <- layMat[df.longlong$value, 2]
  df.longlong$y <- layMat[df.longlong$value, 1]
  
  df.longlong <- df.longlong[order(df.longlong$weight),]
  df.longlong$weight.rank <- as.factor(rep(1:5, each = ceiling(nrow(df.longlong)/5))[1:nrow(df.longlong)])
  
  p1 = ggplot(df.longlong, aes(x, y)) +
    geom_line(aes(group=weight, color = weight, size=weight.rank)) +
    scale_colour_gradientn(colours = cols) +
    guides(size = "none") +
    theme_bw() + 
    scale_size_manual(values = c(0.5, 1, 1.5, 2, 2.5)) +
    labs(x ="", y = "") + 
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.key.height = unit(1.75, "cm"),
          legend.title=element_blank(),
          legend.text=element_text(size=14))
  
  return(p1)
}
