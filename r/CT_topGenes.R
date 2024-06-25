#' @title Display top genes for the cell type of interest given
#' @description
#' For the cell type of interest given, list out genes significant in the corresponding 
#' STANCE individual test with highest proportion of variance explained by this cell type.
#' \'cell_type_of_interest\' indicates the cell type of interest given. 
#' \'top_genes\' lists out all the top genes.
#' \'stacked_bar_plot\' is a stacked bar plot for top genes
#' visualizing the proportion of variance explained by all the cell types and error.
#' @param object STANCE object
#' @param CT_of_interest character string claiming cell type of interest to display.
#' @param numTopGenes (default 20) number of top genes to display.
#' If the number of significant genes is less than this given number,
#' then all significant genes will be displayed.
#' @return return STANCE object.
#' 
#' @import ggplot2
#' @import reshape2
#' @import dplyr
#' 
#' @export


library(ggplot2)
library(reshape2)
library(dplyr)

CT_topGenes <- function(object, CT_of_interest, numTopGenes = 20){
  if(is.null(object@VarComp_estimates)) {
    stop('No estimate of variance components found. Please do \'varcomp_est()\'before this step.')
  }
  n <- nrow(object@gene_expression)
  K <- ncol(object@cell_type_compositions)
  
  VC_prop <- apply(as.matrix(object@VarComp_estimates), MARGIN = 1, 
                   FUN = function(x){
                     return(x / sum(x))
                   })
  VC_prop <- as.data.frame(t(VC_prop))
  
  Test2_result <- object@Test_2[[CT_of_interest]]
  SigGenes.list <- row.names(Test2_result)[Test2_result$p_value < 0.05]
  numTopGenes <- min(length(SigGenes.list), 20)
  
  VC_prop.sorted <- VC_prop[order(VC_prop[,CT_of_interest], decreasing = T),]
  TopGenes.list <- row.names(VC_prop.sorted)[1:numTopGenes]
  
  VC_prop.topGene <- VC_prop.sorted[TopGenes.list,]
  VC_prop.topGene$Gene <- row.names(VC_prop.topGene)
  
  VC_prop.topGene.long <- melt(VC_prop.topGene, id.vars = 'Gene', variable.name = 'Component', value.name = 'Variance')
  
  
  # Use factor with the levels in the order you want to appear in the plot
  VC_prop.topGene.long$Gene <- factor(VC_prop.topGene.long$Gene, levels = TopGenes.list)
  
  # Set color palette
  cols <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set1"))
  myPal <- cols(length(unique(VC_prop.topGene.long$Component)))
  
  stacked.bar.plot <- ggplot(VC_prop.topGene.long, aes(x = Gene, y = Variance, fill = Component)) +
    geom_bar(stat = "identity", position = "stack") +  # Ensure bars are stacked
    theme_minimal() +
    scale_fill_brewer(palette = "Paired") + # Use a color palette, adjust as needed
    labs(x = "Genes", y = "Variance Explained", fill = "Components") +
    scale_fill_manual(values = myPal) +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          axis.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 8, face = "plain"), # Adjust text angle and hjust, vjust for horizontal layout
          axis.text.x = element_text(angle = -45, vjust = 0.5), # Ensure x axis labels (Variance Explained) are readable
          plot.margin = unit(c(1, 1, 1, 1), "cm"),
          axis.title.x = element_text(size = 12, face = "bold"),  # Adjust x-axis title size
          axis.title.y = element_text(size = 12, face = "bold")   # Adjust y-axis title size
          )
  
  output <- list(cell_type_of_interest = CT_of_interest, 
                 top_genes = TopGenes.list, 
                 stacked_bar_plot = stacked.bar.plot)
  
  object@cell_type_top_genes <- output
  return(object)
}
