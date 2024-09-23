#Stack Overflow convert menu item
convertMenuItem <- function(mi,tabName) {
  mi$children[[1]]$attribs['data-toggle']="tab"
  mi$children[[1]]$attribs['data-value'] = tabName
  mi
}

#create GO organisms list
organisms_go <- function() {
  suppressMessages(rbioapi::rba_panther_info(what = "organisms"))
}

#draw counts of a gene
draw_counts <- function(dds_list, gene_name, nrows, ncols) {
  #extract counts of the selected gene
  counts_list <- lapply(dds_list, function(x) {
    plotCounts(x, gene=gene_name, intgroup="condition", returnData=TRUE)
  })
  #create a dataframe of counts
  counts_df <- dplyr::bind_rows(counts_list, .id = "study")
  #create the plot
  plt <- ggplot(counts_df, aes(x=condition, y=count))
  plt <- plt + geom_point(position=position_jitter(w=0.1,h=0))
  plt <- plt + labs(title = paste("Counts of", gene_name, "across studies"),
                    x = "Condition", y = "Normalized counts")
  plt <- plt + scale_y_log10()
  plt <- plt + facet_wrap(~study, nrow = nrows)
  plt
}