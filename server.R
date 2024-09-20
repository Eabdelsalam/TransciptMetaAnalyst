# Define server logic required to draw a histogram
function(input, output, session) {
  # make reactive inputs and examples
  example_counts <- reactive({read.delim("data/example_counts.txt")})
  example_design <- reactive({read.delim("data/example_design.txt")})
  
  #download example counts file
  output$example_counts <- downloadHandler(
    filename = "example_counts.txt",
    content = function(filename){
      write.table(example_counts(), filename, row.names = F, quote = F, sep = "\t")
    },
    contentType = "text/csv"
  )
  
  #download examples experimental design file
  output$example_design <- downloadHandler(
    filename = "example_design.txt",
    content = function(filename){
      write.table(example_design(), filename, row.names = F, quote = F, sep = "\t")
    },
    contentType = "text/csv"
  )
  
  ####### RNA-Seq Analysis #######
  #check against files not uploaded
  #check if "Start Analysis" button was hit before uploading required files
  observeEvent(input$rnaseq_analyze, {
    if (is.null(input$rnaseq_counts)) {
      shinyalert("Missed Counts File!!!", "Please upload the raw counts file first!")
    }
  })
  
  observeEvent(input$rnaseq_analyze, {
    if (is.null(input$rnaseq_design)) {
      shinyalert("Missed Experimental Desgin File!!!", "Please upload the experimental design file first!")
    }
  })
  
  #process counts  file
  rnaseq_counts_input <- eventReactive(input$rnaseq_analyze, {
    uploaded <- input$rnaseq_counts
    if (is.null(uploaded))
      return(NULL)
    #read the data and check it without assigning
    temp <- read.delim(uploaded$datapath, header = T,
                       fill = T, quote = "")
    validate(test_counts(temp))
    return(temp)
  })
  
  #process experimental design file
  rnaseq_design_input <- eventReactive(input$rnaseq_analyze, {
    uploaded <- input$rnaseq_design
    if (is.null(uploaded))
      return(NULL)
    #read the data and check it without assigning
    temp <- read.delim(uploaded$datapath, header = T,
                       fill = T, quote = "")
    validate(test_exp_design(temp))
    return(temp)
  })
  
  # if files uploaded and validated hide info box and show results boxes
  observeEvent(input$rnaseq_analyze, {
    if (!is.null(rnaseq_design_input()) & !is.null(rnaseq_counts_input())) {
      #check for the presence of counts for all samples included in the analysis
      validate(test_match(rnaseq_counts_input(), rnaseq_design_input()))
      
      #provide info about starting the analysis
      shinyalert("Running...", "RNA-Seq meta-analysis is in progress... \nResults will appear soon!",
                 type = "info", closeOnEsc = T, closeOnClickOutside = T)
      #hide information box and show different results boxes
      shinyjs::hide("rnaseq_steps")
      shinyjs::show("general_results_rnaseq")
      shinyjs::show("mainplots")
      shinyjs::show("GO_counts")
      hideSpinner("go_plot_rnaseq")
      hideSpinner("kegg_plot_rnaseq")

      #RNA-Seq analysis progress
      #extract studies info 
      studies <- reactive({unique(rnaseq_design_input()$study)})
      num_studies <- reactive({length(unique(rnaseq_design_input()$study))}) #use for info about number of studies included
      genenames <- reactive({rnaseq_counts_input()$geneid})
      num_genes <- reactive({length(genenames())})
      samples <- reactive({rnaseq_design_input()$sample})
      num_samples <- reactive({length(rnaseq_design_input()$sample)})

      #create a list of counts
      studies_cts <- reactive({
        temp_list <- lapply(studies(),
                            function(x) rnaseq_counts_input()[, colnames(rnaseq_counts_input()) %in% rnaseq_design_input()$sample[rnaseq_design_input()$study == x]])
        temp_list <- lapply(temp_list, function(x) {rownames(x) <- rnaseq_counts_input()$geneid; x}) #add gene names as row names
        names(temp_list) <- studies() #naming the list
        return(temp_list)
      })
      #create a list of designs
      studies_design <- reactive({
        temp_list <- lapply(studies(),
                            function(x) rnaseq_design_input()[rnaseq_design_input()$study == x, ])
        temp_list <- lapply(temp_list, function(x) {rownames(x) <- x$sample; x}) #add sample names as row names
        names(temp_list) <- studies() #naming the list
        return(temp_list)
      })
      #create a list of DESEQ objects
      studies_deseq <- reactive({
        temp_list <- lapply(studies(), function(x) {suppressWarnings(
          DESeqDataSetFromMatrix(countData = studies_cts()[[x]],
                                 colData = studies_design()[[x]],
                                 design = ~ condition))
        })
        temp_list <- lapply(temp_list, function(x) {DESeq(x, quiet = T)})
        names(temp_list) <- studies() #naming the list
        return(temp_list)
      })
      #create list of results
      studies_results <- reactive({
        temp_list <- lapply(studies_deseq(),
                            function(x) {results(x, alpha = 0.05)})
        names(temp_list) <- studies() #naming the list
        return(temp_list)
      })
      #create a list of significant genes
      studies_significant <- reactive({
        temp_list <- lapply(studies_results(),
                            function(x) {subset(x, padj < 0.05 & abs(log2FoldChange) > 1)})
        names(temp_list) <- studies() #naming the list
        return(temp_list)
      })
      
      ############
      ## perform RNA-Seq meta-analysis
      #create a list of raw p-values
      rawpval_rnaseq <- reactive({
        temp_list <- lapply(studies_results(),
                            function(x) {x$pvalue})
        names(temp_list) <- studies()
        return(temp_list)
      })
      
      #create a list of adjusted p-values
      adjpval_rnaseq <- reactive({
        temp_list <- lapply(studies_results(),
                            function(x) {x$padj})
        names(temp_list) <- studies()
        return(temp_list)
      })
      
      #create a list of adjusted p-values
      fold_changes_rnaseq <- reactive({
        temp_list <- lapply(studies_results(),
                            function(x) {x$log2FoldChange})
        names(temp_list) <- studies()
        return(temp_list)
      })
      
      #create data frame of 0-1 scoring of all DEGs across all studies
      DE_rnaseq <- reactive({
        temp_df <- as.data.frame(mapply(
          adjpval_rnaseq(),
          FUN = function(x) ifelse(x <= 0.05, 1, 0)
        ))
        #rename columns and rows of the data frame
        colnames(temp_df) <- studies()
        rownames(temp_df) <- genenames()
        return(temp_df)
      })
      
      #filter DEGs with very low expression values
      #to avoid p-value enrichment near 1
      filtered_rnaseq <- reactive({
        temp_list = lapply(adjpval_rnaseq(), FUN = function(x) {which(is.na(x))})
        names(temp_list) <- studies()
        return(temp_list)
      })
      #apply the filtered values in the raw_pvalues list
      filtered_rawpval_rnaseq <- reactive({
        temp_list = rawpval_rnaseq()
        for (i in studies()) {temp_list[[i]][filtered_rnaseq()[[i]]] = NA}
        return(temp_list)
      })
      
      #Meta analysis using Fisher's method
      fishcomb <- reactive({
        fishercomb(filtered_rawpval_rnaseq(), BHth = 0.05)
      })
      
      #data frame of all DEGs
      DE_results_rnaseq <- reactive({data.frame(DE_rnaseq(),
                                                "DE.fishercomb" = ifelse(fishcomb()$adjpval <= 0.05, 1, 0))})
      
      #filter conflict DEGs
      #identified as up-regulated in one study and down-regulated in another
      signsFC <- reactive({
        temp <- mapply(
          fold_changes_rnaseq(),
          FUN = function(x)
            sign(x))
        return(temp)
        })
      
      sumsigns <- reactive({apply(signsFC(), 1, sum)})
      commonsgFC <- reactive({
        ifelse(abs(sumsigns()) == dim(signsFC())[2], sign(sumsigns()), 0)
      })
      
      #create a dataframe of DEGs with all fold change values
      FC.selectDE <- reactive({
        data.frame(DE_results_rnaseq(),
                   do.call(cbind, fold_changes_rnaseq()),
                   signFC = commonsgFC())
      })
      
      #list of DEGs to be kept or excluded (those with conflict FC)
      keepDE <- reactive({
        FC.selectDE()[which(abs(FC.selectDE()$signFC) == 1), ]
      })
      conflictDE <- reactive({
        FC.selectDE()[which(abs(FC.selectDE()$signFC) == 0), ]
      })
      
      #new identified or excluded DEGs
      fishcomb_de <- reactive({rownames(keepDE())[which(keepDE()[, "DE.fishercomb"] == 1)]})
      indstudy_de <- reactive({
        lapply(studies(),
               function(x) {rownames(keepDE())[which(keepDE()[, x] == 1)]})
      })
      diff <- reactive({IDD.IRR(fishcomb_de(), indstudy_de())})
      
      #create full_counts matrix for PCA
      full_cts_rnaseq <- reactive({
        temp_list <- lapply(studies_deseq(),
                            function(x) {DESeq2::counts(x, normalized = T)})
        temp_matrix <- matrix(unlist(temp_list), ncol = num_samples(), byrow = T)
        colnames(temp_matrix) <- samples()
        rownames(temp_matrix) <- genenames()
        return(temp_matrix)
      })
      
      #create DEGs counts matrix for heatmap
      cts_matrix <- reactive({
        temp_list <- lapply(studies_deseq(),
                            function(x) {DESeq2::counts(x, normalized = T)[fishcomb_de(),]})
        temp_list <- lapply(temp_list,
                            function(x) {t(apply(x, 1, scale))})
        cts_matrix <- matrix(unlist(temp_list), ncol = num_samples(), byrow = T)
        colnames(cts_matrix) <- rnaseq_design_input()$sample
        rownames(cts_matrix) <- fishcomb_de()
        return(cts_matrix)
      })

      #create a list for venn and updset plots
      for_venn <- reactive({
        temp_list <- indstudy_de()
        names(temp_list) <- studies()
        temp_list$Meta <- fishcomb_de()
        return(temp_list)
      })
      
      #create a flag to draw venn only if number of studies below 3
      if (num_studies() > 3) {
        removeTab(inputId = "mainplots", target = "Venn diagram")
      }
      
      #create a full table of final results
      meta_table_rnaseq <- reactive({
        templst <- fishcomb()[2:4]
        tempmatrix <- matrix(unlist(templst), ncol = 3)
        colnames(tempmatrix) <- names(templst)
        rownames(tempmatrix) <- genenames()
        tempdf <- as.data.frame(tempmatrix)
        tempdf$significant <- ifelse(tempdf$adjpval <= 0.05 & rownames(tempdf) %in% rownames(keepDE()), "yes", "no")
        return(tempdf)
      })
      
      ### Rendering outputs ###
      #update the number of the studies
      output$num_studies_rnaseq <- renderInfoBox({
        infoBox(title = "Number of studies", 
                paste(num_studies(), "studies"),
                width = 3,
                icon = icon("list-ol"))
      })
      
      #update the available results to download
      observe({updateSelectInput(session, "download_individual_rnaseq", choices = names(studies_significant()))})
      
      #download individual DEGs list
      output$download_individual_button_rnaseq <- downloadHandler(function() {paste0(input$download_individual_rnaseq,"_DEGs.csv")},
                                                                  content = function(filename) {
                                                                    write.csv(studies_significant()[[as.character(input$download_individual_rnaseq)]], filename, quote = F)
                                                                              }, contentType = "text/csv")
      #update the number of the DEGs
      output$num_degs_rnaseq <- renderInfoBox({degs_rnaseq <- keepDE()[!is.na(keepDE()$DE.fishercomb),]
        degs_rnaseq = degs_rnaseq[degs_rnaseq$DE.fishercomb == 1,]
        degs_rnaseq = dim(degs_rnaseq)[1]
        infobox <-infoBox(title = "Siginficant DEGs (meta-analysis)", 
                          paste(degs_rnaseq, "out of", num_genes(), "genes were identified as DEGs"),
                          paste0("The meta-analysis identified ", diff()[2], "(", diff()[4], "%) new gene(s), while ", diff()[3], "(", diff()[5], "%) are sidetracked."),
                          width = 4,
                          icon = icon("chart-simple"))
        return(infobox)
      })
      
      ####### MAIN PLOTS #######
      #distribution of p-values
      output$p_dist_rnaseq <- renderPlot({
        hist(fishcomb()$rawpval, main = "Fisher Method p-Values Distribution",
             breaks = 100, col = "grey", xlab = "Raw p-values (meta-analysis)")
      })
      #download p-values dist
      output$download_pdist_rnaseq <- downloadHandler("rnaseq_pdist.tiff",
                                                      content = function(file) {
                                                        tiff(file, width = 15, height = 12, units = "cm", res = 300, compression = "lzw")
                                                        hist(fishcomb()$rawpval, main = "Fisher Method p-Values Distribution",
                                                             breaks = 100, col = "grey", xlab = "Raw p-values (meta-analysis)")
                                                        dev.off()
                                                      },
                                                      contentType = "image/png")
      
      #PCA plot
      pca_plt <- reactive({
        full_design <- rnaseq_design_input()
        rownames(full_design) <- full_design$sample
        se <- SummarizedExperiment(log2(full_cts_rnaseq() + 1), colData = full_design)
        return(plotPCA(DESeqTransform(se)))
      })
      output$pca_rnaseq <- renderPlot({
        suppressMessages(pca_plt())
      })
      #download_pca
      output$download_pca_rnaseq <- downloadHandler("rnaseq_pca.tiff",
                                                    content = function(file) {
                                                      #draw to a file
                                                      ggsave(file, pca_plt(), device = "tiff",
                                                             width = 15, height = 12, units = "cm",
                                                             dpi = 300)
                                                    },
                                                    contentType = "image/png")
      
      #Venn diagram
      venn_plt <- reactive({ggvenn(for_venn())})
      output$venn_rnaseq <- renderPlot({
        venn_plt()
      })
      #download venn
      output$download_venn_rnaseq <- downloadHandler("rnaseq_venndiagram.tiff",
                                                     content = function(file) {
                                                       #draw to a file
                                                       ggsave(file, venn_plt(), device = "tiff",
                                                              width = 15, height = 15, units = "cm",
                                                              dpi = 300)
                                                     },
                                                     contentType = "image/png")
      
      #UpSet plot
      upset_plt <- reactive({
        m = make_comb_mat(for_venn())
        return(UpSet(m,
                     top_annotation = upset_top_annotation(m, add_numbers = TRUE),
                     right_annotation = upset_right_annotation(m, add_numbers = TRUE)))
        })
      output$upset_rnaseq <- renderPlot({
        upset_plt()
      })
      #download upset
      output$download_upset_rnaseq <- downloadHandler("rnaseq_upset.tiff",
                                                      content = function(file) {
                                                        #draw to a file
                                                        tiff(file, width = 15, height = 10, units = "cm", res = 300, compression = "lzw")
                                                        ComplexHeatmap::draw(upset_plt())
                                                        dev.off()
                                                      },
                                                      contentType = "image/png")
      
      #heatmap of all genes
      output$heatmap_rnaseq <- renderD3heatmap({
        suppressWarnings(d3heatmap(cts_matrix(),
                                   k_row = num_samples(),
                                   k_col = num_studies()))
      })
      #download button heatmap
      output$download_heatmap_rnaseq <- downloadHandler("rnaseq_heatmap.tiff",
                                                        content = function(filename) {
                                                          annot <- rnaseq_design_input()[, c("condition", "study")]
                                                          rownames(annot) <- rnaseq_design_input()$sample
                                                          pheatmap(cts_matrix(),
                                                                   cutree_rows = num_samples(),
                                                                   cutree_cols = num_studies(),
                                                                   show_rownames = F,
                                                                   color = colorRampPalette(c("yellow", "black", "purple"))(150),
                                                                   annotation_col = annot,
                                                                   filename = filename)
                                                        },
                                                        contentType = "image/png")
      #download heatmap data button
      output$download_heatmap_data_rnaseq <- downloadHandler("rnaseq_heatmap_data.csv",
                                                             content = function(filename) {
                                                               write.csv(cts_matrix(), filename, quote = F)
                                                             }, contentType = "text/csv")
      
      #volcano plot
      for_volcano <- reactive({
        fc_matrix <- matrix(unlist(fold_changes_rnaseq()), ncol = num_studies(), byrow = T)
        temp <- data.frame(
          genes = genenames(),
          fc = rowMeans(fc_matrix, na.rm = T),
          pval = fishcomb()$adjpval
        )
        return(temp)      
      })
      volcano_plt <- reactive({
        plt <- EnhancedVolcano(
          for_volcano(),
          lab = for_volcano()$genes,
          x = 'fc',
          y = 'pval',
          title = NULL,
          subtitle = NULL,
          pCutoff = 0.05,
          border = "full",
          caption = NULL
        )
        return(plt)
      })
      
      output$volcano_rnaseq <- renderPlot({
        suppressWarnings(volcano_plt())
      })
      
      #download volcano plot
      output$download_volcano_rnaseq <- downloadHandler("rnaseq_volcano_plot.tiff",
                                                        content = function(file) {
                                                          #draw to a file
                                                          ggsave(file, volcano_plt(), device = "tiff",
                                                                 width = 22, height = 15, units = "cm",
                                                                 dpi = 300)
                                                        },
                                                        contentType = "image/png")
      
      #download volcano plot data
      output$download_volcano_data_rnaseq <- downloadHandler("rnaseq_volcano_data.csv",
                                                             content = function(filename) {
                                                               write.csv(for_volcano(), filename, quote = F, row.names = F)
                                                             }, contentType = "text/csv")
      
      ## GO analysis
      #create a list of available organisms for analysis
      organisms <- reactive({organisms_go()})
      observe(
        updateSelectizeInput(session, inputId = "organism_go_rnaseq", choices = organisms()$long_name, server = T)
      )
      
      #run the analysis upon hitting the button
      observeEvent(input$run_go_rnaseq, {
        showSpinner("go_plot_rnaseq")
        go_rnaseq_table <- reactive({
          #biological process
          bp <- suppressMessages(rba_panther_enrich(genes = rownames(keepDE()),
                                                    organism = as.numeric(organisms()$taxon_id[organisms()$long_name == input$organism_go_rnaseq]),
                                                    annot_dataset = "GO:0008150", cutoff = 0.05))
          bp_result <- bp$result
          bp_result$category <- "BP"
          #cellular component
          cc <- suppressMessages(rba_panther_enrich(genes = rownames(keepDE()),
                                                    organism = as.numeric(organisms()$taxon_id[organisms()$long_name == input$organism_go_rnaseq]),
                                                    annot_dataset = "GO:0005575", cutoff = 0.05))
          cc_result <- cc$result
          cc_result$category <- "CC"
          #molecular function
          mf <- suppressMessages(rba_panther_enrich(genes = rownames(keepDE()),
                                                    organism = as.numeric(organisms()$taxon_id[organisms()$long_name == input$organism_go_rnaseq]),
                                                    annot_dataset = "GO:0003674", cutoff = 0.05))
          mf_result <- mf$result
          mf_result$category <- "MF"
          #combine in one table
          all_go <- rbind(bp_result, cc_result, mf_result)
          all_go <- all_go[!(is.na(all_go$term.id)), ] #remove unknown categories
          all_go <- all_go[all_go$fold_enrichment > 0, ] #remove categories with 0 genes
          all_go <- all_go[all_go$term.id != "GO:0008150", ] #remove main BP
          all_go <- all_go[all_go$term.id != "GO:0003674", ] #remove main MF
          all_go <- all_go[all_go$term.id != "GO:0005575", ] #remove main CC
          all_go <- all_go %>% arrange(fdr) #sort from the lowerst FDR to the highest FDR
          return(all_go)
        })
        
        #monitor the end of the anlaysis
        observeEvent(!is.null(go_rnaseq_table()), {
          hideSpinner("go_plot_rnaseq")
          #create downloadable table
          output$download_go_table_rnaseq <- downloadHandler("rnaseq_GO_results.txt",
                                                             content = function(filename) {
                                                               write_delim(go_rnaseq_table(), filename, quote = "none", delim = "\t")
                                                             }, contentType = "text/csv")
          #create lollipop plot
          go_plt_rnaseq <- reactive({
            #create data frame to use in constructing the plot
            for_plot <- data.frame(full_term = paste(go_rnaseq_table()$term.id, go_rnaseq_table()$term.label),
                                   fold_enrichment = go_rnaseq_table()$fold_enrichment,
                                   fdr = go_rnaseq_table()$fdr,
                                   genes = go_rnaseq_table()$number_in_list,
                                   plus_minus = go_rnaseq_table()$plus_minus,
                                   category = go_rnaseq_table()$category)
            if (dim(for_plot)[1] > 30) {
              for_plot <- for_plot[1:30, ]
            }
            #building the plot
            plt <- ggplot(data = for_plot, aes(x = reorder(full_term, fold_enrichment),
                                               y = fold_enrichment))
            plt <- plt + geom_segment(aes(x = reorder(full_term, fold_enrichment),
                                          xend = reorder(full_term, fold_enrichment),
                                          y = 0,
                                          yend = fold_enrichment,
                                          colour = fdr),
                                      linewidth = 1.5)
            plt <- plt + scale_color_gradient2(low="blue",
                                               high="red",
                                               mid = "green",
                                               midpoint = mean(for_plot$fdr),
                                               name = bquote(-log[10](FDR)))
            plt <- plt + geom_point(aes(colour = fdr, size = genes))
            plt <- plt + labs(x = "GO terms", y = "Fold enrichment", size = "No. of genes")
            plt <- plt + facet_grid(rows = vars(category),
                                    scales = "free",
                                    space = "free")
            plt <- plt + coord_flip()
            plt <- plt + theme(axis.text.x = element_text(color = "grey20", size = 12),
                               axis.text.y = element_text(color = "grey20", size = 12))
            return(plt)
          })
          #render the plot
          output$go_plot_rnaseq <- renderPlot(
            go_plt_rnaseq()
          )
          #create download plot
          output$download_go_plot_rnaseq <- downloadHandler("rnaseq_go_plot.tiff",
                                                            content = function(file) {
                                                              #draw to a file
                                                              ggsave(file, go_plt_rnaseq(), device = "tiff",
                                                                     width = 30, height = 20, units = "cm",
                                                                     dpi = 300)
                                                            },
                                                            contentType = "image/png")
        })
      })
      
      ## KEGG analysis
      #start the analysis upon hitting the button
      observeEvent(input$run_kegg_rnaseq, {
        showSpinner("kegg_plot_rnaseq")
        kegg_rnaseq_table <- reactive({
          temp <- suppressMessages(enrichKEGG(gene = rownames(keepDE()),
                                              organism = input$organism_kegg_rnaseq,
                                              qvalueCutoff = 0.05,
                                              pAdjustMethod = "BH"))
          result <- as.data.frame(temp)
          return(result)
        })
        
        #check if the analysis completed
        #allow download after analysis
        observeEvent(!(is.null(kegg_rnaseq_table())), {
          hideSpinner("kegg_plot_rnaseq")
          output$download_kegg_table_rnaseq <- downloadHandler("rnaseq_KEGG_results.csv",
                                                               content = function(filename) {
                                                                 write.csv(kegg_rnaseq_table(), filename, quote = F)
                                                               }, contentType = "text/csv")
          #building the plot
          for_plot <- kegg_rnaseq_table()
          if (dim(for_plot)[1] > 30) {
            for_plot <- for_plot[1:30, ]
          }
          for_plot$log_fdr = log10(for_plot$qvalue) * -1
          plt <- ggplot(data = for_plot, aes(x = reorder(Description, FoldEnrichment),
                                             y = FoldEnrichment))
          plt <- plt + geom_segment(aes(x = reorder(Description, FoldEnrichment),
                                        xend = reorder(Description, FoldEnrichment),
                                        y = 0,
                                        yend = FoldEnrichment,
                                        colour = log_fdr),
                                    linewidth = 1.5)
          plt <- plt + scale_color_gradient2(low="blue",
                                             high="red",
                                             mid = "green",
                                             midpoint = mean(for_plot$log_fdr),
                                             name = bquote(-log[10](FDR)))
          plt <- plt + geom_point(aes(colour = log_fdr, size = Count))
          plt <- plt + labs(x = "KEGG pathways", y = "Fold enrichment", size = "No. of genes")
          plt <- plt + facet_grid(rows = vars(subcategory),
                                  scales = "free",
                                  space = "free")
          plt <- plt + coord_flip()
          plt <- plt + theme(axis.text.x = element_text(color = "grey20", size = 12),
                             axis.text.y = element_text(color = "grey20", size = 12))
          #render the plot
          output$kegg_plot_rnaseq <- renderPlot(plt)
          #create download plot
          output$download_kegg_plot_rnaseq <- downloadHandler("rnaseq_kegg_plot.tiff",
                                                              content = function(file) {
                                                                #draw to a file
                                                                ggsave(file, plt, device = "tiff",
                                                                       width = 30, height = 20, units = "cm",
                                                                       dpi = 300)
                                                              },
                                                              contentType = "image/png")
        })
      })
      
      ####### META DEGS #######
      #render data table
      output$rnaseq_DEGs <- renderDataTable(meta_table_rnaseq(), selection = "single")
      #add button to download the data as csv
      output$download_meta_degs_rnaseq <- downloadHandler("rnaseq_meta_analysis_results.csv",
                                                          content = function(filename) {
                                                            write.csv(meta_table_rnaseq(), filename, quote = F)
                                                          }, contentType = "text/csv")
      #write an explanation of the box
      output$counts_message_rnaseq <- renderText("Click on any gene in the above DEGs table to obtain a plot of the gene's normalized counts across the included studies")
      #plot counts of the selected protein from DEGs table
      observeEvent(input$rnaseq_DEGs_rows_selected, {
        #get the selected gene name
        gene_df <- meta_table_rnaseq()[input$rnaseq_DEGs_rows_selected,]
        gene_name <- rownames(gene_df)
        nrows = floor(sqrt(num_studies()))    #define number of rows based on studies
        
        if (!(is.null(gene_df))) { #check if there is a row ro draw counts
          output$counts_message_rnaseq <- renderText("")
          cts_plt_rnaseq <- reactive({draw_counts(studies_deseq(), gene_name, nrows)})
        }
        #render the drawn plot
        output$counts_plots_rnaseq <- renderPlot(cts_plt_rnaseq())
        #allow plot download
        output$download_counts_plot_rnaseq <- downloadHandler(function() {paste(gene_name, "rnaseq_counts_plot.tiff", sep = "_")},
                                                            content = function(file) {
                                                              #draw to a file
                                                              ggsave(file, cts_plt_rnaseq(), device = "tiff",
                                                                     width = 30, height = 20, units = "cm",
                                                                     dpi = 300)
                                                            },
                                                            contentType = "image/png")
      })
      
    } #close rnaseq analysis logic after hitting start analysis
  }) #close RNA-Seq observe statement
  
  
  
  
  ####### Running DEMO #######
  observe(
    if (input$main_tabs == "demo") {
      #provide info about starting the demo analysis
      shinyalert("Demo Running!", "Demo RNA-Seq meta-analysis is in progress... \nResults will appear soon!",
                 type = "info", closeOnEsc = T, closeOnClickOutside = T)
      hideSpinner("go_plot_demo")
      hideSpinner("kegg_plot_demo")
      #Demo analysis progress
      #extract studies info 
      studies <- reactive({unique(example_design()$study)})
      num_studies <- reactive({length(unique(example_design()$study))}) #use for info about number of studies included
      genenames <- reactive({example_counts()$geneid})
      num_genes <- reactive({length(genenames())})
      samples <- reactive({example_design()$sample})
      num_samples <- reactive({length(example_design()$sample)})
      
      #create a list of counts
      studies_cts <- reactive({
        temp_list <- lapply(studies(),
                            function(x) example_counts()[, colnames(example_counts()) %in% example_design()$sample[example_design()$study == x]])
        temp_list <- lapply(temp_list, function(x) {rownames(x) <- example_counts()$geneid; x}) #add gene names as row names
        names(temp_list) <- studies() #naming the list
        return(temp_list)
      })
      #create a list of designs
      studies_design <- reactive({
        temp_list <- lapply(studies(),
                            function(x) example_design()[example_design()$study == x, ])
        temp_list <- lapply(temp_list, function(x) {rownames(x) <- x$sample; x}) #add sample names as row names
        names(temp_list) <- studies() #naming the list
        return(temp_list)
      })
      #create a list of DESEQ objects
      studies_deseq <- reactive({
        temp_list <- lapply(studies(), function(x) {suppressWarnings(
          DESeqDataSetFromMatrix(countData = studies_cts()[[x]],
                                 colData = studies_design()[[x]],
                                 design = ~ condition))
        })
        temp_list <- lapply(temp_list, function(x) {DESeq(x, quiet = T)})
        names(temp_list) <- studies() #naming the list
        return(temp_list)
      })
      #create list of results
      studies_results <- reactive({
        temp_list <- lapply(studies_deseq(),
                            function(x) {results(x, alpha = 0.05)})
        names(temp_list) <- studies() #naming the list
        return(temp_list)
      })
      #create a list of significant genes
      studies_significant <- reactive({
        temp_list <- lapply(studies_results(),
                            function(x) {subset(x, padj < 0.05 & abs(log2FoldChange) > 1)})
        names(temp_list) <- studies() #naming the list
        return(temp_list)
      })
      
      ############
      ## perform RNA-Seq meta-analysis
      #create a list of raw p-values
      rawpval_demo <- reactive({
        temp_list <- lapply(studies_results(),
                            function(x) {x$pvalue})
        names(temp_list) <- studies()
        return(temp_list)
      })
      
      #create a list of adjusted p-values
      adjpval_demo <- reactive({
        temp_list <- lapply(studies_results(),
                            function(x) {x$padj})
        names(temp_list) <- studies()
        return(temp_list)
      })
      
      #create a list of adjusted p-values
      fold_changes_demo <- reactive({
        temp_list <- lapply(studies_results(),
                            function(x) {x$log2FoldChange})
        names(temp_list) <- studies()
        return(temp_list)
      })
      
      #create data frame of 0-1 scoring of all DEGs across all studies
      DE_demo <- reactive({
        temp_df <- as.data.frame(mapply(
          adjpval_demo(),
          FUN = function(x) ifelse(x <= 0.05, 1, 0)
        ))
        #rename columns and rows of the data frame
        colnames(temp_df) <- studies()
        rownames(temp_df) <- genenames()
        return(temp_df)
      })
      
      #filter DEGs with very low expression values
      #to avoid p-value enrichment near 1
      filtered_demo <- reactive({
        temp_list = lapply(adjpval_demo(), FUN = function(x) {which(is.na(x))})
        names(temp_list) <- studies()
        return(temp_list)
      })
      #apply the filtered values in the raw_pvalues list
      filtered_rawpval_demo <- reactive({
        temp_list = rawpval_demo()
        for (i in studies()) {temp_list[[i]][filtered_demo()[[i]]] = NA}
        return(temp_list)
      })
      
      #Meta analysis using Fisher's method
      fishcomb <- reactive({
        fishercomb(filtered_rawpval_demo(), BHth = 0.05)
      })
      
      #data frame of all DEGs
      DE_results_demo <- reactive({data.frame(DE_demo(),
                                              "DE.fishercomb" = ifelse(fishcomb()$adjpval <= 0.05, 1, 0))})
      
      #filter conflict DEGs
      #identified as up-regulated in one study and down-regulated in another
      signsFC <- reactive({
        temp <- mapply(
          fold_changes_demo(),
          FUN = function(x)
            sign(x))
        return(temp)
      })
      
      sumsigns <- reactive({apply(signsFC(), 1, sum)})
      commonsgFC <- reactive({
        ifelse(abs(sumsigns()) == dim(signsFC())[2], sign(sumsigns()), 0)
      })
      
      #create a dataframe of DEGs with all fold change values
      FC.selectDE <- reactive({
        data.frame(DE_results_demo(),
                   do.call(cbind, fold_changes_demo()),
                   signFC = commonsgFC())
      })
      
      #list of DEGs to be kept or excluded (those with conflict FC)
      keepDE <- reactive({
        FC.selectDE()[which(abs(FC.selectDE()$signFC) == 1), ]
      })
      conflictDE <- reactive({
        FC.selectDE()[which(abs(FC.selectDE()$signFC) == 0), ]
      })
      
      #new identified or excluded DEGs
      fishcomb_de <- reactive({rownames(keepDE())[which(keepDE()[, "DE.fishercomb"] == 1)]})
      indstudy_de <- reactive({
        lapply(studies(),
               function(x) {rownames(keepDE())[which(keepDE()[, x] == 1)]})
      })
      diff <- reactive({IDD.IRR(fishcomb_de(), indstudy_de())})
      
      #create full_counts matrix for PCA
      full_cts_demo <- reactive({
        temp_list <- lapply(studies_deseq(),
                            function(x) {DESeq2::counts(x, normalized = T)})
        temp_matrix <- matrix(unlist(temp_list), ncol = num_samples(), byrow = T)
        colnames(temp_matrix) <- samples()
        rownames(temp_matrix) <- genenames()
        return(temp_matrix)
      })
      
      #create DEGs counts matrix for heatmap
      cts_matrix <- reactive({
        temp_list <- lapply(studies_deseq(),
                            function(x) {DESeq2::counts(x, normalized = T)[fishcomb_de(),]})
        temp_list <- lapply(temp_list,
                            function(x) {t(apply(x, 1, scale))})
        cts_matrix <- matrix(unlist(temp_list), ncol = num_samples(), byrow = T)
        colnames(cts_matrix) <- example_design()$sample
        rownames(cts_matrix) <- fishcomb_de()
        return(cts_matrix)
      })
      
      #create a list for venn and updset plots
      for_venn <- reactive({
        temp_list <- indstudy_de()
        names(temp_list) <- studies()
        temp_list$Meta <- fishcomb_de()
        return(temp_list)
      })
      
      #create a full table of final results
      meta_table_demo <- reactive({
        templst <- fishcomb()[2:4]
        tempmatrix <- matrix(unlist(templst), ncol = 3)
        colnames(tempmatrix) <- names(templst)
        rownames(tempmatrix) <- genenames()
        tempdf <- as.data.frame(tempmatrix)
        tempdf$significant <- ifelse(tempdf$adjpval <= 0.05 & rownames(tempdf) %in% rownames(keepDE()), "yes", "no")
        return(tempdf)
      })
      
      ### Rendering outputs ###
      #update the number of the studies
      output$num_studies_demo <- renderInfoBox({
        infoBox(title = "Number of studies", 
                paste(num_studies(), "studies"),
                width = 3,
                icon = icon("list-ol"))
      })
      
      #update the available results to download
      observe({updateSelectInput(session, "download_individual_demo", choices = names(studies_significant()))})
      
      #download individual DEGs list
      output$download_individual_button_demo <- downloadHandler(function() {paste0(input$download_individual_demo,"_DEGs.csv")},
                                                                content = function(filename) {
                                                                  write.csv(studies_significant()[[as.character(input$download_individual_demo)]], filename, quote = F)
                                                                }, contentType = "text/csv")
      #update the number of the DEGs
      output$num_degs_demo <- renderInfoBox({degs_demo <- keepDE()[!is.na(keepDE()$DE.fishercomb),]
      degs_demo = degs_demo[degs_demo$DE.fishercomb == 1,]
      degs_demo = dim(degs_demo)[1]
      infobox <-infoBox(title = "Siginficant DEGs (meta-analysis)", 
                        paste(degs_demo, "out of", num_genes(), "genes were identified as DEGs"),
                        paste0("The meta-analysis identified ", diff()[2], "(", diff()[4], "%) new gene(s), while ", diff()[3], "(", diff()[5], "%) are sidetracked."),
                        width = 4,
                        icon = icon("chart-simple"))
      return(infobox)
      })
      
      ####### MAIN PLOTS #######
      #distribution of p-values
      output$p_dist_demo <- renderPlot({
        hist(fishcomb()$rawpval, main = "Fisher Method p-Values Distribution",
             breaks = 100, col = "grey", xlab = "Raw p-values (meta-analysis)")
      })
      #download p-values dist
      output$download_pdist_demo <- downloadHandler("demo_pdist.tiff",
                                                    content = function(file) {
                                                      tiff(file, width = 15, height = 12, units = "cm", res = 300, compression = "lzw")
                                                      hist(fishcomb()$rawpval, main = "Fisher Method p-Values Distribution",
                                                           breaks = 100, col = "grey", xlab = "Raw p-values (meta-analysis)")
                                                      dev.off()
                                                    },
                                                    contentType = "image/png")
      
      #PCA plot
      pca_plt <- reactive({
        full_design <- example_design()
        rownames(full_design) <- full_design$sample
        se <- SummarizedExperiment(log2(full_cts_demo() + 1), colData = full_design)
        return(plotPCA(DESeqTransform(se)))
      })
      output$pca_demo <- renderPlot({
        suppressMessages(pca_plt())
      })
      #download_pca
      output$download_pca_demo <- downloadHandler("demo_pca.tiff",
                                                  content = function(file) {
                                                    #draw to a file
                                                    ggsave(file, pca_plt(), device = "tiff",
                                                           width = 15, height = 12, units = "cm",
                                                           dpi = 300)
                                                  },
                                                  contentType = "image/png")
      
      #Venn diagram
      venn_plt <- reactive({ggvenn(for_venn())})
      output$venn_demo <- renderPlot({
        venn_plt()
      })
      #download venn
      output$download_venn_demo <- downloadHandler("demo_venndiagram.tiff",
                                                   content = function(file) {
                                                     #draw to a file
                                                     ggsave(file, venn_plt(), device = "tiff",
                                                            width = 15, height = 15, units = "cm",
                                                            dpi = 300)
                                                   },
                                                   contentType = "image/png")
      
      #UpSet plot
      upset_plt <- reactive({
        m = make_comb_mat(for_venn())
        return(UpSet(m,
                     top_annotation = upset_top_annotation(m, add_numbers = TRUE),
                     right_annotation = upset_right_annotation(m, add_numbers = TRUE)))
      })
      output$upset_demo <- renderPlot(upset_plt())
      #download upset
      output$download_upset_demo <- downloadHandler("demo_upset.tiff",
                                                    content = function(file) {
                                                      #draw to a file
                                                      tiff(file, width = 15, height = 10, units = "cm", res = 300, compression = "lzw")
                                                      ComplexHeatmap::draw(upset_plt())
                                                      dev.off()
                                                    },
                                                    contentType = "image/png")
      
      #heatmap of all genes
      output$heatmap_demo <- renderD3heatmap({
        suppressWarnings(d3heatmap(cts_matrix(),
                                   k_row = num_samples(),
                                   k_col = num_studies()))
      })
      #download button heatmap
      output$download_heatmap_demo <- downloadHandler("demo_heatmap.tiff",
                                                      content = function(filename) {
                                                        annot <- example_design()[, c("condition", "study")]
                                                        rownames(annot) <- example_design()$sample
                                                        pheatmap(cts_matrix(),
                                                                 cutree_rows = num_samples(),
                                                                 cutree_cols = num_studies(),
                                                                 show_rownames = F,
                                                                 color = colorRampPalette(c("yellow", "black", "purple"))(150),
                                                                 annotation_col = annot,
                                                                 filename = filename)
                                                      },
                                                      contentType = "image/png")
      #download heatmap data button
      output$download_heatmap_data_demo <- downloadHandler("demo_heatmap_data.csv",
                                                           content = function(filename) {
                                                             write.csv(cts_matrix(), filename, quote = F)
                                                           }, contentType = "text/csv")
      
      #volcano plot
      for_volcano <- reactive({
        fc_matrix <- matrix(unlist(fold_changes_demo()), ncol = num_studies(), byrow = T)
        temp <- data.frame(
          genes = genenames(),
          fc = rowMeans(fc_matrix, na.rm = T),
          pval = fishcomb()$adjpval
        )
        return(temp)      
      })
      volcano_plt <- reactive({
        plt <- EnhancedVolcano(
          for_volcano(),
          lab = for_volcano()$genes,
          x = 'fc',
          y = 'pval',
          title = NULL,
          subtitle = NULL,
          pCutoff = 0.05,
          border = "full",
          caption = NULL
        )
        return(plt)
      })
      
      output$volcano_demo <- renderPlot({
        suppressWarnings(volcano_plt())
      })
      
      #download volcano plot
      output$download_volcano_demo <- downloadHandler("demo_volcanplot.tiff",
                                                      content = function(file) {
                                                        #draw to a file
                                                        ggsave(file, volcano_plt(), device = "tiff",
                                                               width = 22, height = 15, units = "cm",
                                                               dpi = 300)
                                                      },
                                                      contentType = "image/png")
      
      #download volcano plot data
      output$download_volcano_data_demo <- downloadHandler("demo_volcano_data.csv",
                                                           content = function(filename) {
                                                             write.csv(for_volcano(), filename, quote = F, row.names = F)
                                                           }, contentType = "text/csv")
      
      ## GO analysis
      #run the analysis upon hitting the button
      observeEvent(input$run_go_demo, {
        showSpinner("go_plot_demo")
        go_table_demo <- reactive({
          #biological process
          bp <- suppressMessages(rba_panther_enrich(genes = rownames(keepDE()),
                                                    organism = 1111708,
                                                    annot_dataset = "GO:0008150", cutoff = 0.05))
          bp_result <- bp$result
          bp_result$category <- "BP"
          #cellular component
          cc <- suppressMessages(rba_panther_enrich(genes = rownames(keepDE()),
                                                    organism = 1111708,
                                                    annot_dataset = "GO:0005575", cutoff = 0.05))
          cc_result <- cc$result
          cc_result$category <- "CC"
          #molecular function
          mf <- suppressMessages(rba_panther_enrich(genes = rownames(keepDE()),
                                                    organism = 1111708,
                                                    annot_dataset = "GO:0003674", cutoff = 0.05))
          mf_result <- mf$result
          mf_result$category <- "MF"
          #combine in one table
          all_go <- rbind(bp_result, cc_result, mf_result)
          all_go <- all_go[!(is.na(all_go$term.id)), ] #remove unknown categories
          all_go <- all_go[all_go$fold_enrichment > 0, ] #remove categories with 0 genes
          all_go <- all_go[all_go$term.id != "GO:0008150", ] #remove main BP
          all_go <- all_go[all_go$term.id != "GO:0003674", ] #remove main MF
          all_go <- all_go[all_go$term.id != "GO:0005575", ] #remove main CC
          all_go <- all_go %>% arrange(fdr) #sort from the lowerst FDR to the highest FDR
          return(all_go)
        })

        #monitor the end of the anlaysis
        observeEvent(!is.null(go_table_demo()), {
          hideSpinner("go_plot_demo")
          #create downloadable table
          output$download_go_table_demo <- downloadHandler("demo_GO_results.txt",
                                                           content = function(filename) {
                                                             write_delim(go_table_demo(), filename, quote = "none", delim = "\t")
                                                           }, contentType = "text/csv")
          #create lollipop plot
          go_plt_dm <- reactive({
            #create data frame to use in constructing the plot
            for_plot <- data.frame(full_term = paste(go_table_demo()$term.id, go_table_demo()$term.label),
                                   fold_enrichment = go_table_demo()$fold_enrichment,
                                   fdr = go_table_demo()$fdr,
                                   genes = go_table_demo()$number_in_list,
                                   plus_minus = go_table_demo()$plus_minus,
                                   category = go_table_demo()$category)
            if (dim(for_plot)[1] > 30) {
              for_plot <- for_plot[1:30, ]
            }
            #building the plot
            plt <- ggplot(data = for_plot, aes(x = reorder(full_term, fold_enrichment),
                                               y = fold_enrichment))
            plt <- plt + geom_segment(aes(x = reorder(full_term, fold_enrichment),
                                          xend = reorder(full_term, fold_enrichment),
                                          y = 0,
                                          yend = fold_enrichment,
                                          colour = fdr),
                                      linewidth = 1.5)
            plt <- plt + scale_color_gradient2(low="blue",
                                               high="red",
                                               mid = "green",
                                               midpoint = mean(for_plot$fdr),
                                               name = bquote(-log[10](FDR)))
            plt <- plt + geom_point(aes(colour = fdr, size = genes))
            plt <- plt + labs(x = "GO terms", y = "Fold enrichment", size = "No. of genes")
            plt <- plt + facet_grid(rows = vars(category),
                                    scales = "free",
                                    space = "free")
            plt <- plt + coord_flip()
            plt <- plt + theme(axis.text.x = element_text(color = "grey20", size = 12),
                               axis.text.y = element_text(color = "grey20", size = 12))
            return(plt)
          })
          #render the plot
          output$go_plot_demo <- renderPlot(go_plt_dm()) 
          #create download plot
          output$download_go_plot_demo <- downloadHandler("demo_go_plot.tiff",
                                                          content = function(file) {
                                                            #draw to a file
                                                            ggsave(file, go_plt_dm(), device = "tiff",
                                                                   width = 30, height = 20, units = "cm",
                                                                   dpi = 300)
                                                          },
                                                          contentType = "image/png")
        })
      })
      
      ## KEGG analysis
      #start the analysis upon hitting the button
      observeEvent(input$run_kegg_demo, {
        showSpinner("kegg_plot_demo")
        kegg_table_demo <- reactive({
          temp <- suppressMessages(enrichKEGG(gene = rownames(keepDE()),
                                              organism = "syn",
                                              qvalueCutoff = 0.05,
                                              pAdjustMethod = "BH"))
          result <- as.data.frame(temp)
          return(result)
        })
        
        #check if the analysis completed
        #allow download after analysis
        observeEvent(!(is.null(kegg_table_demo())), {
          hideSpinner("kegg_plot_demo")
          output$download_kegg_table_demo <- downloadHandler("demo_KEGG_results.csv",
                                                             content = function(filename) {
                                                               write.csv(kegg_table_demo(), filename, quote = F)
                                                             }, contentType = "text/csv")
          #building the plot
          for_plot <- kegg_table_demo()
          for_plot$log_fdr = log10(for_plot$qvalue) * -1
          plt <- ggplot(data = for_plot, aes(x = reorder(Description, FoldEnrichment),
                                             y = FoldEnrichment))
          plt <- plt + geom_segment(aes(x = reorder(Description, FoldEnrichment),
                                        xend = reorder(Description, FoldEnrichment),
                                        y = 0,
                                        yend = FoldEnrichment,
                                        colour = log_fdr),
                                    linewidth = 1.5)
          plt <- plt + scale_color_gradient2(low="blue",
                                             high="red",
                                             mid = "green",
                                             midpoint = mean(for_plot$log_fdr),
                                             name = bquote(-log[10](FDR)))
          plt <- plt + geom_point(aes(colour = log_fdr, size = Count))
          plt <- plt + labs(x = "KEGG pathways", y = "Fold enrichment", size = "No. of genes")
          plt <- plt + facet_grid(rows = vars(subcategory),
                                  scales = "free",
                                  space = "free")
          plt <- plt + coord_flip()
          plt <- plt + theme(axis.text.x = element_text(color = "grey20", size = 12),
                             axis.text.y = element_text(color = "grey20", size = 12))
          #render the plot
          output$kegg_plot_demo <- renderPlot(plt)
          #create download plot
          output$download_kegg_plot_demo <- downloadHandler("demo_kegg_plot.tiff",
                                                            content = function(file) {
                                                              #draw to a file
                                                              ggsave(file, plt, device = "tiff",
                                                                     width = 30, height = 20, units = "cm",
                                                                     dpi = 300)
                                                            },
                                                            contentType = "image/png")
        })
      })

      ####### META DEGS #######
      #render data table
      output$demo_DEGs <- renderDataTable(meta_table_demo(), selection = "single")
      #add button to download the data as csv
      output$download_meta_degs_demo <- downloadHandler("demo_meta_analysis_results.csv",
                                                        content = function(filename) {
                                                          write.csv(meta_table_demo(), filename, quote = F)
                                                        }, contentType = "text/csv")
      #write an explanation of the box
      output$counts_message_demo <- renderText("Click on any gene in the above DEGs table to obtain a plot of the gene's normalized counts across the included studies")
      #plot counts of the selected protein from DEGs table
      observeEvent(input$demo_DEGs_rows_selected, {
        #get the selected gene name
        gene_df <- meta_table_demo()[input$demo_DEGs_rows_selected,]
        gene_name <- rownames(gene_df)
        nrows = floor(sqrt(num_studies()))    #define number of rows based on studies
        
        if (!(is.null(input$demo_DEGs_rows_selected))) { #check if there is a row ro draw counts
          output$counts_message_demo <- renderText("")
          cts_plt_dm <- reactive({draw_counts(studies_deseq(), gene_name, nrows)})
        }
        #render the drawn plot
        output$counts_plots_demo <- renderPlot(cts_plt_dm())
        #allow plot download
        output$download_counts_plot_demo <- downloadHandler(function() {paste(gene_name, "demo_counts_plot.tiff", sep = "_")},
                                                            content = function(file) {
                                                              #draw to a file
                                                              ggsave(file, cts_plt_dm(), device = "tiff",
                                                                     width = 30, height = 20, units = "cm",
                                                                     dpi = 300)
                                                            },
                                                            contentType = "image/png")
      })
    }
  )
}
