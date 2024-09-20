#package version
VERSION <- "v1.0.0"

#source the supplementary script files
source("R/helper_functions.R")
source("R/test_files.R")

#check if required packages is installed
#list of packages required
required_packages <- c("shiny", "shinydashboard", "shinyjs", "shinyalert", "shinycssloaders", "shiny.info", "DESeq2", "metaRNASeq", "tidyverse", "ggplot2", "ggrepel", "ggvenn", "ComplexHeatmap", "pheatmap", "d3heatmap", "EnhancedVolcano", "DT", "rbioapi", "clusterProfiler", "heatmaply", "gridExtra")

#checking missing packages from list
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

#install missing ones from CRAN
if(length(new_packages)) install.packages(new_packages, dependencies = TRUE)
#install missing ones from BioConductor
if(length(new_packages)) BiocManager::install(new_packages, dependencies = TRUE)

if (!require("devtools")) install.packages("devtools")
if (!require("d3heatmap")) devtools::install_github("talgalili/d3heatmap")

#load required packages
lapply(required_packages, function(x) {
  suppressPackageStartupMessages(suppressMessages(require(x, character.only = T, quietly = T, warn.conflicts = F)))
})