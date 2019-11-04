if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("limma")

setRepositories(ind = 1:3)
devtools::install_github('dylanbeeber/crispRdesignR')

install.packages("tidyverse")

