packages <- c('shiny', 'DT', 'data.table', 'uuid', 'inline', 'dplyr', 'devtools', 'rvest', 'httr', 'magrittr',
              'ggplot2', 'ggbio', 'RCurl', 'Cairo')
install.packages(packages)

library(devtools)
install_github('gokceneraslan/rcppleveldb')

library(BiocInstaller)
biocLite(c('org.Hseg.db', 'GenomicRanges', 'rtracklayer', 'gwascat'))
install.packages('https://cran.r-project.org/src/contrib/Archive/NCBI2R/NCBI2R_1.4.7.tar.gz', repos=NULL)
