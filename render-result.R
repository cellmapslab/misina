
options(shiny.table.class='data table table-hover')

render.result <- function(df) {
  
  df <- unique(df)
  rownames(df) <- NULL
  snp.list <- (split(df, as.factor(df$SNP)))
  rendered.snp.list <- lapply(snp.list, render.SNP)
  
  renderUI({
    fluidPage(
      div(h2(a(span(tags$i(class='fa fa-dot-circle-o'), 
                 "Misina"), 
            href='/', style="font-family: 'Lobster', cursive; color: black; text-decoration: none;")
      )),
      br(),
       tabsetPanel(
         tabPanel('Results',
                  br(),
#                  navlistPanel(
#                    "Header",
#                    tabPanel("First",
                             rendered.snp.list
#                    ),
#                    tabPanel("Second",
#                             h3("This is the second panel")
#                    ),
#                    "-----",
#                    tabPanel("Third",
#                             h3("This is the third panel")
#                    )
#                  )
         ),
        tabPanel('Table view',
                 br(),
                 DT::renderDataTable(df, style='bootstrap', width='100%')
        ),
        tabPanel('Download results', 
                 br(),
                 downloadLink('download.results', 'Click here to download the results.'))))
  })
}

render.SNP <- function(snp) {
  snp.id <- unique(snp$SNP) 
  snp.pos <- unique(snp$SNP.position)
  snp.gene <- unique(snp$gene)
  
  ld.df <- snp[, c('Distance', 'IsProxyOf', 'R2')]
  rendered.ld <- render.LD(ld.df)
  
  mirs <- split(snp, snp$mir)
  rendered.mirs <- lapply(mirs, render.miR)
  #rendered.mirs <- render.miR2(snp)
  rendered.eqtl <- render.eQTL(snp)
  rendered.gwas <- render.gwas(snp)
  
  span(
    fluidRow(column(3, wellPanel(div(h4(strong('SNP: '), snp.id),
                                     h4(strong('Position: '), snp.pos), 
                                     h4(strong('Gene: '), snp.gene), 
                                     h4(strong('LD Information:')),
                                     br(),
                                     rendered.ld, 
                                     rendered.gwas,
                                     style='font-family: monospace;'
    ))),
    column(9, span(rendered.mirs, rendered.eqtl))),
    hr()
  )
}

render.LD <- function(ld) {
  ld <- unique(ld)
  rownames(ld) <- NULL
  ld <- ld[,c('IsProxyOf', 'R2', 'Distance')]
  ld <- ld[order(ld$R2, decreasing = T),]
  colnames(ld)[1] <- 'Original Risk SNP'
  
  
  span(
    renderTable(ld, include.rownames=F),
    style='display: inline-block;vertical-align: top;')
}

render.miR <- function(mir) {
  
  mir.name <- unique(mir$mir)
  mir.target <- unique(mir$gene)
  mir.target.pos <- unique(mir$mir.target.pos)
  mir.accession <- unique(mir$mirbase_acc)
  mir.target.db <- unique(mir$mir.target.db)
  mir.miranda.conserved <- unique(mir$miranda.conserved)
  mir.pred.score <- unique(mir$score)
  mir.seed.category <- unique(mir$seed.category)
  mir.snp.pos <- unique(mir$SNP.position.in.miR)
  
  span(h4(strong(paste0('miRNA: ', mir.name))), 
       tags$table(
         tags$tr(
           tags$td('Target'), tags$td(mir.target)),
         tags$tr(
           tags$td('Target position'), tags$td(mir.target.pos)),
         tags$tr(
           tags$td('Target database'), tags$td(mir.target.db)),
         tags$tr(
           tags$td('miRNA accession'), tags$td(mir.accession)),
         tags$tr(
           tags$td('Target prediction score'), tags$td(mir.pred.score)),
         tags$tr(
           tags$td('miRNA seed category'), tags$td(mir.seed.category)),
         tags$tr(
           tags$td('miRNA SNP position'), tags$td(mir.snp.pos)),
         class='data table table-hover', style='display: inline-block;'),
       style='display: inline-block;')
}

render.miR2 <- function(mir) {
  df <- mir[, c('mir', 'mir.target.pos', 'mirbase_acc', 'mir.target.db',
                'miranda.conserved', 'score', 'seed.category', 'SNP.position.in.miR')]
  df <- unique(df)
  df <- as.data.frame(t(df))
  
  span(
    h4(strong('miRNA')),
    renderTable(df, include.colnames=F),
    style='display: inline-block;vertical-align: top;')
  
}

render.eQTL <- function(snp) {
  df <- snp[, c('eQTL.Gene', 'eQTL.tstat', 'eQTL.pvalue', 'eQTL.Tissue')]
  colnames(df) <- c('Regulated Gene', 't-statistic', 'p-value', 'Tissue')
  df <- unique(df)
  df <- df[!apply(df, 1, function(x)all(is.na(x))),]
  
  if (nrow(df) > 0) {
    df <- as.data.frame(t(df))
    span(
      h4(strong('eQTL Support')),
      renderTable(df, include.colnames=F),
      style='display: inline-block;vertical-align: top;')
  } else {
    span() 
  }
}

render.gwas <- function(snp) {
  df <- snp[, c('DiseaseBySNP', 'DiseaseByGene')]
  colnames(df) <- c('SNP-assoc.', 'Gene-assoc.')
  df <- unique(df)
  df <- df[!apply(df, 1, function(x)all(is.na(x))),]
  
  if (nrow(df) > 0) {
    span(
      br(),
      h4(strong('Disease associations:')),
      renderTable(df, include.rownames=F))
    #style='display: inline-block;vertical-align: top;')
  } else {
    span() 
  }
}