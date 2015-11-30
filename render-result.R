
shiny.table.class <- 'data table-condensed table-hover'
options(shiny.table.class=shiny.table.class)

render.result <- function(df) {
  
  df <- unique(df)
  rownames(df) <- NULL
  
  #use score.snps function to sort the SNP list
  df$snp.priority <- score.snps(df)
  snp.list <- split(df, as.factor(df$SNP))
  snp.ord <- order(sapply(snp.list, function(x)max(x$snp.priority)), decreasing = T)
  rendered.snp.list <- lapply(snp.list, render.SNP)
  rendered.snp.list <- rendered.snp.list[snp.ord]
  df$snp.priority <- NULL
  
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
                 DT::renderDataTable(df, style='bootstrap', width='100%', rownames=F)
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
  score <- max(snp$snp.priority)
  if (score == 0) {
    score.tag <- tags$i(class='fa fa-ban')
  } else if (score > 3) {
    score.tag <- NULL  
  } else {
    score.tag <- span(replicate(score, tags$i(class='fa fa-thumbs-o-up'), simplify = F))
  }
  
  ld.df <- snp[, c('Distance', 'IsProxyOf', 'R2')]
  rendered.ld <- render.LD(ld.df)
  
  mirs <- split(snp, snp$mir)
  rendered.mirs <- lapply(mirs, render.miR)
  #rendered.mirs <- render.miR2(snp)
  
  eqtls <- split(snp, snp$eQTL.Gene)
  rendered.eqtl <- lapply(eqtls, render.eQTL)
  rendered.gwas <- render.gwas(snp)
  
  span(
    fluidRow(column(2, wellPanel(div(h4(strong('SNP: '), snp.id),
                                     h4(strong('Score: '), score.tag), 
                                     h4(strong('Position: '), snp.pos), 
                                     h4(strong('Gene: '), snp.gene), 
                                     h4(strong('LD Information:')),
                                     br(),
                                     rendered.ld, 
                                     rendered.gwas,
                                     style='font-family: monospace;'
    ))),
    column(10, span(rendered.mirs, rendered.eqtl))),
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
  
  s <- tolower(substr(as.character(mir.seed.category), 1, 4))
  if (s == '7mer' | s == '8mer') {
    seed.priority <- tags$i(class='fa fa-check-circle')
  } else {
    seed.priority <- NULL
  }
  
  if (as.integer(mir.snp.pos) < 13)
    snp.pos.priority <- tags$i(class='fa fa-check-circle')
  else
    snp.pos.priority <- NULL
  
  
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
           tags$td('miRNA seed category'), tags$td(span(mir.seed.category, '  ', seed.priority))),
         tags$tr(
           tags$td('miRNA SNP position'), tags$td(span(mir.snp.pos, snp.pos.priority))),
         class=shiny.table.class, style='display: inline-block;'),
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
  
  eqtl.gene <- unique(snp$eQTL.Gene)
  eqtl.tstat <- unique(snp$eQTL.tstat)
  eqtl.pvalue <- unique(snp$eQTL.pvalue)
  eqtl.tissue <- unique(snp$eQTL.Tissue)
  if (unique(snp$eQTL.Gene.Same.as.Target.gene) == T) {
    eqtl.priority <- tags$i(class='fa fa-check-circle')
  } else {
    eqtl.priority <- NULL
  }
  
  span(h4(strong(paste0('eQTL Support'))), 
       tags$table(
         tags$tr(
           tags$td('Gene'), tags$td(span(eqtl.gene, eqtl.priority))),
         tags$tr(
           tags$td('t-statistic'), tags$td(eqtl.tstat)),
         tags$tr(
           tags$td('p-value'), tags$td(eqtl.pvalue)),
         tags$tr(
           tags$td('Tissue'), tags$td(eqtl.tissue)),
         class=shiny.table.class, style='display: inline-block;'),
       style='display: inline-block;vertical-align: top;')
  
}

render.eQTL2 <- function(snp) {
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


score.snps <- function(df) {
 
  #seed category
  df$seed.category[is.na(df$seed.category)] <- ''
  s <- tolower(substr(as.character(df$seed.category), 1, 4))
  score <- as.integer(s == '7mer' | s == '8mer')
  
  #snp position
  score <- score + as.integer(df$SNP.position.in.miR < 13)
  
  #eqtl gene == target gene
  df$eQTL.Gene.Same.as.Target.gene[is.na(df$eQTL.Gene.Same.as.Target.gene)] <- F
  score <- score + as.integer(df$eQTL.Gene.Same.as.Target.gene == T)
  
  score
}