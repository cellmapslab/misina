
shiny.table.class <- 'data table-condensed table-hover'
options(shiny.table.class=shiny.table.class)

render.result <- function(df) {
  
  df <- unique(df)
  rownames(df) <- NULL
  additional.columns <- attr(df, 'additional.columns')
  
  #use score.snps function to sort the SNP list
  df$snp.priority <- score.snps(df)
  snp.list <- split(df, as.factor(df$SNP))
  snp.ord <- order(sapply(snp.list, function(x)max(x$snp.priority)), decreasing = T)
  rendered.snp.list <- lapply(snp.list, render.SNP, additional.columns)
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

render.SNP <- function(snp, additional.columns) {
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
  
  ld.df <- snp[, c('IsProxyOf', 'Distance', 'R2', additional.columns)]
  colnames(ld.df)[1] <- 'Original Risk SNP'
  rendered.ld <- render.LD(ld.df)
  
  mirs <- split(snp, snp$mir)
  rendered.mirs <- lapply(mirs, render.miR)
  #rendered.mirs <- render.miR2(snp)
  
  # TODO: her satir icin ayri eQTL'ler render et!
  eqtls <- split(snp, snp$eQTL.Tissue)
  rendered.eqtl <- lapply(eqtls, render.eQTL)
  rendered.gwas <- render.gwas(snp)
  
  span(
    fluidRow(column(3, wellPanel(div(h4(strong('SNP: '), snp.id),
                                     h4(strong('Score: '), score.tag), 
                                     h4(strong('Position: '), snp.pos), 
                                     h4(strong('Gene: '), snp.gene), 
                                     h4(strong('LD Information:')),
                                     #renderTable(as.data.frame(t(additional.col.df)), include.colnames=F),
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
  ld <- ld[order(ld$R2, decreasing = T),]
  
  span(
    renderTable(as.data.frame(t(ld)), include.colnames=F),
    style='display: inline-block;vertical-align: top;')
}

render.miR <- function(mir) {
  mir <- mir[, c('mir', 'mir.target.pos', 'mirbase_acc', 'mir.target.db', 'miranda.conserved', 'score', 'seed.category', 'SNP.position.in.miR')]
  mir <- unique(mir)
  mir.name <- mir$mir
  mir.target <- mir$gene
  mir.target.pos <- mir$mir.target.pos
  mir.accession <- mir$mirbase_acc
  mir.target.db <- mir$mir.target.db
  mir.miranda.conserved <- mir$miranda.conserved
  mir.pred.score <- mir$score
  mir.seed.category <- mir$seed.category
  mir.seed.category[is.na(mir$seed.category)] <- ''
  mir.snp.pos <- mir$SNP.position.in.miR
  
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
  
  snp <- snp[, c('eQTL.Gene', 'eQTL.tstat', 'eQTL.pvalue', 'eQTL.Tissue', 'eQTL.Gene.Same.as.Target.gene')]
  snp <- unique(snp)
  
  eqtl.gene <- snp$eQTL.Gene
  eqtl.tstat <- snp$eQTL.tstat
  eqtl.pvalue <- snp$eQTL.pvalue
  eqtl.tissue <- snp$eQTL.Tissue
  if (snp$eQTL.Gene.Same.as.Target.gene == T) {
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
  if (all(is.na(s))) 
    score <- 0
  else
    score <- as.integer(s == '7mer' | s == '8mer')
  
  #snp position
  score <- score + as.integer(df$SNP.position.in.miR < 13)
  
  #eqtl gene == target gene
  df$eQTL.Gene.Same.as.Target.gene[is.na(df$eQTL.Gene.Same.as.Target.gene)] <- F
  score <- score + as.integer(df$eQTL.Gene.Same.as.Target.gene == T)
  
  score
}