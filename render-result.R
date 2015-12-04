
shiny.table.class <- 'data table-condensed table-hover'
options(shiny.table.class=shiny.table.class)

library(ggplot2)

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
  
  renderUI({
    fluidPage(
      div(h2(a(span(tags$i(class='fa fa-dot-circle-o'), 
                 "Misina"), 
            href='/', style="font-family: 'Lobster', cursive; color: black; text-decoration: none;")
      )),
      br(),
      fluidRow(
        column(4, render.stats(df)),
        column(6, render.help())
      ),
      br(),
      br(),
      fluidRow(
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
                   render.table.view(df)
          ),
          tabPanel('Download results', 
                   br(),
                   downloadLink('download.results', 'Click here to download the results.')))
        ))
    
  })
}

render.stats <- function(df) { 
 
  original.snp.count <- attr(df, 'original.snp.count')
  total.snp.count    <- attr(df, 'total.snp.count')
  hit.count <- length(unique(df$SNP))
  
  info.table <- tags$table(tags$tr(tags$td(strong('Number of risk SNPs')), tags$td(original.snp.count)),
                           tags$tr(tags$td(strong('Number of total SNPs with LD proxies')), tags$td(total.snp.count)),
                           tags$tr(tags$td(strong('Number of hits'), tags$td(hit.count))),
                           class='data table-condensed table')
  
  tmp <- unique(df[, c('SNP', 'snp.priority')])
  plot <- renderPlot({
    qplot(tmp$snp.priority, xlab='SNP Scores (0-3)') +
      scale_x_discrete(limits=0:3) +
      scale_y_discrete() +
      theme_minimal()
    }, height = 200)
  
 fluidRow(
   column(5, plot),
   column(7, info.table)
 )
  
}

render.help <- function() {
  
  gen.panel <- function(header, body, type='info') {
   div(class=sprintf('panel panel-%s', type),
     div(header, class='panel-heading'),
     div(body, class='panel-body'))
  }
  
  tu <- tags$i(class='fa fa-thumbs-o-up')
  
  main.help.box <- div(
    div('Resulting SNP-miR pairs are scored based on the following criteria:'),
    br(),
    div(span(tu, HTML('&nbsp'), 'miRNA seed type: Scored if seed type is 7- or 8-mer')),
    div(span(tu, HTML('&nbsp'), 'Relative SNP position: Scored if relative position of SNP is within the miRNA binding site (SNPs in 1-12 bp from 5\' end of miRNA)')),
    div(span(tu, HTML('&nbsp'), 'miRNA target - eGene match: Scored if miRNA target gene match eQTL gene (miRNA gene identical to eGene)'))
  )
  
  gen.panel('SNP Scoring', main.help.box)
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
  
  rendered.eqtl <- render.eQTL(snp)
  rendered.gwas <- render.gwas(snp)
  gene.link <- a(snp.gene, 
                 href=paste0('http://www.genecards.org/cgi-bin/carddisp.pl?gene=', snp.gene),
                 target='_blank')
  
  if(substr(snp.id, 1, 2) == 'rs')
    snp.link <- a(snp.id, 
                  href=paste0('http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=', snp.id),
                  target='_blank')
  else
    snp.link <- snp.id
  
  span(
    fluidRow(column(3, wellPanel(div(h4(strong('SNP: '), snp.link),
                                     h4(strong('Score: '), score.tag), 
                                     h4(strong('Position: '), snp.pos), 
                                     h4(strong('Gene: '), gene.link), 
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
  ld$R2 <- round(ld$R2, 2)
  
  res <- lapply(seq_len(nrow(ld)), function(r){
    span(
      renderTable(as.data.frame(t(ld[r,,drop=F])), include.colnames=F),
      hr(),
      style='display: inline-block;vertical-align: top;')  
    })
  return(span(res))
}

render.miR <- function(mir) {
  mir <- mir[, c('mir', 'gene', 'mir.target.pos', 'mirbase_acc', 'mir.target.db', 'miranda.conserved', 'score', 'seed.category', 'SNP.position.in.miR')]
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
    seed.priority <- tags$i(class='fa fa-check-circle', 
                               `data-toggle`="tooltip",
                               `data-placement`="right",
                               title='8mer and 7mer seed types are prioritized')
  } else {
    seed.priority <- NULL
  }
  
  if (as.integer(mir.snp.pos) < 13)
    snp.pos.priority <- tags$i(class='fa fa-check-circle', 
                               `data-toggle`="tooltip",
                               `data-placement`="right",
                               title='SNPs located between 1-12 bp in 5\' end of miRNA are prioritized')
  else
    snp.pos.priority <- NULL
  
  if (!is.na(mir.accession) && mir.accession != '')
    tit <- a(strong(paste0('miRNA: ', mir.name)), 
             href = paste0('http://www.mirbase.org/cgi-bin/mature.pl?mature_acc=', mir.accession),
             target='_blank')
  else 
    tit <- strong(paste0('miRNA: ', mir.name))
  
  span(h4(tit), 
       tags$table(
         #tags$tr(
         #  tags$td('Target'), tags$td(mir.target)),
         tags$tr(
           tags$td(strong('Target position')), tags$td(mir.target.pos)),
         tags$tr(
           tags$td(strong('Target database')), tags$td(mir.target.db)),
         #tags$tr(
         # tags$td('miRNA accession'), tags$td(mir.accession)),
         tags$tr(
           tags$td(strong('Target prediction score')), tags$td(mir.pred.score)),
         tags$tr(
           tags$td(strong('Seed category')), tags$td(span(mir.seed.category, '  ', seed.priority))),
         tags$tr(
           tags$td(strong('miRNA SNP position')), tags$td(span(mir.snp.pos, snp.pos.priority))),
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
  
  snp <- snp[, c('eQTL.Gene', 'eQTL.beta', 'eQTL.tstat', 'eQTL.pvalue', 
                 'eQTL.Tissue', 'eQTL.identical.target')]
  
  snp <- unique(snp)
  no.na <- which(apply(snp, 1, function(x)!all(is.na(x))))
  i <- 1
  
  res <- lapply(no.na, function(row){
    
    s <- snp[row, ,drop=F]
    eqtl.gene <- s$eQTL.Gene
    #eqtl.tstat <- s$eQTL.tstat
    eqtl.effect <- s$eQTL.beta
    eqtl.pvalue <- format.pval(s$eQTL.pvalue, 2)
    eqtl.effect.p.t <- sprintf('%.2f (%s)', eqtl.effect, eqtl.pvalue)
    eqtl.tissue <- s$eQTL.Tissue
    eqtl.sameas <- s$eQTL.identical.target
    
    if (eqtl.sameas == T) {
      eqtl.priority <- tags$i(class='fa fa-check-circle', 
                               `data-toggle`="tooltip",
                               `data-placement`="right",
                               title='SNPs with identical target gene and eGene are prioritized')
    } else {
      eqtl.priority <- NULL
    }
    
    if(i == 1) {
      header <- h4(strong(paste0('eQTL Support')))
      gene.header <- tags$td(strong('Gene'))
      effect.header <- tags$td(strong('Effect size (p-value)'))
      #tstat.header <- tags$td('t-statistic')
      #pvalue.header <- tags$td('p-value')
      tissue.header <- tags$td(strong('Tissue'))
    }
    else {
      header <- h4(strong(HTML('&nbsp;')))
      gene.header <- NULL
      effect.header <- NULL
      #tstat.header <- NULL
      #pvalue.header <- NULL
      tissue.header <- NULL
    }
    
    i <<- i + 1
    span(header, 
         tags$table(
           tags$tr(
             gene.header, tags$td(span(a(eqtl.gene,
                                         href=paste0('http://www.genecards.org/cgi-bin/carddisp.pl?gene=', eqtl.gene),
                                         target='_blank'),
                                       eqtl.priority))),
           tags$tr(
             tissue.header, tags$td(eqtl.tissue)),
           tags$tr(
             effect.header, tags$td(eqtl.effect.p.t)),
           #tags$tr(
           # tstat.header, tags$td(eqtl.tstat)),
           #tags$tr(
           # pvalue.header, tags$td(eqtl.pvalue)),
           class=shiny.table.class, style='display: inline-block;'),
         style='display: inline-block;vertical-align: top;')
  })
  return(span(res))
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
  df$eQTL.identical.target[is.na(df$eQTL.identical.target)] <- F
  score <- score + as.integer(df$eQTL.identical == T)
  
  score
}

render.table.view <-function(df) {
  #colnames(df) <- gsub('.', ' ', colnames(df), fixed = T)
  
  div(#render.table.column.selection(df),
      DT::renderDataTable(df, style='bootstrap', width='100%', rownames=F)
  )
  
}

render.table.column.selection <- function(df) {
  itemfunc <- I("function(item, escape) {
            return '<div>' +
                (item.name ? 'eppek:' + escape(item) : 'duppek') +
            '</div>';}")
  
  selectizeInput('table.view.column.selection', 
                 label='Select columns to be displayed: ',
                 choices=colnames(df),
                 selected=colnames(df),
                 multiple = T,
                 width='100%',
                 options=list(
                   plugins=list('remove_button'),
                   render=list(option=itemfunc))
  )
}