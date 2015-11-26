render.result <- function(df) {
  
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
                 rendered.snp.list
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
  ld.distance <- unique(snp$Distance)
  ld.r2 <- unique(snp$R2)
  orig.snp <- unique(snp$IsProxyOf)
  
  mirs <- split(snp, snp$mir)
  rendered.mirs <- lapply(mirs, render.miR)
  
  span(
    fluidRow(column(3, wellPanel(div(h3(paste0('SNP: ', snp.id)), 
                                     tags$table(
                                       tags$tr(
                                         tags$td('Position: '), tags$td(snp.pos)),
                                       tags$tr(
                                         tags$td('LD Information')
                                       )
#                                        tags$tr(
#                                          tags$td('Distance: '), tags$td(ld.distance)),
#                                        tags$tr(
#                                          tags$td('R2: '), tags$td(ld.r2)),
#                                        tags$tr(
#                                          tags$td('Original risk SNP: '), tags$td(orig.snp))
                                     )
    ))),
    column(9, span(rendered.mirs))),
    hr()
    )
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
  
  span(h4(paste0('miRNA: ', mir.name)), 
       tags$table(
         tags$tr(
           tags$td('Target: '), tags$td(mir.target)),
         tags$tr(
           tags$td('Target position: '), tags$td(mir.target.pos)),
         tags$tr(
           tags$td('Target database: '), tags$td(mir.target.db)),
         tags$tr(
           tags$td('miRNA accession: '), tags$td(mir.accession)),
         tags$tr(
           tags$td('Target prediction score: '), tags$td(mir.pred.score)),
         tags$tr(
           tags$td('miRNA seed category: '), tags$td(mir.seed.category)),
         tags$tr(
           tags$td('miRNA SNP position: '), tags$td(mir.snp.pos)),
         class='table', style='display: inline-block;'),
       style='display: inline-block;')
}