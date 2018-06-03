
render.mir.expression <- function(mir) {
  renderUI({
    fluidPage(
      fluidRow(
        column(6,
          div(h4('miR Tissue Atlas'),
              render.miratlas(mir))
        ),
        column(6,
          div(h4('miRmine'),
              render.mirmine(mir)
          )
        )
      )
    )
  })
}

render.mirmine <- function(x) {
  mirmine <- readRDS('data/processed/mirmine_grouped_2016_april.Rds')
  df <- subset(mirmine, mir == x)
  if(nrow(df) == 0) return(h4('Expression info not found.'))
  
  renderPlot({
    ggplot(df, 
           aes(x=group, y=expression)) + 
      geom_boxplot() + 
      geom_point() + 
      labs(x='Tissue', y=paste0('miRNA Expression (', x, ')')) +
      coord_flip() + 
      theme_minimal()
  }, height = 400)
}

render.miratlas <- function(x) {
  mirtissueatlas <- readRDS('data/processed/mirtissueatlas_grouped_2016_may.Rds')
  df <- subset(mirtissueatlas, mir == x)
  if(nrow(df) == 0) return(h4('Expression info not found.'))
  
  renderPlot({
    ggplot(df, aes(x=tissue, y=expression)) + 
      geom_boxplot() + 
      geom_point() + 
      labs(x='Tissue', y=paste0('miRNA Expression (', x, ')')) +
      coord_flip() + 
      theme_minimal()
  }, height = 800)
}
