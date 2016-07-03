
render.mir.expression <- function(mir) {
  renderUI({
    fluidPage(
      fluidRow(
        column(6,
          div(h2('miR Tissue Atlas'),
              render.miratlas(mir))
        ),
        column(6,
          div(h2('miRmine'),
              render.mirmine(mir)
          )
        )
      )
    )
  })
}

render.mirmine <- function(x) {
  mirmine <- readRDS('data/processed/mirmine_grouped_2016_april.Rds')
  renderPlot({
    ggplot(subset(mirmine, mir == x), 
           aes(x=group, y=expression)) + 
      geom_boxplot() + 
      geom_point() + 
      labs(x='Tissue', y='Expression') +
      coord_flip() + 
      theme_minimal()
  }, height = 700)
}

render.miratlas <- function(x) {
  mirtissueatlas <- readRDS('data/processed/mirtissueatlas_grouped_2016_may.Rds')
  renderPlot({
    ggplot(subset(mirtissueatlas, mir == x), aes(x=tissue, y=expression)) + 
      geom_boxplot() + 
      geom_point() + 
      labs(x='Tissue', y='Expression') +
      coord_flip() + 
      theme_minimal()
  }, height = 700)
}