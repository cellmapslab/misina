
library(shiny)

shinyUI(
  navbarPage("miR-SNP", id='navbar.panel',
             tabPanel("Input", htmlOutput('input.page')),
             
             tabPanel('Results', htmlOutput('result.page')
             ),
             tabPanel('Help', img(src='spinner.gif', id='spinner')
             ),
             singleton(
               tags$head(tags$script(src = "message-handler.js"))
             )
  )
)