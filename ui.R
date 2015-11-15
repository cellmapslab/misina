
library(shiny)

shinyUI(
  navbarPage(a("miR-SNP", href='/'), id='navbar.panel',
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