
library(shiny)

shinyUI(
  #hide top navigation bar
  span(tags$head(tags$style(HTML('.navbar-static-top { display: none; }'))),
       navbarPage("Misina", id='navbar.panel',
                  tabPanel("Input", htmlOutput('input.page')),
                  
                  tabPanel('Results', htmlOutput('result.page')
                  ),
                  tabPanel('Help', img(src='spinner.gif', id='spinner')
                  ),
                  singleton(
                    tags$head(tags$script(src = "message-handler.js"))
                  )
       ))
)