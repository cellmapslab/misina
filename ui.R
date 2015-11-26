
library(shiny)

css.string <- "
  .navbar-static-top { display: none; }
  .row.vdivide [class*='col-']:not(:last-child):after {
    background: #e0e0e0;
    width: 1px;
    content: '';
    display:block;
    position: absolute;
    top:0;
    bottom: 0;
    right: 0;
    min-height: 70px;
  }
"

shinyUI(
  #hide top navigation bar
  span(tags$head(tags$style(HTML(css.string)),
                 tags$link(href="//maxcdn.bootstrapcdn.com/font-awesome/4.3.0/css/font-awesome.min.css", rel="stylesheet"),
                 tags$link(href='https://fonts.googleapis.com/css?family=Lobster', rel='stylesheet', type='text/css')),
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