
library(shiny)

shinyUI(

  navbarPage("miR-SNP", id='navbar.panel',
             tabPanel("Input",

                      fluidPage(

                        # Application title
                        titlePanel("SNP@miR"),
                        br(),

                        fluidRow(
                          column(7,
                                 wellPanel(
                                   h2('Input SNP List'),
                                   fluidRow(
                                     column(6,
                                            h4('Enter one SNP id per line'),
                                            #tags$div(strong("SNP ids")),
                                            tags$textarea(id="snp.id.textarea", rows=8, style="width:100%"),
                                            helpText('or select one of GRASP risk SNPs below'),
                                            selectInput('grasp.selection', "Risk SNPs",
                                                        "Grasp",
                                                        selectize = F)
                                     ),
                                     column(1, br(br(br(br(br(br(h4("or")))))))),
                                     column(5,
                                            h4('Upload a file'),
                                            fileInput('snp.file', 'Upload SNP list in TSV/CSV format',
                                                      accept=c('text/csv', 'text/tsv',
                                                               'text/comma-separated-values',
                                                               'text/tab-separated-values',
                                                               'text/plain', '.csv','.tsv', '.txt'),
                                                      width='100%'),
                                            #htmlOutput('text.below.file')
                                            helpText('then select the column for SNPs below or simply click on a SNP'),
                                            br(),
                                            br(),
                                            #br(),
                                            #br(),
                                            selectInput('snp.col.select', 'SNP column',
                                                        "",
                                                        selectize = F)

                                     )
                                   )
                                 ),

                                 fluidRow(
                                   column(12,
                                          wellPanel(
                                            h2('Configuration'),

                                            fluidRow(
                                              column(6,
                                                     h4('LD Proxy Search'),
                                                     sliderInput('ld.slider', 'LD Cutoff', 0.1, 1, 0.8, 0.1),
                                                     selectInput('ld.population', 'Population', c('African',
                                                                                                  'American',
                                                                                                  'European',
                                                                                                  'East Asia',
                                                                                                  'South Asia'),
                                                                 selected = 'European',
                                                                 selectize=F)
                                              ),


                                              column(4, offset = 1,
                                                     h4('micro-RNA Targets'),
                                                     checkboxGroupInput('mir.target.db', 'Select Target Databases',
                                                                        c('TargetScan','miranda','Starbase'),
                                                                        c('TargetScan','miranda','Starbase')
                                                     )
                                              )
                                            )
                                          ))),

                                 fluidRow(
                                   column(1,
                                          actionButton('submit.button', 'Submit', class='btn-primary')
                                   )
                                 )
                          ),
                          column(5,

                                 fluidRow(
                                   column(12,
                                          DT::dataTableOutput('snp.table')
                                   )
                                 )
                          ),
                          br()
                        )
                      )),

             tabPanel('Results',
                      img(src='spinner.gif', id='spinner'),
                      DT::dataTableOutput('result.table')
             ),
             tabPanel('Help',
                      img(src='spinner.gif', id='spinner')
             ),
             singleton(
             tags$head(tags$script(src = "message-handler.js"))
             )
  )
)