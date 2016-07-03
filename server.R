library(shiny)
library(shinyjs)
library(DT)
library(data.table)
library(dplyr)
library(uuid)
library(inline)
library(grasp2db)

generate.uuid <- function() {
  paste0(strsplit(UUIDgenerate(), '-')[[1]], collapse = '')
}

result.file.from.id <- function(i) {
  if (!dir.exists('jobs'))
    dir.create('jobs')
  
  paste0('jobs/', i, '.rds')
}

error.file.from.id <- function(i) {
  if (!dir.exists('jobs'))
    dir.create('jobs')
  
  paste0('jobs/error-', i, '.rds')
}


get.grasp.cat <- function() {
  ret <- GRASP2() %>% tbl(., 'study') %>% select(PaperPhenotypeCategories) %>% distinct %>% arrange(PaperPhenotypeCategories) %>% as.data.frame %>% `[[`(.,1)
  ret <- sort(unique(trimws(unlist(strsplit(ret, ';')))))
  
  as.list(ret)
}

get.grasp.phenotype <- function(ct) {
  ret <- GRASP2() %>% 
    tbl(., 'study') %>% 
    select(PaperPhenotypeCategories, PaperPhenotypeDescription) %>%
    filter_(.dots=paste0('PaperPhenotypeCategories %like% "%', ct, '%"')) %>%
    select(PaperPhenotypeDescription) %>% 
    distinct %>% 
    as.data.frame %>% `[[`(., 1) %>% 
    sort
  
  #"Beh\xe7et's disease" causes problems
  Encoding(ret) <- 'latin1'
  
  as.list(ret)
}

#set file upload limit to 5M
options(shiny.maxRequestSize=5*1024^2)

shinyServer(function(input, output, session) {
  
  show.button <- T

  source('render-result.R')
  source('render-mirexp.R')
  
  #redirect to result check page if ?job= is given
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query[['job']])) {
      result.check.page(query[['job']])
      updateTabsetPanel(session, 'navbar.panel', selected='Results')
    } else if(!is.null(query[['mir']])) {
      mir.check.page(query[['mir']])
      updateTabsetPanel(session, 'navbar.panel', selected='microRNA Expression')
    }
  })
  
  re.file.data.frame <- reactive({
    snp.file <- input$snp.file
    if(is.null(snp.file))
      return(NULL)
    as.data.frame(fread(snp.file$datapath))
  })
  
  re.selected.column <- reactive({
    #     if (is.null(input$snp.table_columns_selected))
    #       return(NULL)
    #     
    #     cols <- as.integer(input$snp.table_columns_selected) + 1
    
    if (is.null(input$snp.col.select) || input$snp.col.select == '')
      return(NULL)
    which(colnames(re.file.data.frame()) == input$snp.col.select)
  })
  
  re.is.selected.column.valid <- reactive({
    
    if(is.null(re.file.data.frame()) || is.null(re.selected.column()))
      return(FALSE)
    
    snp.df <- re.file.data.frame()
    snp.col <- re.selected.column()
    
    if (any(substr(snp.df[,snp.col], 1, 2) == 'rs'))
      return(TRUE)
    else
      return(FALSE)
  })
  
  re.snp.df <- reactive({
    
    snp.area <- input$snp.id.textarea
    
    if (snp.area != '')
      return(data.frame(SNPs=strsplit(snp.area, '\\s+|,')[[1]], stringsAsFactors = F))
    
    snp.df <- re.file.data.frame()
    snp.col <- re.selected.column()
    snp.col.valid <- re.is.selected.column.valid()
    
    if (!snp.col.valid)
      return(FALSE)
    
    colnames(snp.df)[snp.col] <- 'SNPs'
    return(snp.df)
    
  })
  
  output$snp.table <- DT::renderDataTable({
    if(is.null(re.file.data.frame()))
      return(NULL)
    
    DT::datatable(re.file.data.frame(),
                  rownames = F,
                  style='bootstrap',
                  options = list(pageLength = 19,
                                 ordering=F,
                                 lengthChange=F,
                                 searching=F),
                  #scrollX=T,
                  #scrollY=T),
                  selection = list(mode='single', target='column'))
  })
  
  output$text.below.file <- renderText({
    if(is.null(re.selected.column()))
      return(as.character(h4("Select a SNP from the spreadsheet")))
    
    if(!re.is.selected.column.valid())
      return(as.character(h4("Does not seem like a SNP column...")))
    
    as.character(h4("Perfect!"))
  })
  
  observeEvent(re.file.data.frame(), {
    if(!is.null(re.file.data.frame())) {
      updateSelectInput(session, 'snp.col.select',
                        choices = colnames(re.file.data.frame()),
                        selected = '')
    }
    
  })
  
  observeEvent(input$snp.table_columns_selected, {
    if(!is.null(input$snp.table_columns_selected)) {
      
      cols <- as.integer(input$snp.table_columns_selected) + 1
      updateSelectInput(session, 'snp.col.select',
                        selected = colnames(re.file.data.frame())[cols])
    }
  })
  
  #update grasp phenotypes based on the category selection
  observeEvent(input$grasp.cat, {
    if (!is.null(input$grasp.cat) && input$grasp.cat != '') {
      phs <- get.grasp.phenotype(input$grasp.cat)
      
      #if (!is.null(input$grasp.pheno) && input$grasp.pheno != '')
      #  ext <- as.list(input$grasp.pheno)
      #else
      ext <- list()
      
      updateSelectInput(session, 'grasp.pheno',
                        choices = phs,
                        selected = ext)
    }
  })
  
  #load SNP examples button
  observeEvent(input$loadsnp.button, {
    snps <- 'rs10476052 rs1048920 rs1049633 rs113767110 rs11739062 rs11749762 rs1931895'
    updateTextInput(session, 'snp.id.textarea', value = snps)
  })
  
  # main submit function ----------------------------------------------------
  
  observeEvent(input$submit.button, {
    
    check.res <- check.input()
    if(!check.res[[1]]) {
      session$sendCustomMessage(type = 'testmessage',
                                message =check.res[[2]])
      return(NULL)
      
    }
    
    i <- generate.uuid()
    result.scheduled.page(i)
    process(re.get.all.inputs(), i)
    updateTabsetPanel(session, 'navbar.panel', selected='Results')
  })
  
  check.input <- reactive({
    
    if (input$snp.id.textarea == '' && is.null(re.file.data.frame()) && is.null(input$grasp.pheno)) {
      return(list(FALSE, 'Type SNPs ids, pick a GRASP phenotype or upload a valid file!'))
    }
    
    if (input$snp.id.textarea != '') {
      #check if snp ids are proper
      if (!all(substr(strsplit(input$snp.id.textarea, '\\s+|,')[[1]], 1,2) == 'rs'))
        return(list(FALSE, 'SNP ids are not valid!'))
    }
    
    if (input$snp.id.textarea == '' && is.null(input$grasp.pheno) && !re.is.selected.column.valid()) {
      return(list(FALSE, 'Select a valid SNP column for the given file'))
    }
    
    if(is.null(input$mir.target.db) || length(input$mir.target.db) < 0)
      return(list(FALSE, 'Select at least one target database.'))
    
    return(list(TRUE, ''))
  })
  
  re.get.all.inputs <- reactive({
    
    return(list(snp.df = re.snp.df(),
                grasp.pheno = input$grasp.pheno,
                ld.cutoff = input$ld.slider,
                ld.population = input$ld.population,
                mir.target.db = input$mir.target.db
    ))
  })
  
  process <- function(inputs, i) {
    
    #fork and return the master
    p <- parallel:::mcfork(estranged=F)
    if (!inherits(p, "masterProcess")) {
      return()
    }
    
    # close all file descriptors, apperently closeAllConnections() is
    # not working this way
    ca <- cfunction(body = 'int maxfd=sysconf(_SC_OPEN_MAX); \
                            for(int fd=3;fd<maxfd;fd++)close(fd);return 0;',
                    includes='#include<unistd.h>')
    ca()
    
    if (!dir.exists('jobs'))
      dir.create('jobs')
    
    #something 
    withCallingHandlers({
      #child node executes the rest
      source('helper.R')
      res <- run.pipeline(inputs)
      f <- result.file.from.id(i)
      saveRDS(res, f)
    }, error = function(e){
      t <- sys.calls()
      #humanize output of sys.calls()
      t <- paste0(paste0(rev(seq_along(t)), ': '), rev(t), collapse='\n\n')
      saveRDS(list(error=e, 
                   traceback=t,
                   inputs=inputs), file=error.file.from.id(i))
    })
    stop()
  }
  
  #  Results section --------------------------------------------------------
  
  result.scheduled.page <- function(i){
    output$result.page <- renderUI({
      l <- paste0('?job=', i)
      div(br(),
          span(id='processing_span',
               p('Your analysis is now scheduled. The results will be available ', 
                 a(id='resultlink', 'here', href=l), '.')
               #img(src='spinner.gif', id='spinner')
          ))
    })
  }
  
  output$result.page <- renderUI({h1('...')})
  output$mir.page <- renderUI({h1('...')})
  
  mir.check.page <- function(mir) {
    output$mir.page <- render.mir.expression(mir) 
  }
  
  result.check.page <- function(i) {
    
    f <- result.file.from.id(i)
    if (file.exists(f)) {
      res <- readRDS(f)
      output$result.page <- render.result(res)
    } else { #result not found, show error
      
      if (file.exists(error.file.from.id(i))) {
        err <- readRDS(error.file.from.id(i))
        
        query <- parseQueryString(session$clientData$url_search)
        if (!is.null(query[['debug']])) {
          output$result.page <- renderUI({
            div(p('An error occurred: '),
                pre(as.character(err$error)),
                pre(as.character(err$traceback))
            )
          })
        } else {
          output$result.page <- renderUI({
            div(p('An error occurred: '),
                pre(as.character(err$error))
            )
          })
        }
      } else {
        
        output$result.page <- renderUI({
          div(br(),
              p('Your request is still being processed.'))
        })
      }
    }
  }
  
  #download result handler
  output$download.results <- downloadHandler(
    filename = function() { paste('results-', 
                                  format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), 
                                  '.csv', sep='') },
    content = function(file) {
      query <- parseQueryString(session$clientData$url_search)
      if (!is.null(query[['job']])) {
        i <- query[['job']]
        f <- result.file.from.id(i)
        if (file.exists(f)) {
          res <- readRDS(f)
          write.csv(res, file, row.names = F)
        }
      }
    }
  )
  
  observeEvent(input$conf.button, {
    toggle('configuration.row')
    if (show.button) {
      html(id='conf.button', 'Hide configuration')
      show.button <<- F
    } else {
      html(id='conf.button', 'Show configuration')
      show.button <<- T
    }
  })
  
  output$input.page <- renderUI({
    fluidPage(
      useShinyjs(),
      # Application title
      titlePanel(div(a(span(tags$i(class='fa fa-dot-circle-o'), 
                            "Misina"), 
                       href='/', style="font-family: 'Lobster', cursive; color: black; text-decoration: none;")
      )),
      br(),
      
      fluidRow(
        column(6,
               wellPanel(
                 h4('Input SNP List'),
                 fluidRow(class='vdivide',
                          column(6,
                                 div(
                                   h5('A'),
                                   helpText('Enter SNP ids (separated by space, newline or comma)'),
                                   #tags$div(strong("SNP ids")),
                                   tags$textarea(id="snp.id.textarea", rows=8, style="width:100%"),
                                   div(actionButton('loadsnp.button', 'Load examples', class='btn'), 
                                       style='text-align: right;'),
                                   br(),
                                   style='margin-right: 15px;')
                          ),
                          column(6,
                                 div(
                                   h5('B'),
                                   helpText('or upload a SNP file in TSV/CSV format'),
                                   fileInput('snp.file', '',#'SNP file in TSV/CSV format',
                                             accept=c('text/csv', 'text/tsv',
                                                      'text/comma-separated-values',
                                                      'text/tab-separated-values',
                                                      'text/plain', '.csv','.tsv', '.txt'),
                                             width='100%'),
                                   helpText('then select the column for SNPs below or simply click on a SNP'),
                                   selectInput('snp.col.select', 'SNP column',
                                               "",
                                               selectize = F),
                                   style='margin-left: 15px;')
                          )
                 ),
                 hr(),
                 h5('C'),
                 helpText('or select below a disease to analyze the SNPs associated with it (via ', 
                          a('GRASP', href='http://grasp.nhlbi.nih.gov/', target="_blank"), 
                          ')'),
                 fluidRow(
                   column(5,
                          selectInput('grasp.cat', "Disease/Phenotype Category",
                                      get.grasp.cat(),
                                      multiple = F,
                                      selectize = T)),
                   column(2, 
                          br(),
                          div(
                            span(class="glyphicon glyphicon-chevron-right", `aria-hidden`="true"),
                            style='text-align:center;')),
                   column(5,
                          selectInput('grasp.pheno', "Phenotype",
                                      '',
                                      multiple = T,
                                      selectize = T)
                   )
                 )
               ),
               fluidRow(id='configuration.row', style='display: none;',
                        column(12,
                               wellPanel(
                                 h4('Configuration'),
                                 fluidRow(
                                   column(6,
                                          sliderInput('ld.slider', 
                                                      'LD Proxy Cutoff', 
                                                      0.1, 1, 0.8, 0.1),
                                          selectInput('ld.population', 
                                                      'Population', 
                                                      list('African'    = 'afr',
                                                           'American'   = 'amr',
                                                           'European'   = 'eur',
                                                           'East Asia'  = 'eas',
                                                           'South Asia' = 'sas'),
                                                      selected = 'eur',
                                                      selectize=F)
                                   ),
                                   column(4, offset = 1,
                                          checkboxGroupInput('mir.target.db', 'Select miRNA Target Databases',
                                                             c('TargetScan','miranda','Starbase'),
                                                             c('TargetScan','miranda','Starbase')
                                          )
                                   )
                                 )
                               ))),
               fluidRow(
                 column(12,
                        div(
                          actionButton('submit.button', 'Submit', class='btn-primary'),
                          HTML('&nbsp;'),
                          actionButton('conf.button', 'Show configuration', class='btn')
                        )
                 )
               )
        ),
        column(6,
               fluidRow(
                 column(12,
                        DT::dataTableOutput('snp.table')
                 )
               )
        ),
        br()
      ))
  })
  
})
