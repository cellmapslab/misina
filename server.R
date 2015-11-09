library(shiny)
library(DT)

source('R/parse-SNP.R')
source('R/snipe.R')

library(dplyr)
library(RcppLevelDB)
library(data.table)

# input data files
dbsnp.file <- 'data/processed/dbsnp.leveldb'

starbase.gr <- readRDS('data/processed/starbase.Rds')
targetscan.gr <- readRDS('data/processed/targetscan.Rds')
miranda.gr <- readRDS('data/processed/miranda.Rds')
mir.targets.gr <- merge.granges.aggressively(meta.columns=list(mir.target.db=c('Starbase', 'TargetScan', 'miranda')),
                                             starbase.gr, targetscan.gr, miranda.gr)


#set file upload limit to 35M
options(shiny.maxRequestSize=35*1024^2) 

shinyServer(function(input, output, session) {
  
  re.file.data.frame <- reactive({
    snp.file <- input$snp.file
    if(is.null(snp.file))
      return(NULL)
    
    as.data.frame(fread(snp.file$datapath))
  })
  
  re.selected.column <- reactive({
    if (is.null(input$snp.table_columns_selected))
      return(NULL)
    
    cols <- as.integer(input$snp.table_columns_selected) + 1
    cols
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
                        choices = colnames(re.file.data.frame()))
    }
    
  })
  
  observeEvent(re.selected.column(), {
    if(!is.null(re.selected.column())) {
      updateSelectInput(session, 'snp.col.select',
                        selected = colnames(re.file.data.frame())[re.selected.column()])
    }
  })
  
  observeEvent(input$submit.button, {
    check.res <- check.input()
    if(!check.res[[1]]) {
      session$sendCustomMessage(type = 'testmessage',
                                message =check.res[[2]])
      return(NULL)
      
    }
    
    updateTabsetPanel(session, 'navbar.panel', selected='Results')
  })
  
  check.input <- reactive({
    
    if(input$snp.id.textarea == '' && is.null(re.file.data.frame())){
      return(list(FALSE, 'Type SNPs ids or upload a valid file!'))
    }
    
    if(input$snp.id.textarea != ''){
      #check if snp ids are proper
      if (!all(substr(strsplit(input$snp.id.textarea, ',', fixed = T)[[1]], 1,2) == 'rs'))
        return(list(FALSE, 'SNP ids are not valid!'))
    } 
    
    if (input$snp.id.textarea == '' && !re.is.selected.column.valid()) {
        return(list(FALSE, 'Select a valid SNP column for the given file'))
    }
    
    return(list(TRUE, ''))
  })
  
  #  Results section --------------------------------------------------------
  
#    re.ld.data.frame <- reactive({
#      
#      if (!re.is.selected.column.valid())
#        return(NULL)
#      
#      if(is.null(re.file.data.frame()) || is.null(re.selected.column()))
#        return(NULL)
#      
#      snp.df <- re.file.data.frame()
#      snp.col <- re.selected.column()
#      
#      print('LD SEARCH IN ACTION')
#      colnames(snp.df)[snp.col] <- 'SNPs'
#      extend.with.LD(snp.df)
#    })
   
   output$result.table <- DT::renderDataTable({
     
     #isolate dediysek, hic calismaz demedik
     #results'a gectigin an calisir
     isolate({
       print('isolated code')
       ld.df <- re.ld.data.frame()
       if(is.null(ld.df))
         return(NULL)
       
       DT::datatable(ld.df)
     })
     
   })
  
})