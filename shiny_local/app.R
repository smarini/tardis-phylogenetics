library(shiny)
library(shinyDirectoryInput)

options(shiny.maxRequestSize = 100*1024^2) # Max 100MB per file

ui <- fluidPage(

  includeCSS("tardis.css"),
  
  tags$style(HTML("
        input[type=number] {
              -moz-appearance:textfield;
        }
        input[type=number]::{
              -moz-appearance:textfield;
        }
        input[type=number]::-webkit-outer-spin-button,
        input[type=number]::-webkit-inner-spin-button {
              -webkit-appearance: none;
              margin: 0;
        }
    ")),  
  
  headerPanel("Tardis"),
  
  fluidRow(class = "text-center",
    column( 3, offset = 1,
           numericInput("n.samples", "Sequences per individual", 5, min = 2, max = 500),
           verbatimTextOutput("n.samples"),
           numericInput("gen.size", "Individuals per generation", 100, min = 100, max = 10000),
           verbatimTextOutput("gen.size"),
           numericInput("n.cores", "Number of cores", 1, min = 1, max = 8),
           verbatimTextOutput("n.cores"),
           # numericInput("n.subsamples", "Subsamples to consider", 1, min = 1, max = 10),
           # verbatimTextOutput("n.subsamples")
           ),
    
    column( 3, offeset = 0,
           numericInput("frac.new", "Fraction of random individuals", 0.05, min = 0, max = 1, step = 0.01),
           verbatimTextOutput("frac.new"),
           numericInput("frac.evo", "Fraction of evolved individuals", 0.9, min = 0, max = 1, step = 0.01),
           verbatimTextOutput("frac.evo"),
           numericInput("frac.eli", "Fraction of elite individuals", 0.05, min = 0, max = 0.15, step = 0.01),
           verbatimTextOutput("frac.eli"),
           numericInput("n.gen", "Generations", 10, min = 2, max = 1000),
           verbatimTextOutput("n.gen") ),
    
    column( 3, offset = 0,
           numericInput("w.div", "Genetic diversity weight", 0.5, min = 0, max = 1, step = 0.01),
           verbatimTextOutput("w_div"),  
           numericInput("w.tem", "Temporal diversity weight:", 0.5, min = 0, max = 1, step = 0.01),
           verbatimTextOutput("w_tem"),
           selectInput("dist.opt", "Optimize genetic diversity",
                       c("Max diversity" = "max",
                         "Representative of set mean" = "mean",
                         "Representative of set median" = "median") ),
           actionButton("JC", HTML("Jukes-Cantor <br/> distance"), style="color: #fff; font-weight: 800; background-color: #6495ED;")
           )
    ),
  
  hr(),

  hr(),
  
  fluidRow(  
    column( 3, offset = 1,
            fileInput("distance", "Distance file (csv, rds)",
                      accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv", ".rds", ".RDS"), placeholder = '...') ),
    column( 3, offset = 0,
            fileInput("gen.file", "Sequence file (fasta)",
                      accept = c(".fa", ".fasta"), placeholder = '...') ),
    column( 3, offset = 0,
            fileInput("metadata", "Metadata file (csv)",
                      accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"), placeholder = '...') )
    ),

  fluidRow(  
    column( 4, offset = 1,
            textInput("data.set", "Dataset name", value = "example.dataset", width = NULL)
            ),
    column( 4, offset = 1,
            directoryInput('out.dir', label = 'Select output directory', value = '../shiny_output') )
  ),
  
  fluidRow(  
    column( 9, offset = 1,
            textOutput("to_run") )
    ),
  
  fluidRow(  
    column( 9, offset = 1,
            textOutput("search_space") )
    ),
  
  fluidRow(  
    column( 5, offset = 0),
    column( 3, offset = 0,
            actionButton("Run", "Run Tardis", style="color: #fff; font-weight: 800; background-color: #6495ED;"))
  ),
  
  fluidRow(  
    column( 9, offset = 1,
            h3(textOutput("running")) )
  ),
  
  fluidRow(  
    column( 9, offset = 1,
            h3(textOutput("JC")) )
  ),

  hr()
  
  )
  
server <- function(input, output, session) {

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$out.dir
    },
    handlerExpr = {
      if (input$out.dir > 0) {
        path = choose.dir(default = readDirectoryInput(session, 'out.dir'))
        
        # update the widget value
        updateDirectoryInput(session, 'out.dir', value = path)
        }
    }
  )
    
  output$to_run <- renderText({ paste("Tardis will run on", input$tot.gen.size,
                                      "individuals per generation, handled by", input$n.cores,
                                      "core(s), for", input$n.gen,
                                      "generations, counted as 0 to", (input$n.gen-1), "\n"
                                      )
  })
  
  output$search_space <- reactive({ 
    metadata <- input$metadata$datapath
    
    search_space_message <- "The total combination of possible subsamples is"
    
    if (is.null(metadata)) { return(paste(search_space_message, "not yet calculated (Please load metadata file)")) }
    
    metadata <- read.table(metadata, sep = ',', stringsAsFactors = FALSE, header = TRUE)
    
    search_space = choose(dim(metadata)[1], input$n.samples)
    
    if(search_space == Inf){ return(paste(search_space_message, "too large to be calculated on this machine")) }
    if(search_space == 0) { return(paste(search_space_message, "impossible to calculate. More sequences per subset than sequences in the data set")) }
    return(paste0(search_space_message, " ", search_space, "."))
  })
  
  distance <- eventReactive(input$Run, {
    return(input$distance$datapath)
  })
  
  gen.file <- eventReactive(list(input$Run, input$JC), {
    return(input$gen.file$datapath)
  })
  
  metadata <- eventReactive(input$Run, {
    return(input$metadata$datapath)
  })

  data.set <- eventReactive(input$Run, {
    return(input$data.set)
  })
  
  output.directory <- eventReactive(list(input$Run, input$JC), {
    return(readDirectoryInput(session, 'out.dir'))
  })
  
  output$JC <- eventReactive(input$JC, {
    condition = !is.null(input$gen.file$datapath)
    
    if(condition){
      out.file = paste(output.directory(), "jc.distance.precalc.rds", sep = '/')
      JC.command = paste("Rscript ../bin/JC.pairwise.dist.R",
                         "-i", gen.file(),
                         "-c", input$n.cores,
                         "-d ", out.file)
      message(JC.command)
      system(JC.command)
      if(file.exists(out.file)){
        return(paste("Jukes-Cantor distance calculated and saved in", out.file))
        }else{
          return("Error in processing alignment file.")
        }
      }else{
        return("Please load the sequence aligned fasta file to calculate the Jukes-Cantor distance.")
      }
  })

  output$running <- eventReactive(input$Run, {
    
    conditions = (round(input$frac.new*input$gen.size) +
                    round(input$frac.evo*input$gen.size) +
                    round(input$frac.eli*input$gen.size)) == input$gen.size &
                    !is.null(metadata()) &
                    !is.null(distance()) &
                    !is.null(gen.file()) &
                    !is.null(output.directory())

    if(conditions){
      withProgress(message = 'Executing Tardis...', detail = "Running the GA", value = 0, {
        system(paste('mkdir -p', output.directory()))
  
        incProgress(0, detail = paste("Initializing gen 0"))
        # init
        GA.command = paste("Rscript ../bin/make.gen.R",
                           "--distance", distance(),
                           "--metadata", metadata(),
                           "--gen.size", input$gen.size,
                           "--dist.opt", input$dist.opt,
                           "--n.samples", input$n.samples,
                           "--n.cores", input$n.cores,
#                           "--basedir", '../',
                           "--w.tem", input$w.tem,
                           "--w.div", input$w.div,
                           "--seeds", "../data/seeds.txt",
                           "--out.dir", output.directory(),
                           "--data.set", data.set(),
                           "--generation 0 --frac.new 1 --frac.evo 0 --frac.eli 0")
        message(GA.command)
        system(GA.command)
        message('Init over')
        for (i in 1:(input$n.gen-1)) {
          incProgress(1/(input$n.gen-1), detail = paste("Working on gen", i))
          
          message(i)
          
          GA.command = paste("Rscript ../bin/make.gen.R",
                             "--distance", distance(),
                             "--metadata", metadata(),
                             "--generation", i,
                             "--frac.new", input$frac.new,
                             "--frac.evo", input$frac.evo,
                             "--frac.eli", input$frac.eli,
                             "--gen.size", input$gen.size,
                             "--n.samples", input$n.samples,
                             "--n.cores", input$n.cores,
#                             "--basedir", '../',
                             "--w.tem", input$w.tem,
                             "--w.div", input$w.div,
                             "--seeds", "../data/seeds.txt",
                             "--out.dir", output.directory(),
                             "--data.set", data.set(),
                             "--dist.opt", input$dist.opt)
          message(GA.command)
          system(GA.command)
        }
        
        message('Reading results')
        if(all(file.exists(paste0(output.directory(), "/GA.", data.set(), ".", 0:(input$n.gen-1), ".1.indeces.fitness.csv")),
               file.exists(paste0(output.directory(), "/GA.", data.set(), ".", 0:(input$n.gen-1), ".1.indeces.subsamples.csv"))
               )){

          system(paste("Rscript ../bin/print.results.R",
                       "--distances", distance(),
                       "--n.batches 1",
                       "--ngenerations", input$n.gen,
                       "--outdir", output.directory(),
                       "--data_set", data.set(),
                       "--metadata", metadata()
                         ) )
          message("Reports and results printed")
          command.extractseqs = (paste("python3 ../bin/extractSeqs.py",
                                       output.directory(),
                                       gen.file(),
                                       input$n.gen-1,
                                       1, # n.batches
                                       42,
                                       paste('GA', data.set(), sep = '.'))
                                 )
          system(command.extractseqs)
          
          system(paste0('rm -f ', output.directory(), '/*indeces* ', output.directory(), '/*fitness*'))
          return("Done")
        }else{
          return("Error processing data. Please check input file format and parameter choice, as well as output directory.")
        }
      })
    }else{
      return(paste('Wrong intial conditions. Distance, sequence, and metadata files should all be selected.',
                   'Fractions of random, evolved, and elite generated individuals currently sum to',
                   (input$frac.new + input$frac.evo + input$frac.eli), '(should sum to 1).\n', 
                   'Number of random, evolved, and elite generated individuals currently sum to',
                   (round(input$frac.new*input$tot.gen.size) +
                      round(input$frac.evo*input$tot.gen.size) +
                      round(input$frac.eli*input$tot.gen.size)),
                   paste0('(should sum to ', input$tot.gen.size, ').')
                   ))
    }
  })
}

shinyApp(ui, server)
