taglist <- function(){
      
      tagList(
        h1("This is WheatTest Page!",style="text-align:center"),
        nav_links,


# To add ui part for page6   

useShinyjs(),

column(8,offset = 5, titlePanel("Wheat Alignment output")), # titlePanel

sidebarLayout(

sidebarPanel(

#column(12,textInput("Gene","Your gene name or job ID")),





column(12,wellPanel(div(id='my_textinput1' ,
                   textInput("Gene","Gene name (gene ID is available, such as Traes***) or job ID (for fasta sequence or preferred coordinates) " )))),
                   tags$style(type="text/css", "#my_textinput1 {color: red}","#my_textinput1 {font-size:14px;}"),


column(12,
pickerInput(
  inputId = "Pickgenome", 
  label = "Pick Genome (Please select one!) :", 

  choices = c('','IWGSCNEW','jagger','julius','lancer','landmark','mace','mattis','norin61','stanley','spelta','arinalrfor'),

  selected = c(''), ## by default

  options = list(
    'actions-box' = TRUE, 
    size = 15,
    'selected-text-format' = "count > 1"
  ), 
  multiple = FALSE,
)
),



column(12,
pickerInput(
  inputId = "Chr", 
  label = "Chromosome (Please select one!)", 
  choices = c('','chr1A','chr1B','chr1D','chr2A','chr2B','chr2D','chr3A','chr3B','chr3D','chr4A','chr4B','chr4D','chr5A','chr5B','chr5D','chr6A','chr6B','chr6D','chr7A','chr7B','chr7D'),
  selected = c(''), ## by default
  options = list(
    'actions-box' = TRUE, 
    size = 21,
    'selected-text-format' = "count > 1"
  ), 
  multiple = FALSE,
)
),


column(12,
pickerInput(
  inputId = "Pickformat", 
  label = "CDS or gene (Coding sequence only or full-length gene, We suggest using CDS); OR your fasta sequence; OR your preferred coordinates; OR your large file :", 
  choices = c('CDS','gene','fasta_seq','coordinates','file'),
  selected = c('CDS'), ## by default
  options = list(
    'actions-box' = TRUE, 
    size = 18,
    'selected-text-format' = "count > 1"
  ), 
  multiple = FALSE,
)
),


column(12,textAreaInput("fasta","Your fasta sequence (Please add first line: >query Before pasting your DNA sequence!)",height='100px')),
bsTooltip("fasta", ">query as The first line","right", options = NULL),


column(12,numericRangeInput("Coordinates", "Your preferred genomic coordinates (Size < 10kb) :", min = NA, max = NA, value =c(0,0))),


column(12,fileInput("upload", "Upload large file (chromosome, .gz)", multiple = FALSE)), ## 08/25/22 


column(12,textInput("Upstream","Upstream (bp)",value=0)),
column(12,textInput("Downstream","Downstream (bp)",value=0)),

column(12,
pickerInput(
  inputId = "id", 
  label = "Genomes (Defalt: all genomes selected) :", 

  choices = c('jagger','julius','lancer','landmark','mace','mattis','norin61','stanley','IWGSC','arinalrfor','spelta','query'),

  selected = c('jagger','julius','lancer','landmark','mace','mattis','norin61','stanley','IWGSC','arinalrfor','spelta','query'), ## by default

  options = list(
    'actions-box' = TRUE, 
    size = 15,
    'selected-text-format' = "count > 1"
  ), 
  multiple = TRUE,
)
),


column(12,sliderInput("Distancefilter", "Distance filter between mapped clusters (1kb-50kb) :", min = 1000, max = 50000, value =20000)),
column(12,sliderInput("CDSfilter", "Expected CDS size compared to IWGSC (fold change:0.25-4) :", min = 0.25, max = 4, value =c(0.75,1.25))),

# column(12,sliderInput("Position", "Zoom in (Left to right; bp):", min = 1, max = 50000, value = c(1,50000))),

actionButton("submit", label = "Submit",class = "btn-warning"),
bsTooltip("submit", "Please double check your input (format), and then submit your job","right", options = NULL),




actionButton("submit2", label = "Zoom In",class = "btn-warning"),

actionButton("Haplotypes", "Plot selected haplotypes",style = "background-color:#CCCCFF"),

actionButton("GeneCapture", label = "Involved Genes",style="color: #00F; background-color: #FFFF66; border-color: #c34113; border-radius: 10px; border-width: 2px"),
bsTooltip("GeneCapture", "Show captured genes","right", options = NULL),


  downloadButton('save','save',style = "background-color:#C9DD03"),
  downloadButton('Repeat','Repeat',style = "background-color:#00FFFF"),
  downloadButton('Variation','Variation',style = "background-color:#00FFFF"),
  downloadButton('Haplotype_fa','Haplotype_fa',style = "background-color:#00FFFF"),
  downloadButton('CDS_fa','CDS_fa',style = "background-color:#00FFFF"),



actionButton("Largefile", label = "User's chromosome",style="color: FF66B2; background-color: #FFFF66; border-color: #c34113; border-radius: 10px; border-width: 2px"),
bsTooltip("Largefile", "For large file upload only","right", options = NULL),


actionButton("done", label = "Refresh",style = "background-color:#FF6666"),
bsTooltip("done", "Clean all processed files for your submitted job!","right", options = NULL),


uiOutput("bucket"),

), # sidebarPanel


mainPanel(
fluidRow(
  
     
     column(12,verbatimTextOutput("visits")), ## number of visits ?

     column(12,verbatimTextOutput("info")), ## number of jobs submitted

     # span(textOutput("query_test"), style="color:blue"),
     textOutput('query_test'),
     tags$head(tags$style("#query_test{color: blue;
                                 font-size: 16px;
                                 font-style: italic;
                                 }"
                         )
     ),


     textOutput('coordinates_test'),
     tags$head(tags$style("#coordinates_test{color: blue;
                                 font-size: 16px;
                                 font-style: italic;
                                 }"
                         )
     ),


     # span(textOutput("fasta_test"), style="color:red"),
     textOutput('fasta_test'),
     tags$head(tags$style("#fasta_test{color: red;
                                 font-size: 16px;
                                 font-style: italic;
                                 }"
                         )
     ),

     textOutput('QueryCoordinate_test'),
     tags$head(tags$style("#QueryCoordinate_test{color: green;
                                 font-size: 16px;
                                 font-style: italic;
                                 }"
                         )
     ),


    # column(12,verbatimTextOutput("coordinate_info")),
    # span(textOutput("text"), style="color:green"),
    # span(textOutput("text2"), style="color:red"),
    # span(textOutput("text3"), style="color:black"),
    # span(textOutput("text4"), style="color:blue"),

     column(6, plotOutput("plot",click = "plot_click",dblclick = "plot_dblclick",width = "90%")),   ## Main plot 1, Clustered
     column(6, plotOutput("plot1",click ="plot1_click",dblclick = NULL,width = "90%")),   ## Main plot 1, Tree

     column(9,verbatimTextOutput("info1")),
     column(9,verbatimTextOutput("info2")), 

     column(6, plotOutput("plot2",click = NULL,dblclick = NULL,width = "90%")), # Zoom_in 
     column(6, plotOutput("plot5",click = NULL,dblclick = NULL,width = "90%")), # tree cut
      

      column(12, offset = 0,DT::dataTableOutput("table00"),style='padding-top:5px; padding-bottom:5px'), ## 06/17, cluster information
      
      column(12, offset = 0,DT::dataTableOutput("table0"),style='padding-top:5px; padding-bottom:5px'),

      column(12, offset = 0,DT::dataTableOutput("table2"),style='padding-top:5px; padding-bottom:5px'),

      dataTableOutput("table3"),


column(12, offset = 0,DT::dataTableOutput("table4"),style='padding-top:5px; padding-bottom:5px'),

textOutput('Captured_Gene'),
     tags$head(tags$style("#Captured_Gene{color: black;
                                 font-size: 20px;
                                 font-style: italic;
                                 }"
                         )
     ),


column(12, offset = 0,DT::dataTableOutput("tablecluster"),style='padding-top:5px; padding-bottom:5px')



)
) # mainPanel



) # sidebarLayout





      ) # For tagList

    }



