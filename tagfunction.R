tagfunction <- function(page_title,page_subtitle,Genome_choice,chromosome_choice, default_choice){
      
      tagList(
     #   h1("This is WheatTest Page!",style="text-align:center"),
         h1(page_title,style="text-align:center"),
        nav_links,


# To add ui part for page6   

useShinyjs(),

#column(8,offset = 5, titlePanel("Wheat Alignment output")), # titlePanel
column(8,offset = 5, titlePanel(page_subtitle)), # titlePanel

sidebarLayout(

sidebarPanel(

#column(12,textInput("Gene","Your gene name or job ID")),





column(12,wellPanel(div(id='my_textinput1' ,
                   textInput("Gene","Gene name (gene ID is available) or job ID (for fasta sequence or preferred coordinates) " )))),
                   tags$style(type="text/css", "#my_textinput1 {color: red}","#my_textinput1 {font-size:14px;}"),


column(12,
pickerInput(
  inputId = "Pickgenome", 
  label = "Pick Genome (Please select one!) :", 

 # choices = c('','IWGSCNEW','jagger','julius','lancer','landmark','mace','mattis','norin61','stanley','spelta','arinalrfor'),
  
  choices = Genome_choice,
  selected = c(''), ## by default

  options = list(
    'actions-box' = TRUE, 
    size = 30,
    'selected-text-format' = "count > 1"
  ), 
  multiple = FALSE,
)
),



column(12,
pickerInput(
  inputId = "Chr", 
  label = "Chromosome (Please select one!)", 
 # choices = c('','chr1A','chr1B','chr1D','chr2A','chr2B','chr2D','chr3A','chr3B','chr3D','chr4A','chr4B','chr4D','chr5A','chr5B','chr5D','chr6A','chr6B','chr6D','chr7A','chr7B','chr7D'),
  choices = chromosome_choice,
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
  label = "CDS (Coding sequence); OR your fasta sequence :", 
  choices = c('CDS','fasta_seq'),
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


# column(12,numericRangeInput("Coordinates", "Your preferred genomic coordinates (Size < 10kb) :", min = NA, max = NA, value =c(0,0))),

###############
#column(12,
#pickerInput(
#  inputId = "Pick_Parents", 
#  label = "Haplotype upload (Parent1 or/and Parent2) :", 
#  choices = c('','Parent1','Parent2'),
#  selected = c(''), ## by default
#  options = list(
#    'actions-box' = TRUE, 
#    size = 5,
#    'selected-text-format' = "count > 1"
#  ), 
#  multiple = TRUE,
#)
#),

column(12,fileInput("upload1", "Upload Parent1 (Format: Parent1_chr**.fa.gz)", multiple = FALSE)), ## 08/25/22 
column(12,fileInput("upload2", "Upload Parent2 (Format: Parent2_chr**.fa.gz)", multiple = FALSE)), ## 08/25/22 
###############

column(12,textInput("Upstream","Upstream (bp)",value=0)),
column(12,textInput("Downstream","Downstream (bp)",value=0)),

column(12,
pickerInput(
  inputId = "id", 
  label = "Genomes (Defalt: all genomes selected) :", 

 # choices = c('jagger','julius','lancer','landmark','mace','mattis','norin61','stanley','IWGSC','arinalrfor','spelta','query'),
   choices = default_choice,
 # selected = c('jagger','julius','lancer','landmark','mace','mattis','norin61','stanley','IWGSC','arinalrfor','spelta','query'), ## by default
   selected = default_choice,
  options = list(
    'actions-box' = TRUE, 
    size = 50,
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


uiOutput("Largefile"),

uiOutput("clustertree"),

uiOutput("Haplotypes"),

#actionButton("submit2", label = "Zoom In",class = "btn-warning"),

#actionButton("Haplotypes", "Plot selected haplotypes",style = "background-color:#CCCCFF"),

#uiOutput("GeneCapture"),

#actionButton("GeneCapture", label = "Involved Genes",style="color: #00F; background-color: #FFFF66; border-color: #c34113; border-radius: 10px; border-width: 2px"),
#bsTooltip("GeneCapture", "Show captured genes","right", options = NULL),

#  downloadButton('Repeat','Repeat',style = "background-color:#00FFFF"),
#  downloadButton('Variation','Variation',style = "background-color:#00FFFF"),
#  downloadButton('Haplotype_fa','Haplotype_fa',style = "background-color:#00FFFF"),
#  downloadButton('CDS_fa','CDS_fa',style = "background-color:#00FFFF"),

#actionButton("Largefile", label = "User's chromosome",style="color: FF66B2; background-color: #FFFF66; border-color: #c34113; border-radius: 10px; border-width: 2px"),
#bsTooltip("Largefile", "For large file upload only","right", options = NULL),

#actionButton("done", label = "Refresh",style = "background-color:#FF6666"),
#bsTooltip("done", "Clean all processed files for your submitted job!","right", options = NULL),

uiOutput("bucket"),
uiOutput("submit2"),

#uiOutput("save1"),
#downloadButton('save1','Fig1',style = "background-color:#C9DD03"),

uiOutput("extract_fa"),

uiOutput("Download"),
downloadButton('Save',label = "Save compressed results to ...",style = "background-color:#FFFFFF"),

uiOutput("done"),

), # sidebarPanel


mainPanel(
fluidRow(
  
############### 10/17/22, IP test
# tags$head(
#    tags$script(src="getIP.js")
#  ),
 verbatimTextOutput('IP'),       
############### 10/17/22, IP test

     
     column(12,verbatimTextOutput("visits")), ## number of visits ?

     column(12,verbatimTextOutput("infogene")), ## number of jobs submitted

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

     textOutput('upload1_test'),
     tags$head(tags$style("#upload1_test{color: green;
                                 font-size: 16px;
                                 font-style: italic;
                                 }"
                         )
     ),


     textOutput('upload2_test'),
     tags$head(tags$style("#upload2_test{color: green;
                                 font-size: 16px;
                                 font-style: italic;
                                 }"
                         )
     ),


#### progress in middle ??
tags$head(
    tags$style(
      HTML(".shiny-notification {
              height: 100px;
              width: 800px;
              position:fixed;
              top: calc(50% - 50px);;
              left: calc(50% - 400px);;
            }
           "
      )
    )
  ),
#### progress in middle ??


    # column(12,verbatimTextOutput("coordinate_info")),
    # span(textOutput("text"), style="color:green"),
    # span(textOutput("text2"), style="color:red"),
    # span(textOutput("text3"), style="color:black"),
    # span(textOutput("text4"), style="color:blue"),

     column(6, plotOutput("plot",click = NULL,dblclick = NULL,width = "100%",height = 'auto')),

     column(6, plotOutput("plot2",click = "plot2_click",dblclick = NULL,width = "100%",height = 'auto')),
     
     column(12,verbatimTextOutput("info2")),

     column(12,verbatimTextOutput("info4")),

     column(12,verbatimTextOutput("info3")),


     column(12,verbatimTextOutput("info_TE")),
     column(12, plotOutput("plot3",click = "plot3_click",dblclick = "plot3_dblclick",hover = "plot3_hover",width = "100%",height = 'auto')),


     column(12,verbatimTextOutput("info")),
     column(12, plotOutput("plot4",click = NULL,dblclick = NULL,width = "100%")),
      

     column(12, offset = 0,DT::dataTableOutput("table1"),style='padding-top:5px; padding-bottom:5px'), ## 06/17, cluster information
      
     column(12, offset = 0,DT::dataTableOutput("table2"),style='padding-top:5px; padding-bottom:5px'),

     column(12, offset = 0,DT::dataTableOutput("table3"),style='padding-top:5px; padding-bottom:5px'),

#    dataTableOutput("table3"),


#column(12, offset = 0,DT::dataTableOutput("table4"),style='padding-top:5px; padding-bottom:5px'),

#textOutput('Captured_Gene'),
#     tags$head(tags$style("#Captured_Gene{color: black;
#                                 font-size: 20px;
#                                 font-style: italic;
#                                 }"
#                         )
#     ),


# column(12, offset = 0,DT::dataTableOutput("tablecluster"),style='padding-top:5px; padding-bottom:5px')



)
) # mainPanel



) # sidebarLayout





      ) # For tagList

    }

