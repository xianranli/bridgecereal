############ Page 0

page_0 <- function(){

  page(
    href = "/page0",


    ui <-  function(request){
      tagList(

        h1("Wellcome! This is TestShiny app main page",style="text-align:center"),
        h2("Subheading?"),
        h3("Subheading?"),
   

        nav_links,

# To add ui part for page0


mainPanel(width = 6, helpText("Add some instructions here?"), img(src = 'Wheat.PNG',height="40%", width="40%", align="right")),

verbatimTextOutput("info")


      ) # For tagList
    }, # For ui function of page_0
    

# To add server function part for page0

server <- function(input, output, session){



sometext<-c("References to be added here ???")

output$info <- renderText({ sometext })


    } # server function of Page_0


  ) # page for Page_0

} # Page_0 function