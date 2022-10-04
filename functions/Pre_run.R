
Pre_run <- function() {


## To track login times?
num_visit0 <- as.numeric(read.table("/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny/Visit_Times.txt",header=F));
num_visit1 <- paste("There have been",num_visit0,"visitors on this Wheat page!",sep=" ");
output$visits <- renderText({ num_visit1 });
num_visitNew<- num_visit0+1;
write.table(num_visitNew,"/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny/Visit_Times.txt",row.names=FALSE,col.names=FALSE,quote = FALSE);
## To track login times?



## To count submitted jobs (genes)
num_files0 <- length(list.files("/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/IWGSC/Candidate_genes")); # number of jobs submitted !
num_files1<-paste("There're",num_files0,"genes or jobs have been submitted in Wheat! You can try yours.",sep=' ');
output$info <- renderText({ num_files1 });
##



#### To test existence of gene name ??
  observeEvent(input$Gene, {

    if(input$Gene != ""){

    Test_GeneName0 <- read.table("/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny/Wheat_AllGenomes_GeneName.txt",header=F)
    Test_GeneName1 <- Test_GeneName0[which(Test_GeneName0[,1]==input$Gene),]

    if(length(Test_GeneName1)>=1){
    Test_GeneName2 <- paste("Our database contains gene ",Test_GeneName1,". You can proceed ...",sep='');
    output$query_test <- renderText({ Test_GeneName2 });
    } else if (length(Test_GeneName1)==0){
    Test_GeneName2 <- paste("Wrong gene ID (Please double check! NO space left!). Have you provided your own job ID for fasta sequence or preferred coordinates? You may proceed in this scenario. ",sep=' ');
    output$query_test <- renderText({ Test_GeneName2 });
     }
    
    } 
  

  })
#### To test existence of gene name ??



#### To show gene's coordinates ?? 08/16/22
  observeEvent(input$Chr, {

    if(input$Gene != "" & input$Pickgenome != "" & input$Chr != ""){

    current.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny"
    Genome <- input$Pickgenome

    Genome2<-paste("Triticum.*",Genome,".*","_gene_Working.gff3",sep="") ##
    file0<-list.files(current.folder,pattern = Genome2);
 
    Test_GeneName0 <- read.table(file0,header=F)
    Test_GeneName1 <- Test_GeneName0[which(Test_GeneName0[,9]==input$Gene),c(1,4,5)]

    if(nrow(Test_GeneName1)==1){
    #Test_GeneName2 <- paste0("The query gene is located at ",Test_GeneName1,sep='');
    Test_GeneName2 <- paste(c("The query gene is located at: ", Test_GeneName1), collapse= " ")
    output$coordinates_test <- renderText({ Test_GeneName2 });

    } else if (nrow(Test_GeneName1)==0){
  
    Test_GeneName2 <- paste("No coordinate information for this query gene!",sep='');
    output$coordinates_test <- renderText({ Test_GeneName2 });
     }
    
    } 
  

  })
#### To show gene's coordinates ?? 08/16/22




## enable or disable some of input selections
  observeEvent(input$Pickformat, {
    if(input$Pickformat == "coordinates"){
      shinyjs::enable(id = "Coordinates")


      ### To intersect provided coordinates
    observeEvent(input$Coordinates, {
    current.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny"
    Genome <- input$Pickgenome
    Chromo1<- input$Chr
    Coor1<- input$Coordinates[1]
    Coor2<- input$Coordinates[2]
    Name <- input$Gene
    Genome2<-paste("Triticum.*",Genome,".*","_gene_Working.gff3",sep="") ##
    file0<-list.files(current.folder,pattern = Genome2);
    Tempfile0<-paste(Chromo1,Coor1,Coor2,Name,sep=' ')
    write.table(Tempfile0,"Coordinate_Test.bed",row.names=FALSE,col.names=FALSE,quote = FALSE,sep=" ")
    system("perl -p -i -e 's/ /\t/g' Coordinate_Test.bed")
    system_bedtools0<- paste("bedtools intersect -a ", file0, " -b Coordinate_Test.bed | cut -f9 > Coordinate_Test1.bed ",sep = " ")
    system(system_bedtools0)
    if(file.size("Coordinate_Test1.bed")!=0){
    Test_GeneName0 <- read.table("Coordinate_Test1.bed",header=F)

    Test_GeneName2 <- paste(c("The query coordinates on",Chromo1,"in",Genome,"intersect with",nrow(Test_GeneName0), "annotated genes: ", Test_GeneName0[,1]), collapse = ' ')

    output$QueryCoordinate_test <- renderText({ Test_GeneName2 });
    } else if (file.size("Coordinate_Test1.bed")==0){
    Test_GeneName2 <- paste("The query coordinates intersect with nothing!",sep='');
    output$QueryCoordinate_test <- renderText({ Test_GeneName2 });
    }
    system("rm Coordinate_Test1.bed Coordinate_Test.bed")
    
    })
     ### To intersect provided coordinates
      

      shinyjs::disable(id = "fasta")
      shinyjs::disable(id = "upload")

    } else if (input$Pickformat == "fasta_seq"){
      shinyjs::enable(id = "fasta")

      shinyjs::disable(id = "Coordinates")
      shinyjs::disable(id = "upload")

    } else if (input$Pickformat == "CDS" || input$Pickformat == "gene" ) {
      shinyjs::disable(id = "Coordinates")
      shinyjs::disable(id = "fasta")
      shinyjs::disable(id = "upload")
    } else if (input$Pickformat == "file"){
      shinyjs::enable(id = "upload")

      shinyjs::disable(id = "fasta")
      shinyjs::disable(id = "Coordinates")

      
      shinyjs::disable(id = "Largefile")
      
        ## ip address ?


        observe({
            if (is.null(input$upload)) return()
            file.copy(input$upload$datapath, paste0("/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny/Users/",input$upload$name))

          userfolder0 <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny/Users/";
          
          setwd(userfolder0)
          userfolder01 <- as.character(strsplit(input$upload$name, ".fa.gz"))
          
          userfolder1 <- gsub('_', "", userfolder01, fixed = FALSE)

          dir.create(userfolder1)
          userfolder2<-paste(userfolder0,userfolder1,'/',userfolder1,sep='')
          dir.create(userfolder2)
          move0<-paste('mv',input$upload$name,userfolder2,sep=' ')
          system(move0)
          decompress0<-paste('gunzip ',userfolder2,'/',input$upload$name,sep='')
          system(decompress0)
           
          setwd(userfolder2)
          tempt0<-gsub('_', "", as.character(strsplit(input$upload$name, ".fa.gz")), fixed = FALSE)
          tempt1<-paste(tempt0,'_',input$Chr,'.fa',sep='')
          tempt2<-paste(userfolder01,'.fa',sep='')

          file.rename(tempt2,tempt1)

          blastindex0<-paste("makeblastdb -in ",userfolder2,'/', tempt1 ," -dbtype nucl",sep='')
          system(blastindex0)

          Shiny.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny";
          setwd(Shiny.folder);
          
          updatePickerInput(session = session, inputId = "id",
                        choices = c('jagger','julius','lancer','landmark','mace','mattis','norin61','stanley','IWGSC','arinalrfor','spelta','query',userfolder1),
                        selected = c('jagger','julius','lancer','landmark','mace','mattis','norin61','stanley','IWGSC','arinalrfor','spelta','query',userfolder1))

                

                progress <- Progress$new(session, min=1, max=15)
                on.exit(progress$close())
                progress$set(message = 'Processing your compressed file ...',
                detail = 'Almost there...')
                for (i in 1:10) {
                progress$set(value = i)
                Sys.sleep(0.5)
                }
                


                shinyjs::enable(id = "Largefile")




        }) ## observe



    }

  })
## enable or disable some of input selections

###### To test fasta format
  observeEvent(input$fasta, {

    if(input$fasta != ""){
     
     query0<-input$fasta

   if (startsWith(query0, ">query")){

    Test_fasta2 <- paste("Correct fasta input. You can proceed ...",sep='');
    output$fasta_test <- renderText({ Test_fasta2 });
     
   } else {
    Test_fasta2 <- paste("Incorrect fasta format! Please double check (>query as first line name !!!)",sep='');
    output$fasta_test <- renderText({ Test_fasta2 });
   }

    } 
  

  })
###### To test fasta format





}


