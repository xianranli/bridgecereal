
Pre_run <- function(jobs_folder,AllGenomes_GeneName,default_choice,gff_folder,gff_folder_Species,gpattern,User_folder) {

# hide button Save for now
observe({ toggle(id="Save", condition=!is.null(input$Download))})



## To track login times?
num_visit0 <- as.numeric(read.table("/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny/Visit_Times.txt",header=F));
num_visit1 <- paste("There have been",num_visit0,"visitors!",sep=" ");
output$visits <- renderText({ num_visit1 });
num_visitNew<- num_visit0+1;
write.table(num_visitNew,"/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny/Visit_Times.txt",row.names=FALSE,col.names=FALSE,quote = FALSE);
## To track login times?



## To count submitted jobs (genes)
num_files0 <- length(list.files(jobs_folder)); # number of jobs submitted !

num_files1<-paste("There're",num_files0,"genes or jobs have been submitted on this page! You can try yours.",sep=' ');
output$infogene <- renderText({ num_files1 });
##



#### To test existence of gene name ??
observeEvent(input$Gene, {

    if(input$Gene != ""){

    Test_GeneName0 <- read.table(AllGenomes_GeneName,header=F)
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


   Genome2<-paste(gpattern,input$Pickgenome,".*","_gene_Working.gff3",sep="") ##
   file0<-list.files(path = gff_folder_Species,pattern = Genome2);
   Test_GeneName0 <- read.table(paste(gff_folder_Species,file0,sep=''),header=F)   

    
  # Test_GeneName0 <- read.table('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/gff/Wheat/Triticum_aestivum.IWGSC.53_gene_Working.gff3',header=F)

    Test_GeneName1 <- Test_GeneName0[which(Test_GeneName0[,9]==input$Gene),c(1,4,5,7)]

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

    if (input$Pickformat == "fasta_seq"){
     
      shinyjs::enable(id = "fasta")

    } else if (input$Pickformat == "CDS" ) {
      
      shinyjs::disable(id = "fasta")
      
    } 

    })

###### To test fasta format
observeEvent(input$fasta, {

    if(input$fasta != ""){
     
     query0<-input$fasta

   if (startsWith(query0, ">query")){

    Test_fasta2 <- paste("Correct fasta input. You can proceed ...",sep='');
    output$fasta_test <- renderText({ Test_fasta2 });

              updatePickerInput(session = session, inputId = "id",
                        choices = c(default_choice,'query'),
                        selected = c(default_choice,'query'))
     
   } else {
    Test_fasta2 <- paste("Incorrect fasta format! Please double check (>query as first line name !!!)",sep='');
    output$fasta_test <- renderText({ Test_fasta2 });
   }

    } 
  
  })
###### To test fasta format


######
#observeEvent(input$Pick_Parents, {

#    if (input$Pick_Parents == "" ) {

#      shinyjs::disable(id = "upload1")
#      shinyjs::disable(id = "upload2")
      
      
#    } 

#    if (input$Pick_Parents != "" ){
      
#      shinyjs::enable(id = "upload1")
#      shinyjs::enable(id = "upload2")
      

#    }

#    })
######


observeEvent(input$upload1, {

                progress <- Progress$new(session, min=1, max=10)
                on.exit(progress$close())
                progress$set(message = 'Processing your compressed file Parent1 ...',
                detail = 'Almost there...')
                for (i in 1:5) {
                progress$set(value = i)
                Sys.sleep(0.3)
                }


          if (is.null(input$upload1)) return()

#############
ip_address <- gsub( '\\.', '_', fromJSON(readLines("http://api.hostip.info/get_json.php", warn=F))$ip )
User_Gene<-input$Gene
User_folder0<- paste(ip_address,'_',User_Gene,sep='')

if( file.exists( paste(User_folder, User_folder0 , sep='')) ){

remove_exist_file <- paste('rm -r ',paste(User_folder, User_folder0 , sep=''),sep='')
system(remove_exist_file)

}


Users_folder1<-paste('mkdir ', User_folder , User_folder0 , sep='')
system(Users_folder1)  ## 

Users_folder<-paste(User_folder, User_folder0 , sep='')  ## User's ip_gene !!
#############
       
          file.copy(input$upload1$datapath, paste0(Users_folder,'/',input$upload1$name) )
          
          fileName0<- paste(Users_folder,'/',input$upload1$name,sep='')
          fileName1 <- as.character(strsplit(fileName0, ".fa.gz"))
           
          uncompressed0 <- paste("gzip -d",fileName0, sep=' ') 
          uncompressed1 <- system(uncompressed0)

          fileName2 <- paste(fileName1,'.fa',sep='')
          fileName3 <- paste(fileName1,'.fa.gz',sep='')

          makedb0 <- paste('makeblastdb -in', fileName2 ,'-dbtype nucl',sep=' ') 
          makedb1 <- system(makedb0)
          
          bgzip0 <- paste('bgzip -@ 8 ',fileName2,sep=' ')
          bgzip1 <- system(bgzip0)

          samtools0 <- paste('samtools faidx',fileName3 ,sep=' ')
          samtools1 <- system(samtools0)
          

          ParentName0<- as.character(strsplit(input$upload1$name, ".fa.gz"))
          ParentName1<- gsub('_.*','',ParentName0)

         dir0<- paste(Users_folder,'/','Parent1',sep='')
         dir.create(dir0)
         

         move0<-list.files(Users_folder, pattern='Parent1_')

         for(i in 1:length(move0)){
        
         move1 <- paste(Users_folder,'/',move0[i],sep='')
         move2<- paste('mv',move1,dir0,sep=' ')
         system(move2)

         }


          updatePickerInput(session = session, inputId = "id",
                        choices = c(default_choice,'Parent1'),
                        selected = c(default_choice,'Parent1'))

    Test_upload1 <- paste("Completion of Parent1 upload. You can proceed ...",sep='')
    output$upload1_test <- renderText({ Test_upload1 })

shinyjs::disable(id = "submit")

output$Largefile <- renderUI({

    actionButton("Largefile", label = "Submit",style="color: FF66B2; background-color: #FFFF99; border-color: #c34113; border-radius: 10px; border-width: 2px")

  })               




        } ) ## observe
######

observeEvent(input$upload2, {

                progress <- Progress$new(session, min=1, max=10)
                on.exit(progress$close())
                progress$set(message = 'Processing your compressed file Parent2 ...',
                detail = 'Almost there...')
                for (i in 1:5) {
                progress$set(value = i)
                Sys.sleep(0.3)
                }


          if (is.null(input$upload2)) return()


#############
ip_address <- gsub( '\\.', '_', fromJSON(readLines("http://api.hostip.info/get_json.php", warn=F))$ip )
User_Gene<-input$Gene
User_folder0<- paste(ip_address,'_',User_Gene,sep='')
Users_folder<-paste(User_folder, User_folder0 , sep='')  ## User's ip_gene !!
#############
#############
       
          file.copy(input$upload2$datapath, paste0(Users_folder,'/',input$upload2$name) )
          
          fileName0<- paste(Users_folder,'/',input$upload2$name,sep='')
          fileName1 <- as.character(strsplit(fileName0, ".fa.gz"))
           
          uncompressed0 <- paste("gzip -d",fileName0, sep=' ') 
          uncompressed1 <- system(uncompressed0)

          fileName2 <- paste(fileName1,'.fa',sep='')
          fileName3 <- paste(fileName1,'.fa.gz',sep='')

          makedb0 <- paste('makeblastdb -in', fileName2 ,'-dbtype nucl',sep=' ') 
          makedb1 <- system(makedb0)
          
          bgzip0 <- paste('bgzip -@ 8 ',fileName2,sep=' ')
          bgzip1 <- system(bgzip0)

          samtools0 <- paste('samtools faidx',fileName3 ,sep=' ')
          samtools1 <- system(samtools0)
          

          ParentName0<- as.character(strsplit(input$upload2$name, ".fa.gz"))
          ParentName2<- gsub('_.*','',ParentName0)

         dir0<- paste(Users_folder,'/','Parent2',sep='')
         dir.create(dir0)
         

         move0<-list.files(Users_folder, pattern='Parent2_')

         for(i in 1:length(move0)){
        
         move1 <- paste(Users_folder,'/',move0[i],sep='')
         move2<- paste('mv',move1,dir0,sep=' ')
         system(move2)

         }


          updatePickerInput(session = session, inputId = "id",
                        choices = c(default_choice,'Parent1','Parent2'),
                        selected = c(default_choice,'Parent1','Parent2'))


        Test_upload2 <- paste("Completion of Parent2 upload. You can proceed ...",sep='')
        output$upload2_test <- renderText({ Test_upload2 })

    shinyjs::disable(id = "submit")

    output$Largefile <- renderUI({

    actionButton("Largefile", label = "Submit (For user uploaded haplotype)",style="color: FF66B2; background-color: #FFFF99; border-color: #c34113; border-radius: 10px; border-width: 2px")

  })    

        })


}   # function


