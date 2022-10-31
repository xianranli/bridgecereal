
############ Speciesx
options(shiny.maxRequestSize=300*1024^2) ## Max size of uploaded file (300Mb in this case)


Species <- function(Speciesx,database_folder,gff_folder,script_folder,User_folder){

page_number <- gsub(' ', '/', paste(' ',Speciesx,sep=''))

########################################

Page_species_info_list <- Page_species_info(Speciesx,database_folder,gff_folder,script_folder,User_folder)

page_title <- Page_species_info_list[[1]]
page_subtitle<-Page_species_info_list[[2]]
target_folder0<-Page_species_info_list[[3]]
default_choice <-Page_species_info_list[[4]]
Genome_choice <- Page_species_info_list[[5]]
chromosome_choice <- Page_species_info_list[[6]]
jobs_folder <- Page_species_info_list[[7]]
AllGenomes_GeneName <- Page_species_info_list[[8]]
gpattern <- Page_species_info_list[[9]]
gff_folder_Species <- Page_species_info_list[[10]]
Backup_folder<-Page_species_info_list[[11]]
perlArg0_db_sp <- Page_species_info_list[[12]]
########################################




  page(
#    href = "/page10",

href = page_number ,


ui <- function(request){

   source( paste(script_folder,'tagfunction.R',sep=''), local = TRUE);

      tagfunction(page_title,page_subtitle,Genome_choice,chromosome_choice, default_choice)


    }, # For ui function of page_10



# To add server function part for page10

server <- function(input, output, session){


############### 10/17/22, IP test
#  IP <- reactive({ input$getIP })
#  output$IP <- renderPrint( cat(capture.output(str(IP()), split=TRUE)) )

#  output$IP <- renderPrint( fromJSON(readLines("http://api.hostip.info/get_json.php", warn=F))$ip )
  

############### 10/17/22, IP test


source( paste(script_folder,'Pre_run.R',sep=''), local = TRUE)

Pre_run(jobs_folder,AllGenomes_GeneName,default_choice,gff_folder,gff_folder_Species,gpattern,User_folder)

###### Start submit function! 

observeEvent(input$submit,{


########### with progress ...
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    progress$set(message = 'In progress ...',
                 detail = 'This may take a little while...')
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.1)  ## ??
    }
####################################################
#if(input$Gene==''& input$Chr==''& input$Upstream==''& input$Downstream==''& input$id==''& input$fasta==''){
#  refresh()
#}

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



if (input$fasta=='' ) {


########## 10/20/22 Test CDS
if (input$Pickformat=='CDS'){

Genome<-input$Pickgenome
chromosome<-input$Chr
Selected_gene<-input$Gene ## "TraesJAG1A03G00001180"
query <- "query"


G_gff_pattern <- paste(Genome,".*","_CDS_Working.gff3",sep="") ##
file0<-list.files(gff_folder_Species,pattern = G_gff_pattern);
file1<-read.table(paste(gff_folder_Species,file0,sep=''),header=F)

#file1<-read.table('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/gff/Wheat/Triticum_aestivum.IWGSC.53_CDS_Working.gff3',header=F)


file2<-file1[which(file1[,9]==Selected_gene),c(1,4,5)]

file3<-paste(Users_folder,'/',"positions.txt",sep='')

write.table(file2,file3,row.names=FALSE,col.names=FALSE,quote = FALSE,sep="\t")

#### 10/12/22, strand direction

strand_direction <- unique(file1[which(file1[,9]==Selected_gene),7])

#### 10/12/22, strand direction

setwd(Users_folder)

target_folder1<-paste(target_folder0,Genome,"/",sep = "")
samtools0 <- paste(target_folder1,Genome,"_",chromosome,".fa.gz",sep = "");
samtools1<- paste('perl ', script_folder, 'positions.pl ' ,samtools0 , sep='');
system(samtools1);
system("rm positions.txt");

perlArg4_Users_folder <-Users_folder;


query_fa <- "query.fasta";
dna <- readDNAStringSet("query.fasta");

Gene <- input$Gene  ## from shiny input
cds_ids <- input$Gene ##  ?? from shiny input, for plotting! 

query_length<-length(readDNAStringSet("query.fasta")[[1]]); # target CDS size


system("cat query.fasta > query_COPY.fasta");
system("cat query.fasta > query_COPY2.fasta");

Name_update_1 <- paste(Gene,"_mRNA",sep = "");
Name_update_2 <- paste(Gene,"_CDS",sep = "");

replace_pattern_in(file_contents=query,replace=Name_update_1,file_pattern="query.fasta",file_contents_perl = FALSE);
replace_pattern_in(file_contents=query,replace=Name_update_2,file_pattern="query_COPY.fasta",file_contents_perl = FALSE);

replace_pattern_in(file_contents=query,replace=input$Chr,file_pattern="query_COPY2.fasta",file_contents_perl = FALSE);  ## ?? query_COPY2.fasta

system("cat query.fasta query_COPY.fasta > query_Working.fasta");

Name_update_3 <- paste(Gene,"_ref",sep = "");
file.rename("query_Working.fasta", Name_update_3);

system("rm query.fasta query_COPY.fasta query_COPY2.fasta");  ## ?? query_COPY2.fasta

perlArg1_PickGenome <- input$Pickgenome;
perlArg2_PickGene <- input$Gene;
perlArg3_PickChr <- input$Chr;


perlArg5_PickUp <- input$Upstream;
perlArg6_PickDown <- input$Downstream;



#perlArg10 <- 0;

}
########## 10/20/22 Test CDS

##############
##############

source(paste(administrator_path,"script/cerealBRIDGE_output.R",sep=''), local = TRUE);
cerealBRIDGE_output(User_folder0,perlArg0_db_sp,perlArg1_PickGenome ,perlArg2_PickGene,perlArg3_PickChr,perlArg4_Users_folder,perlArg5_PickUp,perlArg6_PickDown, Backup_folder,strand_direction, database_folder,gff_folder,script_folder,User_folder)

########
########
########
########
########


} else if (input$fasta!='') {   #### input$fasta not null 08/09/22

###
Genome<-input$Pickgenome
chromosome<-input$Chr
Selected_gene<-input$Gene ## "TraesJAG1A03G00001180"
query <- "query";
G_gff_pattern <- paste(Genome,".*","_CDS_Working.gff3",sep="") ##
file0<-list.files(gff_folder_Species,pattern = G_gff_pattern);
file1<-read.table(paste(gff_folder_Species,file0,sep=''),header=F)
file2<-file1[which(file1[,9]==Selected_gene),c(1,4,5)]
file3<-paste(Users_folder,'/',"positions.txt",sep='')
write.table(file2,file3,row.names=FALSE,col.names=FALSE,quote = FALSE,sep="\t")
#### 10/12/22, strand direction
strand_direction <- unique(file1[which(file1[,9]==Selected_gene),7])
#### 10/12/22, strand direction
setwd(Users_folder)
target_folder1<-paste(target_folder0,Genome,"/",sep = "")
samtools0 <- paste(target_folder1,Genome,"_",chromosome,".fa.gz",sep = "");
samtools1<- paste('perl ', script_folder, 'positions.pl ' ,samtools0 , sep='');
system(samtools1);
system("rm positions.txt");
perlArg4_Users_folder <-Users_folder;
query_fa <- "query.fasta";
dna <- readDNAStringSet("query.fasta");
Gene <- input$Gene  ## from shiny input
cds_ids <- input$Gene ##  ?? from shiny input, for plotting! 
query_length<-length(readDNAStringSet("query.fasta")[[1]]); # target CDS size
system("cat query.fasta > query_COPY.fasta");
system("cat query.fasta > query_COPY2.fasta");
Name_update_1 <- paste(Gene,"_mRNA",sep = "");
Name_update_2 <- paste(Gene,"_CDS",sep = "");
replace_pattern_in(file_contents=query,replace=Name_update_1,file_pattern="query.fasta",file_contents_perl = FALSE);
replace_pattern_in(file_contents=query,replace=Name_update_2,file_pattern="query_COPY.fasta",file_contents_perl = FALSE);

replace_pattern_in(file_contents=query,replace=input$Chr,file_pattern="query_COPY2.fasta",file_contents_perl = FALSE);  ## ?? query_COPY2.fasta

system("cat query.fasta query_COPY.fasta > query_Working.fasta");

Name_update_3 <- paste(Gene,"_ref",sep = "");
file.rename("query_Working.fasta", Name_update_3);

system("rm query.fasta query_COPY.fasta query_COPY2.fasta");  ## ?? query_COPY2.fasta
###

###
query0<-input$fasta
tmp <- tempfile(fileext = ".fa")
   if (startsWith(query0, ">")){
     writeLines(query0, tmp)
   } else {
     writeLines(paste0(">Query\n",query0), tmp)
   }
dna <- readDNAStringSet(tmp)
query_fa <- "Pasted_query.fasta"; # paste(Users_folder,'/',"query.fasta",sep='')
writeXStringSet(dna, query_fa, append=FALSE, compress=FALSE, format="fasta");
##########
########## 08/09/22

#writeXStringSet(dna, query_fa, append=FALSE, compress=FALSE, format="fasta");
Name_update_4 <- paste("query","_",input$Chr,".fa",sep = "");
file.rename("Pasted_query.fasta", Name_update_4);
replace_pattern_in(file_contents=query,replace=input$Chr,file_pattern=Name_update_4,file_contents_perl = FALSE);
### ??? 10/27/22 for query ...


   makedb0 <- paste('makeblastdb -in', Name_update_4 ,'-dbtype nucl',sep=' ') 
   makedb1 <- system(makedb0)
          
   bgzip0 <- paste('bgzip -@ 2 ',Name_update_4,sep=' ')
   bgzip1 <- system(bgzip0)

   samtools0 <- paste('samtools faidx',paste(Name_update_4,'.gz',sep='') ,sep=' ')
   samtools1 <- system(samtools0)

   dir0<- paste(Users_folder,'/','query',sep='')
   dir.create(dir0)

   move0<-list.files(Users_folder, pattern='query_');

   for(i in 1:length(move0)){
    
         move1 <- paste(Users_folder,'/',move0[i],sep='')
         move2<- paste('mv',move1,dir0,sep=' ')
         system(move2)

         }


perlArg1_PickGenome <- input$Pickgenome;
perlArg2_PickGene <- input$Gene;
perlArg3_PickChr <- input$Chr;
perlArg4_Users_folder <-Users_folder;


perlArg5_PickUp <- 0;
perlArg6_PickDown <- 0;

##########################################
##########################################

source(paste(administrator_path,"script/cerealBRIDGE_output.R",sep=''), local = TRUE);
cerealBRIDGE_output(User_folder0,perlArg0_db_sp,perlArg1_PickGenome ,perlArg2_PickGene,perlArg3_PickChr,perlArg4_Users_folder,perlArg5_PickUp,perlArg6_PickDown, Backup_folder,strand_direction, database_folder,gff_folder,script_folder,User_folder)

}                       ## input$fasta not null 08/09/22;


}) ## observeEvent submit !!

#################### To test large file upload 08/25/22
####################


observeEvent(input$Largefile,{

########### with progress ...
    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    progress$set(message = 'In progress ...',
                 detail = 'This may take a little while...')
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.2)  ## ??
    }
####################################################
ip_address <- gsub( '\\.', '_', fromJSON(readLines("http://api.hostip.info/get_json.php", warn=F))$ip )
User_Gene<-input$Gene
User_folder0<- paste(ip_address,'_',User_Gene,sep='')

Users_folder<-paste(User_folder, User_folder0 , sep='')  ## User's ip_gene !!
####################################################

Genome<-input$Pickgenome
chromosome<-input$Chr
Selected_gene<-input$Gene ## "TraesJAG1A03G00001180"
query <- "query"

G_gff_pattern <- paste(Genome,".*","_CDS_Working.gff3",sep="") ##
file0<-list.files(gff_folder_Species,pattern = G_gff_pattern);
file1<-read.table(paste(gff_folder_Species,file0,sep=''),header=F)

#G_gff_pattern <- paste(Genome,".*","_CDS_Working.gff3",sep="") ##
#file0<-list.files(gff_folder_Species,pattern = G_gff_pattern);
#file1<-read.table(file0,header=F)
#file1<-read.table('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/gff/Wheat/Triticum_aestivum.IWGSC.53_CDS_Working.gff3',header=F)

file2<-file1[which(file1[,9]==Selected_gene),c(1,4,5)]

file3<-paste(Users_folder,'/',"positions.txt",sep='')

write.table(file2,file3,row.names=FALSE,col.names=FALSE,quote = FALSE,sep="\t")

#### 10/12/22, strand direction

strand_direction <- unique(file1[which(file1[,9]==Selected_gene),7])

#### 10/12/22, strand direction

setwd(Users_folder)

target_folder1<-paste(target_folder0,Genome,"/",sep = "")
samtools0 <- paste(target_folder1,Genome,"_",chromosome,".fa.gz",sep = "");
samtools1<- paste('perl ', script_folder, 'positions.pl ' ,samtools0 , sep='');
system(samtools1);
system("rm positions.txt");


perlArg4_Users_folder <-Users_folder;


query_fa <- "query.fasta";
dna <- readDNAStringSet("query.fasta");

Gene <- input$Gene  ## from shiny input
cds_ids <- input$Gene ##  ?? from shiny input, for plotting! 

query_length<-length(readDNAStringSet("query.fasta")[[1]]); # target CDS size


system("cat query.fasta > query_COPY.fasta");
system("cat query.fasta > query_COPY2.fasta");

Name_update_1 <- paste(Gene,"_mRNA",sep = "");
Name_update_2 <- paste(Gene,"_CDS",sep = "");

replace_pattern_in(file_contents=query,replace=Name_update_1,file_pattern="query.fasta",file_contents_perl = FALSE);
replace_pattern_in(file_contents=query,replace=Name_update_2,file_pattern="query_COPY.fasta",file_contents_perl = FALSE);

replace_pattern_in(file_contents=query,replace=input$Chr,file_pattern="query_COPY2.fasta",file_contents_perl = FALSE);  ## ?? query_COPY2.fasta

system("cat query.fasta query_COPY.fasta > query_Working.fasta");

Name_update_3 <- paste(Gene,"_ref",sep = "");
file.rename("query_Working.fasta", Name_update_3);

system("rm query.fasta query_COPY.fasta query_COPY2.fasta");  ## ?? query_COPY2.fasta

perlArg1_PickGenome <- input$Pickgenome;
perlArg2_PickGene <- input$Gene;
perlArg3_PickChr <- input$Chr;


perlArg5_PickUp <- input$Upstream;
perlArg6_PickDown <- input$Downstream;

#perlArg10 <- 1;

####################

source(paste(administrator_path,"script/cerealBRIDGE_output.R",sep=''), local = TRUE);
cerealBRIDGE_output(User_folder0,perlArg0_db_sp,perlArg1_PickGenome ,perlArg2_PickGene,perlArg3_PickChr,perlArg4_Users_folder,perlArg5_PickUp,perlArg6_PickDown, Backup_folder,strand_direction, database_folder,gff_folder,script_folder,User_folder)

}) # observeEvent input$Largefile

#################### To test large file upload 08/25/22, The End
####################


observeEvent(input$Download,{

                progress <- Progress$new(session, min=1, max=10)
                on.exit(progress$close())
                progress$set(message = 'Prepare your .zip file for downloading ...',
                detail = 'Almost done...')
                for (i in 1:5) {
                progress$set(value = i)
                Sys.sleep(0.3)
                }

ip_address <- gsub( '\\.', '_', fromJSON(readLines("http://api.hostip.info/get_json.php", warn=F))$ip )
User_Gene<-input$Gene
User_folder0<- paste(ip_address,'_',User_Gene,sep='')
Users_folder<-paste(User_folder, User_folder0 , sep='')  ## User's ip_gene !!

## Parent1 or Parent2
if(file.exists(paste(Users_folder,'/','Parent1',sep=''))){

remove_Parent1 <- paste('rm -r ',Users_folder,'/','Parent1',sep='')
system(remove_Parent1)

}
if(file.exists(paste(Users_folder,'/','Parent2',sep=''))){

remove_Parent2 <- paste('rm -r ',Users_folder,'/','Parent2',sep='')
system(remove_Parent2)

}


zip(zipfile=paste(Users_folder,'.zip',sep=''), files=Users_folder)

compress_result0 <- paste(Users_folder,'.zip',sep='')

output$Save <- downloadHandler(
  filename = function() {
   file = paste(User_Gene,'.zip',sep='')
  },
  content = function(file) {
  file.copy(paste(Users_folder,'.zip',sep=''),file)
 }
)


})


####
observeEvent(input$done,{

ip_address <- gsub( '\\.', '_', fromJSON(readLines("http://api.hostip.info/get_json.php", warn=F))$ip )
User_Gene<-input$Gene
User_folder0<- paste(ip_address,'_',User_Gene,sep='')
Users_folder<-paste(User_folder, User_folder0 , sep='')  ## User's ip_gene !!
Gene <- input$Gene  ## from shiny input

if(Gene!=''){

system_clean<-paste('rm -r ',Users_folder,sep=' ')
system(system_clean)

refresh()

}else{refresh()}

}) ## observeEvent DONE, To remove all temp files






    } # server function of Page_10



  ) # page for Page_10
} # Page_10 function
