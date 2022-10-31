
Page_species_info <- function(Speciesx,database_folder,gff_folder,script_folder,User_folder){

if(Speciesx=='Wheat') {

########################################
page_title<-c("This is Wheat Page!");
page_subtitle<-c("Wheat Alignment output");

target_folder0<-paste(database_folder,Speciesx,'/',sep='');
default_choice <-list.files(target_folder0);
Genome_choice <- c('',default_choice);
chromosome_choice <- c('','chr1A','chr1B','chr1D','chr2A','chr2B','chr2D','chr3A','chr3B','chr3D','chr4A','chr4B','chr4D','chr5A','chr5B','chr5D','chr6A','chr6B','chr6D','chr7A','chr7B','chr7D');
jobs_folder <- paste(database_folder,Speciesx,'/','IWGSC','/','Candidate_genes',sep='');
AllGenomes_GeneName <- paste(gff_folder,Speciesx,'/','Wheat_AllGenomes_GeneName.txt',sep='');
gpattern <- "Triticum.*";
gff_folder_Species <- paste(gff_folder,Speciesx,'/',sep='');
Backup_folder<-paste(database_folder,Speciesx,'/',"IWGSC",'/','Candidate_genes','/',sep=''); ## ?? 
perlArg0_db_sp <- paste(database_folder,Speciesx,'/',sep='')
########################################

}

my_list <- list(page_title,page_subtitle,target_folder0, default_choice, Genome_choice,chromosome_choice, jobs_folder,AllGenomes_GeneName,gpattern, gff_folder_Species, Backup_folder, perlArg0_db_sp)
return(my_list)

}
