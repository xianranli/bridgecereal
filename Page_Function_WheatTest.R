
############ Page 6
options(shiny.maxRequestSize=300*1024^2) ## Max size of uploaded file (300Mb in this case)

page_6 <- function(){

  page(
    href = "/page6",


ui <- function(request){

      source('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny/functions/taglist.R',local = TRUE)
      taglist()


    }, # For ui function of page_6



# To add server function part for page6

server <- function(input, output, session){


source('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny/functions/Pre_run.R',local = TRUE)
Pre_run()


###### Start submit function! 

observeEvent(input$submit,{


########### with progress ...
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())

    progress$set(message = 'In progress',
                 detail = 'This may take a while...')
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.5)
    }
####################################################
#if(input$Gene==''& input$Chr==''& input$Upstream==''& input$Downstream==''& input$id==''& input$fasta==''){
#  refresh()
#}


if (input$fasta=='' ) {


########## 08/10/22 Test Gene
if(input$Pickformat=='gene') {

Genome<-input$Pickgenome
chromosome<-input$Chr
Selected_gene<-input$Gene ## "TraesJAG1A03G00001180"
query <- "query"


current.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny"

Genome2<-paste(Genome,".*","_gene_Working.gff3",sep="") ##

file0<-list.files(current.folder,pattern = Genome2);

file00<-paste(current.folder,"/",file0,sep="")

file1<-read.table(file0,header=F)

file2<-file1[which(file1[,9]==Selected_gene),]
write.table(file2,"target0.gff",row.names=FALSE,col.names=FALSE,quote = FALSE,sep="\t")

target_folder0<-"/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/"
target_folder1<-paste(target_folder0,Genome,"/",sep = "")
bedtools0 <- paste(target_folder1,Genome,"_",chromosome,".fa",sep = "");
system_bedtools <- paste("/usr/local/bin/bedtools getfasta -fi", bedtools0 ,"-bed target0.gff -nameOnly -fo query.fasta",sep = " ");
system(system_bedtools);

#replace_pattern_in(file_contents=gene,replace=query,file_pattern="query.fasta",file_contents_perl = FALSE);

system("perl -p -i -e 's/gene/query/g' query.fasta");

query_fa <- "query.fasta";
dna <- readDNAStringSet("query.fasta");
#####

Gene <- input$Gene  ## from shiny input
cds_ids <- input$Gene ## from shiny input, for plotting! 

query_length<-length(readDNAStringSet("query.fasta")[[1]]); # target CDS size

system("cat query.fasta > query_COPY.fasta");

system("cat query.fasta > query_COPY2.fasta");

Name_update_1 <- paste(Gene,"_mRNA",sep = "");
Name_update_2 <- paste(Gene,"_CDS",sep = "");

replace_pattern_in(file_contents=query,replace=Name_update_1,file_pattern="query.fasta",file_contents_perl = FALSE);
replace_pattern_in(file_contents=query,replace=Name_update_2,file_pattern="query_COPY.fasta",file_contents_perl = FALSE);

replace_pattern_in(file_contents=query,replace=input$Chr,file_pattern="query_COPY2.fasta",file_contents_perl = FALSE);

system("cat query.fasta query_COPY.fasta > query_Working.fasta");

Name_update_3 <- paste(Gene,"_ref",sep = "");
file.rename("query_Working.fasta", Name_update_3);

system("rm query.fasta query_COPY.fasta target0.gff");
writeXStringSet(dna, query_fa, append=FALSE, compress=FALSE, format="fasta");

Name_update_4 <- paste("query","_",input$Chr,".fa",sep = "");
file.rename("query_COPY2.fasta", Name_update_4);

system_blast <- paste("/usr/local/bin/makeblastdb -in", Name_update_4 ,"-dbtype nucl",sep = " ");

system(system_blast);

result_folder1<-c('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/');
result_folder2<-paste(result_folder1,"query",sep="");
current.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny"

result_files<-list.files(current.folder, pattern = Name_update_4)
file.copy(result_files, result_folder2)

}
########## 08/10/22 Test Gene

########## 08/11/22 Test CDS
if (input$Pickformat=='CDS'){

Genome<-input$Pickgenome
chromosome<-input$Chr
Selected_gene<-input$Gene ## "TraesJAG1A03G00001180"
query <- "query"


current.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny"

Genome2<-paste(Genome,".*","_CDS_Working.gff3",sep="") ##

file0<-list.files(current.folder,pattern = Genome2);

file00<-paste(current.folder,"/",file0,sep="")

file1<-read.table(file0,header=F)

file2<-file1[which(file1[,9]==Selected_gene),]
write.table(file2,"target0.gff",row.names=FALSE,col.names=FALSE,quote = FALSE,sep="\t")

target_folder0<-"/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/"
target_folder1<-paste(target_folder0,Genome,"/",sep = "")
bedtools0 <- paste(target_folder1,Genome,"_",chromosome,".fa",sep = "");
system_bedtools <- paste("/usr/local/bin/bedtools getfasta -fi", bedtools0 ,"-bed target0.gff -nameOnly -fo query0.fasta",sep = " ");
system(system_bedtools);

#replace_pattern_in(file_contents=gene,replace=query,file_pattern="query.fasta",file_contents_perl = FALSE);

system("bash extract_fasta.sh");

#system("perl -p -i -e 's/gene/query/g' query.fasta");

query_fa <- "query.fasta";
dna <- readDNAStringSet("query.fasta");
#####

Gene <- input$Gene  ## from shiny input
cds_ids <- input$Gene ## from shiny input, for plotting! 

query_length<-length(readDNAStringSet("query.fasta")[[1]]); # target CDS size

system("cat query.fasta > query_COPY.fasta");

system("cat query.fasta > query_COPY2.fasta");

Name_update_1 <- paste(Gene,"_mRNA",sep = "");
Name_update_2 <- paste(Gene,"_CDS",sep = "");

replace_pattern_in(file_contents=query,replace=Name_update_1,file_pattern="query.fasta",file_contents_perl = FALSE);
replace_pattern_in(file_contents=query,replace=Name_update_2,file_pattern="query_COPY.fasta",file_contents_perl = FALSE);

replace_pattern_in(file_contents=query,replace=input$Chr,file_pattern="query_COPY2.fasta",file_contents_perl = FALSE);

system("cat query.fasta query_COPY.fasta > query_Working.fasta");

Name_update_3 <- paste(Gene,"_ref",sep = "");
file.rename("query_Working.fasta", Name_update_3);

system("rm query.fasta query_COPY.fasta target0.gff");
writeXStringSet(dna, query_fa, append=FALSE, compress=FALSE, format="fasta");

Name_update_4 <- paste("query","_",input$Chr,".fa",sep = "");
file.rename("query_COPY2.fasta", Name_update_4);

system_blast <- paste("/usr/local/bin/makeblastdb -in", Name_update_4 ,"-dbtype nucl",sep = " ");

system(system_blast);

result_folder1<-c('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/');
result_folder2<-paste(result_folder1,"query",sep="");
current.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny"

result_files<-list.files(current.folder, pattern = Name_update_4)
file.copy(result_files, result_folder2)



}
########## 08/11/22 Test CDS

########## 08/11/22 Test coordinate
if(input$Pickformat=='coordinates'){

Genome<-input$Pickgenome
chromosome<-input$Chr
Selected_gene<-input$Gene ## ""


coordinate1<-input$Coordinates[1]
coordinate2<-input$Coordinates[2]
file2<-paste(chromosome,coordinate1,coordinate2,'query',sep='_')
write.table(file2,"target0.gff",row.names=FALSE,col.names=FALSE,quote = FALSE,sep="\t")

system("perl -p -i -e 's/_/\t/g' target0.gff");

target_folder0<-"/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/"
target_folder1<-paste(target_folder0,Genome,"/",sep = "")
bedtools0 <- paste(target_folder1,Genome,"_",chromosome,".fa",sep = "");
system_bedtools <- paste("/usr/local/bin/bedtools getfasta -fi", bedtools0 ,"-bed target0.gff -nameOnly -fo query.fasta",sep = " ");
system(system_bedtools);

# replace_pattern_in(file_contents=gene,replace="query",file_pattern="query.fasta",file_contents_perl = FALSE);

system("rm target0.gff");

query_fa <- "query.fasta";
dna <- readDNAStringSet("query.fasta");
#####

Gene <- input$Gene  ## from shiny input
cds_ids <- input$Gene ## from shiny input, for plotting! 

query <- "query"

query_length<-length(readDNAStringSet("query.fasta")[[1]]); # target CDS size

system("cat query.fasta > query_COPY.fasta");

system("cat query.fasta > query_COPY2.fasta");

Name_update_1 <- paste(Gene,"_mRNA",sep = "");
Name_update_2 <- paste(Gene,"_CDS",sep = "");

replace_pattern_in(file_contents=query,replace=Name_update_1,file_pattern="query.fasta",file_contents_perl = FALSE);
replace_pattern_in(file_contents=query,replace=Name_update_2,file_pattern="query_COPY.fasta",file_contents_perl = FALSE);

replace_pattern_in(file_contents=query,replace=input$Chr,file_pattern="query_COPY2.fasta",file_contents_perl = FALSE);

system("cat query.fasta query_COPY.fasta > query_Working.fasta");

Name_update_3 <- paste(Gene,"_ref",sep = "");
file.rename("query_Working.fasta", Name_update_3);

system("rm query.fasta query_COPY.fasta");
writeXStringSet(dna, query_fa, append=FALSE, compress=FALSE, format="fasta");


Name_update_4 <- paste("query","_",input$Chr,".fa",sep = "");
file.rename("query_COPY2.fasta", Name_update_4);


system_blast <- paste("/usr/local/bin/makeblastdb -in", Name_update_4 ,"-dbtype nucl",sep = " ");

system(system_blast);


result_folder1<-c('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/');
result_folder2<-paste(result_folder1,"query",sep="");
current.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny"

result_files<-list.files(current.folder, pattern = Name_update_4)
file.copy(result_files, result_folder2)

}
########## 08/11/22 Test coordinate

##############
##############
source('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny/functions/Main_plot.R',local = TRUE)
Main_plot()

########
########
########
########


} else if (input$fasta!='') {                                           #### input$fasta not null 08/09/22


# output$Testvalue <- renderText({ input$fasta }) 

query0<-input$fasta
tmp <- tempfile(fileext = ".fa")
   if (startsWith(query0, ">")){
     writeLines(query0, tmp)
   } else {
     writeLines(paste0(">Query\n",query0), tmp)
   }
dna <- readDNAStringSet(tmp)
query_fa <- "query.fasta";
writeXStringSet(dna, query_fa, append=FALSE, compress=FALSE, format="fasta");

##########
########## 08/09/22

Gene <- input$Gene  ## from shiny input
cds_ids <- input$Gene ## from shiny input, for plotting! 

query <- "query"

query_length<-length(readDNAStringSet("query.fasta")[[1]]); # target CDS size

system("cat query.fasta > query_COPY.fasta");

system("cat query.fasta > query_COPY2.fasta");

Name_update_1 <- paste(Gene,"_mRNA",sep = "");
Name_update_2 <- paste(Gene,"_CDS",sep = "");

replace_pattern_in(file_contents=query,replace=Name_update_1,file_pattern="query.fasta",file_contents_perl = FALSE);
replace_pattern_in(file_contents=query,replace=Name_update_2,file_pattern="query_COPY.fasta",file_contents_perl = FALSE);

replace_pattern_in(file_contents=query,replace=input$Chr,file_pattern="query_COPY2.fasta",file_contents_perl = FALSE);

system("cat query.fasta query_COPY.fasta > query_Working.fasta");

Name_update_3 <- paste(Gene,"_ref",sep = "");
file.rename("query_Working.fasta", Name_update_3);

system("rm query.fasta query_COPY.fasta");
writeXStringSet(dna, query_fa, append=FALSE, compress=FALSE, format="fasta");


Name_update_4 <- paste("query","_",input$Chr,".fa",sep = "");
file.rename("query_COPY2.fasta", Name_update_4);


system_blast <- paste("/usr/local/bin/makeblastdb -in", Name_update_4 ,"-dbtype nucl",sep = " ");

system(system_blast);


result_folder1<-c('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/');
result_folder2<-paste(result_folder1,"query",sep="");
current.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny"

result_files<-list.files(current.folder, pattern = Name_update_4)
file.copy(result_files, result_folder2)
##########################################
##########################################

### Here is Main_plot() function 
source('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny/functions/Main_plot.R',local = TRUE)
Main_plot()


}                       ## input$fasta not null 08/09/22;







}) ## observeEvent submit !!




#################### To test large file upload 08/25/22
####################

observeEvent(input$Largefile,{

########### with progress ...
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())

    progress$set(message = 'In progress',
                 detail = 'This may take a while...')
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.5)
    }

########## 08/26/22 Test file
########## 08/26/22 Test file
#if (input$Pickformat == "file" & input$upload!=''){


Genome<-input$Pickgenome
chromosome<-input$Chr
Selected_gene<-input$Gene

query <- "query"

current.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny"

Genome2<-paste(Genome,".*","_CDS_Working.gff3",sep="") ##

file0<-list.files(current.folder,pattern = Genome2);

file00<-paste(current.folder,"/",file0,sep="")

file1<-read.table(file0,header=F)

file2<-file1[which(file1[,9]==Selected_gene),]
write.table(file2,"target0.gff",row.names=FALSE,col.names=FALSE,quote = FALSE,sep="\t")

target_folder0<-"/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/"
target_folder1<-paste(target_folder0,Genome,"/",sep = "")
bedtools0 <- paste(target_folder1,Genome,"_",chromosome,".fa",sep = "");
system_bedtools <- paste("/usr/local/bin/bedtools getfasta -fi", bedtools0 ,"-bed target0.gff -nameOnly -fo query0.fasta",sep = " ");
system(system_bedtools);

#replace_pattern_in(file_contents=gene,replace=query,file_pattern="query.fasta",file_contents_perl = FALSE);

system("bash extract_fasta.sh");

#system("perl -p -i -e 's/gene/query/g' query.fasta");

query_fa <- "query.fasta";
dna <- readDNAStringSet("query.fasta");
#####

Gene <- input$Gene  ## from shiny input
cds_ids <- input$Gene ## from shiny input, for plotting! 

query_length<-length(readDNAStringSet("query.fasta")[[1]]); # target CDS size

system("cat query.fasta > query_COPY.fasta");

system("cat query.fasta > query_COPY2.fasta");

Name_update_1 <- paste(Gene,"_mRNA",sep = "");
Name_update_2 <- paste(Gene,"_CDS",sep = "");

replace_pattern_in(file_contents=query,replace=Name_update_1,file_pattern="query.fasta",file_contents_perl = FALSE);
replace_pattern_in(file_contents=query,replace=Name_update_2,file_pattern="query_COPY.fasta",file_contents_perl = FALSE);

replace_pattern_in(file_contents=query,replace=input$Chr,file_pattern="query_COPY2.fasta",file_contents_perl = FALSE);

system("cat query.fasta query_COPY.fasta > query_Working.fasta");

Name_update_3 <- paste(Gene,"_ref",sep = "");
file.rename("query_Working.fasta", Name_update_3);

system("rm query.fasta query_COPY.fasta target0.gff");
writeXStringSet(dna, query_fa, append=FALSE, compress=FALSE, format="fasta");

Name_update_4 <- paste("query","_",input$Chr,".fa",sep = "");
file.rename("query_COPY2.fasta", Name_update_4);

system_blast <- paste("/usr/local/bin/makeblastdb -in", Name_update_4 ,"-dbtype nucl",sep = " ");

system(system_blast);

result_folder1<-c('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/');
result_folder2<-paste(result_folder1,"query",sep="");
current.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny"

result_files<-list.files(current.folder, pattern = Name_update_4)
file.copy(result_files, result_folder2)

#}
########## 08/26/22 Test file
########## 08/26/22 Test file



##############
system("cat extract_syn_fa_Wheat_Test_08_24.pl > extract_syn_fa_Wheat_Test_08_24_Test5_copy_4.pl")


Upstream_update <- paste('p_size=',input$Upstream,sep = "");  ## from shiny input
Downstream_update <- paste('a_size=',input$Downstream,sep = "");  ## from shiny input
replace_pattern_in(file_contents="p_size=10000",replace=Upstream_update,file_pattern="Test5_copy_4",file_contents_perl = TRUE);
replace_pattern_in(file_contents="a_size=1000",replace=Downstream_update,file_pattern="Test5_copy_4",file_contents_perl = TRUE);

Chr_update0 <- paste('target_ch=',input$Chr,sep = "'");  ## from shiny input
symbol<-c("'");
Chr_update <- paste(Chr_update0,symbol,sep = "");
replace_pattern_in(file_contents="target_ch='chr7D'",replace=Chr_update,file_pattern="Test5_copy_4",file_contents_perl = TRUE);

Gene_update0 <- paste('gene=',input$Gene,sep = "'");  ## from shiny input
Gene_update <- paste(Gene_update0,symbol,sep = "");
replace_pattern_in(file_contents="gene='TraesCS7D02G524200'",replace=Gene_update,file_pattern="Test5_copy_4",file_contents_perl = TRUE);



## User folder name
          userfolder01 <- as.character(strsplit(input$upload$name, ".fa.gz"))
          userfolder1 <- gsub('_', "", userfolder01, fixed = FALSE)

replace_pattern_in(file_contents="bszhangTest",replace=userfolder1,file_pattern="Test5_copy_4",file_contents_perl = TRUE);
## User folder name



##############
system("perl extract_syn_fa_Wheat_Test_08_24_Test5_copy_4.pl 1");


######################################################################## To search for the outliers!
result_folder1<-c('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/IWGSC/Candidate_genes/');
result_folder2<-paste(result_folder1,Gene,sep="");
setwd(result_folder2);
## backup Haplotype_syn file
Outlier0<-paste(Gene,'_Haplotype_syn',sep = "");
Outlier00<-paste(Gene,'_Blast_Original',sep = "");
Outlier1<-read.table(Outlier0,header=T); ## _Haplotype_syn
write.table(Outlier1, file = Outlier00,sep= "\t",quote = FALSE,row.names = FALSE); 



########
########
#Outlier5 <- Outlier1[which(Outlier1$Similarity>=input$Similarityfilter),]; ## filter1: Similarity >=96, kept!
#Filtered_output <- list();
#for(i in 1:length(unique(Outlier5[,4]))){
#Outlier2<-Outlier5[Outlier5[,4]==unique(Outlier5[,4])[i],]
#Outlier2_1<-Outlier2[order(-Outlier2$Similarity),]  ## 6/14, rank the Similarity, from high to low!
#Outlier2 <- Outlier2_1 ## 6/14, rank the Similarity, from high to low!
#Outlier3<-Outlier2[,6][as.matrix(dist(Outlier2[,6]))[,1]<=input$Distancefilter]; ## filter2: <= 200 kb, kept
#out_ind<-which(Outlier2[,6] %in% Outlier3);
#Filtered_output[[i]]<-Outlier2[out_ind, ];
#rm(Outlier2,Outlier3,out_ind);
#} 

#Filtered_output2 <- as.data.frame(rbindlist(Filtered_output));

#write.table(Filtered_output2, file = Outlier0,sep= "\t",quote = FALSE,row.names = FALSE);
########
######## 06/21
FilterNew0<-Outlier1



Filter_list0 <- list()

Information_list<-list()

indexNew1<-1

for(indexNew in unique(FilterNew0$Genome)){

FilterNew1<- FilterNew0[which(FilterNew0$Genome==indexNew),]


if(nrow(FilterNew1)>=2){

Trees<-length(unique(cutree(hclust(dist(FilterNew1[,6])), h = input$Distancefilter)))

hclusters_similarity<-list()
hclusters_size<-list()

Information_matrix_Members<-matrix(nrow=Trees,ncol=1)

for(Treecut in 1:Trees) {
hclusters <-FilterNew1[which(cutree(hclust(dist(FilterNew1[,6])), k =Trees, h = input$Distancefilter ) ==Treecut),]
hclusters_similarity[Treecut]<-mean(hclusters[,9])
hclusters_size[Treecut]<-sum(hclusters[,8])

Information_matrix_Members[Treecut,1]<-nrow(hclusters)

}

hcluster_matrix<-matrix(nrow=Trees,ncol=2)

for(Treecut in 1:Trees) {
hcluster_matrix[Treecut,1]<-hclusters_similarity[[Treecut]]
hcluster_matrix[Treecut,2]<-hclusters_size[[Treecut]]
}
hcluster_matrix[,2]<-hcluster_matrix[,2]/query_length
Size_filter1<-which(hcluster_matrix[,2]>=input$CDSfilter[1] & hcluster_matrix[,2]<=input$CDSfilter[2])


if(length(Size_filter1)>1){
Similarity_filter1<-which(max(hcluster_matrix[Size_filter1,][,1])==hcluster_matrix[,1])
}else{

Similarity_filter1<-which(max(hcluster_matrix[,1])==hcluster_matrix[,1])

}

Target_cluster1<-intersect(Size_filter1,Similarity_filter1)

if(length(Target_cluster1)==0){    ###  6/21
Ideal_Size <- 1.0
Target_cluster1<- which(abs(hcluster_matrix[,2] - Ideal_Size) == min(abs(hcluster_matrix[,2] - Ideal_Size)))

if(length(Target_cluster1)>1){
Target_cluster1<-which(max(hcluster_matrix[Target_cluster1,][,1])==hcluster_matrix[,1])
}
}


Filter_list0[[indexNew1]]<-FilterNew1[which(cutree(hclust(dist(FilterNew1[,6])), k =Trees, h = input$Distancefilter ) == Target_cluster1),]


Information_matrix<-matrix(nrow=Trees,ncol=9)
Information_matrix[,1]<-rep(indexNew, Trees)
Information_matrix[1,2]<-nrow(FilterNew1)
Information_matrix[1,3]<-input$Distancefilter/1000
Information_matrix[1,4]<-Trees
Information_matrix[,5]<-1:nrow(hcluster_matrix)
Information_matrix[,6]<-Information_matrix_Members
Information_matrix[,7]<-round(hcluster_matrix[,1], digits = 3)
Information_matrix[,8]<-round(hcluster_matrix[,2], digits = 3)
Information_matrix[Target_cluster1,9]<-c("Selected")

colnames(Information_matrix)<-c("Genomes","Positions","HeightCut (kb)","TotalClusters","ClusterIndex","Members","MeanSimilarity","TotalLength/IWGSC","CandidateCluster")

Information_Table<-as.data.frame(Information_matrix)
Information_list[[indexNew1]]<-Information_Table


}else{

Filter_list0[[indexNew1]]<-FilterNew0[which(FilterNew0$Genome==indexNew),]

}


indexNew1<-indexNew1+1

}

Filtered_output2 <-as.data.frame(rbindlist(Filter_list0))
write.table(Filtered_output2, file= Outlier0,sep= "\t",quote = FALSE,row.names = FALSE);


if(length(Information_list)!=0){
Information_output0 <-as.data.frame(rbindlist(Information_list))
Information_output<-Information_output0[which(Information_output0$Genomes!="IWGSCNEW"),]
}

######## 06/21
########

if(nrow(anti_join(Outlier1,Filtered_output2))!=0){
Outlier4<-anti_join(Outlier1,Filtered_output2); # differences of _Haplotype_syn and _Blast_Original; Not shown in plot!
Outlier4_name <- unique(Outlier4[,4]);
Outlier4_name2 <- intersect(input$id,Outlier4_name);

#output$text <- renderText({ paste('The following haplotypes may have other potential places (Table below) in the same chromosome (at least 200 kb from here):',sep = "") });
#output$text2 <- renderText({Outlier4_name2});
#output$text3 <- renderText({paste0(c('For example, IWGSC contains ',Outlier4[which(Outlier4$Genome=="IWGSC"),]$Size,' bp DNAs mapped onto other places.'), collapse= " ",sep = "")});
}

## 
new.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny";
setwd(new.folder);

#gff3_Coor<-read.table("Triticum_aestivum.IWGSC.53_working_gene_coordinate",header=F);
#gff3_Coor1<-gff3_Coor[which(gff3_Coor$V3==Gene),][,1:2]; # annotation in GFF3 (only gene)

#Mapped0<-data.frame(Filtered_output2[which(Filtered_output2$Genome=="IWGSC"),][,6],Filtered_output2[which(Filtered_output2$Genome=="IWGSC"),][,7]);
#Mapped1<-cbind(min(as.vector(as.matrix(Mapped0))),max(as.vector(as.matrix(Mapped0)))); # Coordinate as shown in plot (all CDS combined)
#output$text4 <- renderText({paste0(c('Submitted gene (IWGSC) annotated between: ', gff3_Coor1 ,'; Plotted CDS (IWGSC) contained between: ',Mapped1), collapse= " ",sep = "")});

######################################################################## To search for the outliers!


system("perl extract_syn_fa_Wheat_Test_08_24_Test5_copy_4.pl 2");

# system("rm extract_syn_fa_Wheat_Test_08_24_Test5_copy_4.pl");


result_folder1<-c('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/IWGSC/Candidate_genes/')
result_folder2<-paste(result_folder1,Gene,sep="")
new.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny"
setwd(result_folder2)
result_files<-list.files(result_folder2, pattern = Gene)
file.copy(result_files, new.folder)


setwd(new.folder)
####################################################
############## Here is plot function ??




var_types <- c('snp', 'ins', 'del');
indel_col <- c("grey", "red", "red");
cds_col <- c("yellowgreen", "brown");
Dir <- paste('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny', '/', sep = '');

CDS_gDNA_blast <- read.table(paste(Dir, Gene, '_ref_mRNA-Haplotype_out_m8', sep = ''), sep = '\t', header = F, stringsAsFactors = F);
gDNAs_blast <- read.table(paste(Dir, Gene, '_Haplotype-Self_out_m8', sep = ''), sep = '\t', header = F, stringsAsFactors = F);

gDNAs_blast <- subset(gDNAs_blast, gDNAs_blast[,4]> 0)    ## >0 ?

N_Gap <- read.table(paste(Dir, Gene, '_Haplotype_N_Gaps', sep = ''), sep = '\t', header = T, stringsAsFactors = F);
Anno <- read.table(paste(Dir, Gene, '_Haplotype_anno', sep = ''), sep = '\t', header = T, stringsAsFactors = F)


variations <- read.table(paste(Dir, Gene, '_Variations', sep = ''), sep = '\t', header = T, stringsAsFactors = F);


if (!file.size(paste(Dir, Gene, '_repMask2', sep = ''))==0) {
repmask <- read.table(paste(Dir, Gene, '_repMask2', sep = ''), header = F, sep = "\t", stringsAsFactors = F);
repmask <- repmask[abs(repmask[,8] - repmask[,7])> 100,]
}

genomes <- input$id ## from shiny input

x_lim <- range(gDNAs_blast[,9:10]) + c(0, 2000);

output_flag = 0
#genomes_r <- genomes
#g_lab <- genomes_r


####################### To check genomes ?
Outlier0<-paste(Gene,'_Haplotype_syn',sep = "");
Outlier1<-read.table(Outlier0,header=T);
Test_Genome0<-unique(Outlier1[,4]);

g_lab <- intersect(genomes,Test_Genome0);
genomes_r <- intersect(genomes,Test_Genome0);
####################### To check genomes ?

output$plot <- renderPlot({


if (!file.size(paste(Dir, Gene, '_repMask2', sep = ''))==0) {
source('Source_Modified.R', local = TRUE) ## To load the Plot_SV_NEW ?? (R) function
}else {
source('Source_Modified_NORepeat.R', local = TRUE) ## To load the Plot_SV_NEW ?? (R) function
}

####################################################
Plot_SV(genomes_r, g_lab, repmask, CDS_gDNA_blast, gDNAs_blast, N_Gap, Anno, variations, output_flag) ## plot in shiny



########################## A table showing clustering results
if(length(Information_list)!=0){

output$table00 <-DT::renderDataTable({
datatable(Information_output,caption = htmltools::tags$caption(
    style = 'caption-side: bottom; text-align: center;',
    'Table 1: ', htmltools::em('Clustering information based on Blast results (At least two places).')
  ), filter = 'top', extensions = 'Buttons',selection = list(target = 'row+column'),
              class="cell-border stripe",
              options = list(dom = "Blfrtip",
                             buttond = list("copy", list(extend = "collection",
                                                         buttons = c("csv"),
                                                         text = "Downloads")), pageLength=10, autoWidth = TRUE,
                             searchHighlight = TRUE, filter = "top")) %>% formatStyle(columns=1:ncol(Information_output), target = c("cell"),backgroundColor = styleEqual(c("Selected"), c('lightblue')))
  }) # DT::renderDataTable

}
########################## A table showing clustering results

########################## A table showing Blast result which is presented in main plot
Filtered_output2_Plotted <- Filtered_output2[which(Filtered_output2$Genome!="IWGSCNEW"),]

output$table0 <-DT::renderDataTable({
datatable(Filtered_output2_Plotted,caption = htmltools::tags$caption(
    style = 'caption-side: bottom; text-align: center;',
    'Table 2: ', htmltools::em('Blast results used for plotting.')
  ), filter = 'top', extensions = 'Buttons',selection = list(target = 'row+column'),
              class="cell-border stripe",
              options = list(dom = "Blfrtip",
                             buttond = list("copy", list(extend = "collection",
                                                         buttons = c("csv"),
                                                         text = "Downloads")), pageLength=10, autoWidth = TRUE,
                             searchHighlight = TRUE, filter = "top")) %>% formatStyle(columns=c(4,8,9), target = c("cell"), backgroundColor = c('yellow'))
  }) # DT::renderDataTable
########################## A table showing Blast result which is presented in main plot

########################## A table showing Blast result which is not presented in main plot ??
Outlier1_0<-read.table(Outlier00,header=T); ## Blast_Original
Outlier1_1 <- Outlier1_0[which(Outlier1_0$Genome!="IWGSCNEW"),]
Outlier1_2 <- Outlier1[which(Outlier1$Genome!="IWGSCNEW"),]

NotShown0 <- anti_join(Outlier1_1,Outlier1_2) # Not shown in plot, other genomic positions.

if(nrow(NotShown0)!=0){

output$table2 <-DT::renderDataTable({
datatable(NotShown0, caption = htmltools::tags$caption(
    style = 'caption-side: bottom; text-align: center;',
    'Table 3: ', htmltools::em('Blast results not shown in plot.')
  ),filter = 'top', extensions = 'Buttons',selection = list(target = 'row+column'),
              class="cell-border stripe",
              options = list(dom = "Blfrtip",
                             buttond = list("copy", list(extend = "collection",
                                                         buttons = c("csv"),
                                                         text = "Downloads")), pageLength=10, autoWidth = TRUE,
                             searchHighlight = TRUE, filter = "top")) %>% formatStyle(columns=c(4,8,9), target = c("cell"), backgroundColor = c('orange'))
  }) # DT::renderDataTable



    df <- reactive({
          NotShown0[input$table2_rows_selected, ];
    })

    output$table3 <- renderDataTable({
      df()
    })


   observeEvent(input$others,{                                                         ## To select and download datatable!
    Other_P <-paste(result_folder2,'/',Outlier0,'_Test',sep='');
    data <- df();
    write.table(data, file = Other_P ,sep= "\t",quote = FALSE,row.names = FALSE);
    })   ## save selected items
#######
#######
####### Try to run SecondRound on selected table items!




#######
#######
####### Try to run SecondRound on selected table items!
}
########################## A table showing Blast result which is not presented in main plot ??


    output$save <- downloadHandler(
    file = "save.png" ,
    content = function(file) {
    png(file = file,width = 1200, height = 1200,res = 240)
    Plot_SV(genomes_r, g_lab, repmask, CDS_gDNA_blast, gDNAs_blast, N_Gap, Anno, variations, output_flag)
    dev.off()
    }) ## save figure from shiny output

    

    output$Repeat <- downloadHandler(
    file = "repmask.csv" ,
    content = function(file) {
      write.csv(repmask, file, row.names = FALSE)
    }
  )   ## save repeat file from shiny output
   
    output$Variation <- downloadHandler(
    file = "Variation.csv" ,
    content = function(file) {
      write.csv(variations, file, row.names = FALSE)
    }
  )  ## save variation file from shiny output

    output$Haplotype_fa <- downloadHandler(
    file = "Haplotype.fa" ,
    content = function(file) {
    file.copy("MOC1_Haplotype.fa", file)
   } 
  )  ## save Haplotype fasta file from shiny output

    output$CDS_fa <- downloadHandler(
    file = "query_CDS.fa" ,
    content = function(file) {
    file.copy("query.fasta", file)
   } 
  )  ## save target's CDS fasta file from shiny output

  }) ## renderPlot Main



############### Updated plot for zoom_in
x_lim1 <- range(gDNAs_blast[,9:10]);
updateSliderInput(session, "Position",value = c(input$Position[1],input$Position[2]), min = 1, max=x_lim1[2]); # sliderInput update ?
# updateSliderInput(session, "Position2", max=x_lim1[2]); # sliderInput update ?

observeEvent(input$submit2,{
x_lim2 <- c(input$Position[1], input$Position[2]); ## New coordinates needed ?
x_lim3 <- (input$Position[1]+input$Position[2])*1/2;
NEWsizes <- (input$Position[2]-input$Position[1])/1000;

if (!file.size(paste(Dir, Gene, '_repMask2', sep = ''))==0) {
repmask <- read.table(paste(Dir, Gene, '_repMask2', sep = ''), header = F, sep = "\t", stringsAsFactors = F);
repmask <- repmask[abs(repmask[,8] - repmask[,7])> 100,]
}

output$plot2 <- renderPlot({
if (!file.size(paste(Dir, Gene, '_repMask2', sep = ''))==0) {
source('Source_Modified2.R', local = TRUE) ## To load the Plot_SV_NEW ?? (R) function
}else {
source('Source_Modified2_NORepeat.R', local = TRUE) ## To load the Plot_SV_NEW ?? (R) function
}
####################################################
Plot_SV2(genomes_r, g_lab, repmask, CDS_gDNA_blast, gDNAs_blast, N_Gap, Anno, variations, output_flag) ## plot in shiny based on x_lim2 ?
  }) ## renderPlot2, another plot ..
}) ## observeEvent submit2 !!
############### Updated plot for zoom_in



############### Another plot for clustering all genomes ?? Figure
observeEvent(input$submit3,{

output$plot3 <- renderPlot({

gDNAs_blast<-gDNAs_blast[which(gDNAs_blast[,1]!='IWGSCNEW' & gDNAs_blast[,2]!='IWGSCNEW'), ]

genomes <- unique(gDNAs_blast[,2]);

n_g <- length(genomes)  
b_matrix <- diag(n_g)
for (g1 in 1:(n_g - 1)) {
 for (g2 in (g1 + 1): n_g ) {
  gDNAs <- subset(gDNAs_blast, gDNAs_blast[,1] == genomes[g1] & gDNAs_blast[,2] == genomes[g2]);
  q_hits <- c(gDNAs[,7], gDNAs[,8]); s_hits <- c(gDNAs[,9,], gDNAs[,10])
  b_matrix[g1, g2] <- lm(q_hits ~ s_hits)$coeff[2]

  gDNAs <- subset(gDNAs_blast, gDNAs_blast[,1] == genomes[g2] & gDNAs_blast[,2] == genomes[g1]);
  q_hits <- c(gDNAs[,7], gDNAs[,8]); s_hits <- c(gDNAs[,9,], gDNAs[,10])
  b_matrix[g2, g1] <- lm(q_hits ~ s_hits)$coeff[2]
 }
}
colnames(b_matrix) <- genomes;
rownames(b_matrix) <- genomes;
h_c <- hclust(dist(b_matrix));

#as.dendrogram(h_c) %>% set("labels_cex", 1.2) %>% sort(type = "nodes") %>% highlight_branches %>% plot(main = "Clustering on all genomes",horiz = TRUE);
##
genomes <- genomes[h_c$order]
hc <- hclust(dist(b_matrix))
cluster_cut <- length(genomes);
hc_r <- cutree(hc, k = cluster_cut)
genomes_r <- c();
for (k in 1:cluster_cut) {
 genomes_r[k] <- names(hc_r)[min(which(hc_r == k))]
}
genomes_r <-  names(sort(hc_r));
g_lab <- genomes_r;
if (!file.size(paste(Dir, Gene, '_repMask2', sep = ''))==0) {
source('Source_Modified.R', local = TRUE) ## To load the Plot_SV_NEW ?? (R) function
}else {
source('Source_Modified_NORepeat.R', local = TRUE) ## To load the Plot_SV_NEW ?? (R) function
}
Plot_SV(genomes_r, g_lab, repmask, CDS_gDNA_blast, gDNAs_blast, N_Gap, Anno, variations, output_flag) ## plot in shiny
##

  }) ## renderPlot3
}) ## observeEvent submit3 !!
############### Another plot for clustering all genomes ?? Figure
############### Another plot for clustering all genomes ?? Tree
observeEvent(input$submit3_2,{

output$plot4 <- renderPlot({

gDNAs_blast<-gDNAs_blast[which(gDNAs_blast[,1]!='IWGSCNEW' & gDNAs_blast[,2]!='IWGSCNEW'), ]

genomes <- unique(gDNAs_blast[,2]);

n_g <- length(genomes)  
b_matrix <- diag(n_g)
for (g1 in 1:(n_g - 1)) {
 for (g2 in (g1 + 1): n_g ) {
  gDNAs <- subset(gDNAs_blast, gDNAs_blast[,1] == genomes[g1] & gDNAs_blast[,2] == genomes[g2]);
  q_hits <- c(gDNAs[,7], gDNAs[,8]); s_hits <- c(gDNAs[,9,], gDNAs[,10])
  b_matrix[g1, g2] <- lm(q_hits ~ s_hits)$coeff[2]

  gDNAs <- subset(gDNAs_blast, gDNAs_blast[,1] == genomes[g2] & gDNAs_blast[,2] == genomes[g1]);
  q_hits <- c(gDNAs[,7], gDNAs[,8]); s_hits <- c(gDNAs[,9,], gDNAs[,10])
  b_matrix[g2, g1] <- lm(q_hits ~ s_hits)$coeff[2]
 }
}
colnames(b_matrix) <- genomes;
rownames(b_matrix) <- genomes;
h_c <- hclust(dist(b_matrix));

as.dendrogram(h_c) %>% set("labels_cex", 0.8) %>% sort(type = "nodes") %>% highlight_branches %>% plot(main = "Clustering on all genomes",horiz = TRUE);

  }) ## renderPlot4
}) ## observeEvent submit3_2 !!
############### Another plot for clustering all genomes ?? Tree







}) # observeEvent input$Largefile

#################### To test large file upload 08/25/22, The End
####################








#### to reveal captured genes! 08/18/22
observeEvent(input$GeneCapture,{


current.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny"
Gene<-input$Gene;


ForMatch_Gene0 <- paste(Gene,"_Blast_Original",sep="")
ForMatch_Gene1 <- read.table(ForMatch_Gene0,header=T)
ForMatch_Gene2 <- ForMatch_Gene1[which(ForMatch_Gene1[,4]!="IWGSC" & ForMatch_Gene1[,4]!="query"),]

for(i in 1:nrow(ForMatch_Gene2)){                   ## switch coordinates for bedtools!

    if(ForMatch_Gene2[i,6]>ForMatch_Gene2[i,7]){
           
           Tem6<-ForMatch_Gene2[i,7]
           Tem7<-ForMatch_Gene2[i,6]
           ForMatch_Gene2[i,6]<-Tem6
           ForMatch_Gene2[i,7]<-Tem7
    }

}

ForMatch_Gene3 <- ForMatch_Gene2[,c(5,6,7)] 
ForMatch_Gene4 <- data.frame(ForMatch_Gene2[,4:9])
ForMatch_Gene3$Newone<- with(ForMatch_Gene4, paste(ForMatch_Gene4$Genome, ForMatch_Gene4$ch, ForMatch_Gene4$sbj_St, ForMatch_Gene4$sbj_E,ForMatch_Gene4$Size,ForMatch_Gene4$Similarity,sep='_'))
ForMatch_Gene3$Newone1<- ForMatch_Gene4$Genome


for(uniqGenome in unique(ForMatch_Gene3$Newone1) ){

ForMatch_Gene5 <- subset(ForMatch_Gene3, ForMatch_Gene3$Newone1==uniqGenome)[,1:4]
write.table(ForMatch_Gene5,"Gene_captured1.bed",row.names=FALSE,col.names=FALSE,quote = FALSE,sep=" ")
system("perl -p -i -e 's/ /\t/g' Gene_captured1.bed")
Genome2<-paste("Triticum.*",uniqGenome,".*","_gene_Working.gff3",sep="") ##
file0<-list.files(current.folder,pattern = Genome2);
system_bedtools0<- paste("bedtools intersect -a ", file0, " -b Gene_captured1.bed >> Gene_captured_All.txt ",sep = " ")
system(system_bedtools0)
system("rm Gene_captured1.bed")

}


    if(file.size("Gene_captured_All.txt")!=0){
   
    Test_GeneName0 <- read.table("Gene_captured_All.txt",header=F)
    colnames(Test_GeneName0)<-c("chr","column2","gene","coor1","coor2","column6","column7","Size","Gene_Captured")
    Test_GeneName0[,8]<-Test_GeneName0[,5]-Test_GeneName0[,4]+1


output$table4 <-DT::renderDataTable({
datatable(Test_GeneName0,caption = htmltools::tags$caption(
    style = 'caption-side: bottom; text-align: center;',
    'Table 4: ', htmltools::em('Captured genes (overlap of blast results with annotated genes) in all genomes for Both Table2 And Table3.')
  ), filter = 'top', extensions = 'Buttons',selection = list(target = 'row+column'),
              class="cell-border stripe",
              options = list(dom = "Blfrtip",
                             buttond = list("copy", list(extend = "collection",
                                                         buttons = c("csv"),
                                                         text = "Downloads")), pageLength=50, autoWidth = TRUE,
                             searchHighlight = TRUE, filter = "top")) %>% formatStyle(columns=c(8,9), target = c("cell"), backgroundColor = c('gold'))
  }) # DT::renderDataTable


    } else if (file.size("Gene_captured_All.txt")==0){
    Test_GeneName2 <- paste("No gene captured in all wheat genomes!",sep='');
    output$Captured_Gene <- renderText({ Test_GeneName2 });   ## 
    }

system("rm Gene_captured_All.txt")


})
#### to reveal captured genes! 08/18/22


####
observeEvent(input$done,{

Gene <- input$Gene  ## from shiny input

if(Gene!=''){
new.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny"
result_files<-list.files(new.folder, pattern = Gene)
file.remove(result_files)
file.remove("query.fasta","grasrep.fa.nhr","grasrep.fa.nin","grasrep.fa.nsq")

file.remove("extract_syn_fa_Wheat_1_05_17_Test4_copy_3.pl")

file.remove("gDNAs_blast_NEW","Selected_lines_Coordinates.bed")

Name_update_4 <- paste("query","_",input$Chr,".fa",sep = "");
result_files2<-list.files(new.folder, pattern = Name_update_4)
file.remove(result_files2)

if(input$fasta!=''){
query.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/query"
setwd(query.folder)
result_files3<-list.files(query.folder, pattern = Name_update_4)
file.remove(result_files3)
setwd(new.folder)
}


if(input$Pickgenome!=''){
query.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/query"
setwd(query.folder)
result_files3<-list.files(query.folder, pattern = Name_update_4)
file.remove(result_files3)
setwd(new.folder)
}


refresh()

}else{refresh()}

}) ## observeEvent DONE, To remove all temp files






    } # server function of Page_6



  ) # page for Page_6
} # Page_6 function
