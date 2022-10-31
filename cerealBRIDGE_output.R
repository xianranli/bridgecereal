cerealBRIDGE_output <- function(User_folder0,perlArg0_db_sp,perlArg1_PickGenome ,perlArg2_PickGene,perlArg3_PickChr,perlArg4_Users_folder,perlArg5_PickUp,perlArg6_PickDown, Backup_folder,strand_direction,database_folder,gff_folder,script_folder,User_folder) {

setwd(perlArg4_Users_folder);

system1 <- paste("perl", paste(script_folder,'extract_syn_fa.pl',sep=''),perlArg0_db_sp,perlArg1_PickGenome ,perlArg2_PickGene,perlArg3_PickChr,perlArg4_Users_folder,perlArg5_PickUp,perlArg6_PickDown, 1, 0,sep=' ');
system2 <- paste("perl", paste(script_folder,'extract_syn_fa.pl',sep=''),perlArg0_db_sp,perlArg1_PickGenome ,perlArg2_PickGene,perlArg3_PickChr,perlArg4_Users_folder,perlArg5_PickUp,perlArg6_PickDown, 0, 2,sep=' ');


Dir <- paste(User_folder, User_folder0,'/',sep = '');




##############
system(system1);


######################################################################## To search for the outliers!

Backup_folder_Gene<-paste(Backup_folder,Gene,sep="");

setwd(Backup_folder_Gene);

## backup Haplotype_syn file
BlastSyn<-paste(Gene,'_Haplotype_syn',sep = "");
Blast_Ori<-paste(Gene,'_Blast_Original',sep = "");
BlastSynWorking<-read.table(BlastSyn,header=T); ## _Haplotype_syn
write.table(BlastSynWorking, file = Blast_Ori,sep= "\t",quote = FALSE,row.names = FALSE); 


#### 10/24/22 syn file column 6-7 start-end formation
for(i in 2:nrow(BlastSynWorking)) {           
if(BlastSynWorking[i,6]>BlastSynWorking[i,7]){  
   temp6<- BlastSynWorking[i,6]
   temp7<- BlastSynWorking[i,7]
   BlastSynWorking[i,7]<- temp6
   BlastSynWorking[i,6]<- temp7
}
}
#### 10/24/22 

########
########
######## 06/21
FilterNew0<-BlastSynWorking

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

Filtered_HaplotypeSyn <-as.data.frame(rbindlist(Filter_list0))

write.table(Filtered_HaplotypeSyn, file= BlastSyn,sep= "\t",quote = FALSE,row.names = FALSE);


if(length(Information_list)!=0){
Information_output0 <-as.data.frame(rbindlist(Information_list))
Information_output<-Information_output0[which(Information_output0$Genomes!=''),]
}

######## 06/21


########

if(nrow(anti_join(BlastSynWorking,Filtered_HaplotypeSyn))!=0){
Outlier4<-anti_join(BlastSynWorking,Filtered_HaplotypeSyn); # differences of _Haplotype_syn (NEW or filtered) and _Blast_Original; Not shown in plot!
Outlier4_name <- unique(Outlier4[,4]);
Outlier4_name2 <- intersect(input$id,Outlier4_name);

}

## 

setwd(perlArg4_Users_folder);

######################################################################## To search for the outliers!


system(system2);

Backup_folder_Gene<-paste(Backup_folder,Gene,sep="")

setwd(Backup_folder_Gene)

result_files<-list.files(Backup_folder_Gene, pattern = Gene)

file.copy(result_files, perlArg4_Users_folder)

setwd(perlArg4_Users_folder)


####################################################
####################### To check genomes ?
#BlastSyn<-paste(Gene,'_Haplotype_syn',sep = "");
#BlastSynWorking<-read.table(BlastSyn,header=T);
#Test_Genome0<-unique(BlastSynWorking[,4]);

#g_lab <- intersect(genomes,Test_Genome0);
#genomes_r <- intersect(genomes,Test_Genome0);
####################### To check genomes ?

############ 10/06/22 filter to remove crossover
Gene <- input$Gene
#Dir <- paste('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/User', '/', User_folder0,'/',sep = '');

CrossFilter0 <- read.table(paste(Dir, Gene, '_Haplotype-Self_out_m8', sep = ''), sep = '\t', header = F, stringsAsFactors = F);

combine_num <- ncol(combn(unique(CrossFilter0$V1),2))

for(index in 1:combine_num){

name1<- combn(unique(CrossFilter0$V1),2)[,index][1]
name2<- combn(unique(CrossFilter0$V1),2)[,index][2]

#print(index,name1,name2)

Genome1 <- CrossFilter0[which(CrossFilter0$V1 == name1 & CrossFilter0$V2 == name2),]
Genome1_1 <- CrossFilter0[which(CrossFilter0$V1 == name1 & CrossFilter0$V2 == name1),]
Genome2 <- CrossFilter0[which(CrossFilter0$V1 == name2 & CrossFilter0$V2 == name1),]
Genome2_2 <- CrossFilter0[which(CrossFilter0$V1 == name2 & CrossFilter0$V2 == name2),]

Genome_1_2 <- rbind(Genome1,Genome2, Genome1_1,Genome2_2)  ## V1==V2 or V1 !=V2

Genome_1_3 <- rbind(Genome1,Genome2) ##  Only V1 != V2

Genome_1_4 <- rbind(Genome1_1,Genome2_2) ##  Only V1 == V2

Self_candidate1<- Genome1_1[which(Genome1_1$V3 !=100), ][,7:8]
Self_candidate2<- Genome2_2[which(Genome2_2$V3 !=100), ][,7:8]

Genome_1_3[,13]<- c(row.names(Genome_1_3))

Genome_1_3_1 <-Genome_1_3[,c(13,7,8)]

colnames(Genome_1_3_1)<-c("ID","V7","V8")

Self_candidate1[,3]<- c(row.names(Self_candidate1))
Self_candidate2[,3]<- c(row.names(Self_candidate2))

Self_candidate1_1 <- Self_candidate1[,c(3,1,2)]
Self_candidate2_1 <- Self_candidate2[,c(3,1,2)]

colnames(Self_candidate1_1)<-c("ID","V7","V8")
colnames(Self_candidate2_1)<-c("ID","V7","V8")

Setdiff_ID1 <- (Genome_1_3_1  %>% semi_join(Self_candidate1_1,by = c("V7","V8"))) $ID
Setdiff_ID2 <- (Genome_1_3_1  %>% semi_join(Self_candidate2_1,by = c("V7","V8"))) $ID

Setdiff_ID <-c(Setdiff_ID1,Setdiff_ID2)

#CrossFilter1 <-CrossFilter0[!(row.names(CrossFilter0) %in% Setdiff_ID), ]

CrossFilter0 <-CrossFilter0[!(row.names(CrossFilter0) %in% Setdiff_ID), ]

# Genome_1_2[!(row.names(Genome_1_2) %in% Setdiff_ID), ]

}

write.table(CrossFilter0, file=paste(Gene, '_Haplotype-Self_out_m8_Crossover', sep = ''), sep="\t", quote = FALSE,row.names = FALSE,col.names = FALSE)
############ 10/06/22 filter to remove crossover


output$plot <- renderPlot({

output$plot2 <- NULL #09/28/22

output$plot3 <- NULL #09/28/22
output$plot4 <- NULL #09/28/22
output$info <- NULL #09/28/22
output$submit2 <-NULL #09/28/22

output$info2 <- NULL #09/28/22 
output$bucket <- NULL #09/28/22

output$Haplotypes <- NULL #09/28/22


Gene <- input$Gene  ## from shiny input
cds_ids <- input$Gene ## from shiny input, for plotting! 
Ref_genome <- input$Pickgenome

var_types <- c('snp', 'ins', 'del');
indel_col <- c("grey", "red", "red");
cds_col <- c("yellowgreen", "brown");

#Dir <- paste('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/User', '/', User_folder0, '/', sep = '');

CDS_gDNA_blast <- read.table(paste(Dir, Gene, '_ref_mRNA-Haplotype_out_m8', sep = ''), sep = '\t', header = F, stringsAsFactors = F);

gDNAs_blast <- read.table(paste(Dir, Gene, '_Haplotype-Self_out_m8', sep = ''), sep = '\t', header = F, stringsAsFactors = F);
gDNAs_blast <- subset(gDNAs_blast, gDNAs_blast[,4]> 0)

N_Gap <- read.table(paste(Dir, Gene, '_Haplotype_N_Gaps', sep = ''), sep = '\t', header = T, stringsAsFactors = F);
Anno <- read.table(paste(Dir, Gene, '_Haplotype_anno', sep = ''), sep = '\t', header = T, stringsAsFactors = F)
# variations <- read.table(paste(Dir, Gene, '_Variations', sep = ''), sep = '\t', header = T, stringsAsFactors = F);

if (!file.size(paste(Dir, Gene, '_repMask2', sep = ''))==0) {
repmask <- read.table(paste(Dir, Gene, '_repMask2', sep = ''), header = F, sep = "\t", stringsAsFactors = F);
repmask <- repmask[abs(repmask[,8] - repmask[,7])> 100,]
}

x_lim <- range(gDNAs_blast[,9:10]) + c(0, 2000)
output_flag = 0
##
##

#genomes <- sort(unique(gDNAs_blast[,2])); ## 09/26/22

genomes <- input$id ## from shiny input

####################### To check genomes ?
BlastSyn<-paste(Gene,'_Haplotype_syn',sep = "");
BlastSynWorking<-read.table(BlastSyn,header=T);
Test_Genome0<-unique(BlastSynWorking[,4]);

Test_Genome1 <- intersect(genomes,Test_Genome0);

if('query' %in% Test_Genome1){
 genomes <- c(sort(Test_Genome1[which(Test_Genome1 != 'query')]), 'query')  ## 10/05/22
} else {
 genomes <- Test_Genome1
}
####################### To check genomes ?

genomes_r <- genomes
g_lab <- genomes_r;

#b_matrix_groups2 <- 0
haplotypes <- 0


if (!file.size(paste(Dir, Gene, '_repMask2', sep = ''))==0) {
repeats<-1 ## To load the Plot_SV_NEW ?? (R) function
}else {
repeats<-0 ## To load the Plot_SV_NEW ?? (R) function
}

#source('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/script/Plot_SV.R', local = TRUE)
source(paste(script_folder, 'Plot_SV.R', sep=''), local = TRUE)
Plot_SV(genomes_r, g_lab, repmask, CDS_gDNA_blast, gDNAs_blast, N_Gap, Anno, output_flag,Gene,Ref_genome,haplotypes,repeats,strand_direction) ## plot in shiny
## 08/30/22
## 08/30/22


output$Download <- renderUI({

actionButton("Download", label = "Prepare To Download",style = "background-color:#FFFFFF")

})



#output$Save <- renderUI({
#downloadButton('Save',label = "Save compressed results to ...",style = "background-color:#C9DD03")
#  })  


 }, height = function() {length(input$id)*10+400} )  ## ## renderPlot 09/26/22


#output$GeneCapture <- renderUI({
#    actionButton("GeneCapture", label = "Involved Genes",style="color: #00F; background-color: #FFFF66; border-color: #c34113; border-radius: 10px; border-width: 2px")
#  })


output$done <- renderUI({
    actionButton("done", label = "DONE",style = "background-color:#FF6666")
  })


output$clustertree <- renderUI({
    actionButton("clustertree", label = "TREE",,style = "background-color:#3399FF")
  })



observeEvent(input$clustertree,{


#color_option <- c("blue","black","red","orange","grey","green")


val <- reactiveValues(clickx = NULL, clicky = NULL)
  
  observe({

    input$plot2_click

    isolate({
      val$clickx = c(val$clickx, input$plot2_click$x)
      val$clicky = c(val$clicky, input$plot2_click$y) 
      
    })

  }) #adding clicks to list


output$plot3 <- NULL #09/28/22
output$plot4 <- NULL #09/28/22
output$info <- NULL #09/28/22
output$submit2 <-NULL #09/28/22
output$info2 <- NULL
output$plot2 <- NULL  ## 09/29/22 
output$Haplotypes <- NULL #09/28/22
output$bucket <- NULL #09/28/22


output$plot2 <- renderPlot({

Gene <- input$Gene  ## from shiny input

gDNAs_blast <- read.table(paste(Dir, Gene, '_Haplotype-Self_out_m8_Crossover', sep = ''), sep = '\t', header = F, stringsAsFactors = F);   ## 10/10/22 crossover removed !
gDNAs_blast <- subset(gDNAs_blast, gDNAs_blast[,4]> 0)

genomes <- input$id ## from shiny input

####################### To check genomes ?
BlastSyn<-paste(Gene,'_Haplotype-Self_out_m8',sep = "");
BlastSynWorking<-read.table(BlastSyn,header=F);
Test_Genome0<-unique(BlastSynWorking$V2);

Test_Genome1 <- intersect(genomes,Test_Genome0);

if('query' %in% Test_Genome1){
 genomes <- Test_Genome1[which(Test_Genome1 != 'query')]  ## 10/05/22
} else {
 genomes <- Test_Genome1
}
####################### To check genomes ?
## 09/26
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


Title<-paste("Clustering on all haplotypes", "You can do tree cut on Height (y-axis) using single_click",sep='\n')

#as.dendrogram(h_c) %>% set("labels_cex", 0.9) %>% sort(type = "nodes") %>% highlight_branches %>% plot(main = "Clustering on all haplotypes",ylab = "Height",horiz = FALSE);
as.dendrogram(h_c) %>% set("labels_cex", 0.9) %>% sort(type = "nodes") %>% highlight_branches %>% plot(main = Title,ylab = "Height",horiz = FALSE);


color_option <- c("blue","black","red","orange","grey","green")
abline(a=NULL, b=NULL, val$clicky, col=color_option[ceiling(as.numeric(val$clicky)) %% 6 +1],lty = 2,lwd=3 );


observeEvent(input$plot2_click,{

output$plot3 <- NULL #09/28/22
output$plot4 <- NULL #09/28/22
output$info <- NULL #09/28/22
output$submit2 <-NULL #09/28/22


output$info2 <- NULL

output$info4 <- NULL
output$info3 <- NULL


    memb<-cutree(as.dendrogram(h_c), h=input$plot2_click$y)

    b_matrix0 <- cbind(b_matrix, cluster =as.data.frame(memb)) 

    b_matrix_groups <- b_matrix0[,'memb',drop=FALSE]
    b_matrix_groups1<-b_matrix_groups[order(b_matrix_groups$memb), , drop = FALSE]

 #   output$tablecluster <-DT::renderDataTable({
 #   datatable(b_matrix_groups1, extensions = 'Buttons',
 #             class="cell-border stripe",
 #             options = list(dom = "Blfrtip",
 #                            buttond = list("copy", list(extend = "collection",
 #                                                        buttons = c("csv"),
 #                                                        text = "Downloads")), pageLength=50, autoWidth = TRUE,
 #                            searchHighlight = TRUE, filter = "top"))
 #          }) # DT::renderDataTable   

   
   b_matrix_groups3 <- cbind(rownames(b_matrix_groups1), data.frame(b_matrix_groups1, row.names=NULL))

   colnames(b_matrix_groups3)<-c('genomes','memb')

  uniq_genome<-list()
  for(i in unique(b_matrix_groups3$memb)){
  uniq_genome[[i]]<-b_matrix_groups3[which(b_matrix_groups3$memb==i),][1,]
  }
  b_matrix_groups2 <- as.data.frame(rbindlist(uniq_genome))



  memb_count<-b_matrix_groups1 %>% group_by(memb) %>% summarize(count=n()) ## 09/26/22
  b_matrix_groups2<- cbind(b_matrix_groups2,memb_count$count) ## 09/26/22
  colnames(b_matrix_groups2)<-c('genomes_rep','main_clusters','haplotypes_rep') ## 09/26/22
  write.table(b_matrix_groups2,"b_matrix_groups2.txt",row.names=FALSE,col.names=TRUE,quote = FALSE,sep="\t") ## 09/26/22
#    output$tablecluster2 <-DT::renderDataTable({
#    datatable(b_matrix_groups2, extensions = 'Buttons',
#              class="cell-border stripe",
#              options = list(dom = "Blfrtip",
#                             buttond = list("copy", list(extend = "collection",
#                                                         buttons = c("csv"),
#                                                         text = "Downloads")), pageLength=50, autoWidth = TRUE,
#                             searchHighlight = TRUE, filter = "top"))
#           }) # DT::renderDataTable  
#info2_text<-paste('Cut tree based on your selected height ~', round(input$plot2_click$y,1),',', 'You can proceed ...',sep=' ')

info2_text<- paste0(c('Cut tree based on your selected height ~', round(input$plot2_click$y,1),',', 'with color:',color_option[ceiling(input$plot2_click$y) %% 6 +1],',', 'Click on "Plot selected haplotypes" to view haplotypes ...'), collapse= " ")

output$info2 <- renderText({info2_text})


  output$bucket <- renderUI({
    
    bucket_list(
      header = "Candidate haplotypes for plotting",
      group_name = "bucket_list_group",
      orientation = "horizontal",

    add_rank_list(
      text = "Drag from here",
      labels = b_matrix_groups2$genomes_rep,
    ),

      add_rank_list(text = "Order of plot",
          #          labels = b_matrix_groups2$genomes_rep, 
                    input_id = "list_2")
    )  
  })


output$Haplotypes <- renderUI({
    actionButton("Haplotypes", "Plot selected haplotypes",style = "background-color:#CCCCFF")
  })



 }) # input$plot2_click

  

  }, height = function() {length(input$id)*10+400}) ## renderplot2 09/26/22



 })   ## input$clustertree


observeEvent(input$Haplotypes,{

output$info <- NULL
output$plot4 <- NULL
output$submit2 <- NULL


info4_text<- paste('Left single_click and right double_click on top of figure to select preferred coordinates for Trimming ... ',sep='')
output$info4 <- renderText({info4_text})



color_option <- c("blue","black","red","orange","grey","green")

val1 <- reactiveValues(clickx = NULL, clicky = NULL)
val2 <- reactiveValues(clickx = NULL, clicky = NULL)

#vals1 <- reactiveValues(count = 1)
#vals2 <- reactiveValues(count = 1)

 observe({

    input$plot3_click
    input$plot3_dblclick

    isolate({

      val1$clickx = c(val1$clickx, input$plot3_click$x)
      val1$clicky = c(val1$clicky, input$plot3_click$y) 
      
      val2$clickx = c(val2$clickx, input$plot3_dblclick$x)
      val2$clicky = c(val2$clicky, input$plot3_dblclick$y)
      
      
   #  vals1$count =  1
   #   vals2$count =  length(c(input$plot3_dblclick$x))

   })

  }) #adding clicks to list


#recttext <- function(xl, yb, xr, yt, text, rectArgs = NULL, textArgs = NULL) {
#  center <- c(mean(c(xl, xr)), mean(c(yb, yt)))
#  do.call('rect', c(list(xleft = xl, ybottom = yb, xright = xr, ytop = yt), rectArgs))
#  do.call('text', c(list(x = center[1], y = center[2], labels = text), textArgs))
#}

##
##

Gene <- input$Gene  ## from shiny input
cds_ids <- input$Gene ## from shiny input
Ref_genome <- input$Pickgenome

var_types <- c('snp', 'ins', 'del');
indel_col <- c("grey", "red", "red");
cds_col <- c("yellowgreen", "brown");

#Dir <- paste('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/User', '/', User_folder0,'/', sep = '');

CDS_gDNA_blast <- read.table(paste(Dir, Gene, '_ref_mRNA-Haplotype_out_m8', sep = ''), sep = '\t', header = F, stringsAsFactors = F);
gDNAs_blast <- read.table(paste(Dir, Gene, '_Haplotype-Self_out_m8', sep = ''), sep = '\t', header = F, stringsAsFactors = F);
gDNAs_blast <- subset(gDNAs_blast, gDNAs_blast[,4]> 0)
N_Gap <- read.table(paste(Dir, Gene, '_Haplotype_N_Gaps', sep = ''), sep = '\t', header = T, stringsAsFactors = F);
Anno <- read.table(paste(Dir, Gene, '_Haplotype_anno', sep = ''), sep = '\t', header = T, stringsAsFactors = F)
#variations <- read.table(paste(Dir, Gene, '_Variations', sep = ''), sep = '\t', header = T, stringsAsFactors = F);

if (!file.size(paste(Dir, Gene, '_repMask2', sep = ''))==0) {
repmask <- read.table(paste(Dir, Gene, '_repMask2', sep = ''), header = F, sep = "\t", stringsAsFactors = F);
repmask <- repmask[abs(repmask[,8] - repmask[,7])> 100,]
}



output_flag = 0
##
##

x_lim <- range(gDNAs_blast[,9:10]) + c(0, 2000)



output$plot3 <- renderPlot({
##
Genome_order <- input$list_2 ## New haplotypes's order
genomes_r <- Genome_order
genomes <- Genome_order
n_g <- length(genomes)
g_lab <- genomes_r;
##

# b_matrix_groups2 <- read.table("b_matrix_groups2.txt",header=T)  # 09/26/22
#b_matrix_groups2 <- 1
haplotypes <- 1

if (!file.size(paste(Dir, Gene, '_repMask2', sep = ''))==0) {
repeats<-1 ## To load the Plot_SV_NEW ?? (R) function
}else {
repeats<-0 ## To load the Plot_SV_NEW ?? (R) function
}
####################################################
source(paste(script_folder, 'Plot_SV.R', sep=''), local = TRUE)
Plot_SV(genomes_r, g_lab, repmask, CDS_gDNA_blast, gDNAs_blast, N_Gap, Anno, output_flag,Gene,Ref_genome,haplotypes,repeats,strand_direction) ## plot in shiny

arrows(as.numeric(val1$clickx), as.numeric(val1$clicky)+0.5, as.numeric(val1$clickx), as.numeric(val1$clicky)+0.1,length = 0.25, lwd=3,col=color_option[1])
arrows(as.numeric(val2$clickx), as.numeric(val2$clicky)+0.5, as.numeric(val2$clickx), as.numeric(val2$clicky)+0.1,length = 0.25, lwd=3,col=color_option[3])

recttext <- function(xl, yb, text, rectArgs = NULL, textArgs = NULL) {
  center<-c(xl,yb)
  do.call('text', c(list(x = center[1], y = center[2], labels = text), textArgs))
}

recttext(as.numeric(val1$clickx), as.numeric(val1$clicky)+0.7, 'left', textArgs = list(col = 'blue', cex = 1.5))
recttext(as.numeric(val2$clickx), as.numeric(val2$clicky)+0.7, 'Right',textArgs = list(col = 'red', cex = 1.5))


}, height = function() {length(input$id)*10+300}) ## ## renderplot3 09/26/22




## 09/28/22
observeEvent(input$plot3_click,{


Genome_order <- input$list_2 ## New haplotypes's order
genomes_r <- Genome_order
genomes <- Genome_order

x_p_s <- c();
x_p_s[1] <- floor(input$plot3_click$x[1]);
n_g <- length(genomes_r);
for (g1 in 1:(n_g - 1)) {
 gDNAs <- subset(gDNAs_blast, gDNAs_blast[,1] == genomes_r[g1] & gDNAs_blast[,2] == genomes_r[g1 + 1])
 gDNAs <- subset(gDNAs, gDNAs[,7] <= x_p_s[g1] & gDNAs[, 8] >= x_p_s[g1]);
 x_p_s[g1 + 1] <- gDNAs[1,9] - gDNAs[1,7] + x_p_s[g1];
 # segments(x_p_s[g1], n_g - g1 - 0.01, x_p_s[g1 + 1], n_g - g1 - 0.9 - 0.1) ??
}
# output$info <- renderText({
#   print(x_p_s)
#  })

# output$table1 <- renderTable({
# x_p_s})

observeEvent(input$plot3_dblclick,{

x_p_ss <- c();
x_p_ss[1] <- floor(input$plot3_dblclick$x[1]);
n_g <- length(genomes_r);
for (g1 in 1:(n_g - 1)) {
 gDNAs <- subset(gDNAs_blast, gDNAs_blast[,1] == genomes_r[g1] & gDNAs_blast[,2] == genomes_r[g1 + 1])
 gDNAs <- subset(gDNAs, gDNAs[,7] <= x_p_ss[g1] & gDNAs[, 8] >= x_p_ss[g1]);
 x_p_ss[g1 + 1] <- gDNAs[1,9] - gDNAs[1,7] + x_p_ss[g1];
 # segments(x_p_s[g1], n_g - g1 - 0.01, x_p_s[g1 + 1], n_g - g1 - 0.9 - 0.1) ??
}
# output$info <- renderText({
#   print(x_p_s)
#  })

Combined<-cbind(x_p_s,x_p_ss);
colnames(Combined)<-c("start","end");
rownames(Combined)<-genomes;


Combined<-rbind(Combined, average =colMeans(Combined, na.rm=FALSE))  # 08/31/22


#output$table2 <-DT::renderDataTable({
#datatable(Combined, extensions = 'Buttons',
#              class="cell-border stripe",
#              options = list(dom = "Blfrtip",
#                             buttond = list("copy", list(extend = "collection",
#                                                         buttons = c("csv"),
#                                                         text = "Downloads")), pageLength=25, autoWidth = TRUE,
#                             searchHighlight = TRUE, filter = "top"))
#  
#  }) # DT::renderDataTable


write.table(Combined, file="Selected_lines_Coordinates.bed", sep="\t")

Combined <-read.table("Selected_lines_Coordinates.bed",header=T)

if(is.na(tail(Combined,1)[,1])) {

info_text<-paste('Please click/double click again ...',',','Coordinates not complete!',sep=' ')

} else if (is.na(tail(Combined,1)[,2])){

info_text<-paste('Please click/double click again ...',',','Coordinates not complete!',sep=' ')

} else {

info_text<-paste('left coordinate ~', Combined[1,1],';', 'right coordinate ~', Combined[1,2], ',','You can click on Trim Button ...',sep=' ')

}

output$info <- renderText({ info_text })


output$submit2 <- renderUI({


        if(!is.na(tail(Combined,1)[,1]) & !is.na(tail(Combined,1)[,2])) {
           
            if(Combined[1,1] < Combined[1,2])   {

              actionButton("submit2", label = "Trim",style = "background-color:#66FF66")
                     
                     }

                }       

    })



output$extract_fa <- renderUI({


        if(!is.na(tail(Combined,1)[,1]) & !is.na(tail(Combined,1)[,2])) {
           
            if(Combined[1,1] < Combined[1,2])   {

              actionButton("extract_fa", label = "Extract trimmed fasta",style = "background-color:#66FF66")
                     
                     }

                }       

    })

info3_text<- paste('You may reset arrows (selected coordinates) using button: "Plot selected haplotypes" ', sep="")
output$info3 <- renderText({info3_text})



})  ## plot3_dblclick

}) ## plot3_click



#Genome_order <- input$list_2
#g_rep_y<-list()

#for (g in 1:length(Genome_order)) {

# g_rep_y[[g]] <- length(Genome_order) - g - 0.025 + 0.01 ;

#}


#observeEvent(input$plot3_hover,{

#if (!file.size(paste(Dir, Gene, '_repMask2', sep = ''))==0) {

#Genome_order <- input$list_2
#repmask_hover <- read.table(paste(Dir, Gene, '_repMask2', sep = ''), header = F, sep = "\t", stringsAsFactors = F);
#repmask_hover <- repmask_hover[which(repmask_hover$V1 %in% Genome_order),c(1,2,7,8)]



# output$info_TE <- renderText({ paste(input$plot3_hover$x,input$plot3_hover$y) })




#  filtered_data <- reactive({
#    filter(data, name == input$x)
#  })




#  displayed_text <- reactive({
#    req(input$plot3_hover)
#    hover <- input$plot3_hover

#    dist <- sqrt((hover$x - filtered_data()$x)^2 + (hover$y - filtered_data()$y)^2)


#  })

#  output$hover_info <- renderPrint({
#    req(displayed_text())

#      cat("Name\n")
#      displayed_text()
#  })



#}   ## repmask exists ..


##  }) ## plot3_hover






})  ## input$Haplotypes



########## Start To Trim

observeEvent(input$submit2,{


Gene <- input$Gene ## from shiny input
cds_ids <- input$Gene ## from shiny input
Ref_genome <- input$Pickgenome

var_types <- c('snp', 'ins', 'del');
indel_col <- c("grey", "red", "red");
cds_col <- c("yellowgreen", "brown");

##
#Dir <- paste('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/User', '/', User_folder0,'/', sep = '');

# CoordinateFilter1<-read.table("MOC1_Haplotype-Self_out_m8",header=F)

#Name0<- input$id

CoordinateFilter0<-read.table("Selected_lines_Coordinates.bed",header=T) ## New selected coordinates!
CoordinateFilter0<-round(CoordinateFilter0,0)
CoordinateFilter0<-CoordinateFilter0[which(row.names(CoordinateFilter0)!='average'),] # remove average values

#Name0<- rownames(CoordinateFilter0)  ## 09/23/22

Name0<- input$list_2  ## 09/23/22

CDS_gDNA_blast <- read.table(paste(Dir, Gene, '_ref_mRNA-Haplotype_out_m8', sep = ''), sep = '\t', header = F, stringsAsFactors = F);
#### 10/10/22 filter for CDS_gDNA_blast
CoordinateFilter1<-CDS_gDNA_blast

CoordinateFilter1 <- CoordinateFilter1[which(CoordinateFilter1$V2 %in% Name0),]

CoordinateFilter2 <- CoordinateFilter1


CDS_gDNA_blast_New1 <- list()
index_coor1 <- 1

for(gname2 in unique(CoordinateFilter2$V2)){

CoordinateFilter3<- CoordinateFilter2[which(CoordinateFilter2$V2==gname2),]

#
for(i in 1:nrow(CoordinateFilter3)) {           
if(CoordinateFilter3[i,9]>CoordinateFilter3[i,10]){  
   temp10<- CoordinateFilter3[i,9]
   temp9<- CoordinateFilter3[i,10]
   CoordinateFilter3[i,10]<- temp10
   CoordinateFilter3[i,9]<- temp9
}
}
#

Target1<-CoordinateFilter0[which(rownames(CoordinateFilter0)==unique(CoordinateFilter3$V2)),]


CoordinateFilter3_1 <- CoordinateFilter3[which((CoordinateFilter3$V9>=Target1$start & CoordinateFilter3$V9<=Target1$end) | (CoordinateFilter3$V10>=Target1$start & CoordinateFilter3$V10<=Target1$end ) | (CoordinateFilter3$V9<=Target1$start & CoordinateFilter3$V10>=Target1$end ) ),]  ### filter

if(nrow(CoordinateFilter3_1)==0){

CoordinateFilter2_2 <-CoordinateFilter2

for(i in 1:nrow(CoordinateFilter2_2)){

CoordinateFilter2_2[i,9]<-0
CoordinateFilter2_2[i,10]<-1

}


CDS_gDNA_blast_New1[[index_coor1]] <- CoordinateFilter2_2

index_coor1 <- index_coor1+1

CDS_gDNA_blast <- as.data.frame(rbindlist(CDS_gDNA_blast_New1)) ## new input files


} else {

for(i in 1:nrow(CoordinateFilter3_1)){

 if(CoordinateFilter3_1[i,]$V9>=Target1$start & CoordinateFilter3_1[i,]$V9<=Target1$end & CoordinateFilter3_1[i,]$V10>=Target1$end) {
           
     CoordinateFilter3_1[i,]$V10 <- Target1$end

  }

  if(CoordinateFilter3_1[i,]$V10>=Target1$start & CoordinateFilter3_1[i,]$V10<=Target1$end & CoordinateFilter3_1[i,]$V9<=Target1$start) {
           
     CoordinateFilter3_1[i,]$V9 <- Target1$start

  }


  if(CoordinateFilter3_1[i,]$V9>=Target1$start & CoordinateFilter3_1[i,]$V10<=Target1$end) {
           
     Temp9 <- CoordinateFilter3_1[i,]$V9
     Temp10 <- CoordinateFilter3_1[i,]$V10
     
     CoordinateFilter3_1[i,]$V9 <- Temp9
     CoordinateFilter3_1[i,]$V10 <- Temp10

  }



    if(CoordinateFilter3_1[i,]$V9<=Target1$start & CoordinateFilter3_1[i,]$V10>=Target1$end) {
           
     Temp9 <- Target1$start
     Temp10 <- Target1$end
     
     CoordinateFilter3_1[i,]$V9 <- Temp9
     CoordinateFilter3_1[i,]$V10 <- Temp10

  }

}

CDS_gDNA_blast_New1[[index_coor1]] <- CoordinateFilter3_1

index_coor1 <- index_coor1+1

CDS_gDNA_blast <- as.data.frame(rbindlist(CDS_gDNA_blast_New1)) ## new input files

} # else
##



}

write.table(CDS_gDNA_blast, file=paste(Dir, Gene, '_ref_mRNA-Haplotype_out_m8_New', sep = ''), sep=" ", quote = FALSE,row.names = FALSE,col.names = FALSE)
#### 10/10/22 filter for CDS_gDNA_blast


gDNAs_blast <- read.table(paste(Dir, Gene, '_Haplotype-Self_out_m8', sep = ''), sep = '\t', header = F, stringsAsFactors = F);
gDNAs_blast <- subset(gDNAs_blast, gDNAs_blast[,4]> 0)
#### 10/10/22 filter for gDNAs_blast
CoordinateFilter0<-read.table("Selected_lines_Coordinates.bed",header=T) ## New selected coordinates!
CoordinateFilter0<-round(CoordinateFilter0,0)
CoordinateFilter0<-CoordinateFilter0[which(row.names(CoordinateFilter0)!='average'),] # remove average values

CoordinateFilter1 <- gDNAs_blast

#Name0<-rownames(CoordinateFilter0)
Name0<- input$list_2

CoordinateFilter1 <- CoordinateFilter1[which(CoordinateFilter1$V1 %in% Name0 & CoordinateFilter1$V2 %in% Name0),]

gDNAs_blast_New0 <- list()

index_coor <- 1

for(gname in rownames(CoordinateFilter0)) {

CoordinateFilter2 <- CoordinateFilter1[which(CoordinateFilter1$V1==gname & CoordinateFilter1$V2!=gname ),]


gDNAs_blast_New1 <- list()
index_coor1 <- 1

for(gname2 in unique(CoordinateFilter2$V2)){

CoordinateFilter3<- CoordinateFilter2[which(CoordinateFilter2$V2==gname2),]

#
for(i in 1:nrow(CoordinateFilter3)) {           
if(CoordinateFilter3[i,9]>CoordinateFilter3[i,10]){  
   temp10<- CoordinateFilter3[i,9]
   temp9<- CoordinateFilter3[i,10]
   CoordinateFilter3[i,10]<- temp10
   CoordinateFilter3[i,9]<- temp9
}
}
#

Target1<-CoordinateFilter0[which(rownames(CoordinateFilter0)==unique(CoordinateFilter3$V1)),]
Target2<-CoordinateFilter0[which(rownames(CoordinateFilter0)==unique(CoordinateFilter3$V2)),]

###
CoordinateFilter4 <- CoordinateFilter1[which(CoordinateFilter1$V1==gname & CoordinateFilter1$V2==gname ),]
CoordinateFilter4$V7<-Target1$start
CoordinateFilter4$V8<-Target1$end
CoordinateFilter4$V9<-Target1$start
CoordinateFilter4$V10<-Target1$end
###

CoordinateFilter3_1 <- CoordinateFilter3[which((CoordinateFilter3$V7>=Target1$start & CoordinateFilter3$V7<=Target1$end) | (CoordinateFilter3$V8>=Target1$start & CoordinateFilter3$V8<=Target1$end ) | (CoordinateFilter3$V7<=Target1$start & CoordinateFilter3$V8>=Target1$end ) ),]  ### filter

CoordinateFilter3_1 <- CoordinateFilter3_1[which((CoordinateFilter3_1$V9>=Target2$start & CoordinateFilter3_1$V9<=Target2$end) | (CoordinateFilter3_1$V10>=Target2$start & CoordinateFilter3_1$V10<=Target2$end ) | (CoordinateFilter3_1$V9<=Target2$start & CoordinateFilter3_1$V10>=Target2$end ) ),]  ### filter

##
for(i in 1:nrow(CoordinateFilter3_1)){

  if(CoordinateFilter3_1[i,]$V7>=Target1$start & CoordinateFilter3_1[i,]$V7<=Target1$end & CoordinateFilter3_1[i,]$V8>=Target1$end) {
           
     CoordinateFilter3_1[i,]$V8 <- Target1$end

  }

  if(CoordinateFilter3_1[i,]$V8>=Target1$start & CoordinateFilter3_1[i,]$V8<=Target1$end & CoordinateFilter3_1[i,]$V7<=Target1$start) {
           
     CoordinateFilter3_1[i,]$V7 <- Target1$start

  }
  
  if(CoordinateFilter3_1[i,]$V7>=Target1$start & CoordinateFilter3_1[i,]$V8<=Target1$end) {
           
     Temp7 <- CoordinateFilter3_1[i,]$V7
     Temp8 <- CoordinateFilter3_1[i,]$V8
     
     CoordinateFilter3_1[i,]$V7 <- Temp7
     CoordinateFilter3_1[i,]$V8 <- Temp8

  }


  if(CoordinateFilter3_1[i,]$V7<=Target1$start & CoordinateFilter3_1[i,]$V8>=Target1$end) {
           
     Temp7 <- Target1$start
     Temp8 <- Target1$end
     
     CoordinateFilter3_1[i,]$V7 <- Temp7
     CoordinateFilter3_1[i,]$V8 <- Temp8

  }

##

  if(CoordinateFilter3_1[i,]$V9>=Target2$start & CoordinateFilter3_1[i,]$V9<=Target2$end & CoordinateFilter3_1[i,]$V10>=Target2$end) {
           
     CoordinateFilter3_1[i,]$V10 <- Target2$end

  }

  if(CoordinateFilter3_1[i,]$V10>=Target2$start & CoordinateFilter3_1[i,]$V10<=Target2$end & CoordinateFilter3_1[i,]$V9<=Target2$start) {
           
     CoordinateFilter3_1[i,]$V9 <- Target2$start

  }


  if(CoordinateFilter3_1[i,]$V9>=Target2$start & CoordinateFilter3_1[i,]$V10<=Target2$end) {
           
     Temp9 <- CoordinateFilter3_1[i,]$V9
     Temp10 <- CoordinateFilter3_1[i,]$V10
     
     CoordinateFilter3_1[i,]$V9 <- Temp9
     CoordinateFilter3_1[i,]$V10 <- Temp10

  }



    if(CoordinateFilter3_1[i,]$V9<=Target2$start & CoordinateFilter3_1[i,]$V10>=Target2$end) {
           
     Temp9 <- Target2$start
     Temp10 <- Target2$end
     
     CoordinateFilter3_1[i,]$V9 <- Temp9
     CoordinateFilter3_1[i,]$V10 <- Temp10

  }



}
##

gDNAs_blast_New1[[index_coor1]] <- CoordinateFilter3_1

index_coor1 <- index_coor1+1

}

CoordinateFilter4_1 <- as.data.frame(rbindlist(gDNAs_blast_New1))



CoordinateFilter5_1<-rbind(CoordinateFilter4,CoordinateFilter4_1)

gDNAs_blast_New0[[index_coor]]<- CoordinateFilter5_1

index_coor <- index_coor+1

}

gDNAs_blast <- as.data.frame(rbindlist(gDNAs_blast_New0)) ## new input files
##

write.table(gDNAs_blast, file=paste(Dir, Gene, '_Haplotype-Self_out_m8_New', sep = ''), sep=" ", quote = FALSE,row.names = FALSE,col.names = FALSE)

#### 10/10/22 filter for gDNAs_blast



N_Gap <- read.table(paste(Dir, Gene, '_Haplotype_N_Gaps', sep = ''), sep = '\t', header = T, stringsAsFactors = F);
#### 10/10/22 N_Gap filter based on Selected_lines_Coordinates.bed
CoordinateFilter1<-N_Gap

CoordinateFilter0<-read.table("Selected_lines_Coordinates.bed",header=T) ## New selected coordinates!
CoordinateFilter0<-round(CoordinateFilter0,0)
CoordinateFilter0<-CoordinateFilter0[which(row.names(CoordinateFilter0)!='average'),] # remove average values

#Name0<-rownames(CoordinateFilter0)
Name0<- input$list_2


if(length(intersect(Name0,CoordinateFilter1$Genome)) == 0){

write.table(N_Gap, file=paste(Dir, Gene, '_Haplotype_N_Gaps_New', sep = ''), sep="\t", quote = FALSE,row.names = FALSE,col.names = TRUE)

N_Gap <- read.table(paste(Dir, Gene, '_Haplotype_N_Gaps_New', sep = ''), sep = '\t', header = T, stringsAsFactors = F);

}

if(length(intersect(Name0,CoordinateFilter1$Genome)) != 0){


CoordinateFilter1 <- CoordinateFilter1[which(CoordinateFilter1$Genome %in% Name0),]
CoordinateFilter2 <- CoordinateFilter1


N_Gap_New1 <- list()
index_coor1 <- 1

for(gname2 in unique(CoordinateFilter2$Genome)){

CoordinateFilter3<- CoordinateFilter2[which(CoordinateFilter2$Genome==gname2),]

for(i in 1:nrow(CoordinateFilter3)) {           
if(CoordinateFilter3[i,2]>CoordinateFilter3[i,3]){  
   temp10<- CoordinateFilter3[i,2]
   temp9<- CoordinateFilter3[i,3]
   CoordinateFilter3[i,3]<- temp10
   CoordinateFilter3[i,2]<- temp9
}
}


Target1<-CoordinateFilter0[which(rownames(CoordinateFilter0)==unique(CoordinateFilter3$Genome)),]


N_Gaps1 <- matrix(nrow=nrow(CoordinateFilter3), ncol=3, byrow=TRUE)

for(i in 1:nrow(CoordinateFilter3)){

if(CoordinateFilter3[i,]$GAP_Start>Target1$end | CoordinateFilter3[i,]$GAP_End<Target1$start){

N_Gaps1[i,1]<-gname2
N_Gaps1[i,2]<-0
N_Gaps1[i,3]<-1

}

if(CoordinateFilter3[i,]$GAP_Start<Target1$start & CoordinateFilter3[i,]$GAP_End>Target1$end){

N_Gaps1[i,1]<-gname2
N_Gaps1[i,2]<-Target1$start
N_Gaps1[i,3]<-Target1$end

}

if(CoordinateFilter3[i,]$GAP_Start>Target1$start & CoordinateFilter3[i,]$GAP_End<Target1$end){

N_Gaps1[i,1]<-gname2
N_Gaps1[i,2]<-CoordinateFilter3[i,]$GAP_Start
N_Gaps1[i,3]<-CoordinateFilter3[i,]$GAP_End

}

if(CoordinateFilter3[i,]$GAP_Start>Target1$start & CoordinateFilter3[i,]$GAP_Start<Target1$end & CoordinateFilter3[i,]$GAP_End>Target1$end ){

N_Gaps1[i,1]<-gname2
N_Gaps1[i,2]<-CoordinateFilter3[i,]$GAP_Start
N_Gaps1[i,3]<-Target1$end

}


if(CoordinateFilter3[i,]$GAP_End>Target1$start & CoordinateFilter3[i,]$GAP_End<Target1$end & CoordinateFilter3[i,]$GAP_Start<Target1$start){

N_Gaps1[i,1]<-gname2
N_Gaps1[i,2]<-Target1$start
N_Gaps1[i,3]<-CoordinateFilter3[i,]$GAP_End

}


}

colnames(N_Gaps1)<-c("Genome","GAP_Start","GAP_End")


N_Gap2 <- as.data.frame(N_Gaps1) ## new input files


N_Gap_New1[[index_coor1]] <- N_Gap2

index_coor1 <- index_coor1+1

}

N_Gap_New2 <- as.data.frame(rbindlist(N_Gap_New1)) ## new input files
write.table(N_Gap_New2, file=paste(Dir, Gene, '_Haplotype_N_Gaps_New', sep = ''), sep="\t", quote = FALSE,row.names = FALSE,col.names = TRUE)

N_Gap <- read.table(paste(Dir, Gene, '_Haplotype_N_Gaps_New', sep = ''), sep = '\t', header = T, stringsAsFactors = F);

}
#### 10/10/22 N_Gap filter based on Selected_lines_Coordinates.bed



#### repeats

if (!file.size(paste(Dir, Gene, '_repMask2', sep = ''))==0) {
repmask <- read.table(paste(Dir, Gene, '_repMask2', sep = ''), header = F, sep = "\t", stringsAsFactors = F);
repmask <- repmask[abs(repmask[,8] - repmask[,7])> 100,]


CoordinateFilter1<-repmask


CoordinateFilter0<-read.table("Selected_lines_Coordinates.bed",header=T) ## New selected coordinates!
CoordinateFilter0<-round(CoordinateFilter0,0)
CoordinateFilter0<-CoordinateFilter0[which(row.names(CoordinateFilter0)!='average'),] # remove average values

#Name0<-rownames(CoordinateFilter0)
Name0<- input$list_2

CoordinateFilter1 <- CoordinateFilter1[which(CoordinateFilter1$V1 %in% Name0),]

CoordinateFilter2 <- CoordinateFilter1



repmask_New1 <- list()
index_coor1 <- 1

for(gname2 in unique(CoordinateFilter2$V1)){

CoordinateFilter3<- CoordinateFilter2[which(CoordinateFilter2$V1==gname2),]

#
for(i in 1:nrow(CoordinateFilter3)) {           
if(CoordinateFilter3[i,7]>CoordinateFilter3[i,8]){  
   temp10<- CoordinateFilter3[i,7]
   temp9<- CoordinateFilter3[i,8]
   CoordinateFilter3[i,8]<- temp10
   CoordinateFilter3[i,7]<- temp9
}
}
#

Target1<-CoordinateFilter0[which(rownames(CoordinateFilter0)==unique(CoordinateFilter3$V1)),]


CoordinateFilter3_1 <- CoordinateFilter3[which((CoordinateFilter3$V7>=Target1$start & CoordinateFilter3$V7<=Target1$end) | (CoordinateFilter3$V8>=Target1$start & CoordinateFilter3$V8<=Target1$end ) | (CoordinateFilter3$V7<=Target1$start & CoordinateFilter3$V8>=Target1$end ) ),]  ### filter

if(nrow(CoordinateFilter3_1)==0){

CoordinateFilter2_2 <-CoordinateFilter2

for(i in 1:nrow(CoordinateFilter2_2)){

CoordinateFilter2_2[i,7]<-0
CoordinateFilter2_2[i,8]<-1

}


repmask_New1[[index_coor1]] <- CoordinateFilter2_2

index_coor1 <- index_coor1+1

repmask <- as.data.frame(rbindlist(repmask_New1)) ## new input files


} else {

for(i in 1:nrow(CoordinateFilter3_1)){

 if(CoordinateFilter3_1[i,]$V7>=Target1$start & CoordinateFilter3_1[i,]$V7<=Target1$end & CoordinateFilter3_1[i,]$V8>=Target1$end) {
           
     CoordinateFilter3_1[i,]$V8 <- Target1$end

  }

  if(CoordinateFilter3_1[i,]$V8>=Target1$start & CoordinateFilter3_1[i,]$V8<=Target1$end & CoordinateFilter3_1[i,]$V7<=Target1$start) {
           
     CoordinateFilter3_1[i,]$V7 <- Target1$start

  }


  if(CoordinateFilter3_1[i,]$V7>=Target1$start & CoordinateFilter3_1[i,]$V8<=Target1$end) {
           
     Temp7 <- CoordinateFilter3_1[i,]$V7
     Temp8 <- CoordinateFilter3_1[i,]$V8
     
     CoordinateFilter3_1[i,]$V7 <- Temp7
     CoordinateFilter3_1[i,]$V8 <- Temp8

  }



    if(CoordinateFilter3_1[i,]$V7<=Target1$start & CoordinateFilter3_1[i,]$V8>=Target1$end) {
           
     Temp7 <- Target1$start
     Temp8 <- Target1$end
     
     CoordinateFilter3_1[i,]$V7 <- Temp7
     CoordinateFilter3_1[i,]$V8 <- Temp8

  }

}

repmask_New1[[index_coor1]] <- CoordinateFilter3_1

index_coor1 <- index_coor1+1

repmask <- as.data.frame(rbindlist(repmask_New1)) ## new input files

} # else
##



}

##

write.table(repmask, file=paste(Dir, Gene, '_repMask2_New', sep = ''), sep=" ", quote = FALSE,row.names = FALSE,col.names = FALSE)

#### 10/10/22 repeat filter based on Selected_lines_Coordinates.bed

}


####
####
Anno <- read.table(paste(Dir, Gene, '_Haplotype_anno', sep = ''), sep = '\t', header = T, stringsAsFactors = F)
#variations <- read.table(paste(Dir, Gene, '_Variations', sep = ''), sep = '\t', header = T, stringsAsFactors = F);

N_Gap <- read.table(paste(Dir, Gene, '_Haplotype_N_Gaps_New', sep = ''), sep = '\t', header = T, stringsAsFactors = F);

gDNAs_blast <- read.table(paste(Dir, Gene, '_Haplotype-Self_out_m8_New', sep = ''),header=F)

output_flag = 0

b_matrix_groups2 <- read.table("b_matrix_groups2.txt",header=T)  # 09/26/22

genomes <- input$list_2

genomes_r <- genomes

n_g <- length(genomes)

g_lab <- genomes_r;
##################### 09/07/22

x_lim <- range(CoordinateFilter0)+c(0,2000)

haplotypes<-1

output$plot4 <- renderPlot({

if (!file.size(paste(Dir, Gene, '_repMask2', sep = ''))==0) {

repeats<-1 ## To load the Plot_SV_NEW ?? (R) function
}else {
repeats<-0 ## To load the Plot_SV_NEW ?? (R) function
}
####################################################

#source('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/script/Plot_SV.R', local = TRUE)
source(paste(script_folder, 'Plot_SV.R', sep=''), local = TRUE)
Plot_SV(genomes_r, g_lab, repmask, CDS_gDNA_blast, gDNAs_blast, N_Gap, Anno, output_flag,Gene,Ref_genome,haplotypes,repeats,strand_direction) ## plot in shiny



}) 


})  ## input$submit2  Trim


########################## A table showing clustering results
if(length(Information_list)!=0){

output$table1 <-DT::renderDataTable({
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
Filtered_HaplotypeSyn_Plotted <- Filtered_HaplotypeSyn[which(Filtered_HaplotypeSyn$Genome!=''),]

output$table2 <-DT::renderDataTable({
datatable(Filtered_HaplotypeSyn_Plotted,caption = htmltools::tags$caption(
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
BlastSynWorking_0<-read.table(Blast_Ori,header=T); ## Blast_Original
BlastSynWorking_1 <- BlastSynWorking_0[which(BlastSynWorking_0$Genome!=''),]
BlastSynWorking_2 <- BlastSynWorking[which(BlastSynWorking$Genome!=''),]

NotShown0 <- anti_join(BlastSynWorking_1,BlastSynWorking_2) # Not shown in plot, other genomic positions.

if(nrow(NotShown0)!=0){

output$table3 <-DT::renderDataTable({
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
}
########################## A table showing Blast result which is not presented in main plot ??




## Based on selected haplotypes and region
observeEvent(input$extract_fa,{

    progress <- Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    progress$set(message = 'Extract trimmed fasta ...',
                 detail = 'Almost done...')
    for (i in 1:10) {
      progress$set(value = i)
      Sys.sleep(0.1)  ## ??
    }

Gene <- input$Gene ## from shiny input

CoordinateFilter0<-read.table("Selected_lines_Coordinates.bed",header=T) ## New selected coordinates!
CoordinateFilter0<-round(CoordinateFilter0,0)
Sel_Hap_Coor <- CoordinateFilter0[which(row.names(CoordinateFilter0)!='average'),] # remove average values

query_extract_fa<-list();

dna_Haplotype_fa <- readDNAStringSet(paste(Dir, Gene, '_Haplotype.fa', sep = ''));


for (x_lines in 1:nrow(Sel_Hap_Coor)) {
  query_name <- row.names(Sel_Hap_Coor[x_lines,]); 
  start_1 <- Sel_Hap_Coor[x_lines,][,1]; 
  end_1 <- Sel_Hap_Coor[x_lines,][,2];
  query_extract <- dna_Haplotype_fa[grepl(query_name, dna_Haplotype_fa@ranges@NAMES)];
  query_extract_fa[[x_lines]] <- subseq(query_extract, start=start_1, end=end_1);
  writeXStringSet(query_extract_fa[[x_lines]], file=paste0(x_lines, '_Selected_New.fa'));
  }

system_fasta0<-paste("cat *_Selected_New.fa >", paste(Gene,'_User_Selected.fa',sep=''), sep=' ')
system(system_fasta0);
system("rm *_Selected_New.fa");

  }) # observeEvent input$extract_fa

## Based on selected haplotypes and region



####### hover function to get the coordinates

# output$showing_TE <- renderText({

#    xy_str <- function(e) {
#      if(is.null(e)) return("NULL\n")
     # paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
#       paste0("Seq coordinate=", round(e$x, 0), " y=", NA, "\n")
#    }

#    paste0(
#      
#      "hover: ", xy_str(input$plot_hover)
#      
#    )
#  })












}  ## function


