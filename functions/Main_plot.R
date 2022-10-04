Main_plot <- function() {

system("cat extract_syn_fa_Wheat_1_05_17.pl > extract_syn_fa_Wheat_1_05_17_Test4_copy_3.pl")


Upstream_update <- paste('p_size=',input$Upstream,sep = "");  ## from shiny input
Downstream_update <- paste('a_size=',input$Downstream,sep = "");  ## from shiny input
replace_pattern_in(file_contents="p_size=10000",replace=Upstream_update,file_pattern="Test4_copy_3",file_contents_perl = TRUE);
replace_pattern_in(file_contents="a_size=1000",replace=Downstream_update,file_pattern="Test4_copy_3",file_contents_perl = TRUE);

Chr_update0 <- paste('target_ch=',input$Chr,sep = "'");  ## from shiny input
symbol<-c("'");
Chr_update <- paste(Chr_update0,symbol,sep = "");
replace_pattern_in(file_contents="target_ch='chr7D'",replace=Chr_update,file_pattern="Test4_copy_3",file_contents_perl = TRUE);

Gene_update0 <- paste('gene=',input$Gene,sep = "'");  ## from shiny input
Gene_update <- paste(Gene_update0,symbol,sep = "");
replace_pattern_in(file_contents="gene='TraesCS7D02G524200'",replace=Gene_update,file_pattern="Test4_copy_3",file_contents_perl = TRUE);


##############
system("perl extract_syn_fa_Wheat_1_05_17_Test4_copy_3.pl 1");


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


system("perl extract_syn_fa_Wheat_1_05_17_Test4_copy_3.pl 2");

# system("rm extract_syn_fa_Wheat_1_05_17_Test4_copy_3.pl");


result_folder1<-c('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Wheat/IWGSC/Candidate_genes/')
result_folder2<-paste(result_folder1,Gene,sep="")
new.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny"
setwd(result_folder2)
result_files<-list.files(result_folder2, pattern = Gene)
file.copy(result_files, new.folder)


setwd(new.folder)


####################################################


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


## 08/30/22
## 08/30/22
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
## 08/30/22
## 08/30/22



## 08/30/22 PM

### Click to get start and end coordinates ??

observeEvent(input$plot_click,{
x_p_s <- c();
x_p_s[1] <- floor(input$plot_click$x[1]);
n_g <- length(genomes_r);
for (g1 in 1:(n_g - 1)) {
 gDNAs <- subset(gDNAs_blast, gDNAs_blast[,1] == genomes_r[g1] & gDNAs_blast[,2] == genomes_r[g1 + 1])
 gDNAs <- subset(gDNAs, gDNAs[,7] <= x_p_s[g1] & gDNAs[, 8] >= x_p_s[g1]);
 x_p_s[g1 + 1] <- gDNAs[1,9] - gDNAs[1,7] + x_p_s[g1];

 # segments(x_p_s[g1], n_g - g1 - 0.01, x_p_s[g1 + 1], n_g - g1 - 0.9 - 0.1) ## ?? 08/30/22 PM
}
# output$info <- renderText({
#   print(x_p_s)
#  })

# output$table1 <- renderTable({
# x_p_s})

observeEvent(input$plot_dblclick,{
x_p_ss <- c();
x_p_ss[1] <- floor(input$plot_dblclick$x[1]);
n_g <- length(genomes_r);
for (g1 in 1:(n_g - 1)) {
 gDNAs <- subset(gDNAs_blast, gDNAs_blast[,1] == genomes_r[g1] & gDNAs_blast[,2] == genomes_r[g1 + 1])
 gDNAs <- subset(gDNAs, gDNAs[,7] <= x_p_ss[g1] & gDNAs[, 8] >= x_p_ss[g1]);
 x_p_ss[g1 + 1] <- gDNAs[1,9] - gDNAs[1,7] + x_p_ss[g1];

 # segments(x_p_s[g1], n_g - g1 - 0.01, x_p_s[g1 + 1], n_g - g1 - 0.9 - 0.1) ## ?? 08/30/22 PM
}
# output$info <- renderText({
#   print(x_p_s)
#  })

Combined<-cbind(x_p_s,x_p_ss);
colnames(Combined)<-c("start","end");
rownames(Combined)<-genomes;

Combined<-rbind(Combined, average =colMeans(Combined, na.rm=FALSE))

#output$tableclick <-DT::renderDataTable({
#datatable(Combined, extensions = 'Buttons',
#              class="cell-border stripe",
#              options = list(dom = "Blfrtip",
#                             buttond = list("copy", list(extend = "collection",
#                                                         buttons = c("csv"),
#                                                         text = "Downloads")), pageLength=25, autoWidth = TRUE,
#                             searchHighlight = TRUE, filter = "top")) 
#  }) # DT::renderDataTable

write.table(Combined, file="Selected_lines_Coordinates.bed", sep="\t")

Combined <-read.table("Selected_lines_Coordinates.bed",header=T)

if(is.na(tail(Combined,1)[,1])) {

info_text<-paste('Please click/double click again ...',',','Coordinates not complete!',sep=' ')

} else if (is.na(tail(Combined,1)[,2])){

info_text<-paste('Please click/double click again ...',',','Coordinates not complete!',sep=' ')

} else {

info_text<-paste('left coordinate ~', Combined[1,1],';', 'right coordinate ~', Combined[1,2], ',','You can click on Zoom in ...',sep=' ')

}

output$info1 <- renderText({ info_text })

})  ## plot_dblclick


}) ## plot_click

### Click to get start and end coordinates ??

## 08/30/22 PM





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


observeEvent(input$submit2,{


CoordinateFilter0<-read.table("Selected_lines_Coordinates.bed",header=T) ## New selected coordinates!


CoordinateFilter0<-round(CoordinateFilter0,0)
CoordinateFilter0<-CoordinateFilter0[which(row.names(CoordinateFilter0)!='average'),] # remove average values

##
CoordinateFilter1<-read.table(paste(Dir, Gene, '_Haplotype-Self_out_m8', sep = ''), sep = '\t', header = F, stringsAsFactors = F);



Name0<- input$id

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

CoordinateFilter3_1 <- CoordinateFilter3[which(CoordinateFilter3$V8>=Target1$start | CoordinateFilter3$V10>=Target2$start ),]  ### filter

##
for(i in 1:nrow(CoordinateFilter3_1)){

  if(CoordinateFilter3_1[i,]$V7<Target1$start) {
           
     temp_0 <- Target1$start
     CoordinateFilter3_1[i,]$V7 <- temp_0

  }
  if(CoordinateFilter3_1[i,]$V8>Target1$end) {
           
     temp_1 <- Target1$end
     CoordinateFilter3_1[i,]$V8 <- temp_1

  }
  
  if(CoordinateFilter3_1[i,]$V9<Target2$start) {
           
     temp_00 <- Target2$start
     CoordinateFilter3_1[i,]$V9 <- temp_00

  }
  if(CoordinateFilter3_1[i,]$V10>Target2$end) {
           
     temp_11 <- Target2$end
     CoordinateFilter3_1[i,]$V10 <- temp_11

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

write.table(gDNAs_blast, file="gDNAs_blast_NEW", sep=" ", quote = FALSE,row.names = FALSE,col.names = FALSE)

##################### 09/07/22
genomes <- unique(gDNAs_blast[,2]);

n_g <- length(genomes)


g_lab <- genomes_r;
##################### 09/07/22

x_lim <- range(CoordinateFilter0)

output$plot2 <- renderPlot({




if (!file.size(paste(Dir, Gene, '_repMask2', sep = ''))==0) {
source('Source_Modified.R', local = TRUE) ## To load the Plot_SV_NEW ?? (R) function
}else {
source('Source_Modified_NORepeat.R', local = TRUE) ## To load the Plot_SV_NEW ?? (R) function
}
Plot_SV(genomes_r, g_lab, repmask, CDS_gDNA_blast, gDNAs_blast, N_Gap, Anno, variations, output_flag) ## plot in shiny


})


})
####### 09/01/22
####### 09/01/22
####### 09/01/22
####### 09/01/22

#observeEvent(input$submit3_2,{

output$plot1 <- renderPlot({

gDNAs_blast<- read.table(paste(Dir, Gene, '_Haplotype-Self_out_m8', sep = ''), sep = '\t', header = F, stringsAsFactors = F);
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

b_matrix<-b_matrix*10 ## ???? 08/30/22 !!!

h_c <- hclust(dist(b_matrix));


if(length(get_nodes_attr(as.dendrogram(h_c), "height"))==length(which(get_nodes_attr(as.dendrogram(h_c), "height")==0))) {

as.dendrogram(h_c) %>% set("labels_cex", 0.8) %>% sort(type = "nodes") %>% plot(main = "Clustering on all haplotypes",horiz = TRUE);

} else {

as.dendrogram(h_c) %>% set("labels_cex", 0.8) %>% sort(type = "nodes") %>% highlight_branches %>% plot(main = "Clustering on all haplotypes",horiz = TRUE);

}
######## 08/30/22


observeEvent(input$plot1_click,{


    memb<-cutree(as.dendrogram(h_c), h=input$plot1_click$x)

    b_matrix0 <- cbind(b_matrix, cluster =as.data.frame(memb)) 

    b_matrix_groups <- b_matrix0[,'memb',drop=FALSE]
    b_matrix_groups1<-b_matrix_groups[order(b_matrix_groups$memb), , drop = FALSE]

    #output$tablecluster <-DT::renderDataTable({
    #datatable(b_matrix_groups1, extensions = 'Buttons',
    #          class="cell-border stripe",
    #          options = list(dom = "Blfrtip",
    #                         buttond = list("copy", list(extend = "collection",
    #                                                     buttons = c("csv"),
    #                                                     text = "Downloads")), pageLength=50, autoWidth = TRUE,
    #                         searchHighlight = TRUE, filter = "top"))
    #       }) # DT::renderDataTable   

   
   b_matrix_groups3 <- cbind(rownames(b_matrix_groups1), data.frame(b_matrix_groups1, row.names=NULL))
   colnames(b_matrix_groups3)<-c('genomes','memb')

  uniq_genome<-list()
  for(i in unique(b_matrix_groups3$memb)){
  uniq_genome[[i]]<-b_matrix_groups3[which(b_matrix_groups3$memb==i),][1,]
  }
  b_matrix_groups2 <- as.data.frame(rbindlist(uniq_genome))

  #  output$tablecluster2 <-DT::renderDataTable({
  #  datatable(b_matrix_groups2, extensions = 'Buttons',
  #            class="cell-border stripe",
  #            options = list(dom = "Blfrtip",
  #                           buttond = list("copy", list(extend = "collection",
  #                                                       buttons = c("csv"),
  #                                                       text = "Downloads")), pageLength=50, autoWidth = TRUE,
  #                           searchHighlight = TRUE, filter = "top"))
  #         }) # DT::renderDataTable  


info2_text<-paste('Cut tree based on your selected height ~', round(input$plot1_click$x,1),',', 'You can proceed ...',sep=' ')

output$info2 <- renderText({ info2_text })



  output$bucket <- renderUI({
    bucket_list(
      header = "Candidate haplotypes for plotting",
      group_name = "bucket_list_group",
      orientation = "horizontal",
      add_rank_list(text = "Drag from here",
                    labels = b_matrix_groups2$genomes, 
                    input_id = "list_1"),

      add_rank_list(text = "to here (as order of plot)",
                   # labels = NULL,
                    labels = b_matrix_groups2$genomes, 
                    input_id = "list_2")
    )  
  })



 }) # input$plot1_click
######## 08/30/22

  
  }) ## renderPlot1



observeEvent(input$Haplotypes,{

#CoordinateFilter0<-read.table("Selected_lines_Coordinates.bed",header=T) ## New selected coordinates!
#CoordinateFilter0<-round(CoordinateFilter0,0)
#CoordinateFilter0<-CoordinateFilter0[which(row.names(CoordinateFilter0)!='average'),] # remove average values
#gDNAs_blast <- read.table("gDNAs_blast_NEW",header=F)


gDNAs_blast<- read.table(paste(Dir, Gene, '_Haplotype-Self_out_m8', sep = ''), sep = '\t', header = F, stringsAsFactors = F);
gDNAs_blast<-gDNAs_blast[which(gDNAs_blast[,1]!='IWGSCNEW' & gDNAs_blast[,2]!='IWGSCNEW'), ]

genomes <- unique(gDNAs_blast[,2]);
x_lim <- range(gDNAs_blast[,9:10]) + c(0, 2000);

#x_lim <- range(CoordinateFilter0)


Genome_order <- input$list_2 ## New haplotypes's order


genomes_r <- Genome_order
genomes <- Genome_order
n_g <- length(genomes)
g_lab <- genomes_r;


output$plot5 <- renderPlot({


if (!file.size(paste(Dir, Gene, '_repMask2', sep = ''))==0) {
source('Source_Modified.R', local = TRUE) ## To load the Plot_SV_NEW ?? (R) function
}else {
source('Source_Modified_NORepeat.R', local = TRUE) ## To load the Plot_SV_NEW ?? (R) function
}
Plot_SV(genomes_r, g_lab, repmask, CDS_gDNA_blast, gDNAs_blast, N_Gap, Anno, variations, output_flag) ## plot in shiny
})


})








  }) ## render Main Plot. output$plot




}


