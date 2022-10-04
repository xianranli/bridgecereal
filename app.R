library(shiny)
library(Biostrings)
library(hutils)
library(shinyjs)
library(brochure)
library(shinyWidgets)
library(DT)
library(data.table)
library(dplyr)
library(msa)
library(seqinr)
library(ape)
library(dendextend)

library(shinyBS)
library(sortable)

############################################################ Creating a navlink 

nav_links <- tags$ul(

flowLayout(

  tags$li(
    tags$a(href = "/TestShiny/page0", "Main"), 
  ),
  tags$li(
    tags$a(href = "/TestShiny", "Wheat"), 
  ),
  tags$li(
    tags$a(href = "/TestShiny/page2", "Barley"), 
  ),
   tags$li(
    tags$a(href = "/TestShiny/page3", "Sorghum"), 
  ),
   tags$li(
    tags$a(href = "/TestShiny/page4", "Maize"), 
  ),
   tags$li(
    tags$a(href = "/TestShiny/page5", "Rice"), 
  ),
   tags$li(
    tags$a(href = "/TestShiny/page6", "WheatTest"),

  ),

tags$style(
"li a {font-size: 20px;font-weight: bold;}",
)

)


)


############ Page 1

page_1 <- function(){

  page(
    href = "/",


    ui <- function(request){
      tagList(
        h1("This is wheat page",style="text-align:center"),
        nav_links,


# To add ui part for page1   

useShinyjs(),

column(8,offset = 5, titlePanel("Wheat Alignment output")), # titlePanel

sidebarLayout(

sidebarPanel(

column(12,textInput("Gene","Gene Name (TraesCS4A02G464500 or others)")),

# column(6,textInput("Gene_ID","Sobic")),

column(12,textInput("Chr","Chromosome (chr4A or others)")),

column(12,textInput("Upstream","Upstream (bp)")),
column(12,textInput("Downstream","Downstream (bp)")),

column(12,
pickerInput(
  inputId = "id", 
  label = "Genomes (Select IWGSC and others) :", 
  choices = c('jagger','julius','lancer','landmark','mace','mattis','norin61','stanley','IWGSC','spelta','arinalrfor'),

  selected = c('jagger','julius','lancer','landmark','mace','mattis','norin61','stanley','IWGSC','spelta','arinalrfor'), ## by default

  options = list(
    'actions-box' = TRUE, 
    size = 18,
    'selected-text-format' = "count > 2"
  ), 
  multiple = TRUE,
)
),

# fileInput("file1", "Upload a fasta file"),
#column(12,sliderInput("Similarityfilter", "Similarity filter (96-100) :", min = 96, max = 100, value = 96,step=0.5 )),
column(12,sliderInput("Distancefilter", "Distance filter between mapped clusters (1kb-50kb) :", min = 1000, max = 50000, value =20000)),
column(12,sliderInput("CDSfilter", "Expected CDS size compared to IWGSC (fold change:0.25-4) :", min = 0.25, max = 4, value =c(0.75,1.25))),

column(12,sliderInput("Position", "Zoom in (Left to right; bp):", min = 1, max = 50000, value = c(1,50000))),

actionButton("submit", label = "Submit",class = "btn-warning"),

  downloadButton('save','save',style = "background-color:#C9DD03"),
  downloadButton('Repeat','Repeat',style = "background-color:#00FFFF"),
  downloadButton('Variation','Variation',style = "background-color:#00FFFF"),
  downloadButton('Haplotype_fa','Haplotype_fa',style = "background-color:#00FFFF"),
  downloadButton('CDS_fa','CDS_fa',style = "background-color:#00FFFF"),

#sliderInput("Position2", "Right (bp):", min = 1, max = 50000, value = 500),
actionButton("submit2", label = "Zoom_In",class = "btn-warning"),

actionButton("submit3", label = "Clustering genomes",style="color: #00F; background-color: #FFCCCC; border-color: #c34113; border-radius: 10px; border-width: 2px"),
actionButton("submit3_2",label = "Tree",style="color: #00F; background-color: #FFCCCC; border-color: #c34113; border-radius: 10px; border-width: 2px"),

# actionButton('others',  label ='OthersSelected',style = "color: #333"),
# actionButton("submit4", label = "SecondRound", style = "color: #333"),


actionButton("done", label = "Refresh",style = "background-color:#FF6666"),

# actionButton("refresh", "Refresh page",style = "background-color:#99FF99"),
# actionButton("Second_round", "Zoom_In",style = "background-color:#CCCCFF"),


), # sidebarPanel


mainPanel(
fluidRow(
  
  # column(12, plotOutput("plot",click = "plot_click",dblclick = "plot_dblclick",hover = "plot_hover",brush = "plot_brush",width = "100%")),
  #  textOutput("text"), ## Outlier4
     span(textOutput("text"), style="color:green"),
     span(textOutput("text2"), style="color:red"),
     span(textOutput("text3"), style="color:black"),
     span(textOutput("text4"), style="color:blue"),

     column(12, plotOutput("plot",click = "plot_click",dblclick = "plot_dblclick",width = "100%")),

     column(9, offset = 1,plotOutput("plot2",click = "plot_click",dblclick = "plot_dblclick",width = "90%"),style='padding-left:0px; padding-right:0px; padding-top:5px; padding-bottom:5px'), # Zoom_in

     column(6, offset = 0,plotOutput("plot3",click = "plot_click",dblclick = "plot_dblclick",width = "80%"),style='padding-top:5px; padding-bottom:7px'), # clustering 1
     column(6, offset = 0,plotOutput("plot4",click = "plot_click",dblclick = "plot_dblclick",width = "90%"),style='padding-top:5px; padding-bottom:5px'), # clustering 1


  #   column(4, plotOutput("plot_tree")),
      column(12, offset = 0,DT::dataTableOutput("table00"),style='padding-top:5px; padding-bottom:5px'), ## 06/17, cluster information
      
      column(12, offset = 0,DT::dataTableOutput("table0"),style='padding-top:5px; padding-bottom:5px'),
    # DT::dataTableOutput("table0"),
      column(12, offset = 0,DT::dataTableOutput("table2"),style='padding-top:5px; padding-bottom:5px'),
    # DT::dataTableOutput("table2"),
  #   verbatimTextOutput("info"),
      dataTableOutput("table3")

)
) # mainPanel



) # sidebarLayout



      ) # For tagList
    }, # For ui function of page_1



# To add server function part for page1

server <- function(input, output, session){




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
if(input$Gene==''|| input$Chr==''||input$Upstream==''||input$Downstream==''||input$id==''){
  refresh()

  }else{


Gene <- input$Gene  ## from shiny input
cds_ids <- input$Gene  ## from shiny input

query <- input$Gene ## from shiny input
dna <- readDNAStringSet("Triticum_aestivum.IWGSC.cds_1.fa"); ## Try CDS.1
query_hit<-dna[grepl(query, dna@ranges@NAMES)];
query_fa <- "query.fasta";
writeXStringSet(query_hit, query_fa, append=FALSE, compress=FALSE, format="fasta");

query_length<-length(readDNAStringSet("query.fasta")[[1]]); # target CDS size

system("cat query.fasta > query_COPY.fasta");

system("perl -p -i -e 's/.[0-9]$//g' query.fasta")
system("perl -p -i -e 's/.[0-9]$//g' query_COPY.fasta")

Name_update_1 <- paste(Gene,"_mRNA",sep = "");
Name_update_2 <- paste(Gene,"_CDS",sep = "");

replace_pattern_in(file_contents=query,replace=Name_update_1,file_pattern="query.fasta",file_contents_perl = FALSE);
replace_pattern_in(file_contents=query,replace=Name_update_2,file_pattern="query_COPY.fasta",file_contents_perl = FALSE);

system("cat query.fasta query_COPY.fasta > query_Working.fasta");
Name_update_3 <- paste(Gene,"_ref",sep = "");
file.rename("query_Working.fasta", Name_update_3)
system("rm query.fasta query_COPY.fasta");
writeXStringSet(query_hit, query_fa, append=FALSE, compress=FALSE, format="fasta");

##############
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

output$text <- renderText({ paste('The following haplotypes may have other potential places (Table below) in the same chromosome (at least 200 kb from here):',sep = "") });
output$text2 <- renderText({Outlier4_name2});
output$text3 <- renderText({paste0(c('For example, IWGSC contains ',Outlier4[which(Outlier4$Genome=="IWGSC"),]$Size,' bp DNAs mapped onto other places.'), collapse= " ",sep = "")});
}

## 
new.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny";
setwd(new.folder);

gff3_Coor<-read.table("Triticum_aestivum.IWGSC.53_working_gene_coordinate",header=F);
gff3_Coor1<-gff3_Coor[which(gff3_Coor$V3==Gene),][,1:2]; # annotation in GFF3 (only gene)

Mapped0<-data.frame(Filtered_output2[which(Filtered_output2$Genome=="IWGSC"),][,6],Filtered_output2[which(Filtered_output2$Genome=="IWGSC"),][,7]);
Mapped1<-cbind(min(as.vector(as.matrix(Mapped0))),max(as.vector(as.matrix(Mapped0)))); # Coordinate as shown in plot (all CDS combined)
output$text4 <- renderText({paste0(c('Submitted gene (IWGSC) annotated between: ', gff3_Coor1 ,'; Plotted CDS (IWGSC) contained between: ',Mapped1), collapse= " ",sep = "")});

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







} ## To check NULL input after clicking on submit button with a refresh()


}) ## observeEvent submit !!



observeEvent(input$done,{

Gene <- input$Gene  ## from shiny input

if(Gene!=''){
new.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny"
result_files<-list.files(new.folder, pattern = Gene)
file.remove(result_files)
file.remove("query.fasta","grasrep.fa.nhr","grasrep.fa.nin","grasrep.fa.nsq")

file.remove("extract_syn_fa_Wheat_1_05_17_Test4_copy_3.pl")

refresh()

}else{refresh()}

}) ## observeEvent DONE, To remove all temp files


# observeEvent(input$refresh, {
#      refresh()
#    }) ## observeEvent refresh



######################################################## another sumbit ?? to view other parts of blast result in the same chromosome!


######################################################## another sumbit ?? to view other parts of blast result in the same chromosome!





    } # server function of Page_1



  ) # page for Page_1
} # Page_1 function











############ Page 0

source('Page_Function_Main.R', local = TRUE); ## To load Main page function!

############ Page 2

source('Page_Function_Barley.R', local = TRUE); ## To load Page2 function for Barley!

############ Page 3

source('Page_Function_Sorghum.R', local = TRUE); ## To load Page3 function for Shorghum!

############ Page 4

source('Page_Function_Maize.R', local = TRUE); ## To load Page4 function for Maize!

############ Page 5

source('Page_Function_Rice.R', local = TRUE); ## To load Page5 function for Rice!

############ Page 6

source('Page_Function_WheatTest.R', local = TRUE); ## To load Page6 function for WheatTest!

############ To combine pages together

 brochureApp(
  page_0(),
  page_1(),
  page_2(),
  page_3(),
  page_4(),
  page_5(),
  page_6()
# To add many other pages

)



