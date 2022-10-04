############ Page 4

page_4 <- function(){

  page(
    href = "/page4",


    ui <- function(request){
      tagList(
        h1("This is Maize genome page",style="text-align:center"),
        nav_links,


# To add ui part for page4   

useShinyjs(),

column(8,offset = 5, titlePanel("Maize Alignment output")), # titlePanel

sidebarLayout(

sidebarPanel(

column(12,textInput("Gene","Gene Name (Zm00001eb093920 or others)")),

# column(6,textInput("Gene_ID","Sobic")),

column(12,textInput("Chr","Chromosome (chr2 or others)")),

column(12,textInput("Upstream","Upstream (bp)")),
column(12,textInput("Downstream","Downstream (bp)")),

column(12,
pickerInput(
  inputId = "id", 
  label = "Genomes (Select B73 and others) :", 
  choices = c('B73', 'B97', 'CML52', 'CML69', 'CML103', 'CML228',  'CML247',  'CML277',  'CML322',  'CML333', 'HP301', 'Il14H', 'Ki3', 'Ki11', 'Ky21', 'M37W', 'M162W', 'Mo17', 'Mo18W', 'Ms71', 'Ky21', 'NC350', 'NC358', 'Oh7B', 'Oh43', 'P39','Tx303','Tzi8'),
  options = list(
    'actions-box' = TRUE, 
    size = 27,
    'selected-text-format' = "count > 2"
  ), 
  multiple = TRUE,
)
),

# fileInput("file1", "Upload a fasta file"),

column(12,sliderInput("Position", "Left to right (bp):", min = 1, max = 50000, value = c(1,50000))),

actionButton("submit", label = "Submit",class = "btn-warning"),

  downloadButton('save','save',style = "background-color:#C9DD03"),
  downloadButton('Repeat','Repeat',style = "background-color:#00FFFF"),
  downloadButton('Variation','Variation',style = "background-color:#00FFFF"),
  downloadButton('Haplotype_fa','Haplotype_fa',style = "background-color:#00FFFF"),
  downloadButton('CDS_fa','CDS_fa',style = "background-color:#00FFFF"),

#sliderInput("Position2", "Right (bp):", min = 1, max = 50000, value = 500),

actionButton("submit2", label = "Zoom_In",class = "btn-warning"),


actionButton("done", label = "DONE",style = "background-color:#FF6666"),
actionButton("refresh", "Refresh page",style = "background-color:#99FF99"),
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
     column(9, offset = 1,plotOutput("plot2",click = "plot_click",dblclick = "plot_dblclick",width = "90%")),
  #   column(4, plotOutput("plot_tree")),

      DT::dataTableOutput("table2")
   #  DT::dataTableOutput("table3")
   # verbatimTextOutput("info")

)
) # mainPanel



) # sidebarLayout



      ) # For tagList
    }, # For ui function of page_4



# To add server function part for page4

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
dna <- readDNAStringSet("Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.cds.fa"); ## Try CDS.1
query_hit<-dna[grepl(query, dna@ranges@NAMES)];
query_fa <- "query.fasta";
writeXStringSet(query_hit, query_fa, append=FALSE, compress=FALSE, format="fasta");

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
system("cat extract_syn_fa_Maize_1_06_01.pl > extract_syn_fa_Maize_1_06_01_Test4_copy_3.pl")


Upstream_update <- paste('p_size=',input$Upstream,sep = "");  ## from shiny input
Downstream_update <- paste('a_size=',input$Downstream,sep = "");  ## from shiny input
replace_pattern_in(file_contents="p_size=10000",replace=Upstream_update,file_pattern="Test4_copy_3",file_contents_perl = TRUE);
replace_pattern_in(file_contents="a_size=1000",replace=Downstream_update,file_pattern="Test4_copy_3",file_contents_perl = TRUE);

Chr_update0 <- paste('target_ch=',input$Chr,sep = "'");  ## from shiny input
symbol<-c("'");
Chr_update <- paste(Chr_update0,symbol,sep = "");
replace_pattern_in(file_contents="target_ch='chr2'",replace=Chr_update,file_pattern="Test4_copy_3",file_contents_perl = TRUE);

Gene_update0 <- paste('gene=',input$Gene,sep = "'");  ## from shiny input
Gene_update <- paste(Gene_update0,symbol,sep = "");
replace_pattern_in(file_contents="gene='Zm00001eb093920'",replace=Gene_update,file_pattern="Test4_copy_3",file_contents_perl = TRUE);


##############
system("perl extract_syn_fa_Maize_1_06_01_Test4_copy_3.pl 1");


######################################################################## To search for the outliers!
result_folder1<-c('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Maize/B73/Candidate_genes/');
result_folder2<-paste(result_folder1,Gene,sep="");
setwd(result_folder2);
## backup Haplotype_syn file
Outlier0<-paste(Gene,'_Haplotype_syn',sep = "");
Outlier00<-paste(Gene,'_Blast_Original',sep = "");
Outlier1<-read.table(Outlier0,header=T); ## _Haplotype_syn
write.table(Outlier1, file = Outlier00,sep= "\t",quote = FALSE,row.names = FALSE);


Outlier5 <- Outlier1[which(Outlier1$Similarity>=96),]; ## filter1: Similarity >=96, kept


Filtered_output <- list();
for(i in 1:length(unique(Outlier5[,4]))){
Outlier2<-Outlier5[Outlier5[,4]==unique(Outlier5[,4])[i],]
Outlier3<-sort(Outlier2[,6])[as.matrix(dist(sort(Outlier2[,6])))[,1]<=200000]; ## filter2: <= 200 kb, kept
out_ind<-which(Outlier2[,6] %in% Outlier3);
Filtered_output[[i]]<-Outlier2[out_ind, ];
rm(Outlier2,Outlier3,out_ind);
} 
Filtered_output2 <- as.data.frame(rbindlist(Filtered_output));
write.table(Filtered_output2, file = Outlier0,sep= "\t",quote = FALSE,row.names = FALSE);
##

if(nrow(anti_join(Outlier1,Filtered_output2))!=0){
Outlier4<-anti_join(Outlier1,Filtered_output2); # differences of _Haplotype_syn and _Blast_Original; Not shown in plot!
Outlier4_name <- unique(Outlier4[,4]);
Outlier4_name2 <- intersect(input$id,Outlier4_name);

output$text <- renderText({ paste('The following haplotypes may have other potential places in the same chromosome (at least 200 kb from here):',sep = "") });
output$text2 <- renderText({Outlier4_name2});
output$text3 <- renderText({paste0(c('For example, B73 contains ',Outlier4[which(Outlier4$Genome=="B73"),]$Size,' bp DNAs mapped onto other places.'), collapse= " ",sep = "")});
}

## 
new.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny";
setwd(new.folder);

gff3_Coor<-read.table("Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_working_gene_coordinate",header=F);
gff3_Coor1<-gff3_Coor[which(gff3_Coor$V3==Gene),][,1:2]; # annotation in GFF3 (only gene)

Mapped0<-data.frame(Filtered_output2[which(Filtered_output2$Genome=="B73"),][,6],Filtered_output2[which(Filtered_output2$Genome=="B73"),][,7]);
Mapped1<-cbind(min(as.vector(as.matrix(Mapped0))),max(as.vector(as.matrix(Mapped0)))); # Coordinate as shown in plot (all CDS combined)
output$text4 <- renderText({paste0(c('Submitted gene (B73) annotated between: ', gff3_Coor1 ,'; Plotted CDS (B73) contained between: ',Mapped1), collapse= " ",sep = "")});

######################################################################## To search for the outliers!


system("perl extract_syn_fa_Maize_1_06_01_Test4_copy_3.pl 2");
system("rm extract_syn_fa_Maize_1_06_01_Test4_copy_3.pl");

result_folder1<-c('/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/Maize/B73/Candidate_genes/')
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
gDNAs_blast <- subset(gDNAs_blast, gDNAs_blast[,4]> 300)
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


} ## To check NULL input after clicking on submit button with a refresh()


}) ## observeEvent submit !!



observeEvent(input$done,{

Gene <- input$Gene  ## from shiny input

if(Gene!=''){
new.folder <- "/mnt/946c1663-fcbd-4a78-8887-c55f23c5b496/bszhang/ShinyApps/TestShiny"
result_files<-list.files(new.folder, pattern = Gene)
file.remove(result_files)
file.remove("query.fasta","grasrep.fa.nhr","grasrep.fa.nin","grasrep.fa.nsq")
}else{refresh()}

}) ## observeEvent DONE, To remove all temp files


observeEvent(input$refresh, {
      refresh()
    }) ## observeEvent refresh



######################################################## another sumbit ?? to view other parts of blast result in the same chromosome!




######################################################## another sumbit ?? to view other parts of blast result in the same chromosome!






    } # server function of Page_4
  ) # page for Page_4
} # Page_4 function

