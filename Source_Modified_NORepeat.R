



Plot_SV <- function(genomes, g_lab, repmask, CDS_gDNA_blast, gDNAs_blast, N_Gap, Anno, variations,output_flag ) {
# if (output_flag == 1) {png(file = paste(Dir, Gene, '_NAM.png', sep = ''), width= 19 * .5, height= 12 * .75 , pointsize= 10 , units = "in", res = 600)};
 if (output_flag == 1) {png(file = paste(Dir, Gene, '_NAM.png', sep = ''), width= 8, height= 11 , pointsize= 10 , units = "in", res = 600)};
 
 par(mar = c(1.2, 1.0, 0, 0) , mgp = c(1, 0.1, 0), tck = -0.01, cex.axis = .9, cex.lab = 1, family = "mono"); ## ??
 
 plot(-100, -100, xlim = x_lim, ylim = c(-0.5,length(genomes) + 0), xlab = '', ylab = '', xaxt = "n", yaxt = "n", bty = "n");

 
 for (g in 1:length(genomes)) {
 

## To add arrows and axis ##
  if (g == length(genomes)) { 
   arrows(range(gDNAs_blast[,9:10])[1], length(genomes)-0.5, range(gDNAs_blast[,9:10])[2], length(genomes)-0.5,lwd=2.5);
     }
  if (g == length(genomes)) {
   axis(1, at = c(seq(from=range(gDNAs_blast[,9:10])[1],to=range(gDNAs_blast[,9:10])[2],by=(range(gDNAs_blast[,9:10])[2]-range(gDNAs_blast[,9:10])[1])/(10-1))),labels=paste0(c(ceiling(seq(from=range(gDNAs_blast[,9:10])[1],to=range(gDNAs_blast[,9:10])[2],by=(range(gDNAs_blast[,9:10])[2]-range(gDNAs_blast[,9:10])[1])/(10-1)))),'bp'),cex.axis=0.8,cex.lab=0.8,tick = TRUE,col = "blue", lty = 1, lwd = 2.5, lwd.ticks=1.5,tck =0.02,col.ticks = "black",col.axis = "black");
     }
## To add arrows and axis ##



  if (g < length(genomes)) {
   gDNAs <- subset(gDNAs_blast, gDNAs_blast[,1] == genomes[g] & gDNAs_blast[,2] == genomes[g + 1]);
   for (k in 1:nrow(gDNAs)) {
    polygon(c(gDNAs[k,7:8], gDNAs[k,c(10,9)]), length(genomes) - g - c(0.3, .3, 0.9, 0.9), col = adjustcolor( "gray", alpha.f = 0.5), border = "NA");
   }
  }
 
  ## for mRNA
  for (cds_i in 1:length(cds_ids) ) {
    CDSs <- subset(CDS_gDNA_blast, CDS_gDNA_blast[,1] == paste(cds_ids[cds_i], '_mRNA', sep = '') & CDS_gDNA_blast[,2] == genomes[g] );
    arrow_code <- 1;
    if (CDSs[1,9] > CDSs[1, 10]) {arrow_code <- 2}
    if (nrow(CDSs) > 0) {
 #    arrows(max(CDSs[,9:10]) + 100, length(genomes) - g - 0.175, min(CDSs[,9:10]) - 100, length(genomes) - g - 0.175, code = arrow_code, angle = 15, length = .05);
     rect(CDSs[,9], length(genomes) - g - 0.25, CDSs[,10] , length(genomes) - g - 0.1, col = "burlywood",  border = "NA");
    }
   }
  ### for CDS
  for (cds_i in 1:length(cds_ids) ) {
    CDSs <- subset(CDS_gDNA_blast, CDS_gDNA_blast[,1] == paste(cds_ids[cds_i], '_CDS', sep = '') & CDS_gDNA_blast[,2] == genomes[g] );
    arrow_code <- 1;
    if (CDSs[1,9] > CDSs[1, 10]) {arrow_code <- 2}
    if (nrow(CDSs) > 0) {
 #    arrows(max(CDSs[,9:10]) + 100, length(genomes) - g - 0.175, min(CDSs[,9:10]) - 100, length(genomes) - g - 0.175, code = arrow_code, angle = 15, length = .05);
     rect(CDSs[,9], length(genomes) - g - 0.25, CDSs[,10] , length(genomes) - g - 0.1, col = cds_col[cds_i], border ="NA");
    }
   }
   
    snp_v <- subset(variations, variations$Inbred == genomes[g] & variations$Type == 'snp');
   if (nrow(snp_v) >= 1) {
    segments(snp_v$Bp, rep(length(genomes) - g - 0.25, nrow(snp_v)), snp_v$Bp, rep(length(genomes) - g - 0.1, nrow(snp_v)), lwd = .01, col = "blue");
   }
   ins_v <- subset(variations, variations$Inbred == genomes[g] & variations$Type == 'ins');
   if (nrow(ins_v) >= 1) {
    points(ins_v$Bp, rep(length(genomes) - g - 0.1 , nrow(ins_v)) , pch = 25, bg = indel_col[ins_v$Size %%3 + 1], col = 1, cex = 0.6);
    text(ins_v$Bp, rep(length(genomes) - g - 0.1 + 0.3, nrow(ins_v)), paste(ins_v$Size, 'bp', sep = ''), srt = 90, cex = .4)
   }
   ins_v <- subset(variations, variations$Inbred == genomes[g] & variations$Type == 'del');
   if (nrow(ins_v) >= 1) {
    points(ins_v$Bp, rep(length(genomes) - g - 0.25, nrow(ins_v)), pch = 24, bg = indel_col[ins_v$Size %%3 + 1], col = 1, cex = 0.6);
    text(ins_v$Bp, rep(length(genomes) - g - 0.1 - 0.4, nrow(ins_v)), paste(ins_v$Size, 'bp', sep = ''), srt = 90, cex = .4)
   }

  self <- subset(gDNAs_blast,  gDNAs_blast[,1] == genomes[g] & gDNAs_blast[,2] == genomes[g ])
  
  sizes <- round((max(self[,7:8]) - min(self[,7:8]) )/ 1000, 2);
 # legend(max(self[,7:8]),length(genomes) - g + 0.3, c(g_lab[g], paste('(', sizes, 'kb)', sep = '')), bty = "n", adj = c(0, 0), cex = .8 )
  legend(max(self[,7:8]),length(genomes) - g + 0.3, c(g_lab[g], paste('(', sizes, 'kb)', sep = '')), bty = "n", adj = c(0, 0), cex = .8 )
 
  gap <- subset(N_Gap, N_Gap[,1] == genomes[g]);
  if (nrow(gap) > 0) {
   segments(gap[,2], rep(length(genomes) - g - 0.025, nrow(gap)), gap[,3], rep(length(genomes) - g - 0.025, nrow(gap)), lty = 2, col = "black")
  }
  
   for (k in 1:nrow(self)) {
   if (self[k,7] == self[k,9] & self[k,8] == self[k,10]) {
    rect(self[k,7] - 100, length(genomes) - g, self[k,8] + 100, length(genomes) - g - 0.05, col = "darksalmon", border = "NA");
   } 
 #  else {
 #   mid <- mean(t(self[k,8:9]));
 #   polygon(c(self[k, 7], mid, self[k, 10], self[k,9], mid, self[k,8]), length(genomes) - g - c(0, -0.5, 0, 0, -0.45, 0), col = "purple", border = "NA")
 #  }
  }
   # for annotation
   gIDs <- subset(Anno, Anno$Genome == genomes[g] & Anno$Type == 'gene');
   gCDS <- subset(Anno, Anno$Genome == genomes[g] & Anno$Type == 'CDS');
   if (nrow(gIDs) > 0) {
    arrows(gIDs$Start, rep(length(genomes) - g + 0.05, nrow(gIDs)), gIDs$End, rep(length(genomes) - g + 0.05, nrow(gIDs)), code = gIDs$Strand, angle = 15, length = .05, lwd = .5 )
 #   text(gIDs$Start, rep(length(genomes) - g + 0.075, nrow(gIDs)), gIDs$refID, cex = .7 )
    }
    
   if (nrow(gCDS) > 0) {
    
    for (k in 1:nrow(gCDS)) {
     rect(gCDS$Start[k], length(genomes) - g + 0.02, gCDS$End[k], length(genomes) - g - 0.05 - 0.02, col = "darksalmon", lwd = .5)
    }
   };
 
 # primers <- subset(primer_blast, primer_blast[,2] == genomes[g]);
 # for (k in 1:nrow(primers)) {
 #  arrows(primers[k,9], length(genomes) - g - 0.175, primers[k,10], length(genomes) - g - 0.175, col = 2);
 #  text(mean(primers[k,9], primers[k,10]), length(genomes) - g - 0.175, primers[k,1])
 # }

## Repeats
#   if (nrow(repmask) > 0) {
#   g_repmask <- subset(repmask, repmask[,1] == genomes[g]);
#   if (nrow(g_repmask) > 0) {
#    g_rep_y <- length(genomes) - g - 0.025 + 0.01 ;
#    rect(g_repmask[,7], g_rep_y + 0.1, g_repmask[,8], g_rep_y +  0.2, col = "black", border = "NA");
#    g_repmask_L <- g_repmask[abs(g_repmask[,8] - g_repmask[,7]) > 200,]
#   #g_repmask_Harb <- subset(g_repmask, g_repmask[,7] == 'DNA/Harbinger')
# #   if(nrow(g_repmask_Harb) > 0) {text(rowMeans(g_repmask_Harb[,2:3]), rep(g_rep_y + 0.3, nrow(g_repmask_Harb)), g_repmask_Harb[,4], cex = .6, col = "red")}
#    if(nrow(g_repmask_L) > 0) {text(rowMeans(g_repmask_L[,7:8]), rep(g_rep_y + 0.3, nrow(g_repmask_L)), g_repmask_L[,2], cex = .6)}
#   }
#  }

 }

 if (output_flag == 1) { dev.off()}
}


