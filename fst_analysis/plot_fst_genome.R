#calculate chromosome coordinates
chrom.coords <- function(scafL,chromNames,gap = 9000000) {
  chromosome = vector()
  chromLengths  = vector()
  chromStarts = vector()
  chromEnds = vector()
  chromMid = vector()
  chrom = 1
  endLast = 0
  scafCurrent <- subset(scafL, chromosome == chromNames[1])
  chromosome[chrom] <- chrom
  chromLengths[chrom] <- sum(scafCurrent$length)
  chromStarts[chrom] <- endLast + 1
  chromEnds[chrom] <- endLast + chromLengths[chrom]
  chromMid[chrom] <- endLast + chromLengths[chrom]/2
  endLast = chromEnds[chrom]
  chrom = chrom + 1
  for (i in 2:length(chromNames)) {
    chromosome[chrom] <- chrom
    scafCurrent <- subset(scafL, chromosome == chromNames[i])
    chromLengths[chrom] <- sum(scafCurrent$length)
    chromStarts[chrom] <- endLast + gap + 1
    chromEnds[chrom] <- endLast + gap + chromLengths[chrom]
    chromMid[chrom] <- endLast + gap + chromLengths[chrom]/2
    endLast = chromEnds[chrom]
    chrom = chrom + 1
  }
  table <- as.data.frame(cbind(chromosome,chromLengths,chromStarts,chromEnds,chromMid))
  return(table)
}

#calculate scaffold coordinates
scaf.coords <- function(scafL,gap = 0) {
  scaffold = vector()
  scafStarts = vector()
  scafEnds = vector()
  chrom = 1
  endLast = 0
  for (e in 1:nrow(scafL)) {
    scaffold[chrom] <- levels(scafL$scaffold)[e]
    scafStarts[chrom] <- endLast + gap + 1
    scafEnds[chrom] <- endLast + gap + scafL$length[e]
    endLast = scafEnds[chrom]
    chrom = chrom + 1
  }
  table <- as.data.frame(cbind(scaffold,scafStarts,scafEnds))
  return(table)
}

#scafL <- read.table("Pogonus_LG_lengths_dud.txt", h=T)
scafL <- read.table("Pogonus_LG_lengths_nie.txt", h=T)
chromNames <-c(1:11)
chrom_coords <- chrom.coords(scafL, chromNames)
# scaf_coords <- scaf.coords(scafL)

scaf_coords <- merge(scafL, chrom_coords, by="chromosome", all.x=TRUE)

## read files
# pattern_BE <- "Pogonus_reseqALL_dudPrim\\.chr_(1|2|3|4|5|6|7|8|9|10)\\.stats_B_NIE_T_B_DUD_S_w50000_s50000_eggStats\\.stats"
# pattern_FR <- "Pogonus_reseqALL_dudPrim\\.chr_(1|2|3|4|5|6|7|8|9|10)\\.stats_F_GUE_T_F_GUE_S_w50000_s50000_eggStats\\.stats"
# pattern_PO <- "Pogonus_reseqALL_dudPrim\\.chr_(1|2|3|4|5|6|7|8|9|10)\\.stats_P_AVE_T_P_AVE_S_w50000_s50000_eggStats\\.stats"
# pattern_SP <- "Pogonus_reseqALL_dudPrim\\.chr_(1|2|3|4|5|6|7|8|9|10)\\.stats_S_HUE_T_S_COT_S_w50000_s50000_eggStats\\.stats"
# pattern_UM <- "Pogonus_reseqALL_dudPrim\\.chr_(1|2|3|4|5|6|7|8|9|10)\\.stats_E_SEV_T_F_CAM_S_w50000_s50000_eggStats\\.stats"
# pattern_HE <- "Pogonus_reseqALL_dudPrim\\.chr_(1|2|3|4|5|6|7|8|9|10)\\.stats_B_HEI_2000_B_HEI_2018_w50000_s50000_eggStats\\.stats"

pattern_BE <- "Pogonus_reseqALL_nieuPrim\\.chr_(1|2|3|4|5|6|7|8|9|10)\\.stats_B_NIE_T_B_DUD_S_w50000_s50000_eggStats\\.stats"
pattern_FR <- "Pogonus_reseqALL_nieuPrim\\.chr_(1|2|3|4|5|6|7|8|9|10)\\.stats_F_GUE_T_F_GUE_S_w50000_s50000_eggStats\\.stats"
pattern_PO <- "Pogonus_reseqALL_nieuPrim\\.chr_(1|2|3|4|5|6|7|8|9|10)\\.stats_P_AVE_T_P_AVE_S_w50000_s50000_eggStats\\.stats"
pattern_SP <- "Pogonus_reseqALL_nieuPrim\\.chr_(1|2|3|4|5|6|7|8|9|10)\\.stats_S_HUE_T_S_COT_S_w50000_s50000_eggStats\\.stats"
pattern_UM <- "Pogonus_reseqALL_nieuPrim\\.chr_(1|2|3|4|5|6|7|8|9|10)\\.stats_E_SEV_T_F_CAM_S_w50000_s50000_eggStats\\.stats"
pattern_HE <- "Pogonus_reseqALL_nieuPrim\\.chr_(1|2|3|4|5|6|7|8|9|10)\\.stats_B_HEI_2000_B_HEI_2018_w50000_s50000_eggStats\\.stats"


file_list_BE <- list.files(pattern = pattern_BE, full.names = TRUE)
file_list_FR <- list.files(pattern = pattern_FR, full.names = TRUE)
file_list_PO <- list.files(pattern = pattern_PO, full.names = TRUE)
file_list_SP <- list.files(pattern = pattern_SP, full.names = TRUE)
file_list_UM <- list.files(pattern = pattern_UM, full.names = TRUE)
file_list_HE <- list.files(pattern = pattern_HE, full.names = TRUE)

library(data.table)
stats_BE <- rbindlist(lapply(file_list_BE, fread))
stats_FR <- rbindlist(lapply(file_list_FR, fread))
stats_PO <- rbindlist(lapply(file_list_PO, fread))
stats_SP <- rbindlist(lapply(file_list_SP, fread))
stats_UM <- rbindlist(lapply(file_list_UM, fread))
stats_HE <- rbindlist(lapply(file_list_HE, fread))

comp=list(stats_BE,stats_FR,stats_PO,stats_SP,stats_UM,stats_HE)
names=c("Belgium T/S","France T/S","Portugal T/S", "Spain T/S", "UK (T) / FR (S)", "Heist 2000/2018")
col=c("black","black","black","black","black","black","black","black", "black")

par(mai=c(0.05,0.4,0.05,0.4), oma=c(2,0,1,0)+0)

layout(matrix(c(1:8), nrow=8, byrow=TRUE), height = c(1,1,1,1,1,1,0.5,2))
layout.show(n=8)

begin = chrom_coords$chromStarts[1]/1000000
end = chrom_coords$chromEnds[10]/1000000

top = 1
bot = 0
bot2 = 0

###
# plot Fst

colfuncR <- colorRampPalette(c("black", "red"))
colL <- colfuncR(100)

for (i in 1:length(comp)){
# for (i in 1:3){
  comb <- merge(comp[[i]], scaf_coords, by = 'scaffold', all.x=TRUE)
  comb <- as.data.frame(comb)
  comb$chromPos <- comp[[i]]$mid + comb$chromStarts
  
  comb <- na.omit(comb)
  
  comb[,9][comb[,9] < 0] <- 0
  

  
  comb$col <- colfuncR(100)[as.integer((comb[,9]/0.8)*100)+1]
  
  
  plot(0, pch = "",xlim = c(begin,end), ylim = c(bot,top), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE)
  rect(chrom_coords$chromStarts/1000000,rep(bot,length(chrom_coords$chromStarts)),chrom_coords$chromEnds/1000000,rep(top,length(chrom_coords$chromStarts)), col = c("gray95","gray90"), lwd = 0, border = c("gray95","gray90"))
  
  par(new=TRUE)
  plot(comb$chromPos/1000000, comb[,9], type="p", pch=19, cex=0.7,col=adjustcolor(comb$col, alpha=1),xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n")
  
  axis(2,cex.axis = 1, line = -1.5)
  mtext(side = 2, text = expression(paste("Fst")), cex=0.5,line = 0.8)
  mtext(side = 4, text = names[i], cex=0.6)
    
}

axis(1, at=chrom_coords[,5][1:10]/1000000, labels=(1:10),lwd=0, lwd.ticks=0)
segments(x0=5, y0=0.9, x1=15, y1=0.9, lwd = 1)
text(10,0.8,labels = "10Mb", cex = 1)

# plot scaffolds
#agp <- read.table('sorted_prim_dud.agp')
agp <- read.table('sorted_prim_nieu.agp')
agp <- subset(agp, agp$V5 == 'W')[,c(1:3)]

rows_with_pattern <- apply(agp, 1, function(row) any(grepl('CM.*RagTag', row)))
agp <- agp[rows_with_pattern, ]

colnames(agp) <- c('scaffold', 'start', 'end')

agpM <- merge(agp, scaf_coords, by = 'scaffold', all.x=TRUE)
agpM$chromPosStart <- agpM$start + agpM$chromStarts -1
agpM$chromPosEnd <- agpM$end + agpM$chromStarts -1

plot(0, pch = "",xlim = c(begin,end), ylim = c(bot,top), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE)

rect(agpM$chromPosStart/1000000,
     rep(bot2,length(agpM$chromPosStart)),
     agpM$chromPosEnd/1000000,
     rep(top,length(agpM$chromPosEnd)), 
     col = c("#fbb4ae", "#b3cde3"), lwd = 0, border =c("#fbb4ae", "#b3cde3"))
