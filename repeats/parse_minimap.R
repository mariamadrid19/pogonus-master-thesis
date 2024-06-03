#Script written by Steven Van Belleghem 2024

mini <- read.table("alignment_nie_to_dud.paf", h=F, sep='\t', fill=T)

colnames(mini) <- c('Qname', 'Qlen', 'Qstart', 'Qend', 'ori', 
                    'Tname', 'Tlen', 'Tstart', 'Tend', 'match_bases',
                    'match_len', 'mapQ', 'Tag1', 'Tag2', 'Tag3', 'Tag4', 'Tag5')

mini <- mini[order(mini$Qname, mini$Qstart),]
head(mini)
nrow(mini)

# filter for alignments of at least 1000 bp and an identity score of >0.7 (~blast dentity score)
miniS <- subset(mini, mini$match_bases > 500 & mini$match_bases/mini$match_len > 0.5)
nrow(miniS)
head(miniS)

# loop to find for each target alignment interval overlapping intervals. If interval falls within another, add one to that row

miniS$seqCov <- 1 # every alignment starts with a coverage of 1 (also needed for the summing)

for(e in 1:nrow(miniS)){ # loop through each line
  scafE  <- miniS[e,]$Qname
  startE <- miniS[e,]$Qstart
  endE   <- miniS[e,]$Qend
  countE <- miniS[e,]$seqCov
  
  miniSubset <- subset(miniS, miniS$Qname == scafE)
  
  for(i in 1:nrow(miniSubset)){ # for each line, make smaller subse from the original table which includes lines from the same scaffold
    
    scafI  <- miniSubset[i,]$Qname
    startI <- miniSubset[i,]$Qstart
    endI   <- miniSubset[i,]$Qend
    countI <- miniSubset[i,]$seqCov
    
    if(startE >= startI & endE <= endI){ # check if interval falls within previous interval
      if(!identical(miniS[e,], miniSubset[i,])){ # check if we're not just looking at the exactly same line
        if(miniS[e,]$Tstart != miniSubset[i,]$Tstart & 
           miniS[e,]$Tend != miniSubset[i,]$Tend ){ # make sure the alignment is on a different region in the target
          
          miniS[e,]$seqCov <- countE + 1
          countE <- miniS[e,]$seqCov
        }
      }
    }
  }
}

# You can start here
write.table(miniS, file = "minimap2_repeatscore.txt", quote = F, row.names = F)


#plots
miniS <- read.table("minimap2_repeatscore_matchL05.txt", h=T, sep=' ', fill=T)
repeats <- subset(miniS, miniS$seqCov >20)

# plotting
test <- subset(miniS, miniS$Qname == 'CM008230.1_RagTag')
# Plot the intervals
plot(0, 0, type = "n", xlim = c(0, max(test$Qend)), ylim = c(0, max(test$seqCov)), xlab = "Position", ylab = "Repeat score")
# Add segments
segments(test$Qstart, test$seqCov, test$Qend, test$seqCov, col = "blue", lwd = 2)

#subset only the linkage groups
LGs <- subset(miniS, grepl("CM", Qname))

#convert scores > 2 into 2 
LGs$seqCov <- ifelse(LGs$seqCov > 2, 2, LGs$seqCov)

# Create a new color column based on seqCov values
LGs$color <- ifelse(LGs$seqCov == 1, "orange", ifelse(LGs$seqCov == 2, "orange3", NA))

library(ggplot2)
# Plot the data
ggplot(LGs, aes(x = Qend, y = seqCov, color = color)) +
  geom_point(size = 0.1) +  # Adjust the size for smaller points
  scale_color_identity() +  # Use the color column directly
  facet_grid(Qname ~ .) +  # This will create a different row for each chromosome
  labs(x = "Position", y = "Repeat Score", title = "Repeat Score Across Chromosomes") +
  ylim(0, 3) + xlim(0,125000000) +# Set the y-axis limit
  theme_minimal()
