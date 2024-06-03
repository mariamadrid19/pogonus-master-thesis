#Script written by Steven Van Belleghem 2024

library(GenomicRanges)
library(rtracklayer)

# Read the GFF file into a data frame
gff_file <- "sorted_prim_dud.fasta.out.gff.gz"

gff_data <- read.table(gff_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(gff_data) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# Create GRanges object
gr <- GRanges(seqnames = gff_data$seqname,
              ranges = IRanges(start = gff_data$start, end = gff_data$end),
              strand = gff_data$strand)

# Get unique scaffolds
scaffolds <- unique(seqnames(gr))
scaffolds <- unique(scaffolds[grep("RagTag", scaffolds)])[c(1:10)]

result_list <- lapply(scaffolds, function(scaffold) {
  # Subset the GRanges object for the current scaffold
  gr_scaffold <- gr[seqnames(gr) == scaffold]
  
  # Determine the length of the scaffold
  scaffold_length <- max(end(gr_scaffold))
  
  # Create sliding windows for the scaffold
  windows <- slidingWindows(GRanges(seqnames = scaffold, ranges = IRanges(start = 1, end = scaffold_length)), 
                            width = 100000, step = 100000)[[1]]
  
  # Calculate the amount covered by features in each window
  coverage <- sapply(seq_along(windows), function(i) {
    window <- windows[i]
    overlaps <- pintersect(rep(window, length(gr_scaffold)), gr_scaffold)
    sum(width(overlaps))
  })
  
  # Create a data frame for the results
  data.frame(
    seqnames = scaffold,
    start = start(windows),
    end = end(windows),
    coverage = coverage
  )
})

# Combine the results for all scaffolds
result <- do.call(rbind, result_list)

# Print the result
print(result)

result$mid <- (result$start+result$end)/2
result$perc <- result$coverage/100000

scaffolds <- c("CM008230.1_RagTag", "CM008231.1_RagTag", "CM008233.1_RagTag", "CM008234.1_RagTag", "CM008235.1_RagTag", "CM008236.1_RagTag", "CM008237.1_RagTag", "CM008238.1_RagTag","CM008239.1_RagTag","CM008240.1_RagTag")

for (scaffold in scaffolds) {
  result_scaffolds <- subset(result, result$seqnames == scaffold)
  
  plot(result_scaffolds$start, result_scaffolds$perc, pch=19, cex = 0.3, col="orange4",main = scaffold, xlab = "Chromosome Position",ylab="Repeat Content")
}
