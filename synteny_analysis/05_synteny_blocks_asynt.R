#Uses functions from asynt.R algorithm from Simon Martin (https://github.com/simonhmartin/asynt)

setwd("/Users/mariamadrid/Documents/thesis/asynt-scaffolds")

# Define the list of scaffolds
scaffolds <- c("A", "B", "C", "D", "E", "F", "G", "H")

# Define the directory containing your files
file_dir <- "/Users/mariamadrid/Documents/thesis/asynt-scaffolds/"

par(mfrow=c(4, 2))  # Adjust margins as needed

# Loop through each scaffold
for (scaffold in scaffolds) {
  # Construct file paths
  paf_file <- paste0(file_dir, "scaffold_", scaffold, ".paf.gz")
  ref_fai_file <- paste0(file_dir, "scaffold_", scaffold, "_dud.fasta.fai")
  query_fai_file <- paste0(file_dir, "scaffold_", scaffold, "_nieu.fasta.fai")
  
  # Import alignments
  alignments <- import.paf(paf_file)
  alignments_plotter <- alignments[,1:7]
  alignments_plotter$ref_sp <- paste0(scaffold, "_dud")
  alignments_plotter$query_sp <- paste0(scaffold, "_nieu")
  
  # Import scaffold length data
  ref_data <- import.genome(fai_file = ref_fai_file)
  query_data <- import.genome(fai_file = query_fai_file)
  
  # Remove short alignments
  alignments <- subset(alignments, Rlen >= 20000 & Qlen >= 20000)
  
  # Visualize the alignment
  plot.alignments.multi(alignments, reference_lens = ref_data$seq_len, query_lens = query_data$seq_len, sigmoid = TRUE)
}



#import alignments
alignments <- import.paf("/Users/mariamadrid/Documents/thesis/asynt-scaffolds/s6_ref.paf.gz")
alignments_plotter <- alignments[,1:7]
alignments_plotter$ref_sp <- "old_reference"
alignments_plotter$query_sp <- "s6"

#Import scaffold length data which is necessary for some functions of asynt
ref_data <- import.genome(fai_file="/Users/mariamadrid/Documents/thesis/asynt-scaffolds/10_old_ref.fasta.fai")
query_data <- import.genome(fai_file="/Users/mariamadrid/Documents/thesis/asynt-scaffolds/s6_dh2.fasta.fai")

# Remove short alignments
# Rlen and Qlen give the number of bases on the reference and query that are included in an alignment tract
alignments <- subset(alignments, Rlen >= 20000 & Qlen >= 20000)

# Only care about query scaffolds that share a large proportion of aligned sequence with the reference scaffold
# Find these by looking at the total alignemnt length for each scaffold
query_aln_len <- get.query.aln.len(alignments)
#barplot(query_aln_len, las=2)

#a similar approach can be used with aligned proportion, which may be more appropriate if query scaffolds are very variable in length
query_aln_prop <- get.query.aln.prop(alignments, query_lens = query_data$seq_len)
#barplot(query_aln_prop, las=2)

#subset the alignments for only those contigs with large proportion aligned to the reference contig of interest
alignments <- subset(alignments, query_aln_prop[query] > 0.1)


#visualise the alignment
synblocks <- get.synteny.blocks.multi(alignments, min_subblock_size=200)
synblocks <- get.synteny.blocks.multi(synblocks, min_subblock_size=2000)
synblocks <- get.synteny.blocks.multi(synblocks, min_subblock_size=20000)
plot.alignments.multi(synblocks, reference_lens=ref_data$seq_len, query_lens=query_data$seq_len, sigmoid=T)

alignments <- subset(alignments, reference == "CM008231.1")
synblocks <- get.synteny.blocks.multi(alignments, min_subblock_size=200)
synblocks <- get.synteny.blocks.multi(synblocks, min_subblock_size=2000)
synblocks <- get.synteny.blocks.multi(synblocks, min_subblock_size=20000)
plot.alignments.multi(synblocks, reference_lens=ref_data$seq_len, query_lens=query_data$seq_len, sigmoid=T)
