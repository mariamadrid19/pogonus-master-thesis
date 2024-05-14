library(ggplot2)
library(ggfortify)
library(dplyr)

setwd("/Users/mariamadrid/Documents/thesis/coverage")

# Custom function to rename scaffolds
rename_scaffold <- function(scaffold) {
  scaffold_number <- as.numeric(sub("scaffold_", "", scaffold))
  ifelse(scaffold_number < 10, 
         sprintf("scaffold_0%d", scaffold_number), 
         sprintf("scaffold_%d", scaffold_number))
}


#belgian individuals, mapped against Dudzele assembly (chromosomes)
#import coverage data
GC129388_dudPrim <- read.delim("~/Documents/thesis/coverage/GC129388.dudPrim.pcov", header=FALSE, comment.char="#")
GC129394_dudPrim <- read.delim("~/Documents/thesis/coverage/GC129394.dudPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC129388_dudPrim) <- c("region", "coverage", "length")
names(GC129394_dudPrim) <- c("region", "coverage", "length")


#INDIVIDUAL FROM NIEUWPOORT (TIDAL, SHORT WINGS)
scaf_coords <- read.delim("~/Documents/thesis/coverage/scaf_coords.tsv")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC129388_dudPrim <- GC129388_dudPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC129388_dudPrim <- GC129388_dudPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC129388_dudPrim <- GC129388_dudPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC129388_dudPrim <- GC129388_dudPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC129388_dudPrim <- GC129388_dudPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC129388_dudPrim_chromosomes <- GC129388_dudPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC129388_dudPrim_chromosomes <- GC129388_dudPrim_chromosomes %>%
  mutate(scaffold = case_when(
    grepl("CM008233", scaffold) ~ "scaffold_02",
    grepl("CM008230", scaffold) ~ "scaffold_01",
    grepl("CM008234", scaffold) ~ "scaffold_03",
    grepl("CM008235", scaffold) ~ "scaffold_04",
    grepl("CM008236", scaffold) ~ "scaffold_05",
    grepl("CM008237", scaffold) ~ "scaffold_06",
    grepl("CM008238", scaffold) ~ "scaffold_07",
    grepl("CM008239", scaffold) ~ "scaffold_08",
    grepl("CM008240", scaffold) ~ "scaffold_09",
    grepl("CM008231", scaffold) ~ "scaffold_10",
    TRUE ~ scaffold # Keep the original scaffold if no pattern matches
  ))

#Merge data with scaffold lengths
GC129388_dudPrim_chromosomes <- merge(GC129388_dudPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC129388_dudPrim_chromosomes <- GC129388_dudPrim_chromosomes[, !(names(GC129388_dudPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC129388_dudPrim_chromosomes)[names(GC129388_dudPrim_chromosomes) == "length.x"] <- "window_length"

#Calculate genomic position based on position and chromosome start position
GC129388_dudPrim_chromosomes$genomic_position <- (GC129388_dudPrim_chromosomes$position + GC129388_dudPrim_chromosomes$chromStarts)-1

#Add new columns
GC129388_dudPrim_chromosomes$individual <- "GC129388"
GC129388_dudPrim_chromosomes$country <- "Belgium"
GC129388_dudPrim_chromosomes$ecotype <- "Tidal (SW)"

#Calculate normalized coverage per individual
average_coverage_88D <- mean(GC129388_dudPrim_chromosomes$coverage)
GC129388_dudPrim_chromosomes <- transform(GC129388_dudPrim_chromosomes, normalized_coverage = coverage / average_coverage_88D)


#INDIVIDUAL FROM DUDZELE (SEASONAL, LONG WINGS)
#import scaffold coordinates
scaf_coords_N <- read.delim("~/Documents/thesis/coverage/scaf_coords_N.tsv")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC129394_dudPrim <- GC129394_dudPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC129394_dudPrim <- GC129394_dudPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC129394_dudPrim <- GC129394_dudPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC129394_dudPrim <- GC129394_dudPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC129394_dudPrim <- GC129394_dudPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC129394_dudPrim_chromosomes <- GC129394_dudPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC129394_dudPrim_chromosomes <- GC129394_dudPrim_chromosomes %>%
  mutate(scaffold = case_when(
    grepl("CM008233", scaffold) ~ "scaffold_02",
    grepl("CM008230", scaffold) ~ "scaffold_01",
    grepl("CM008234", scaffold) ~ "scaffold_03",
    grepl("CM008235", scaffold) ~ "scaffold_04",
    grepl("CM008236", scaffold) ~ "scaffold_05",
    grepl("CM008237", scaffold) ~ "scaffold_06",
    grepl("CM008238", scaffold) ~ "scaffold_07",
    grepl("CM008239", scaffold) ~ "scaffold_08",
    grepl("CM008240", scaffold) ~ "scaffold_09",
    grepl("CM008231", scaffold) ~ "scaffold_10",
    TRUE ~ scaffold # Keep the original scaffold if no pattern matches
  ))

#Merge data with scaffold lengths
GC129394_dudPrim_chromosomes <- merge(GC129394_dudPrim_chromosomes, scaf_coords_N, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC129394_dudPrim_chromosomes <- GC129394_dudPrim_chromosomes[, !(names(GC129394_dudPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC129394_dudPrim_chromosomes)[names(GC129394_dudPrim_chromosomes) == "length.x"] <- "window_length"

#Calculate genomic position based on position and chromosome start position
GC129394_dudPrim_chromosomes$genomic_position <- (GC129394_dudPrim_chromosomes$position + GC129394_dudPrim_chromosomes$chromStarts)-1

#Add new columns
GC129394_dudPrim_chromosomes$individual <- "GC129388"
GC129394_dudPrim_chromosomes$country <- "Belgium"
GC129394_dudPrim_chromosomes$ecotype <- "Seasonal (LW)"

#Calculate normalized coverage per individual
average_coverage_94D <- mean(GC129394_dudPrim_chromosomes$coverage)
GC129394_dudPrim_chromosomes <- transform(GC129394_dudPrim_chromosomes, normalized_coverage = coverage / average_coverage_94D)

#Merge datasets
dudzele_chromosomes <- rbind(GC129388_dudPrim_chromosomes, GC129394_dudPrim_chromosomes)
#save dataset
write.table(dudzele_chromosomes, "dudzele_chromosomes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Plot
plot_dudPrim <- ggplot(data=dudzele_chromosomes)+geom_hline(yintercept=1,alpha=0.5,color="chartreuse4",linewidth=0.5)+geom_hline(yintercept=0.5,alpha=0.5,color="chartreuse3",linewidth=0.5)+geom_hline(yintercept=0.00001,alpha=0.5,color="chartreuse2",linewidth=0.5)
plot_dudPrim <- plot_dudPrim+geom_line(aes(x=position,y=normalized_coverage,group=individual,color=ecotype))
plot_dudPrim <- plot_dudPrim+facet_grid(chromosome~.)+ylim(0,2)+xlim(0,60000000)
plot_dudPrim <- plot_dudPrim+ggtitle("Dudzele assembly, Belgian re-sequenced individuals")+scale_color_manual(values = c("Seasonal (LW)" = "red", "Tidal (SW)" = "blue"))+theme_classic()
plot_dudPrim

#belgian individuals, mapped against Nieuwpoort assembly (chromosomes)
#import coverage data
GC129388_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC129388.nieuPrim.pcov", header=FALSE, comment.char="#")
GC129394_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC129394.nieuPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC129388_nieuPrim) <- c("region", "coverage", "length")
names(GC129394_nieuPrim) <- c("region", "coverage", "length")


#INDIVIDUAL FROM NIEUWPOORT (TIDAL, SHORT WINGS)
scaf_coords <- read.delim("~/Documents/thesis/coverage/scaf_coords.tsv")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC129388_nieuPrim <- GC129388_nieuPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC129388_nieuPrim <- GC129388_nieuPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC129388_nieuPrim <- GC129388_nieuPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC129388_nieuPrim <- GC129388_nieuPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC129388_nieuPrim <- GC129388_nieuPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC129388_nieuPrim_chromosomes <- GC129388_nieuPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC129388_nieuPrim_chromosomes <- GC129388_nieuPrim_chromosomes %>%
  mutate(scaffold = case_when(
    grepl("CM008233", scaffold) ~ "scaffold_02",
    grepl("CM008230", scaffold) ~ "scaffold_01",
    grepl("CM008234", scaffold) ~ "scaffold_03",
    grepl("CM008235", scaffold) ~ "scaffold_04",
    grepl("CM008236", scaffold) ~ "scaffold_05",
    grepl("CM008237", scaffold) ~ "scaffold_06",
    grepl("CM008238", scaffold) ~ "scaffold_07",
    grepl("CM008239", scaffold) ~ "scaffold_08",
    grepl("CM008240", scaffold) ~ "scaffold_09",
    grepl("CM008231", scaffold) ~ "scaffold_10",
    TRUE ~ scaffold # Keep the original scaffold if no pattern matches
  ))

#Merge data with scaffold lengths
GC129388_nieuPrim_chromosomes <- merge(GC129388_nieuPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC129388_nieuPrim_chromosomes <- GC129388_nieuPrim_chromosomes[, !(names(GC129388_nieuPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC129388_nieuPrim_chromosomes)[names(GC129388_nieuPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC129388_nieuPrim_chromosomes$individual <- "GC129388"
GC129388_nieuPrim_chromosomes$country <- "Belgium"
GC129388_nieuPrim_chromosomes$ecotype <- "Tidal (SW)"

#Calculate normalized coverage per individual
average_coverage_88D <- mean(GC129388_nieuPrim_chromosomes$coverage)
GC129388_nieuPrim_chromosomes <- transform(GC129388_nieuPrim_chromosomes, normalized_coverage = coverage / average_coverage_88D)


#INDIVIDUAL FROM DUDZELE (SEASONAL, LONG WINGS)
#import scaffold coordinates
scaf_coords_N <- read.delim("~/Documents/thesis/coverage/scaf_coords_N.tsv")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC129394_nieuPrim <- GC129394_nieuPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC129394_nieuPrim <- GC129394_nieuPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC129394_nieuPrim <- GC129394_nieuPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC129394_nieuPrim <- GC129394_nieuPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC129394_nieuPrim <- GC129394_nieuPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC129394_nieuPrim_chromosomes <- GC129394_nieuPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC129394_nieuPrim_chromosomes <- GC129394_nieuPrim_chromosomes %>%
  mutate(scaffold = case_when(
    grepl("CM008233", scaffold) ~ "scaffold_02",
    grepl("CM008230", scaffold) ~ "scaffold_01",
    grepl("CM008234", scaffold) ~ "scaffold_03",
    grepl("CM008235", scaffold) ~ "scaffold_04",
    grepl("CM008236", scaffold) ~ "scaffold_05",
    grepl("CM008237", scaffold) ~ "scaffold_06",
    grepl("CM008238", scaffold) ~ "scaffold_07",
    grepl("CM008239", scaffold) ~ "scaffold_08",
    grepl("CM008240", scaffold) ~ "scaffold_09",
    grepl("CM008231", scaffold) ~ "scaffold_10",
    TRUE ~ scaffold # Keep the original scaffold if no pattern matches
  ))

#Merge data with scaffold lengths
GC129394_nieuPrim_chromosomes <- merge(GC129394_nieuPrim_chromosomes, scaf_coords_N, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC129394_nieuPrim_chromosomes <- GC129394_nieuPrim_chromosomes[, !(names(GC129394_nieuPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC129394_nieuPrim_chromosomes)[names(GC129394_nieuPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC129394_nieuPrim_chromosomes$individual <- "GC129388"
GC129394_nieuPrim_chromosomes$country <- "Belgium"
GC129394_nieuPrim_chromosomes$ecotype <- "Seasonal (LW)"

#Calculate normalized coverage per individual
average_coverage_94D <- mean(GC129394_nieuPrim_chromosomes$coverage)
GC129394_nieuPrim_chromosomes <- transform(GC129394_nieuPrim_chromosomes, normalized_coverage = coverage / average_coverage_94D)

#Merge datasets
nieuwpoort_chromosomes <- rbind(GC129388_nieuPrim_chromosomes, GC129394_nieuPrim_chromosomes)
#save dataset
write.table(nieuwpoort_chromosomes, "nieuwpoort_chromosomes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Plot
plot_nieuPrim <- ggplot(data=nieuwpoort_chromosomes)+geom_hline(yintercept=1,alpha=0.5,color="chartreuse4",linewidth=0.5)+geom_hline(yintercept=0.5,alpha=0.5,color="chartreuse3",linewidth=0.5)+geom_hline(yintercept=0.00001,alpha=0.5,color="chartreuse2",linewidth=0.5)
plot_nieuPrim <- plot_nieuPrim+geom_line(aes(x=position,y=normalized_coverage,group=individual,color=ecotype))
plot_nieuPrim <- plot_nieuPrim+facet_grid(chromosome~.)+ylim(0,4)
plot_nieuPrim <- plot_nieuPrim+ggtitle("Nieuwpoort assembly, Belgian re-sequenced individuals")+scale_color_manual(values = c("Seasonal (LW)" = "red", "Tidal (SW)" = "blue"))+theme_classic()
plot_nieuPrim

