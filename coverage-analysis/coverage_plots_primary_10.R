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

#primary assemblies
GC129388_dud10 <- read.delim("~/Documents/thesis/coverage/GC129388_dud10.raw.pcov", header=FALSE, comment.char="#")
GC129389_dud10 <- read.delim("~/Documents/thesis/coverage/GC129389_dud10.raw.pcov", header=FALSE, comment.char="#")
GC129394_dud10 <- read.delim("~/Documents/thesis/coverage/GC129394_dud10.raw.pcov", header=FALSE, comment.char="#")
GC129395_dud10 <- read.delim("~/Documents/thesis/coverage/GC129395_dud10.raw.pcov", header=FALSE, comment.char="#")

GC129388_nieu10 <- read.delim("~/Documents/thesis/coverage/GC129388_nieu10.raw.pcov", header=FALSE, comment.char="#")
GC129389_nieu10 <- read.delim("~/Documents/thesis/coverage/GC129389_nieu10.raw.pcov", header=FALSE, comment.char="#")
GC129394_nieu10 <- read.delim("~/Documents/thesis/coverage/GC129394_nieu10.raw.pcov", header=FALSE, comment.char="#")
GC129395_nieu10 <- read.delim("~/Documents/thesis/coverage/GC129395_nieu10.raw.pcov", header=FALSE, comment.char="#")

names(GC129388_dud10) <- c("region", "coverage", "length")
names(GC129389_dud10) <- c("region", "coverage", "length")
names(GC129394_dud10) <- c("region", "coverage", "length")
names(GC129395_dud10) <- c("region", "coverage", "length")

names(GC129388_nieu10) <- c("region", "coverage", "length")
names(GC129389_nieu10) <- c("region", "coverage", "length")
names(GC129394_nieu10) <- c("region", "coverage", "length")
names(GC129395_nieu10) <- c("region", "coverage", "length")

scaf_coords <- read.delim("~/Documents/thesis/coverage/scaf_coords.tsv")
scaf_coords_N <- read.delim("~/Documents/thesis/coverage/scaf_coords_N.tsv")

##GC129388_dud10
# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC129388_dud10 <- GC129388_dud10 %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC129388_dud10 <- GC129388_dud10 %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC129388_dud10 <- GC129388_dud10 %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC129388_dud10 <- GC129388_dud10 %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC129388_dud10 <- GC129388_dud10 %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Apply the custom function to the "scaffold" column
GC129388_dud10$scaffold <- rename_scaffold(GC129388_dud10$scaffold)

# Sort the data frame by the "scaffold" column
GC129388_dud10 <- GC129388_dud10 %>%
  arrange(scaffold) %>%
  group_by(scaffold)

# Save the dataframe as a TSV file
#write.table(GC129388_dud10, "GC129388_dud10.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Merge data with scaffold lengths
merged_GC129388_dud10 <- merge(GC129388_dud10, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the position column from the merged data frame
merged_GC129388_dud10 <- merged_GC129388_dud10[, !(names(merged_GC129388_dud10) %in% c("length.y"))]

# Rename the position column to match the original column name
names(merged_GC129388_dud10)[names(merged_GC129388_dud10) == "length.x"] <- "window_length"

#Change names
merged_GC129388_dud10$genomic_position <- (merged_GC129388_dud10$position + merged_GC129388_dud10$chromStarts)-1

#Add new columns
merged_GC129388_dud10$individual <- "GC129388"
merged_GC129388_dud10$population <- "Nieuwpoort"

average_coverage <- mean(merged_GC129388_dud10$coverage)
merged_GC129388_dud10 <- transform(merged_GC129388_dud10, normalized_coverage = coverage / average_coverage)

ggplot(merged_GC129388_dud10, aes(x=genomic_position, y = normalized_coverage))+ geom_point()+ theme_minimal() + ggtitle("Mapped against PRIM_DUD (10)")

##GC129389_dud10
# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC129389_dud10 <- GC129389_dud10 %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC129389_dud10 <- GC129389_dud10 %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC129389_dud10 <- GC129389_dud10 %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC129389_dud10 <- GC129389_dud10 %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC129389_dud10 <- GC129389_dud10 %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Custom function to rename scaffolds
rename_scaffold <- function(scaffold) {
  scaffold_number <- as.numeric(sub("scaffold_", "", scaffold))
  ifelse(scaffold_number < 10, 
         sprintf("scaffold_0%d", scaffold_number), 
         sprintf("scaffold_%d", scaffold_number))
}

# Apply the custom function to the "scaffold" column
GC129389_dud10$scaffold <- rename_scaffold(GC129389_dud10$scaffold)

# Sort the data frame by the "scaffold" column
GC129389_dud10 <- GC129389_dud10 %>%
  arrange(scaffold) %>%
  group_by(scaffold)

#Merge data with scaffold lengths
merged_GC129389_dud10 <- merge(GC129389_dud10, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the position column from the merged data frame
merged_GC129389_dud10 <- merged_GC129389_dud10[, !(names(merged_GC129389_dud10) %in% c("length.y"))]

# Rename the position column to match the original column name
names(merged_GC129389_dud10)[names(merged_GC129389_dud10) == "length.x"] <- "window_length"

#Change names
merged_GC129389_dud10$genomic_position <- (merged_GC129389_dud10$position + merged_GC129389_dud10$chromStarts)-1

#Add new columns
merged_GC129389_dud10$individual <- "GC129389"
merged_GC129389_dud10$population <- "Nieuwpoort"

average_coverage <- mean(merged_GC129389_dud10$coverage)
merged_GC129389_dud10 <- transform(merged_GC129389_dud10, normalized_coverage = coverage / average_coverage)

ggplot(merged_GC129389_dud10, aes(x=genomic_position, y = normalized_coverage))+ geom_point()+ theme_minimal() + ggtitle("Mapped against PRIM_DUD (10)")

##GC129394_dud10
# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC129394_dud10 <- GC129394_dud10 %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC129394_dud10 <- GC129394_dud10 %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC129394_dud10 <- GC129394_dud10 %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC129394_dud10 <- GC129394_dud10 %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC129394_dud10 <- GC129394_dud10 %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Custom function to rename scaffolds
rename_scaffold <- function(scaffold) {
  scaffold_number <- as.numeric(sub("scaffold_", "", scaffold))
  ifelse(scaffold_number < 10, 
         sprintf("scaffold_0%d", scaffold_number), 
         sprintf("scaffold_%d", scaffold_number))
}

# Apply the custom function to the "scaffold" column
GC129394_dud10$scaffold <- rename_scaffold(GC129394_dud10$scaffold)

# Sort the data frame by the "scaffold" column
GC129394_dud10 <- GC129394_dud10 %>%
  arrange(scaffold) %>%
  group_by(scaffold)

#Merge data with scaffold lengths
merged_GC129394_dud10 <- merge(GC129394_dud10, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the position column from the merged data frame
merged_GC129394_dud10 <- merged_GC129394_dud10[, !(names(merged_GC129394_dud10) %in% c("length.y"))]

# Rename the position column to match the original column name
names(merged_GC129394_dud10)[names(merged_GC129394_dud10) == "length.x"] <- "window_length"

#Change names
merged_GC129394_dud10$genomic_position <- (merged_GC129394_dud10$position + merged_GC129394_dud10$chromStarts)-1

#Add new columns
merged_GC129394_dud10$individual <- "GC129394"
merged_GC129394_dud10$population <- "Dudzele"

average_coverage <- mean(merged_GC129394_dud10$coverage)
merged_GC129394_dud10 <- transform(merged_GC129394_dud10, normalized_coverage = coverage / average_coverage)

ggplot(merged_GC129394_dud10, aes(x=genomic_position, y = normalized_coverage))+ geom_point()+ theme_minimal() + ggtitle("Mapped against PRIM_DUD (10)")

##GC129395_dud10
# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC129395_dud10 <- GC129395_dud10 %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC129395_dud10 <- GC129395_dud10 %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC129395_dud10 <- GC129395_dud10 %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC129395_dud10 <- GC129395_dud10 %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC129395_dud10 <- GC129395_dud10 %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Custom function to rename scaffolds
rename_scaffold <- function(scaffold) {
  scaffold_number <- as.numeric(sub("scaffold_", "", scaffold))
  ifelse(scaffold_number < 10, 
         sprintf("scaffold_0%d", scaffold_number), 
         sprintf("scaffold_%d", scaffold_number))
}

# Apply the custom function to the "scaffold" column
GC129395_dud10$scaffold <- rename_scaffold(GC129395_dud10$scaffold)

# Sort the data frame by the "scaffold" column
GC129395_dud10 <- GC129395_dud10 %>%
  arrange(scaffold) %>%
  group_by(scaffold)

#Merge data with scaffold lengths
merged_GC129395_dud10 <- merge(GC129395_dud10, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the position column from the merged data frame
merged_GC129395_dud10 <- merged_GC129395_dud10[, !(names(merged_GC129395_dud10) %in% c("length.y"))]

# Rename the position column to match the original column name
names(merged_GC129395_dud10)[names(merged_GC129395_dud10) == "length.x"] <- "window_length"

#Change names
merged_GC129395_dud10$genomic_position <- (merged_GC129395_dud10$position + merged_GC129395_dud10$chromStarts)-1

#Add new columns
merged_GC129395_dud10$individual <- "GC129395"
merged_GC129395_dud10$population <- "Dudzele"

average_coverage <- mean(merged_GC129395_dud10$coverage)
merged_GC129395_dud10 <- transform(merged_GC129395_dud10, normalized_coverage = coverage / average_coverage)

ggplot(merged_GC129395_dud10, aes(x=genomic_position, y = normalized_coverage))+ geom_point()+ theme_minimal() + ggtitle("Mapped against PRIM_DUD (10)")



#Combine everything for DUD10
DUD10 <- rbind(merged_GC129388_dud10, merged_GC129389_dud10, merged_GC129394_dud10, merged_GC129395_dud10)
write.table(DUD10, "DUD10.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
ggplot(DUD10, aes(x=genomic_position, y = normalized_coverage, color=individual))+ geom_smooth(method = "loess", span=0.1)+ theme_minimal() + ggtitle("Mapped against PRIM_DUD (10)")

ggplot(DUD10, aes(x=genomic_position, y = normalized_coverage, color=population))+geom_smooth(method = "loess", span=0.1)+theme_minimal() + ggtitle("Mapped against PRIM_DUD (10)")

DUD10 <- DUD10 %>% arrange(genomic_position)
ggplot(DUD10, aes(x = genomic_position, y = normalized_coverage, color=population)) +
  geom_line() +
  labs(x = "Genomic Position", y = "Coverage") +
  ggtitle("Coverage vs. Genomic Position")

DUD10_scaffold_01 <- subset(DUD10, scaffold=="scaffold_01")
ggplot(DUD10_scaffold_01, aes(x=genomic_position, y = normalized_coverage, color=population))+geom_smooth(method = "loess", span=0.05)+theme_minimal() + ggtitle("Scaffold 1, Mapped against PRIM_DUD (10)")

DUD10_scaffold_05 <- subset(DUD10, scaffold=="scaffold_05")
ggplot(DUD10_scaffold_05, aes(x=genomic_position, y = normalized_coverage, color=population))+geom_smooth(method = "loess", span=0.05)+theme_minimal() + ggtitle("Scaffold 5, Mapped against PRIM_DUD (10)")


##GC129388_nieu10
# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC129388_nieu10 <- GC129388_nieu10 %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC129388_nieu10 <- GC129388_nieu10 %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC129388_nieu10 <- GC129388_nieu10 %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC129388_nieu10 <- GC129388_nieu10 %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC129388_nieu10 <- GC129388_nieu10 %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Apply the custom function to the "scaffold" column
GC129388_nieu10$scaffold <- rename_scaffold(GC129388_nieu10$scaffold)

# Sort the data frame by the "scaffold" column
GC129388_nieu10 <- GC129388_nieu10 %>%
  arrange(scaffold) %>%
  group_by(scaffold)

#Merge data with scaffold lengths
merged_GC129388_nieu10 <- merge(GC129388_nieu10, scaf_coords_N, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the position column from the merged data frame
merged_GC129388_nieu10 <- merged_GC129388_nieu10[, !(names(merged_GC129388_nieu10) %in% c("length.y"))]

# Rename the position column to match the original column name
names(merged_GC129388_nieu10)[names(merged_GC129388_nieu10) == "length.x"] <- "window_length"

#Add new column, genomic position
merged_GC129388_nieu10$genomic_position <- (merged_GC129388_nieu10$position + merged_GC129388_nieu10$chromStarts)-1

#Add new columns
merged_GC129388_nieu10$individual <- "GC129388"
merged_GC129388_nieu10$population <- "Nieuwpoort"

average_coverage <- mean(merged_GC129388_nieu10$coverage)
merged_GC129388_nieu10 <- transform(merged_GC129388_nieu10, normalized_coverage = coverage / average_coverage)

ggplot(merged_GC129388_nieu10, aes(x=genomic_position, y = normalized_coverage))+ geom_point()+ theme_minimal() + ggtitle("Mapped against PRIM_NIEU (10)")


##GC129389_nieu10
# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC129389_nieu10 <- GC129389_nieu10 %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC129389_nieu10 <- GC129389_nieu10 %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC129389_nieu10 <- GC129389_nieu10 %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC129389_nieu10 <- GC129389_nieu10 %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC129389_nieu10 <- GC129389_nieu10 %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Apply the custom function to the "scaffold" column
GC129389_nieu10$scaffold <- rename_scaffold(GC129389_nieu10$scaffold)

# Sort the data frame by the "scaffold" column
GC129389_nieu10 <- GC129389_nieu10 %>%
  arrange(scaffold) %>%
  group_by(scaffold)

#Merge data with scaffold lengths
merged_GC129389_nieu10 <- merge(GC129389_nieu10, scaf_coords_N, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the position column from the merged data frame
merged_GC129389_nieu10 <- merged_GC129389_nieu10[, !(names(merged_GC129389_nieu10) %in% c("length.y"))]

# Rename the position column to match the original column name
names(merged_GC129389_nieu10)[names(merged_GC129389_nieu10) == "length.x"] <- "window_length"

#Add new column, genomic position
merged_GC129389_nieu10$genomic_position <- (merged_GC129389_nieu10$position + merged_GC129389_nieu10$chromStarts)-1

#Add new columns
merged_GC129389_nieu10$individual <- "GC129389"
merged_GC129389_nieu10$population <- "Nieuwpoort"

average_coverage <- mean(merged_GC129389_nieu10$coverage)
merged_GC129389_nieu10 <- transform(merged_GC129389_nieu10, normalized_coverage = coverage / average_coverage)

ggplot(merged_GC129389_nieu10, aes(x=genomic_position, y = normalized_coverage))+ geom_point()+ theme_minimal() + ggtitle("Mapped against PRIM_NIEU (10)")

##GC129394_nieu10
# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC129394_nieu10 <- GC129394_nieu10 %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC129394_nieu10 <- GC129394_nieu10 %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC129394_nieu10 <- GC129394_nieu10 %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC129394_nieu10 <- GC129394_nieu10 %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC129394_nieu10 <- GC129394_nieu10 %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Apply the custom function to the "scaffold" column
GC129394_nieu10$scaffold <- rename_scaffold(GC129394_nieu10$scaffold)

# Sort the data frame by the "scaffold" column
GC129394_nieu10 <- GC129394_nieu10 %>%
  arrange(scaffold) %>%
  group_by(scaffold)

#Merge data with scaffold lengths
merged_GC129394_nieu10 <- merge(GC129394_nieu10, scaf_coords_N, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the position column from the merged data frame
merged_GC129394_nieu10 <- merged_GC129394_nieu10[, !(names(merged_GC129394_nieu10) %in% c("length.y"))]

# Rename the position column to match the original column name
names(merged_GC129394_nieu10)[names(merged_GC129394_nieu10) == "length.x"] <- "window_length"

#Add new column, genomic position
merged_GC129394_nieu10$genomic_position <- (merged_GC129394_nieu10$position + merged_GC129394_nieu10$chromStarts)-1

#Add new columns
merged_GC129394_nieu10$individual <- "GC129394"
merged_GC129394_nieu10$population <- "Dudzele"

average_coverage <- mean(merged_GC129394_nieu10$coverage)
merged_GC129394_nieu10 <- transform(merged_GC129394_nieu10, normalized_coverage = coverage / average_coverage)

ggplot(merged_GC129394_nieu10, aes(x=genomic_position, y = normalized_coverage))+ geom_point()+ theme_minimal() + ggtitle("Mapped against PRIM_NIEU (10)")


##GC129395_nieu10
# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC129395_nieu10 <- GC129395_nieu10 %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC129395_nieu10 <- GC129395_nieu10 %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC129395_nieu10 <- GC129395_nieu10 %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC129395_nieu10 <- GC129395_nieu10 %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC129395_nieu10 <- GC129395_nieu10 %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Apply the custom function to the "scaffold" column
GC129395_nieu10$scaffold <- rename_scaffold(GC129395_nieu10$scaffold)

# Sort the data frame by the "scaffold" column
GC129395_nieu10 <- GC129395_nieu10 %>%
  arrange(scaffold) %>%
  group_by(scaffold)

#Merge data with scaffold lengths
merged_GC129395_nieu10 <- merge(GC129395_nieu10, scaf_coords_N, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the position column from the merged data frame
merged_GC129395_nieu10 <- merged_GC129395_nieu10[, !(names(merged_GC129395_nieu10) %in% c("length.y"))]

# Rename the position column to match the original column name
names(merged_GC129395_nieu10)[names(merged_GC129395_nieu10) == "length.x"] <- "window_length"

#Add new column, genomic position
merged_GC129395_nieu10$genomic_position <- (merged_GC129395_nieu10$position + merged_GC129395_nieu10$chromStarts)-1

#Add new columns
merged_GC129395_nieu10$individual <- "GC129395"
merged_GC129395_nieu10$population <- "Dudzele"

average_coverage <- mean(merged_GC129395_nieu10$coverage)
merged_GC129395_nieu10 <- transform(merged_GC129395_nieu10, normalized_coverage = coverage / average_coverage)

ggplot(merged_GC129395_nieu10, aes(x=genomic_position, y = normalized_coverage))+ geom_point()+ theme_minimal() + ggtitle("Mapped against PRIM_NIEU (10)")


#Combine everything for NIEU10
NIEU10 <- rbind(merged_GC129388_nieu10, merged_GC129389_nieu10, merged_GC129394_nieu10, merged_GC129395_nieu10)
write.table(NIEU10, "NIEU10.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
ggplot(NIEU10, aes(x=genomic_position, y = normalized_coverage, color=individual))+ geom_smooth(method = "loess", span=0.1)+ theme_minimal() + ggtitle("Mapped against PRIM_NIEU (10)")

ggplot(NIEU10, aes(x=genomic_position, y = normalized_coverage, color=population))+geom_smooth(method = "loess", span=0.1)+theme_minimal() + ggtitle("Mapped against PRIM_NIEU (10)")

NIEU10 <- NIEU10 %>% arrange(genomic_position)
ggplot(NIEU10, aes(x = genomic_position, y = normalized_coverage, color=population, group=individual)) +
  geom_line() +
  labs(x = "Genomic Position", y = "Normalized Coverage") +
  ggtitle("Coverage vs. Genomic Position")

