library(ggplot2)
library(ggfortify)
library(dplyr)

setwd("/Users/mariamadrid/Documents/thesis/coverage")

scaf_coords <- read.delim("~/Documents/thesis/coverage/scaf_coords.tsv")
scaf_coords_N <- read.delim("~/Documents/thesis/coverage/scaf_coords_N.tsv")

#Spanish individuals, mapped against Dudzele assembly (chromosomes)
#import coverage data, individual GC136108
GC136108_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136108.dudPrim.pcov", header=FALSE, comment.char="#")
GC136108_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136108.dudPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136108_dudPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136108_dudPrim <- GC136108_dudPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136108_dudPrim <- GC136108_dudPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136108_dudPrim <- GC136108_dudPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136108_dudPrim <- GC136108_dudPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136108_dudPrim <- GC136108_dudPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136108_dudPrim_chromosomes <- GC136108_dudPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136108_dudPrim_chromosomes <- GC136108_dudPrim_chromosomes %>%
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
GC136108_dudPrim_chromosomes <- merge(GC136108_dudPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136108_dudPrim_chromosomes <- GC136108_dudPrim_chromosomes[, !(names(GC136108_dudPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136108_dudPrim_chromosomes)[names(GC136108_dudPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136108_dudPrim_chromosomes$individual <- "GC136108"
GC136108_dudPrim_chromosomes$country <- "Spain"
GC136108_dudPrim_chromosomes$ecotype <- "Tidal (SW)"

#Calculate normalized coverage per individual
average_coverage_08D <- mean(GC136108_dudPrim_chromosomes$coverage)
GC136108_dudPrim_chromosomes <- transform(GC136108_dudPrim_chromosomes, normalized_coverage = coverage / average_coverage_08D)

#import coverage data, individual GC136107
GC136107_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136107.dudPrim.pcov", header=FALSE, comment.char="#")
GC136107_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136107.dudPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136107_dudPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136107_dudPrim <- GC136107_dudPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136107_dudPrim <- GC136107_dudPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136107_dudPrim <- GC136107_dudPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136107_dudPrim <- GC136107_dudPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136107_dudPrim <- GC136107_dudPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136107_dudPrim_chromosomes <- GC136107_dudPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136107_dudPrim_chromosomes <- GC136107_dudPrim_chromosomes %>%
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
GC136107_dudPrim_chromosomes <- merge(GC136107_dudPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136107_dudPrim_chromosomes <- GC136107_dudPrim_chromosomes[, !(names(GC136107_dudPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136107_dudPrim_chromosomes)[names(GC136107_dudPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136107_dudPrim_chromosomes$individual <- "GC136107"
GC136107_dudPrim_chromosomes$country <- "Spain"
GC136107_dudPrim_chromosomes$ecotype <- "Tidal (SW)"

#Calculate normalized coverage per individual
average_coverage_07D <- mean(GC136107_dudPrim_chromosomes$coverage)
GC136107_dudPrim_chromosomes <- transform(GC136107_dudPrim_chromosomes, normalized_coverage = coverage / average_coverage_07D)



#import coverage data, individual GC136109
GC136109_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136109.dudPrim.pcov", header=FALSE, comment.char="#")
GC136109_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136109.dudPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136109_dudPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136109_dudPrim <- GC136109_dudPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136109_dudPrim <- GC136109_dudPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136109_dudPrim <- GC136109_dudPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136109_dudPrim <- GC136109_dudPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136109_dudPrim <- GC136109_dudPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136109_dudPrim_chromosomes <- GC136109_dudPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136109_dudPrim_chromosomes <- GC136109_dudPrim_chromosomes %>%
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
GC136109_dudPrim_chromosomes <- merge(GC136109_dudPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136109_dudPrim_chromosomes <- GC136109_dudPrim_chromosomes[, !(names(GC136109_dudPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136109_dudPrim_chromosomes)[names(GC136109_dudPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136109_dudPrim_chromosomes$individual <- "GC136109"
GC136109_dudPrim_chromosomes$country <- "Spain"
GC136109_dudPrim_chromosomes$ecotype <- "Seasonal (LW)"

#Calculate normalized coverage per individual
average_coverage_09D <- mean(GC136109_dudPrim_chromosomes$coverage)
GC136109_dudPrim_chromosomes <- transform(GC136109_dudPrim_chromosomes, normalized_coverage = coverage / average_coverage_09D)

#import coverage data, individual GC136110
GC136110_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136110.dudPrim.pcov", header=FALSE, comment.char="#")
GC136110_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136110.dudPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136110_dudPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136110_dudPrim <- GC136110_dudPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136110_dudPrim <- GC136110_dudPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136110_dudPrim <- GC136110_dudPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136110_dudPrim <- GC136110_dudPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136110_dudPrim <- GC136110_dudPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136110_dudPrim_chromosomes <- GC136110_dudPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136110_dudPrim_chromosomes <- GC136110_dudPrim_chromosomes %>%
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
GC136110_dudPrim_chromosomes <- merge(GC136110_dudPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136110_dudPrim_chromosomes <- GC136110_dudPrim_chromosomes[, !(names(GC136110_dudPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136110_dudPrim_chromosomes)[names(GC136110_dudPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136110_dudPrim_chromosomes$individual <- "GC136110"
GC136110_dudPrim_chromosomes$country <- "Spain"
GC136110_dudPrim_chromosomes$ecotype <- "Seasonal (LW)"

#Calculate normalized coverage per individual
average_coverage_10D <- mean(GC136110_dudPrim_chromosomes$coverage)
GC136110_dudPrim_chromosomes <- transform(GC136110_dudPrim_chromosomes, normalized_coverage = coverage / average_coverage_10D)

#import coverage data, individual GC136111
GC136111_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136111.dudPrim.pcov", header=FALSE, comment.char="#")
GC136111_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136111.dudPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136111_dudPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136111_dudPrim <- GC136111_dudPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136111_dudPrim <- GC136111_dudPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136111_dudPrim <- GC136111_dudPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136111_dudPrim <- GC136111_dudPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136111_dudPrim <- GC136111_dudPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136111_dudPrim_chromosomes <- GC136111_dudPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136111_dudPrim_chromosomes <- GC136111_dudPrim_chromosomes %>%
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
GC136111_dudPrim_chromosomes <- merge(GC136111_dudPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136111_dudPrim_chromosomes <- GC136111_dudPrim_chromosomes[, !(names(GC136111_dudPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136111_dudPrim_chromosomes)[names(GC136111_dudPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136111_dudPrim_chromosomes$individual <- "GC136111"
GC136111_dudPrim_chromosomes$country <- "Spain"
GC136111_dudPrim_chromosomes$ecotype <- "Seasonal (LW)"

#Calculate normalized coverage per individual
average_coverage_11D <- mean(GC136111_dudPrim_chromosomes$coverage)
GC136111_dudPrim_chromosomes <- transform(GC136111_dudPrim_chromosomes, normalized_coverage = coverage / average_coverage_11D)

#import coverage data, individual GC136112
GC136112_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136112.dudPrim.pcov", header=FALSE, comment.char="#")
GC136112_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136112.dudPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136112_dudPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136112_dudPrim <- GC136112_dudPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136112_dudPrim <- GC136112_dudPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136112_dudPrim <- GC136112_dudPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136112_dudPrim <- GC136112_dudPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136112_dudPrim <- GC136112_dudPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136112_dudPrim_chromosomes <- GC136112_dudPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136112_dudPrim_chromosomes <- GC136112_dudPrim_chromosomes %>%
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
GC136112_dudPrim_chromosomes <- merge(GC136112_dudPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136112_dudPrim_chromosomes <- GC136112_dudPrim_chromosomes[, !(names(GC136112_dudPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136112_dudPrim_chromosomes)[names(GC136112_dudPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136112_dudPrim_chromosomes$individual <- "GC136112"
GC136112_dudPrim_chromosomes$country <- "Spain"
GC136112_dudPrim_chromosomes$ecotype <- "Seasonal (LW)"

#Calculate normalized coverage per individual
average_coverage_12D <- mean(GC136112_dudPrim_chromosomes$coverage)
GC136112_dudPrim_chromosomes <- transform(GC136112_dudPrim_chromosomes, normalized_coverage = coverage / average_coverage_12D)

#import coverage data, individual GC136113
GC136113_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136113.dudPrim.pcov", header=FALSE, comment.char="#")
GC136113_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136113.dudPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136113_dudPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136113_dudPrim <- GC136113_dudPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136113_dudPrim <- GC136113_dudPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136113_dudPrim <- GC136113_dudPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136113_dudPrim <- GC136113_dudPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136113_dudPrim <- GC136113_dudPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136113_dudPrim_chromosomes <- GC136113_dudPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136113_dudPrim_chromosomes <- GC136113_dudPrim_chromosomes %>%
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
GC136113_dudPrim_chromosomes <- merge(GC136113_dudPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136113_dudPrim_chromosomes <- GC136113_dudPrim_chromosomes[, !(names(GC136113_dudPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136113_dudPrim_chromosomes)[names(GC136113_dudPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136113_dudPrim_chromosomes$individual <- "GC136113"
GC136113_dudPrim_chromosomes$country <- "Spain"
GC136113_dudPrim_chromosomes$ecotype <- "Seasonal (LW)"

#Calculate normalized coverage per individual
average_coverage_13D <- mean(GC136113_dudPrim_chromosomes$coverage)
GC136113_dudPrim_chromosomes <- transform(GC136113_dudPrim_chromosomes, normalized_coverage = coverage / average_coverage_13D)

#import coverage data, individual GC136114
GC136114_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136114.dudPrim.pcov", header=FALSE, comment.char="#")
GC136114_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136114.dudPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136114_dudPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136114_dudPrim <- GC136114_dudPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136114_dudPrim <- GC136114_dudPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136114_dudPrim <- GC136114_dudPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136114_dudPrim <- GC136114_dudPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136114_dudPrim <- GC136114_dudPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136114_dudPrim_chromosomes <- GC136114_dudPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136114_dudPrim_chromosomes <- GC136114_dudPrim_chromosomes %>%
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
GC136114_dudPrim_chromosomes <- merge(GC136114_dudPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136114_dudPrim_chromosomes <- GC136114_dudPrim_chromosomes[, !(names(GC136114_dudPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136114_dudPrim_chromosomes)[names(GC136114_dudPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136114_dudPrim_chromosomes$individual <- "GC136114"
GC136114_dudPrim_chromosomes$country <- "Spain"
GC136114_dudPrim_chromosomes$ecotype <- "Tidal (SW)"

#Calculate normalized coverage per individual
average_coverage_14D <- mean(GC136114_dudPrim_chromosomes$coverage)
GC136114_dudPrim_chromosomes <- transform(GC136114_dudPrim_chromosomes, normalized_coverage = coverage / average_coverage_14D)

#import coverage data, individual GC136115
GC136115_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136115.dudPrim.pcov", header=FALSE, comment.char="#")
GC136115_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136115.dudPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136115_dudPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136115_dudPrim <- GC136115_dudPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136115_dudPrim <- GC136115_dudPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136115_dudPrim <- GC136115_dudPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136115_dudPrim <- GC136115_dudPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136115_dudPrim <- GC136115_dudPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136115_dudPrim_chromosomes <- GC136115_dudPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136115_dudPrim_chromosomes <- GC136115_dudPrim_chromosomes %>%
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
GC136115_dudPrim_chromosomes <- merge(GC136115_dudPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136115_dudPrim_chromosomes <- GC136115_dudPrim_chromosomes[, !(names(GC136115_dudPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136115_dudPrim_chromosomes)[names(GC136115_dudPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136115_dudPrim_chromosomes$individual <- "GC136115"
GC136115_dudPrim_chromosomes$country <- "Spain"
GC136115_dudPrim_chromosomes$ecotype <- "Tidal (SW)"

#Calculate normalized coverage per individual
average_coverage_15D <- mean(GC136115_dudPrim_chromosomes$coverage)
GC136115_dudPrim_chromosomes <- transform(GC136115_dudPrim_chromosomes, normalized_coverage = coverage / average_coverage_15D)


#import coverage data, individual GC136116
GC136116_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136116.dudPrim.pcov", header=FALSE, comment.char="#")
GC136116_dudPrim <- read.delim("~/Documents/thesis/coverage/GC136116.dudPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136116_dudPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136116_dudPrim <- GC136116_dudPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136116_dudPrim <- GC136116_dudPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136116_dudPrim <- GC136116_dudPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136116_dudPrim <- GC136116_dudPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136116_dudPrim <- GC136116_dudPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136116_dudPrim_chromosomes <- GC136116_dudPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136116_dudPrim_chromosomes <- GC136116_dudPrim_chromosomes %>%
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
GC136116_dudPrim_chromosomes <- merge(GC136116_dudPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136116_dudPrim_chromosomes <- GC136116_dudPrim_chromosomes[, !(names(GC136116_dudPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136116_dudPrim_chromosomes)[names(GC136116_dudPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136116_dudPrim_chromosomes$individual <- "GC136116"
GC136116_dudPrim_chromosomes$country <- "Spain"
GC136116_dudPrim_chromosomes$ecotype <- "Tidal (SW)"

#Calculate normalized coverage per individual
average_coverage_16D <- mean(GC136116_dudPrim_chromosomes$coverage)
GC136116_dudPrim_chromosomes <- transform(GC136116_dudPrim_chromosomes, normalized_coverage = coverage / average_coverage_16D)


#Merge datasets
spanish_primary_dudzele <- rbind(GC136107_dudPrim_chromosomes, GC136108_dudPrim_chromosomes, GC136109_dudPrim_chromosomes,GC136110_dudPrim_chromosomes,GC136111_dudPrim_chromosomes,GC136112_dudPrim_chromosomes,GC136113_dudPrim_chromosomes, GC136114_dudPrim_chromosomes,GC136115_dudPrim_chromosomes, GC136116_dudPrim_chromosomes)
#save dataset
write.table(spanish_primary_dudzele, "spanish_primary_dudzele.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Plot
plot_spainDud <- ggplot(data=spanish_primary_dudzele)+geom_hline(yintercept=1,alpha=0.5,color="chartreuse4",linewidth=0.5)+geom_hline(yintercept=0.5,alpha=0.5,color="chartreuse3",linewidth=0.5)+geom_hline(yintercept=0.00001,alpha=0.5,color="chartreuse2",linewidth=0.5)
plot_spainDud <- plot_spainDud+geom_line(aes(x=position,y=normalized_coverage,group=individual,color=ecotype))
plot_spainDud <- plot_spainDud+facet_grid(chromosome~.)+ylim(0,4)
plot_spainDud <- plot_spainDud+ggtitle("Dudzele assembly, Spanish re-sequenced individuals")+scale_color_manual(values = c("Seasonal (LW)" = "red", "Tidal (SW)" = "blue"))+theme_classic()
plot_spainDud