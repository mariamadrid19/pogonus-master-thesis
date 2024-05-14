library(ggplot2)
library(ggfortify)
library(dplyr)

setwd("/Users/mariamadrid/Documents/thesis/coverage")

scaf_coords <- read.delim("~/Documents/thesis/coverage/scaf_coords.tsv")
scaf_coords_N <- read.delim("~/Documents/thesis/coverage/scaf_coords_N.tsv")

#Spanish individuals, mapped against Nieuwpoort assembly (chromosomes)
#import coverage data, individual GC136108
GC136108_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136108.nieuPrim.pcov", header=FALSE, comment.char="#")
GC136108_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136108.nieuPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136108_nieuPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136108_nieuPrim <- GC136108_nieuPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136108_nieuPrim <- GC136108_nieuPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136108_nieuPrim <- GC136108_nieuPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136108_nieuPrim <- GC136108_nieuPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136108_nieuPrim <- GC136108_nieuPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136108_nieuPrim_chromosomes <- GC136108_nieuPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136108_nieuPrim_chromosomes <- GC136108_nieuPrim_chromosomes %>%
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
GC136108_nieuPrim_chromosomes <- merge(GC136108_nieuPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136108_nieuPrim_chromosomes <- GC136108_nieuPrim_chromosomes[, !(names(GC136108_nieuPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136108_nieuPrim_chromosomes)[names(GC136108_nieuPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136108_nieuPrim_chromosomes$individual <- "GC136108"
GC136108_nieuPrim_chromosomes$country <- "Spain"
GC136108_nieuPrim_chromosomes$ecotype <- "Tidal (SW)"

#Calculate normalized coverage per individual
average_coverage_08D <- mean(GC136108_nieuPrim_chromosomes$coverage)
GC136108_nieuPrim_chromosomes <- transform(GC136108_nieuPrim_chromosomes, normalized_coverage = coverage / average_coverage_08D)

#import coverage data, individual GC136107
GC136107_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136107.nieuPrim.pcov", header=FALSE, comment.char="#")
GC136107_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136107.nieuPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136107_nieuPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136107_nieuPrim <- GC136107_nieuPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136107_nieuPrim <- GC136107_nieuPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136107_nieuPrim <- GC136107_nieuPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136107_nieuPrim <- GC136107_nieuPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136107_nieuPrim <- GC136107_nieuPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136107_nieuPrim_chromosomes <- GC136107_nieuPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136107_nieuPrim_chromosomes <- GC136107_nieuPrim_chromosomes %>%
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
GC136107_nieuPrim_chromosomes <- merge(GC136107_nieuPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136107_nieuPrim_chromosomes <- GC136107_nieuPrim_chromosomes[, !(names(GC136107_nieuPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136107_nieuPrim_chromosomes)[names(GC136107_nieuPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136107_nieuPrim_chromosomes$individual <- "GC136107"
GC136107_nieuPrim_chromosomes$country <- "Spain"
GC136107_nieuPrim_chromosomes$ecotype <- "Tidal (SW)"

#Calculate normalized coverage per individual
average_coverage_07D <- mean(GC136107_nieuPrim_chromosomes$coverage)
GC136107_nieuPrim_chromosomes <- transform(GC136107_nieuPrim_chromosomes, normalized_coverage = coverage / average_coverage_07D)



#import coverage data, individual GC136109
GC136109_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136109.nieuPrim.pcov", header=FALSE, comment.char="#")
GC136109_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136109.nieuPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136109_nieuPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136109_nieuPrim <- GC136109_nieuPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136109_nieuPrim <- GC136109_nieuPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136109_nieuPrim <- GC136109_nieuPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136109_nieuPrim <- GC136109_nieuPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136109_nieuPrim <- GC136109_nieuPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136109_nieuPrim_chromosomes <- GC136109_nieuPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136109_nieuPrim_chromosomes <- GC136109_nieuPrim_chromosomes %>%
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
GC136109_nieuPrim_chromosomes <- merge(GC136109_nieuPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136109_nieuPrim_chromosomes <- GC136109_nieuPrim_chromosomes[, !(names(GC136109_nieuPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136109_nieuPrim_chromosomes)[names(GC136109_nieuPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136109_nieuPrim_chromosomes$individual <- "GC136109"
GC136109_nieuPrim_chromosomes$country <- "Spain"
GC136109_nieuPrim_chromosomes$ecotype <- "Seasonal (LW)"

#Calculate normalized coverage per individual
average_coverage_09D <- mean(GC136109_nieuPrim_chromosomes$coverage)
GC136109_nieuPrim_chromosomes <- transform(GC136109_nieuPrim_chromosomes, normalized_coverage = coverage / average_coverage_09D)

#import coverage data, individual GC136110
GC136110_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136110.nieuPrim.pcov", header=FALSE, comment.char="#")
GC136110_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136110.nieuPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136110_nieuPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136110_nieuPrim <- GC136110_nieuPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136110_nieuPrim <- GC136110_nieuPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136110_nieuPrim <- GC136110_nieuPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136110_nieuPrim <- GC136110_nieuPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136110_nieuPrim <- GC136110_nieuPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136110_nieuPrim_chromosomes <- GC136110_nieuPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136110_nieuPrim_chromosomes <- GC136110_nieuPrim_chromosomes %>%
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
GC136110_nieuPrim_chromosomes <- merge(GC136110_nieuPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136110_nieuPrim_chromosomes <- GC136110_nieuPrim_chromosomes[, !(names(GC136110_nieuPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136110_nieuPrim_chromosomes)[names(GC136110_nieuPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136110_nieuPrim_chromosomes$individual <- "GC136110"
GC136110_nieuPrim_chromosomes$country <- "Spain"
GC136110_nieuPrim_chromosomes$ecotype <- "Seasonal (LW)"

#Calculate normalized coverage per individual
average_coverage_10D <- mean(GC136110_nieuPrim_chromosomes$coverage)
GC136110_nieuPrim_chromosomes <- transform(GC136110_nieuPrim_chromosomes, normalized_coverage = coverage / average_coverage_10D)

#import coverage data, individual GC136111
GC136111_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136111.nieuPrim.pcov", header=FALSE, comment.char="#")
GC136111_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136111.nieuPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136111_nieuPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136111_nieuPrim <- GC136111_nieuPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136111_nieuPrim <- GC136111_nieuPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136111_nieuPrim <- GC136111_nieuPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136111_nieuPrim <- GC136111_nieuPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136111_nieuPrim <- GC136111_nieuPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136111_nieuPrim_chromosomes <- GC136111_nieuPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136111_nieuPrim_chromosomes <- GC136111_nieuPrim_chromosomes %>%
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
GC136111_nieuPrim_chromosomes <- merge(GC136111_nieuPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136111_nieuPrim_chromosomes <- GC136111_nieuPrim_chromosomes[, !(names(GC136111_nieuPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136111_nieuPrim_chromosomes)[names(GC136111_nieuPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136111_nieuPrim_chromosomes$individual <- "GC136111"
GC136111_nieuPrim_chromosomes$country <- "Spain"
GC136111_nieuPrim_chromosomes$ecotype <- "Seasonal (LW)"

#Calculate normalized coverage per individual
average_coverage_11D <- mean(GC136111_nieuPrim_chromosomes$coverage)
GC136111_nieuPrim_chromosomes <- transform(GC136111_nieuPrim_chromosomes, normalized_coverage = coverage / average_coverage_11D)

#import coverage data, individual GC136112
GC136112_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136112.nieuPrim.pcov", header=FALSE, comment.char="#")
GC136112_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136112.nieuPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136112_nieuPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136112_nieuPrim <- GC136112_nieuPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136112_nieuPrim <- GC136112_nieuPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136112_nieuPrim <- GC136112_nieuPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136112_nieuPrim <- GC136112_nieuPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136112_nieuPrim <- GC136112_nieuPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136112_nieuPrim_chromosomes <- GC136112_nieuPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136112_nieuPrim_chromosomes <- GC136112_nieuPrim_chromosomes %>%
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
GC136112_nieuPrim_chromosomes <- merge(GC136112_nieuPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136112_nieuPrim_chromosomes <- GC136112_nieuPrim_chromosomes[, !(names(GC136112_nieuPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136112_nieuPrim_chromosomes)[names(GC136112_nieuPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136112_nieuPrim_chromosomes$individual <- "GC136112"
GC136112_nieuPrim_chromosomes$country <- "Spain"
GC136112_nieuPrim_chromosomes$ecotype <- "Seasonal (LW)"

#Calculate normalized coverage per individual
average_coverage_12D <- mean(GC136112_nieuPrim_chromosomes$coverage)
GC136112_nieuPrim_chromosomes <- transform(GC136112_nieuPrim_chromosomes, normalized_coverage = coverage / average_coverage_12D)

#import coverage data, individual GC136113
GC136113_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136113.nieuPrim.pcov", header=FALSE, comment.char="#")
GC136113_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136113.nieuPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136113_nieuPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136113_nieuPrim <- GC136113_nieuPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136113_nieuPrim <- GC136113_nieuPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136113_nieuPrim <- GC136113_nieuPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136113_nieuPrim <- GC136113_nieuPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136113_nieuPrim <- GC136113_nieuPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136113_nieuPrim_chromosomes <- GC136113_nieuPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136113_nieuPrim_chromosomes <- GC136113_nieuPrim_chromosomes %>%
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
GC136113_nieuPrim_chromosomes <- merge(GC136113_nieuPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136113_nieuPrim_chromosomes <- GC136113_nieuPrim_chromosomes[, !(names(GC136113_nieuPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136113_nieuPrim_chromosomes)[names(GC136113_nieuPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136113_nieuPrim_chromosomes$individual <- "GC136113"
GC136113_nieuPrim_chromosomes$country <- "Spain"
GC136113_nieuPrim_chromosomes$ecotype <- "Seasonal (LW)"

#Calculate normalized coverage per individual
average_coverage_13D <- mean(GC136113_nieuPrim_chromosomes$coverage)
GC136113_nieuPrim_chromosomes <- transform(GC136113_nieuPrim_chromosomes, normalized_coverage = coverage / average_coverage_13D)

#import coverage data, individual GC136114
GC136114_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136114.nieuPrim.pcov", header=FALSE, comment.char="#")
GC136114_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136114.nieuPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136114_nieuPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136114_nieuPrim <- GC136114_nieuPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136114_nieuPrim <- GC136114_nieuPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136114_nieuPrim <- GC136114_nieuPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136114_nieuPrim <- GC136114_nieuPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136114_nieuPrim <- GC136114_nieuPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136114_nieuPrim_chromosomes <- GC136114_nieuPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136114_nieuPrim_chromosomes <- GC136114_nieuPrim_chromosomes %>%
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
GC136114_nieuPrim_chromosomes <- merge(GC136114_nieuPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136114_nieuPrim_chromosomes <- GC136114_nieuPrim_chromosomes[, !(names(GC136114_nieuPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136114_nieuPrim_chromosomes)[names(GC136114_nieuPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136114_nieuPrim_chromosomes$individual <- "GC136114"
GC136114_nieuPrim_chromosomes$country <- "Spain"
GC136114_nieuPrim_chromosomes$ecotype <- "Tidal (SW)"

#Calculate normalized coverage per individual
average_coverage_14D <- mean(GC136114_nieuPrim_chromosomes$coverage)
GC136114_nieuPrim_chromosomes <- transform(GC136114_nieuPrim_chromosomes, normalized_coverage = coverage / average_coverage_14D)

#import coverage data, individual GC136115
GC136115_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136115.nieuPrim.pcov", header=FALSE, comment.char="#")
GC136115_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136115.nieuPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136115_nieuPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136115_nieuPrim <- GC136115_nieuPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136115_nieuPrim <- GC136115_nieuPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136115_nieuPrim <- GC136115_nieuPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136115_nieuPrim <- GC136115_nieuPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136115_nieuPrim <- GC136115_nieuPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136115_nieuPrim_chromosomes <- GC136115_nieuPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136115_nieuPrim_chromosomes <- GC136115_nieuPrim_chromosomes %>%
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
GC136115_nieuPrim_chromosomes <- merge(GC136115_nieuPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136115_nieuPrim_chromosomes <- GC136115_nieuPrim_chromosomes[, !(names(GC136115_nieuPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136115_nieuPrim_chromosomes)[names(GC136115_nieuPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136115_nieuPrim_chromosomes$individual <- "GC136115"
GC136115_nieuPrim_chromosomes$country <- "Spain"
GC136115_nieuPrim_chromosomes$ecotype <- "Tidal (SW)"

#Calculate normalized coverage per individual
average_coverage_15D <- mean(GC136115_nieuPrim_chromosomes$coverage)
GC136115_nieuPrim_chromosomes <- transform(GC136115_nieuPrim_chromosomes, normalized_coverage = coverage / average_coverage_15D)


#import coverage data, individual GC136116
GC136116_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136116.nieuPrim.pcov", header=FALSE, comment.char="#")
GC136116_nieuPrim <- read.delim("~/Documents/thesis/coverage/GC136116.nieuPrim.pcov", header=FALSE, comment.char="#")

#change column names
names(GC136116_nieuPrim) <- c("region", "coverage", "length")

# Add a new column "scaffold" by extracting scaffold names from the "region" column
GC136116_nieuPrim <- GC136116_nieuPrim %>%
  mutate(scaffold = gsub("_[0-9]+$", "", region))

# Filter out rows without ">" symbol in the "region" column
GC136116_nieuPrim <- GC136116_nieuPrim %>%
  filter(grepl(">", region))

# Remove everything before the last underscore in the "region" column
GC136116_nieuPrim <- GC136116_nieuPrim %>%
  mutate(region = gsub(".*_", "", region))

# Remove ">" symbol from the "scaffold" column
GC136116_nieuPrim <- GC136116_nieuPrim %>%
  mutate(scaffold = gsub(">", "", scaffold))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC136116_nieuPrim <- GC136116_nieuPrim %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Extract rows where the scaffold column contains "CM" and create a new dataframe
GC136116_nieuPrim_chromosomes <- GC136116_nieuPrim %>% 
  filter(grepl("CM", scaffold))

#Rename the scaffolds
GC136116_nieuPrim_chromosomes <- GC136116_nieuPrim_chromosomes %>%
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
GC136116_nieuPrim_chromosomes <- merge(GC136116_nieuPrim_chromosomes, scaf_coords, by.x = "scaffold", by.y = "scaffold", all.x = TRUE)

# Remove the length.y column from the merged data frame
GC136116_nieuPrim_chromosomes <- GC136116_nieuPrim_chromosomes[, !(names(GC136116_nieuPrim_chromosomes) %in% c("length.y"))]

# Rename the length.x column to window_length (from pcov)
names(GC136116_nieuPrim_chromosomes)[names(GC136116_nieuPrim_chromosomes) == "length.x"] <- "window_length"

#Add new columns
GC136116_nieuPrim_chromosomes$individual <- "GC136116"
GC136116_nieuPrim_chromosomes$country <- "Spain"
GC136116_nieuPrim_chromosomes$ecotype <- "Tidal (SW)"

#Calculate normalized coverage per individual
average_coverage_16D <- mean(GC136116_nieuPrim_chromosomes$coverage)
GC136116_nieuPrim_chromosomes <- transform(GC136116_nieuPrim_chromosomes, normalized_coverage = coverage / average_coverage_16D)


#Merge datasets
spanish_primary_nieuwpoort <- rbind(GC136107_nieuPrim_chromosomes, GC136108_nieuPrim_chromosomes, GC136109_nieuPrim_chromosomes,GC136110_nieuPrim_chromosomes,GC136111_nieuPrim_chromosomes,GC136112_nieuPrim_chromosomes,GC136113_nieuPrim_chromosomes, GC136114_nieuPrim_chromosomes,GC136115_nieuPrim_chromosomes, GC136116_nieuPrim_chromosomes)
#save dataset
write.table(spanish_primary_nieuwpoort, "spanish_primary_nieuwpoort.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Plot
plot_spainNieu <- ggplot(data=spanish_primary_nieuwpoort)+geom_hline(yintercept=1,alpha=0.5,color="chartreuse4",linewidth=0.5)+geom_hline(yintercept=0.5,alpha=0.5,color="chartreuse3",linewidth=0.5)+geom_hline(yintercept=0.00001,alpha=0.5,color="chartreuse2",linewidth=0.5)
plot_spainNieu <- plot_spainNieu+geom_line(aes(x=position,y=normalized_coverage,group=individual,color=ecotype))
plot_spainNieu <- plot_spainNieu+facet_grid(chromosome~.)+ylim(0,4)
plot_spainNieu <- plot_spainNieu+ggtitle("Nieuwpoort assembly, Spanish re-sequenced individuals")+scale_color_manual(values = c("Seasonal (LW)" = "red", "Tidal (SW)" = "blue"))+theme_classic()
plot_spainNieu