# Define the folder containing raw files
folder_path <- "P:/Projects/WFB_SIA_2024_Jaitovich_LongCOVID/Lipidomics/RawFile"  # Adjust as needed

# Recursively list all files inside batch folders
raw_files <- list.files(folder_path, full.names = TRUE, recursive = TRUE)  # Get full paths
raw_files_names <- basename(raw_files)  # Extract just the filenames (without folder path)

# Get timestamps for actual raw files, not folders
timestamps <- file.info(raw_files)$mtime  # Modification timestamps

# Create a dataframe with filenames and timestamps
file_info <- data.frame(
  rawfile_name = raw_files_names,  # Just filenames (no path)
  timestamp = timestamps,
  stringsAsFactors = FALSE
)

# Convert timestamp to proper date-time format
file_info$timestamp <- as.POSIXct(file_info$timestamp, origin="1970-01-01")

# Read Metadata dataframe
metadata_file <- "D:/LongCOVID/Metadata/RunOrderMetadata_v2.csv"
Metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

# Merge Metadata with timestamps based on rawfile_name_R
Metadata <- merge(Metadata, file_info, by = "rawfile_name", all.x = TRUE)

# Print updated Metadata with timestamps
print(head(Metadata))

# Save the updated Metadata with timestamps
write.csv(Metadata, file.path("D:/LongCOVID/Metadata", "Updated_Metadata.csv"), row.names = FALSE)

cat("Metadata with timestamps saved successfully!\n")
