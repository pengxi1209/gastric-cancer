#!/usr/bin/env Rscript
library(gGnome)

dir.create("ggnome/result", showWarnings = FALSE, recursive = TRUE)

jabba_data = readRDS("ggnome/rds_dir/${sample}_bedpe_jabba.rds")
genome_data = gG(jabba = jabba_data)
genome_data = events(genome_data, verbose = TRUE)

event_counts = table(genome_data$meta$event$type)

counts_df = as.data.frame(t(as.matrix(event_counts)))

result_df = data.frame(
  sample = "${sample}_bedpe",
  counts_df,
  row.names = NULL
)


write.table(result_df, file = "ggnome/result/${sample}_bedpe_event_counts.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(gGnome)
})

analyze_complex_sv <- function(genome_bedpe, sample_name, file_type, output_dir) {
  # Extract events table
  events_table <- genome_bedpe$meta$event

  if (is.null(events_table) || nrow(events_table) == 0) {
    warning(paste("No events found for sample:", sample_name, "with file type:", file_type))
    return(NULL)
  }

  all_sv_details <- data.frame(
    Sample = character(),
    File_Type = character(),
    Event_Type = character(),
    Event_ID = integer(),
    Chromosome = character(),
    Start = numeric(),
    End = numeric(),
    Length = numeric(),
    Strand = character(),  
    stringsAsFactors = FALSE
  )

  for (event_type in unique(events_table$type)) {
    current_events <- events_table[events_table$type == event_type,]

    for (i in 1:nrow(current_events)) {
      footprint <- current_events$footprint[i]
      event_id <- current_events$ev.id[i]

      if (is.na(footprint) || footprint == "") {
        next
      }

      if (grepl(",", footprint)) {
        # For events with multiple coordinates (like tyfonas)
        coord_parts <- unlist(strsplit(footprint, ","))

        for (coord in coord_parts) {
          strand <- ifelse(grepl("\\+$", coord), "+",
                          ifelse(grepl("-$", coord), "-", "."))

          coord <- gsub("[+\\-]$", "", coord)

          if (!grepl(":", coord) || !grepl("-", coord)) {
            warning(paste("Invalid coordinate format:", coord, "for sample:", sample_name))
            next
          }

          chrom <- gsub(":.*", "", coord)
          range <- gsub(".*:", "", coord)
          start_pos <- as.numeric(gsub("-.*", "", range))
          end_pos <- as.numeric(gsub(".*-", "", range))

          all_sv_details <- rbind(all_sv_details, data.frame(
            Sample = sample_name,
            File_Type = file_type,
            Event_Type = event_type,
            Event_ID = event_id,
            Chromosome = chrom,
            Start = start_pos,
            End = end_pos,
            Length = end_pos - start_pos + 1,
            Strand = strand
          ))
        }
      } else {
        strand <- ifelse(grepl("\\+$", footprint), "+",
                        ifelse(grepl("-$", footprint), "-", "."))

        footprint <- gsub("[+\\-]$", "", footprint)
        if (!grepl(":", footprint) || !grepl("-", footprint)) {
          warning(paste("Invalid coordinate format:", footprint, "for sample:", sample_name))
          next
        }

        chrom <- gsub(":.*", "", footprint)
        range <- gsub(".*:", "", footprint)
        start_pos <- as.numeric(gsub("-.*", "", range))
        end_pos <- as.numeric(gsub(".*-", "", range))

        all_sv_details <- rbind(all_sv_details, data.frame(
          Sample = sample_name,
          File_Type = file_type,
          Event_Type = event_type,
          Event_ID = event_id,
          Chromosome = chrom,
          Start = start_pos,
          End = end_pos,
          Length = end_pos - start_pos + 1,
          Strand = strand
        ))
      }
    }
  }

  if (nrow(all_sv_details) == 0) {
    warning(paste("No valid coordinates found for sample:", sample_name, "with file type:", file_type))
    return(NULL)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  output_file <- file.path(output_dir, paste0(sample_name, "_", file_type, "_sv_coordinates.txt"))

  write.table(all_sv_details, file = output_file, row.names = FALSE, sep = "\t", quote = FALSE)

  cat("Processed", sample_name, file_type, "- Results written to:", output_file, "\n")

  return(all_sv_details)
}


rds_file <- "ggnome/rds_dir/${sample}_bedpe_jabba.rds"
output_dir <- "ggnome/coordinate/result/bedpe_coordinates"

if (!file.exists(rds_file)) {
  stop("File not found: ", rds_file)
}

jabba_bedpe <- readRDS(rds_file)

genome_bedpe <- gG(jabba = jabba_bedpe)
genome_bedpe <- events(genome_bedpe, verbose = FALSE)

analyze_complex_sv(genome_bedpe, "${sample}", "bedpe", output_dir)

cat("Processing completed successfully.\n")

