# EDA - Import RDS file as dataframe
df <- readRDS("dat.rds")

# Load required packages
if (!require(dplyr)) install.packages("dplyr")
library(dplyr)
if (!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)

# Display initial structure
cat("Original dataset dimensions:", dim(df), "\n")
print(paste("Total columns:", ncol(df)))

# Identify metadata columns vs measurement columns
metadata_cols <- grep("^Metadata_", colnames(df), value = TRUE)
measurement_cols <- setdiff(colnames(df), metadata_cols)

cat("Metadata columns:", length(metadata_cols), "\n")
cat("Measurement columns:", length(measurement_cols), "\n")

# COMPOUND ANALYSIS
cat("\n=== COMPOUND CELL COUNT ANALYSIS ===\n")

# MISSING DATA ANALYSIS
cat("\n=== MISSING DATA ANALYSIS ===\n")

# Count missing values by column
missing_counts <- sapply(df, function(x) sum(is.na(x)))
missing_percentages <- sapply(df, function(x) round(sum(is.na(x))/length(x) * 100, 2))

# Create missing data summary
missing_summary <- data.frame(
  Column = names(missing_counts),
  Missing_Count = missing_counts,
  Missing_Percentage = missing_percentages,
  Total_Rows = nrow(df),
  stringsAsFactors = FALSE
)

# Sort by missing percentage (descending)
missing_summary <- missing_summary[order(missing_summary$Missing_Percentage, decreasing = TRUE), ]

# Overall missing data statistics
total_cells <- nrow(df) * ncol(df)
total_missing <- sum(missing_counts)
overall_missing_pct <- round(total_missing / total_cells * 100, 2)

cat("Overall missing data statistics:\n")
cat("  Total cells:", total_cells, "\n")
cat("  Total missing:", total_missing, "\n")
cat("  Overall missing percentage:", overall_missing_pct, "%\n\n")

# Show columns with missing data
cols_with_missing <- missing_summary[missing_summary$Missing_Count > 0, ]
if(nrow(cols_with_missing) > 0) {
  cat("Columns with missing data:\n")
  print(cols_with_missing, row.names = FALSE)
  
  # Separate analysis for metadata vs measurement columns
  metadata_missing <- cols_with_missing[grepl("^Metadata_", cols_with_missing$Column), ]
  measurement_missing <- cols_with_missing[!grepl("^Metadata_", cols_with_missing$Column), ]
  
  cat("\nMetadata columns with missing data:", nrow(metadata_missing), "\n")
  if(nrow(metadata_missing) > 0) {
    print(metadata_missing, row.names = FALSE)
  }
  
  cat("\nMeasurement columns with missing data:", nrow(measurement_missing), "\n")
  if(nrow(measurement_missing) > 0) {
    cat("Top 10 measurement columns with highest missing percentages:\n")
    print(head(measurement_missing, 10), row.names = FALSE)
  }
  
} else {
  cat("No missing data found in the dataset\n")
}

# Missing data patterns analysis
cat("\n=== MISSING DATA PATTERNS ===\n")

# Check for completely missing rows
complete_rows <- complete.cases(df)
missing_rows_count <- sum(!complete_rows)
missing_rows_pct <- round(missing_rows_count / nrow(df) * 100, 2)

cat("Rows with any missing data:", missing_rows_count, "(", missing_rows_pct, "%)\n")
cat("Complete rows (no missing data):", sum(complete_rows), "\n")

# Identify columns with high missing data (>50%)
high_missing_cols <- missing_summary[missing_summary$Missing_Percentage > 50, ]
if(nrow(high_missing_cols) > 0) {
  cat("\nColumns with >50% missing data (candidates for removal):\n")
  print(high_missing_cols, row.names = FALSE)
}

# Check for missing data correlation patterns
if(sum(missing_counts) > 0) {
  # Create missing indicator matrix
  missing_matrix <- is.na(df)
  
  # Find columns that are always missing together
  if(ncol(missing_matrix) > 1) {
    missing_correlations <- cor(missing_matrix * 1, use = "pairwise.complete.obs")
    perfect_missing_cor <- which(missing_correlations == 1 & row(missing_correlations) != col(missing_correlations), arr.ind = TRUE)
    
    if(length(perfect_missing_cor) > 0) {
      cat("\nColumns with perfectly correlated missing patterns:\n")
      for(i in 1:nrow(perfect_missing_cor)) {
        col1 <- rownames(missing_correlations)[perfect_missing_cor[i,1]]
        col2 <- colnames(missing_correlations)[perfect_missing_cor[i,2]]
        cat(sprintf("  %s <-> %s\n", col1, col2))
      }
    }
  }
}

# COMPOUND CELL COUNT RANKING
cat("\n=== COMPOUND RANKING BY CELL COUNT ===\n")

# Use Metadata_EOS as compound identifier
compound_col <- "Metadata_EOS"

if(compound_col %in% colnames(df) && "Metadata_Object_Count" %in% colnames(df)) {
  # Count cells per compound using Metadata_Object_Count
  compound_counts <- df %>%
    group_by(!!sym(compound_col)) %>%
    summarise(
      cell_count = mean(Metadata_Object_Count, na.rm = TRUE),
      cell_count_sd = sd(Metadata_Object_Count, na.rm = TRUE),
      cell_count_min = min(Metadata_Object_Count, na.rm = TRUE),
      cell_count_max = max(Metadata_Object_Count, na.rm = TRUE),
      total_observations = n(),
      .groups = 'drop'
    ) %>%
    arrange(cell_count)
  
  cat("\nCompound cell count summary:\n")
  cat("Total unique compounds:", nrow(compound_counts), "\n")
  cat("Min cells per compound:", min(compound_counts$cell_count), "\n")
  cat("Max cells per compound:", max(compound_counts$cell_count), "\n")
  cat("Mean cells per compound:", round(mean(compound_counts$cell_count), 2), "\n")
  cat("Median cells per compound:", median(compound_counts$cell_count), "\n")

  # Create summary by compound type and store reference lines
  compound_type_summary <- NULL
  dmso_avg <- NA
  nocodazole_avg <- NA  
  tetrandrine_avg <- NA
  dmso_sd <- NA
  nocodazole_sd <- NA
  tetrandrine_sd <- NA
  
  if("Metadata_Compound_type" %in% colnames(df)) {
    cat("\n=== COMPOUND TYPE ANALYSIS ===\n")
    
    compound_type_summary <- df %>%
      group_by(Metadata_Compound_type) %>%
      summarise(
        observations = n(),
        min_cell_count = min(Metadata_Object_Count, na.rm = TRUE),
        max_cell_count = max(Metadata_Object_Count, na.rm = TRUE),
        avg_cell_count = round(mean(Metadata_Object_Count, na.rm = TRUE), 2),
        sd_cell_count = round(sd(Metadata_Object_Count, na.rm = TRUE), 2),
        .groups = 'drop'
      ) %>%
      arrange(desc(observations))
    
    cat("Summary by compound type:\n")
    print(compound_type_summary)
    
    # Extract values for vertical lines
    if("DMSO" %in% compound_type_summary$Metadata_Compound_type) {
      dmso_avg <- compound_type_summary$avg_cell_count[compound_type_summary$Metadata_Compound_type == "DMSO"]
      dmso_sd <- compound_type_summary$sd_cell_count[compound_type_summary$Metadata_Compound_type == "DMSO"]
    }
    if("Nocodazole" %in% compound_type_summary$Metadata_Compound_type) {
      nocodazole_avg <- compound_type_summary$avg_cell_count[compound_type_summary$Metadata_Compound_type == "Nocodazole"]
      nocodazole_sd <- compound_type_summary$sd_cell_count[compound_type_summary$Metadata_Compound_type == "Nocodazole"]
    }
    if("Tetrandrine" %in% compound_type_summary$Metadata_Compound_type) {
      tetrandrine_avg <- compound_type_summary$avg_cell_count[compound_type_summary$Metadata_Compound_type == "Tetrandrine"]
      tetrandrine_sd <- compound_type_summary$sd_cell_count[compound_type_summary$Metadata_Compound_type == "Tetrandrine"]
    }
  }

  # Show compounds with least cells
  cat("\nTop 20 compounds with LEAST cells:\n")
  print(head(compound_counts, 20))
  
  # Create histogram of cell count distribution
  cat("\nCreating histogram of compound cell count distribution...\n")
  
  histogram_plot <- ggplot(compound_counts, aes(x = cell_count)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
    labs(
      title = "Distribution of Cell Counts per Compound",
      x = "Average Number of Cells per Compound",
      y = "Frequency (Number of Compounds)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # Add vertical lines and print color legend to console
  reference_lines <- c()
  
  if(!is.na(dmso_avg)) {
    histogram_plot <- histogram_plot + 
      geom_vline(xintercept = dmso_avg, color = "green", linetype = "dotted", linewidth = 1)
    reference_lines <- c(reference_lines, "Green dotted line: DMSO")
  }
  if(!is.na(nocodazole_avg)) {
    histogram_plot <- histogram_plot + 
      geom_vline(xintercept = nocodazole_avg, color = "purple", linetype = "dotted", linewidth = 1)
    reference_lines <- c(reference_lines, "Purple dotted line: Nocodazole")
  }
  if(!is.na(tetrandrine_avg)) {
    histogram_plot <- histogram_plot + 
      geom_vline(xintercept = tetrandrine_avg, color = "red", linetype = "dotted", linewidth = 1)
    reference_lines <- c(reference_lines, "Red dotted line: Tetrandrine")
  }
  
  # Display the plot
  print(histogram_plot)
  
  # Create histogram of cell count standard deviation distribution
  cat("\nCreating histogram of compound cell count standard deviation distribution...\n")
  
  sd_plot <- ggplot(compound_counts, aes(x = cell_count_sd)) +
    geom_histogram(bins = 30, fill = "darkred", color = "black", alpha = 0.7) +
    labs(
      title = "Distribution of Cell Count Standard Deviation per Compound",
      x = "Standard Deviation of Cell Counts per Compound",
      y = "Frequency (Number of Compounds)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # Add vertical lines for standard deviations and print color legend to console
  sd_reference_lines <- c()
  
  if(!is.na(dmso_sd)) {
    sd_plot <- sd_plot + 
      geom_vline(xintercept = dmso_sd, color = "green", linetype = "dotted", linewidth = 1)
    sd_reference_lines <- c(sd_reference_lines, "Green dotted line: DMSO")
  }
  if(!is.na(nocodazole_sd)) {
    sd_plot <- sd_plot + 
      geom_vline(xintercept = nocodazole_sd, color = "purple", linetype = "dotted", linewidth = 1)
    sd_reference_lines <- c(sd_reference_lines, "Purple dotted line: Nocodazole")
  }
  if(!is.na(tetrandrine_sd)) {
    sd_plot <- sd_plot + 
      geom_vline(xintercept = tetrandrine_sd, color = "red", linetype = "dotted", linewidth = 1)
    sd_reference_lines <- c(sd_reference_lines, "Red dotted line: Tetrandrine")
  }
  
  # Print combined reference line legend to console
  if(length(reference_lines) > 0 || length(sd_reference_lines) > 0) {
    cat("\nReference lines on histograms:\n")
    cat("  Green dotted lines: DMSO\n")
    cat("  Purple dotted lines: Nocodazole\n")
    cat("  Red dotted lines: Tetrandrine\n")
  }
  
  # Display the standard deviation plot
  print(sd_plot)
  
  # Create compound analysis separated by site
  if("Metadata_Site" %in% colnames(df)) {
    cat("\nCreating site-specific compound analysis...\n")
    
    compound_counts_by_site <- df %>%
      group_by(!!sym(compound_col), Metadata_Site) %>%
      summarise(
        cell_count = mean(Metadata_Object_Count, na.rm = TRUE),
        cell_count_sd = sd(Metadata_Object_Count, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Create faceted histogram for mean cell counts by site
    mean_by_site_plot <- ggplot(compound_counts_by_site, aes(x = cell_count)) +
      geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
      facet_wrap(~ Metadata_Site, scales = "free_y", labeller = label_both) +
      labs(
        title = "Distribution of Mean Cell Counts per Compound by Site",
        x = "Average Number of Cells per Compound",
        y = "Frequency (Number of Compounds)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 11, face = "bold")
      )
    
    # Calculate site-specific compound type averages
    site_compound_averages <- df %>%
      filter(Metadata_Compound_type %in% c("DMSO", "Nocodazole", "Tetrandrine")) %>%
      group_by(Metadata_Site, Metadata_Compound_type) %>%
      summarise(
        avg_cell_count = mean(Metadata_Object_Count, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Add vertical lines for site-specific means
    mean_by_site_plot <- mean_by_site_plot + 
      geom_vline(data = site_compound_averages, 
                 aes(xintercept = avg_cell_count, color = Metadata_Compound_type), 
                 linetype = "dotted", linewidth = 1) +
      scale_color_manual(
        values = c("DMSO" = "green", "Nocodazole" = "purple", "Tetrandrine" = "red")
      ) +
      theme(legend.position = "none")
    
    
    # Create faceted histogram for standard deviation by site
    sd_by_site_plot <- ggplot(compound_counts_by_site, aes(x = cell_count_sd)) +
      geom_histogram(bins = 30, fill = "darkred", color = "black", alpha = 0.7) +
      facet_wrap(~ Metadata_Site, scales = "free_y", labeller = label_both) +
      labs(
        title = "Distribution of Cell Count Standard Deviation per Compound by Site",
        x = "Standard Deviation of Cell Counts per Compound",
        y = "Frequency (Number of Compounds)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 11, face = "bold")
      )
    
    # Calculate site-specific compound type standard deviations
    site_compound_sds <- df %>%
      filter(Metadata_Compound_type %in% c("DMSO", "Nocodazole", "Tetrandrine")) %>%
      group_by(Metadata_Site, Metadata_Compound_type) %>%
      summarise(
        sd_cell_count = sd(Metadata_Object_Count, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Add vertical lines for site-specific SDs
    sd_by_site_plot <- sd_by_site_plot + 
      geom_vline(data = site_compound_sds, 
                 aes(xintercept = sd_cell_count, color = Metadata_Compound_type), 
                 linetype = "dotted", linewidth = 1) +
      scale_color_manual(
        values = c("DMSO" = "green", "Nocodazole" = "purple", "Tetrandrine" = "red")
      ) +
      theme(legend.position = "none")
    
    # Print combined reference line legend for site plots
    cat("\nReference lines on site-specific histograms:\n")
    cat("  Green dotted lines: DMSO (site-specific values)\n")
    cat("  Purple dotted lines: Nocodazole (site-specific values)\n")
    cat("  Red dotted lines: Tetrandrine (site-specific values)\n")
    
    # Display the site plots
    print(mean_by_site_plot)
    print(sd_by_site_plot)
    
  } else {
    cat("Warning: Metadata_Site column not found - skipping site analysis\n")
  }
  
  # CREATE WELL POSITION GRID ANALYSIS
  if("Metadata_Row" %in% colnames(df) && "Metadata_Col" %in% colnames(df)) {
    cat("\nCreating well position grid analysis...\n")
    
    # Create unique well identifier
    df_wells <- df %>%
      mutate(Well_ID = paste0(Metadata_Row, sprintf("%02d", as.numeric(Metadata_Col)))) %>%
      mutate(Well_Number = match(Well_ID, sort(unique(Well_ID))))
    
    # Calculate statistics per well
    well_stats <- df_wells %>%
      group_by(Metadata_Row, Metadata_Col, Well_ID, Well_Number) %>%
      summarise(
        avg_cell_count = mean(Metadata_Object_Count, na.rm = TRUE),
        sd_cell_count = sd(Metadata_Object_Count, na.rm = TRUE),
        observations = n(),
        .groups = 'drop'
      ) %>%
      arrange(Well_Number)
    
    # Display well mapping
    cat("\nFirst 20 well mappings:\n")
    print(head(well_stats[c("Well_ID", "Well_Number", "avg_cell_count")], 20))
    
    # Create 16x24 grid (assuming standard 384-well plate: A-P rows, 1-24 columns)
    # Convert row letters to numbers for plotting
    well_stats <- well_stats %>%
      mutate(
        Row_Num = match(Metadata_Row, LETTERS[1:16]),
        Col_Num = as.numeric(Metadata_Col)
      )
    
    # Create the heatmap separated by site
    if("Metadata_Site" %in% colnames(df)) {
      cat("\nCreating 16x24 grid heatmap of average cell counts by site...\n")
      
      # Calculate statistics per well per site
      well_stats_by_site <- df %>%
        mutate(Well_ID = paste0(Metadata_Row, sprintf("%02d", as.numeric(Metadata_Col)))) %>%
        group_by(Metadata_Site, Metadata_Row, Metadata_Col, Well_ID) %>%
        summarise(
          avg_cell_count = mean(Metadata_Object_Count, na.rm = TRUE),
          sd_cell_count = sd(Metadata_Object_Count, na.rm = TRUE),
          observations = n(),
          .groups = 'drop'
        ) %>%
        mutate(
          Row_Num = match(Metadata_Row, LETTERS[1:16]),
          Col_Num = as.numeric(Metadata_Col)
        )
      
      well_heatmap <- ggplot(well_stats_by_site, aes(x = Col_Num, y = Row_Num, fill = avg_cell_count)) +
        geom_tile(color = "white", size = 0.1) +
        facet_wrap(~ Metadata_Site, labeller = label_both) +
        scale_fill_gradient(low = "blue", high = "red", name = "Average\nCell Count") +
        scale_x_continuous(breaks = seq(1, 24, by = 4), expand = c(0,0)) +
        scale_y_continuous(breaks = seq(1, 16, by = 2), labels = LETTERS[seq(1, 16, by = 2)], expand = c(0,0), trans = "reverse") +
        labs(
          title = "Average Cell Count by Well Position (16x24 Grid) by Site",
          x = "Column",
          y = "Row"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          panel.grid = element_blank(),
          strip.text = element_text(size = 11, face = "bold")
        ) +
        coord_fixed()
      
      # Display the heatmap
      print(well_heatmap)
      
      # Create heatmap for standard deviation by site
      cat("\nCreating 16x24 grid heatmap of cell count standard deviation by site...\n")
      
      well_sd_heatmap <- ggplot(well_stats_by_site, aes(x = Col_Num, y = Row_Num, fill = sd_cell_count)) +
        geom_tile(color = "white", size = 0.1) +
        facet_wrap(~ Metadata_Site, labeller = label_both) +
        scale_fill_gradient(low = "blue", high = "red", name = "Cell Count\nStd Dev") +
        scale_x_continuous(breaks = seq(1, 24, by = 4), expand = c(0,0)) +
        scale_y_continuous(breaks = seq(1, 16, by = 2), labels = LETTERS[seq(1, 16, by = 2)], expand = c(0,0), trans = "reverse") +
        labs(
          title = "Cell Count Standard Deviation by Well Position (16x24 Grid) by Site",
          x = "Column",
          y = "Row"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          panel.grid = element_blank(),
          strip.text = element_text(size = 11, face = "bold")
        ) +
        coord_fixed()
      
      # Display the SD heatmap
      print(well_sd_heatmap)
      
      # Print site-specific summary statistics
      site_well_summary <- well_stats_by_site %>%
        group_by(Metadata_Site) %>%
        summarise(
          unique_wells = n(),
          min_cell_count = round(min(avg_cell_count, na.rm = TRUE), 2),
          max_cell_count = round(max(avg_cell_count, na.rm = TRUE), 2),
          mean_cell_count = round(mean(avg_cell_count, na.rm = TRUE), 2),
          min_sd = round(min(sd_cell_count, na.rm = TRUE), 2),
          max_sd = round(max(sd_cell_count, na.rm = TRUE), 2),
          mean_sd = round(mean(sd_cell_count, na.rm = TRUE), 2),
          .groups = 'drop'
        )
      
      cat("\nWell position summary by site:\n")
      print(site_well_summary)
      
    } else {
      cat("\nCreating single 16x24 grid heatmap of average cell counts...\n")
      
      well_heatmap <- ggplot(well_stats, aes(x = Col_Num, y = Row_Num, fill = avg_cell_count)) +
        geom_tile(color = "white", size = 0.1) +
        scale_fill_gradient(low = "blue", high = "red", name = "Average\nCell Count") +
        scale_x_continuous(breaks = 1:24, expand = c(0,0)) +
        scale_y_continuous(breaks = 1:16, labels = LETTERS[1:16], expand = c(0,0), trans = "reverse") +
        labs(
          title = "Average Cell Count by Well Position (16x24 Grid)",
          x = "Column",
          y = "Row"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          panel.grid = element_blank(),
          aspect.ratio = 16/24
        ) +
        coord_fixed()
      
      # Display the heatmap
      print(well_heatmap)
    }
    
    # Print summary statistics
    cat("\nWell position summary:\n")
    cat("Total unique wells:", nrow(well_stats), "\n")
    cat("Well number range:", min(well_stats$Well_Number), "to", max(well_stats$Well_Number), "\n")
    cat("Average cell count range:", round(min(well_stats$avg_cell_count), 2), "to", round(max(well_stats$avg_cell_count), 2), "\n")
    
  } else {
    cat("Warning: Metadata_Row and/or Metadata_Col columns not found - skipping well analysis\n")
  }
} else {
  if(!"Metadata_EOS" %in% colnames(df)) {
    cat("Error: Metadata_EOS column not found in dataset\n")
  }
  if(!"Metadata_Object_Count" %in% colnames(df)) {
    cat("Error: Metadata_Object_Count column not found in dataset\n")
  }
}

# ------- 7 plates separately; raw, mean, sd, standardize by site ------     
{
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(ggplot2); library(stringr)
  })
  
  # =========================
  # Config
  # =========================
  # ============================================================================
  site_to_plot    <- "FMP"     # <-- change to the site you want
  # ============================================================================
  plates_per_site <- 7
  reps_per_site   <- 4
  rep_levels      <- paste0("R", 1:reps_per_site)
  
  controls <- c("DMSO", "Nocodazole", "Tetrandrine")
  control_color_map <- c("DMSO"="green","Nocodazole"="purple","Tetrandrine"="red",
                         "Other"="grey80","Missing"="black")
  missing_fill_color <- "grey30"
  tile_border_width  <- 0.25
  
  # =========================
  # Checks
  # =========================
  stopifnot(exists("df"))
  req_cols <- c("Metadata_Site","Metadata_Row","Metadata_Col",
                "Metadata_Object_Count","Metadata_Compound_type")
  miss <- setdiff(req_cols, names(df))
  if (length(miss)) stop("Missing cols: ", paste(miss, collapse=", "))
  
  # Pick your plate column (edit if needed)
  plate_col <- c("Metadata_Unique_plate","Metadata_Plate_num","Metadata_Plate",
                 "Plate","PlateID","Barcode","Metadata_Barcode")
  plate_col <- plate_col[plate_col %in% names(df)][1]
  if (is.na(plate_col)) stop("Set `plate_col` to your plate ID column.")
  
  # =========================
  # Parse Plate & Replicate from plate string
  # =========================
  df_parsed <- df %>%
    mutate(
      Metadata_Row = as.character(Metadata_Row),
      Metadata_Col = suppressWarnings(as.numeric(Metadata_Col)),
      PlateRaw     = as.character(.data[[plate_col]]),
      .rep_digit   = str_match(PlateRaw, "(?i)(?:^|[_-])R\\s*([1-4])$")[,2],
      Replicate    = ifelse(!is.na(.rep_digit), paste0("R", str_trim(.rep_digit)), NA_character_),
      Plate        = ifelse(!is.na(Replicate),
                            str_replace(PlateRaw, "(?i)[_-]R\\s*[1-4]$", ""),
                            PlateRaw)
    ) %>%
    select(-.rep_digit) %>%
    filter(!is.na(Metadata_Col), Metadata_Col>=1, Metadata_Col<=24,
           Metadata_Row %in% LETTERS[1:16])
  
  to_well_coords <- function(d) d %>% mutate(Row_Num = match(Metadata_Row, LETTERS[1:16]),
                                             Col_Num = Metadata_Col)
  
  # Plate order helper (within THIS site)
  plate_order_for_site <- function(site_name) {
    p <- df_parsed %>% filter(Metadata_Site==site_name) %>% distinct(Plate) %>% pull(Plate)
    if (!length(p)) return(character(0))
    num <- suppressWarnings(as.numeric(str_extract(p, "\\d+")))
    p[order(is.na(num), num, p)]
  }
  
  # =========================
  # Build complete **4 × 7** layout for ONE site
  # =========================
  build_site_grid <- function(site_name) {
    wells <- expand_grid(Metadata_Row=LETTERS[1:16], Metadata_Col=1:24)
    
    # pick exactly 7 plates (pad with placeholders if fewer exist)
    p_order <- plate_order_for_site(site_name)
    sel_plates <- head(p_order, plates_per_site)
    if (length(sel_plates) < plates_per_site) {
      sel_plates <- c(sel_plates,
                      paste0("MISSING_PLATE_", seq_len(plates_per_site - length(sel_plates))))
    }
    
    # skeleton: 7 plates × 4 reps × wells
    skel <- expand_grid(
      Metadata_Site = site_name,
      Plate         = sel_plates,
      Replicate     = rep_levels
    ) %>%
      crossing(wells) %>%
      to_well_coords() %>%
      mutate(
        PlateColIdx = match(Plate, sel_plates),          # columns = plates
        RepRowIdx   = match(Replicate, rep_levels),      # rows    = replicates
        PlateFacet  = factor(sprintf("Plate %d\n%s", PlateColIdx, Plate),
                             levels = sprintf("Plate %d\n%s", 1:plates_per_site, sel_plates)),
        RepFacet    = factor(sprintf("Rep %d\n%s", RepRowIdx, Replicate),
                             levels = sprintf("Rep %d\n%s", 1:reps_per_site, rep_levels))
      )
    
    # actual averages and dominant control (filtered to this site)
    dat <- df_parsed %>%
      filter(Metadata_Site==site_name) %>%
      group_by(Metadata_Site, Plate, Replicate, Metadata_Row, Metadata_Col) %>%
      summarise(avg_cell_count=mean(Metadata_Object_Count, na.rm=TRUE), .groups="drop")
    
    ctrl <- df_parsed %>%
      filter(Metadata_Site==site_name) %>%
      group_by(Metadata_Site, Plate, Replicate, Metadata_Row, Metadata_Col, Metadata_Compound_type) %>%
      summarise(n=dplyr::n(), .groups="drop") %>%
      group_by(Metadata_Site, Plate, Replicate, Metadata_Row, Metadata_Col) %>%
      slice_max(order_by=n, n=1, with_ties=FALSE) %>%
      ungroup() %>%
      mutate(control_label = ifelse(Metadata_Compound_type %in% controls,
                                    Metadata_Compound_type, "Other")) %>%
      select(Metadata_Site, Plate, Replicate, Metadata_Row, Metadata_Col, control_label)
    
    merged <- skel %>%
      left_join(dat,  by=c("Metadata_Site","Plate","Replicate","Metadata_Row","Metadata_Col")) %>%
      left_join(ctrl, by=c("Metadata_Site","Plate","Replicate","Metadata_Row","Metadata_Col"))
    
    if (!"control_label" %in% names(merged)) merged$control_label <- NA_character_
    
    merged %>%
      mutate(
        is_missing    = is.na(avg_cell_count),
        control_label = ifelse(is_missing, "Missing", control_label)
        # keep avg_cell_count = NA so fill uses `missing_fill_color`
      )
  }
  
  site_df <- build_site_grid(site_to_plot)
  
  # =========================
  # 1) RAW (single site) — **4 × 7** (rows = replicates, cols = plates)
  # =========================
  p_raw <- ggplot(site_df, aes(x=Col_Num, y=Row_Num)) +
    geom_tile(aes(fill=avg_cell_count, color=control_label), linewidth=tile_border_width) +
    facet_grid(RepFacet ~ PlateFacet, switch="y", drop=FALSE) +   # <<< 4 × 7 layout
    scale_fill_gradient(low="blue", high="red", name="Avg cell count",
                        na.value=missing_fill_color) +
    scale_color_manual(values=control_color_map, drop=FALSE, name="Border") +
    scale_x_continuous(breaks=seq(1,24,4), expand=c(0,0)) +
    scale_y_continuous(breaks=seq(1,16,2), labels=LETTERS[seq(1,16,2)],
                       trans="reverse", expand=c(0,0)) +
    labs(
      title=paste0("RAW — Site: ", site_to_plot, "  (4 Replicates × 7 Plates)"),
      subtitle="Missing plate–replicate wells are shaded in a distinct color; borders show control type",
      x="Column", y="Row"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid=element_blank(),
          strip.text=element_text(face="bold")) +
    coord_fixed()
  print(p_raw)
  
  # =========================
  # 2) AVG across replicates (single site) — **1 × 7** (one row, 7 columns)
  # =========================
  avg_site <- site_df %>%
    group_by(Metadata_Site, Plate, PlateFacet, Metadata_Row, Metadata_Col, Row_Num, Col_Num) %>%
    summarise(avg_over_reps = mean(avg_cell_count, na.rm=TRUE), .groups="drop")
  
  p_avg <- ggplot(avg_site, aes(x=Col_Num, y=Row_Num, fill=avg_over_reps)) +
    geom_tile(color="white", linewidth=tile_border_width) +
    facet_grid(. ~ PlateFacet, switch="y", drop=FALSE) +          # <<< 1 × 7
    scale_fill_gradient(low="blue", high="red", name="Avg over reps",
                        na.value=missing_fill_color) +
    scale_x_continuous(breaks=seq(1,24,4), expand=c(0,0)) +
    scale_y_continuous(breaks=seq(1,16,2), labels=LETTERS[seq(1,16,2)],
                       trans="reverse", expand=c(0,0)) +
    labs(
      title=paste0("AVG across replicates — Site: ", site_to_plot),
      subtitle="One row × 7 plates; grey indicates missing",
      x="Column", y="Row"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          panel.grid=element_blank(),
          strip.text=element_text(face="bold")) +
    coord_fixed()
  print(p_avg)
  
  # =========================
  # 3) NORMALIZED (per-site z-score of AVG) — **1 × 7**
  # =========================
  # -------- NORMALIZED (per-site z-score of AVG) — robust against SD == 0/NA -----
  # avg_site is the 7×1 dataset built just above (per-plate, per-well averages)
  
  # Helper: safe z-score on a numeric vector
  safe_z <- function(x) {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) {
      # fall back to mean-centering so we still see structure
      return(x - m)
    }
    (x - m) / s
  }
  
  # Compute z per site (single site selected upstream)
  avg_norm <- avg_site %>%
    mutate(z_over_reps = safe_z(avg_over_reps))
  
  # Optional: sanity printout
  site_stats <- avg_site %>%
    summarise(
      n_non_na = sum(!is.na(avg_over_reps)),
      site_mean = mean(avg_over_reps, na.rm = TRUE),
      site_sd   = sd  (avg_over_reps, na.rm = TRUE)
    )
  print(site_stats)
  
  # Plot (unchanged layout: 1 × 7)
  p_norm <- ggplot(avg_norm, aes(x = Col_Num, y = Row_Num, fill = z_over_reps)) +
    geom_tile(color = "white", linewidth = tile_border_width) +
    facet_grid(. ~ PlateFacet, switch = "y", drop = FALSE) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      name = "Z-score (AVG)", na.value = missing_fill_color
    ) +
    scale_x_continuous(breaks = seq(1, 24, 4), expand = c(0, 0)) +
    scale_y_continuous(
      breaks = seq(1, 16, 2),
      labels = LETTERS[seq(1, 16, 2)],
      trans  = "reverse", expand = c(0, 0)
    ) +
    labs(
      title = paste0("NORMALIZED (Z-score of AVG) — Site: ", site_to_plot),
      subtitle = "Per-site z = (AVG − mean)/sd; falls back to mean-centering if sd ≤ 0 or NA",
      x = "Column", y = "Row"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    coord_fixed()
  
  print(p_norm)
  
  
}




