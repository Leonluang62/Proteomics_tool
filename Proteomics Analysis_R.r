rm(list = ls())

# Install packages
library(data.table)
library(ggplot2)
library(limma)

AddMissingSample <- function(df) {
    dt_long <- copy(df)
    
    # Generate completed combination of Protein-Condition-Replicate
    complete_combinations <- CJ(PG.ProteinAccessions = unique(dt_long$PG.ProteinAccessions),
                               R.Condition = unique(dt_long$R.Condition),
                               R.Replicate = unique(dt_long$R.Replicate))
    
    # Merge to the original long table
    dt_complete <- merge(complete_combinations, dt_long,
                         by = c("PG.ProteinAccessions", "R.Condition", "R.Replicate"),
                         all.x = TRUE)
    
    # Set the value column with same type in case it was changed during transformation
    dt_complete[, R.Replicate := as(R.Replicate, class(dt_long$R.Replicate))]
    
    # Replace NA in "nPep" by 0 (if nPep column exists)
    if ("nPep" %in% names(dt_complete)) {
      dt_complete[is.na(nPep), nPep := 0]
    }
    
    # Calculate nRep: total number of replicates per condition (theoretical maximum)
    total_reps_per_condition <- dt_complete[, .(nRep = length(unique(R.Replicate))), by = R.Condition]
    
    # Add nRep column
    dt_complete <- merge(dt_complete, total_reps_per_condition, 
                         by = "R.Condition", all.x = TRUE)
    
    # Calculate nObsProt: number of non-NA PG.Quantity per protein per condition
    obs_counts <- dt_complete[, .(nObsProt = sum(!is.na(PG.Quantity))), 
                              by = .(PG.ProteinAccessions, R.Condition)]
    
    # Add nObsProt column
    dt_complete <- merge(dt_complete, obs_counts,
                         by = c("PG.ProteinAccessions", "R.Condition"), all.x = TRUE)
    
    # Calculate completeness
    dt_complete[, completeness := nObsProt / nRep]
        
    return(dt_complete)
}

MultiGeneAnnotation <- function(ProteinAccessions, fasta_map) {
  # Create lookup_tables
  lookup_tables <- list(
    genes = setNames(ifelse(is.na(fasta_map$Genes), "NA", fasta_map$Genes), 
                    fasta_map$UniProtAccession),
    species = setNames(ifelse(is.na(fasta_map$Species), "NA", fasta_map$Species), 
                      fasta_map$UniProtAccession),
    entry_names = setNames(ifelse(is.na(fasta_map$Entry_name), "NA", fasta_map$Entry_name), 
                          fasta_map$UniProtAccession)
  )
  
  # Batch annotation
  results <- lapply(ProteinAccessions, function(pg_string) {
    protein_ids <- unlist(strsplit(pg_string, ";"))
    
    list(
      Genes = paste(lookup_tables$genes[protein_ids], collapse = ";"),
      Species = paste(lookup_tables$species[protein_ids], collapse = ";"),
      Entry_names = paste(lookup_tables$entry_names[protein_ids], collapse = ";")
    )
  })
  
  return(results)
}

GetProteinMatrix <- function(df, fasta) {
    # read fasta file
    fasta_map <- fasta
    
    # Make a copy of rawdata as the protein table
    dt_prot <- copy(df)

    # Percentage of multi-proteins group by rows
    semicolon_count <- sum(grepl(";", dt_prot$PG.ProteinAccessions))
    total_count <- nrow(dt_prot)
    cat("Row contains semicolon:", semicolon_count, "/", total_count, "(", round(semicolon_count/total_count*100, 2), "%)\n")

    # calculate peptide number of each protein
    dt_prot[, nPep := uniqueN(EG.ModifiedPeptide), by = c("PG.ProteinAccessions", "R.Condition", "R.Replicate")]
    
    # Extract protein table removing FG groups
    dt_prot <- unique(dt_prot[, .(PG.ProteinAccessions, 
                                  PG.Quantity, 
                                  R.Condition, 
                                  R.Replicate,
                                  nPep
                                 )])
    
    # Add the missing samples to long table. Long table doesn't contains samples without protein quantity
    dt_prot = AddMissingSample(dt_prot)

    # Map gene with species
    dt_prot <- merge(dt_prot,
                     fasta_map[, .(UniProtAccession, Species, Genes, Entry_name)],
                     by.x = "PG.ProteinAccessions",
                     by.y = "UniProtAccession",
                     all.x = TRUE)

    
    # Annotate multi ProteinAssessions
    unmapped_rows <- dt_prot[is.na(Species), which = TRUE]
    if (length(unmapped_rows) > 0) {
          results <- MultiGeneAnnotation(
                dt_prot$PG.ProteinAccessions[unmapped_rows], 
                fasta_map)
        
          # Extract annotation results for unmapped rows
          species_results <- sapply(results, function(x) x$Species)
          genes_results <- sapply(results, function(x) x$Genes)
          entry_results <- sapply(results, function(x) x$Entry_names)
                                  
          # Update table
          dt_prot[unmapped_rows, Species := species_results]
          dt_prot[unmapped_rows, Genes := genes_results]
          dt_prot[unmapped_rows, Entry_name := entry_results]
    }

    # Print number of proteotypic groups (single protein)
    pg_total_count <- length(unique(dt_prot$PG.ProteinAccessions))
    pg_proteotypic_values <- unique(dt_prot[!is.na(Genes)]$PG.ProteinAccessions)
    pg_proteotypic_count <- length(pg_proteotypic_values)
    cat("Protein group proteolytic (single protein):", pg_proteotypic_count, "/", pg_total_count, "(", round(pg_proteotypic_count/pg_total_count*100, 2), "%)\n")
    
    # Print number of multi-proteins groups (pg containing ";")
    pg_semicolon_values <- unique(dt_prot[grepl(";", PG.ProteinAccessions)]$PG.ProteinAccessions)
    pg_semicolon_count <- length(pg_semicolon_values)
    cat("Protein group contains semicolon:", pg_semicolon_count, "/", pg_total_count, "(", round(pg_semicolon_count/pg_total_count*100, 2), "%)\n")
    #pg_semicolon_values
    
    # Print number of unmapped groups
    pg_unmapped_values <- unique(dt_prot[is.na(Genes) & !grepl(";", PG.ProteinAccessions)]$PG.ProteinAccessions)
    pg_unmapped_count <- length(pg_unmapped_values)
    cat("Protein group ummapped:", pg_unmapped_count, "/", pg_total_count, "(", round(pg_unmapped_count/pg_total_count*100, 2), "%)\n")
    #unmapped_values

    return(dt_prot)
}

MissingValueFilter <- function(df, nSample, nCondition, nUnipep) {
    # Assign the cutoff value according to the given method 
    sample_n_cutoff = nSample
    condition_n_cutoff = nCondition
    unipep_n_cutoff = nUnipep

    # Filter the missing value based on sample_n_cutoff in each condition
    dt_filtered <- copy(df)

    # Print number of rows before filter
    total_count <- nrow(dt_filtered)
    cat("Number of rows before filter:", total_count, "\n")
    
    # Filter the protein-condition with given number of observed replicates
    condition_valid <- dt_filtered[nObsProt >= sample_n_cutoff, .(PG.ProteinAccessions, R.Condition)]
    condition_valid_count <- condition_valid[, .(valid_count = length(unique(R.Condition))), by = PG.ProteinAccessions]

    # Filter the missing value based on the condition_n_cutoff
    proteins_keep <- condition_valid_count[valid_count >= condition_n_cutoff, PG.ProteinAccessions]  # filter proteins found in at least condition_n_cutoff conditions
    dt_filtered <- dt_filtered[PG.ProteinAccessions %in% proteins_keep]

    # Filter proteins with given number of unique peptides
    if (nUnipep != 0) {
        dt_filtered <- dt_filtered[!is.na(Genes)]
        protein_valid <- dt_filtered[nPep >= unipep_n_cutoff, .(nPep_valid = TRUE), by = PG.ProteinAccessions]
        dt_filtered <- dt_filtered[PG.ProteinAccessions %in% protein_valid$PG.ProteinAccessions]
    }
    
    # Print number of rows after filter
    total_count_filtered <- nrow(dt_filtered)
    cat("Number of rows after filter:", total_count_filtered, "\n")

    return(dt_filtered)
}

ProteinNormalization <- function(df, method, species) {
    if (method == "Median") { 
        # Check if datatable contains the normalizer proteins
        if (nrow(df[Species == species]) == 0) {
            stop("No ", species, " Protein found. The normalization process is aborted")
        }

        # Calculate median of each replicate and normalization factor
        medians <- df[Species == species, .(species_median = median(PG.Quantity, na.rm = TRUE)), by = c("R.Condition", "R.Replicate")]
        reference_median <- medians[1]$species_median
        medians[, normalizer := species_median / reference_median]

        # Normalize the protein quantity by meadian of each replicate
        dt_normalized <- copy(df[medians[, c("R.Condition", "R.Replicate", "normalizer")], on = c("R.Condition", "R.Replicate")])
        dt_normalized <- dt_normalized[, PG.Quantity.Normlized := PG.Quantity / normalizer]

        # Print nomalizer for each condition
        dt_normalized[, .(Normalizer = unique(normalizer)), by = .(R.Condition, R.Replicate)][
            , sprintf("Condition: %-20s\tReplicate: %-10s\tNormalizer: %s", R.Condition, R.Replicate, Normalizer)
        ] |> cat(sep = "\n")

        return(dt_normalized)
    } else {
        cat("No such normalization method. Please use the following methods 'Median'.")
    }
}

Log2AndImputation <- function(df, method, par) {
    # Input validation
    valid_methods <- c("Min", "NormDist_Fixed", "NormDist_Quantile")
    if (!method %in% valid_methods) {
        stop("Invalid method. Please use: ", paste(valid_methods, collapse = ", "))
    }
    
    # Parameter validation
    if (!is.numeric(par)) {
        stop("Parameter 'par' must be numeric")
    }
    if (method == "Min" && (par <= 0 || par >= 1)) {
        stop("For Min method, parameter must be between 0 and 1")
    }
    if (method == "NormDist_Quantile" && (par <= 0 || par >= 1)) {
        stop("For NormDist_Quantile method, parameter must be between 0 and 1")
    }
    
    # Determine quantity column
    value_col <- if ("PG.Quantity.Normlized" %in% names(df)) {
        "PG.Quantity.Normlized"
    } else if ("PG.Quantity" %in% names(df)) {
        "PG.Quantity"
    } else {
        stop("Neither 'PG.Quantity' nor 'PG.Quantity.Normalized' found")
    }
    
    # Create copy and mark imputed values
    dt_imputed <- copy(df)
    dt_imputed[, Imputed := is.na(get(value_col))]
    
    # Define imputation function for min-based methods
    min_based_imputation <- function(data, factor) {
        data[, log2Quantity := {
            values <- .SD[[value_col]]   # Get values of given column to be imputed
            if (any(is.na(values))) {   # Impute only when NA exists to save calculation time
                min_val <- min(values, na.rm = TRUE)
                if (min_val <= 0) {
                    warning("Non-positive values found before log2 transformation")
                }
                values[is.na(values)] <- min_val * factor
            }
            log2(values)
        }, by = PG.ProteinAccessions, .SDcols = value_col]  # .SDcols = value_col 使用 'PG.Quantity' or 'PG.Quantity.Normalized' 做计算
    }
    
    # Define function for normal distribution based imputation
    norm_dist_imputation <- function(data, mean_method, factor) {
        # First log2 transform
        data[, log2Quantity := {
            if (any(.SD[[value_col]] <= 0, na.rm = TRUE)) {
                warning("Non-positive values found before log2 transformation")
            }
            log2(.SD[[value_col]])
        }]
        
        # Then perform normal distribution imputation
        data[, log2Quantity := {
            values <- .SD[["log2Quantity"]]   # Get values of given column to be imputed
            if (any(is.na(values))) {    # Impute only when NA exists to save calculation time
                mean_imputed <- if(mean_method == "fixed") {
                    mean(values, na.rm = TRUE) - factor
                } else if(mean_method == "quantile") {
                    quantile(values, probs = factor, na.rm = TRUE)
                }
                sd_imputed <- sd(values, na.rm = TRUE) * 0.3
                n_na <- sum(is.na(values))
                values[is.na(values)] <- rnorm(n_na, mean_imputed, sd_imputed)
            }
            values
        }, by = PG.ProteinAccessions, .SDcols = "log2Quantity"]   # .SDcols = "log2Quantity" 使用 log2 quantity 做计算
    }
    
    # Perform imputation based on method
    if (method == "Min") {
        min_based_imputation(dt_imputed, par)
    } else if (method == "NormDist_Fixed") {
        norm_dist_imputation(dt_imputed, mean_method = "fixed", par)
    } else if (method == "NormDist_Quantile") {
        norm_dist_imputation(dt_imputed, mean_method = "quantile", par)
    }
    
    return(dt_imputed)
}

TwoSampleTest <- function(df, workflow_text, equal_var = TRUE) {
    # Generate all comparisons
    comparisons <- CJ(Cond1 = unique(df$R.Condition),
                      Cond2 = unique(df$R.Condition))[Cond1 < Cond2]
    
    # Analyze each comparison
    all_comparison_results <- lapply(1:nrow(comparisons), function(comp_idx) {
        # Get 2 conditions name of the comparison
        Cond1 <- comparisons[comp_idx, Cond1]
        Cond2 <- comparisons[comp_idx, Cond2]
        
        # Filter data for these two conditions
        subset_df <- df[R.Condition %in% c(Cond1, Cond2)]
        
        # === Prepare data in one step ===
        wide_data <- dcast(subset_df,
                          PG.ProteinAccessions ~ R.Condition + R.Replicate,
                          value.var = "log2Quantity")
        
        # Convert to matrix
        protein_matrix <- as.matrix(wide_data[, -1])
        rownames(protein_matrix) <- wide_data$PG.ProteinAccessions
        
        # Create condition mapping - efficient version
        sample_names <- colnames(protein_matrix)
        conditions <- substring(sample_names, 1, regexpr("_", sample_names) - 1)
        cond1_cols <- which(conditions == Cond1)
        cond2_cols <- which(conditions == Cond2)
        
        # === Fully vectorized statistical calculations ===
        n1 <- length(cond1_cols)
        n2 <- length(cond2_cols)
        
        # Batch calculate means and standard deviations
        mean_cond1 <- rowMeans(protein_matrix[, cond1_cols, drop = FALSE], na.rm = TRUE)
        mean_cond2 <- rowMeans(protein_matrix[, cond2_cols, drop = FALSE], na.rm = TRUE)
        
        # Calculate standard deviations - using more efficient method
        var_cond1 <- rowMeans((protein_matrix[, cond1_cols, drop = FALSE] - mean_cond1)^2, na.rm = TRUE) * n1 / (n1 - 1)
        var_cond2 <- rowMeans((protein_matrix[, cond2_cols, drop = FALSE] - mean_cond2)^2, na.rm = TRUE) * n2 / (n2 - 1)
        sd_cond1 <- sqrt(var_cond1)
        sd_cond2 <- sqrt(var_cond2)
        
        log2FC <- mean_cond1 - mean_cond2
        
        # === Efficient t-test calculations ===
        if (equal_var) {
            # Equal variance t-test (Student's t-test)
            pooled_var <- ((n1-1) * var_cond1 + (n2-1) * var_cond2) / (n1 + n2 - 2)
            se_diff <- sqrt(pooled_var * (1/n1 + 1/n2))
            t_stat <- (mean_cond1 - mean_cond2) / se_diff
            df_t <- n1 + n2 - 2
        } else {
            # Unequal variance t-test (Welch's t-test)
            se_diff <- sqrt(var_cond1/n1 + var_cond2/n2)
            t_stat <- (mean_cond1 - mean_cond2) / se_diff
            
            # Welch-Satterthwaite degrees of freedom
            df_t <- (var_cond1/n1 + var_cond2/n2)^2 / 
                    ((var_cond1/n1)^2/(n1-1) + (var_cond2/n2)^2/(n2-1))
        }
        
        # Batch calculate p-values
        t_p_value <- 2 * pt(abs(t_stat), df_t, lower.tail = FALSE)
        
        # === Limma analysis ===
        # Create design matrix
        group <- factor(ifelse(conditions == Cond1, "A", "B"), levels = c("A", "B"))
        design <- model.matrix(~0 + group)
        colnames(design) <- c("A", "B")
        
        # Perform limma analysis
        fit <- lmFit(protein_matrix, design)
        my.contrasts <- makeContrasts("A-B", levels = design)
        contrast_fit <- contrasts.fit(fit, my.contrasts)
        contrast_fit <- eBayes(contrast_fit)
        
        # Get limma results
        limma_results <- topTable(contrast_fit, adjust = "BH", number = Inf, sort.by = "none")
        
        # === Create final results directly - zero-copy optimization ===
        data.table(
            PG.ProteinAccessions = wide_data$PG.ProteinAccessions,
            cond1 = Cond1,
            cond2 = Cond2,
            t_p_value = t_p_value,
            t_adj_p_value = p.adjust(t_p_value, method = "BH"),
            limma_p_value = limma_results$P.Value,
            limma_adj_p_value = limma_results$adj.P.Val,
            mean_cond1 = mean_cond1,
            mean_cond2 = mean_cond2,
            sd_cond1 = sd_cond1,
            sd_cond2 = sd_cond2,
            log2FC = log2FC,
            Comparison = paste0(Cond1, "_v_", Cond2),
            Workflow = workflow_text
        )
    })
    
    # Efficiently merge results
    return(rbindlist(all_comparison_results))
}

# Function to calculate counts at different p-value thresholds
roc_pvalue_counts <- function(df, p_thresholds, fixed_fc) {
    counts_pvalue <- lapply(p_thresholds, function(p_thresh) {
        true_positive_count <- df[log2FC >= fixed_fc & t_adj_p_value <= p_thresh & grepl("^(HUMAN;?)+$", Species), .N]
        false_positive_count <- df[log2FC >= fixed_fc & t_adj_p_value <= p_thresh & !grepl("^(HUMAN;?)+$", Species), .N]
        data.table(
            p_threshold = p_thresh,
            p_true_positive = true_positive_count,
            p_false_positive = false_positive_count
        )
    })
    rbindlist(counts_pvalue)
}

# Function to calculate counts at different log2FC thresholds
roc_log2fc_counts <- function(df, fc_thresholds, fixed_p) {
    counts_fc <- lapply(fc_thresholds, function(fc_thresh) {
        true_positive_count <- df[log2FC >= fc_thresh & t_adj_p_value <= fixed_p & grepl("^(HUMAN;?)+$", Species), .N]
        false_positive_count <- df[log2FC >= fc_thresh & t_adj_p_value <= fixed_p & !grepl("^(HUMAN;?)+$", Species), .N]
        data.table(
            fc_threshold = fc_thresh,
            fc_true_positive = true_positive_count,
            fc_false_positive = false_positive_count
        )
    })
    rbindlist(counts_fc)
}

CalculationROC <- function(df, fixed_fc, fixed_p) {
    # Create p-value & log2FC thresholds
    down_lim <- min(df$t_adj_p_value, na.rm = TRUE) * 0.1
    p_thresholds <- seq(down_lim, 1, length.out = 100)
    up_lim <- max(df$log2FC, na.rm = TRUE)
    fc_thresholds <- seq(0, up_lim, length.out = 100)
    
    # Calculate counts for each workflow
    roc_pvalue <- rbindlist(lapply(unique(df$Workflow), function(wf) {
        counts <- roc_pvalue_counts(df[Workflow == wf], p_thresholds, fixed_fc)
        counts[, Workflow := wf]
    }))
    
    roc_log2fc <- rbindlist(lapply(unique(df$Workflow), function(wf) {
        counts <- roc_log2fc_counts(df[Workflow == wf], fc_thresholds, fixed_p)
        counts[, Workflow := wf]
    }))
    
    # Calculate TPR and FPR
    roc_pvalue[, `:=`(
        TPR = p_true_positive / max(p_true_positive),
        FPR = p_false_positive / max(p_false_positive)
    ), by = Workflow]
    
    roc_log2fc[, `:=`(
        TPR = fc_true_positive / max(fc_true_positive),
        FPR = fc_false_positive / max(fc_false_positive)
    ), by = Workflow]
    

    return(list(roc_pvalue = roc_pvalue, roc_log2fc = roc_log2fc))
}


# -------------------------------input--------------------------------
# Define input and output path
input_dir <- "C:/Users/LAUX/Downloads"

# Create output directory if it doesn't exist
if (!dir.exists("C:/Users/LAUX/Downloads")) {
  dir.create("C:/Users/LAUX/Downloads")
}

output_dir <- "C:/Users/LAUX/Downloads"

# Read fasta and Expected_FC
fasta_map <- fread("../../databases/fasta_map_4species.csv")

# pre-define the processing methods and parameters
process_workflow <- data.table(
    workflow = c("w1", "w2", "w3", "w4", "w5"),
    nSample = c(3, 3, 3, 3, 3),
    nCondition = c(2, 1, 2, 2, 2),
    nUnipep = c(0, 0, 2, 0, 0),
    normlizer = c("Median", "Median", "Median", "Median", "Median"),
    species = c("HUMAN", "HUMAN", "HUMAN", "HUMAN", "HUMAN"),
    imputbase = c("Min", "Min", "Min", "Min", "NormDist_Fixed"),
    imputfactor = c(0.5, 0.5, 0.5, 0.8, 1.8)
    )

process_workflow

# Choose the process workflow
wk = "w1"

# Get all CSV files in the input folder
input_files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)

# Check if there are CSV files to process
if (length(input_files) == 0) {
    stop("No CSV files found in ./input directory")}

# Read processing parameters from workflow configuration
nSample <- process_workflow[workflow == wk, nSample]
nCondition <- process_workflow[workflow == wk, nCondition]
nUnipep <- process_workflow[workflow == wk, nUnipep]
normalizer <- process_workflow[workflow == wk, normlizer]
species <- process_workflow[workflow == wk, species]
imputbase <- process_workflow[workflow == wk, imputbase]
imputfactor <- process_workflow[workflow == wk, imputfactor]

# Process each file according to given process workflow
for (file_path in input_files) {
    # Get filename without path and extension for output naming
    file_name <- tools::file_path_sans_ext(basename(file_path))
    cat("Processing file:", file_name, "\n")
    
    # Read the data & process
    dt <- fread(file_path)  
    dt_prot <- GetProteinMatrix(dt, fasta_map)
    dt_filtered <- MissingValueFilter(dt_prot, nSample, nCondition, nUnipep)
    dt_normalized <- ProteinNormalization(dt_filtered, normalizer, species)
    dt_imputed <- Log2AndImputation(dt_normalized, imputbase, imputfactor)
    
    # Add workflow label and save the protein data matrix
    dt_imputed[, Workflow := wk]
    pg_report <- file.path(output_dir, paste0(file_name, "_PG report.csv"))
    write.csv(dt_imputed, pg_report, row.names = FALSE)

    dt_imputed_wide <- dcast(dt_imputed,
                          PG.ProteinAccessions + Genes + Workflow ~ R.Condition + R.Replicate,
                          value.var = "log2Quantity")
    pg_matrix <- file.path(output_dir, paste0(file_name, "_PG matrix.csv"))
    write.csv(dt_imputed_wide, pg_matrix, row.names = FALSE)

    # Perform statistic test
    dt_dep <- TwoSampleTest(dt_imputed, wk, TRUE)

    # Add Genes and Species column to "dep" table
    dt_dep <- fasta_map[, .(UniProtAccession, Species, Genes)][dt_dep, on = c("UniProtAccession" = "PG.ProteinAccessions")]
    rows_to_process <- dt_dep[is.na(Species), which = TRUE]
    if (length(rows_to_process) > 0) {
      results <- MultiGeneAnnotation(
        dt_dep$UniProtAccession[rows_to_process], 
        fasta_map
      )
      # Extract annotation results for unmapped rows
      species_results <- sapply(results, function(x) x$Species)
      genes_results <- sapply(results, function(x) x$Genes)
      
      # Update table
      dt_dep[rows_to_process, Species := species_results]
      dt_dep[rows_to_process, Genes := genes_results]
    }

    # Save "dep" matrix
    dep_matrix <- file.path(output_dir, paste0(file_name, "_DEP matrix.csv"))
    write.csv(dt_dep, dep_matrix, row.names = FALSE)
    
}
        

# PG report example:
dt_imputed

# PG matrix example:
dt_imputed_wide

# DEP matrix example:
dt_dep
