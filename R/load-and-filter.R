# Load data from MQ and remove contaminants

int_column_pattern <- "Intensity\\.H\\."


#' Load and filter MaxQuant proteinGroups.txt
#'
#' Load the results and remove contaminants.
#'
#' @param proteinGroups_path Path to proteinGroups.txt
#' @examples
#' proteinGroups_path <-
#' system.file("extdata", "Conde_9508_sub_2.txt", package = "pumbaR")
#' pg <- load_MQ(proteinGroups_path)
#' @export
load_MQ <- function(proteinGroups_path){

  # read proteinGroups.txt
  pg <- utils::read.table(proteinGroups_path, quote="\"", row.names=NULL,
             header=TRUE, sep="\t", fill=TRUE, na.strings=c("Non Num\303\251rique"))

  # remove contaminants and keep only columns of interest
  is_contaminant <- PAFcontaminants::contaminants_MQ(pg)
  keep_columns_names <- c("Majority.protein.IDs", "Gene.names", "Fasta.headers", "Peptides", "Score", "Only.identified.by.site",
                          "Reverse", "Potential.contaminant", "id", "Peptide.IDs", "Mol..weight..kDa.")
  keep_columns <- colnames(pg) %in% keep_columns_names
  keep_columns[grep("Intensity\\.H\\.", colnames(pg))] <- TRUE
  pg_flt <- pg[! is_contaminant, keep_columns]

  pg_flt
}

#' Get intensities
#'
#' Extract the Intensity.H values and put slice numbers as the column names.
#'
#' @param pg ProteinGroups data.frame.
#' @examples
#' proteinGroups_path <-
#' system.file("extdata", "Conde_9508_sub.txt", package = "pumbaR")
#' pg <- load_MQ(proteinGroups_path)
#' ints <- get_intensities(pg)
#' @export
get_intensities <- function(pg){
  # we extract the Intensity.H columns
  col_ids <- grep(int_column_pattern, colnames(pg))
  col_names <- colnames(pg)[col_ids]
  slice_nrs <- get_slice_numbers(col_names)

  # create a data.frame with the intensities
  res <- pg[, col_ids]
  colnames(res) <- slice_nrs
  res
}

#' Get the slice numbers
#'
#' Extract the slice numbers from the column names
#'
#' @param col_names Column names.
#' @examples
#' proteinGroups_path <-
#' system.file("extdata", "Conde_9508_sub.txt", package = "pumbaR")
#' pg <- load_MQ(proteinGroups_path)
#' get_slice_numbers(colnames(pq))
#' @export
get_slice_numbers <- function(col_names){
  col_names_flt <- grep(int_column_pattern, col_names, value=TRUE)
  match_pattern <- "Intensity[\\.|A-Z|a-z|_]+"
  as.numeric(sub(match_pattern, "", col_names_flt, perl=TRUE))
}



