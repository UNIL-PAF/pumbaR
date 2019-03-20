# Load data from MQ and remove contaminants

int_column_pattern <- "Intensity\\.H\\."
alt_int_column_pattern <- "Intensity\\."


#' Load and filter MaxQuant proteinGroups.txt
#'
#' Load the results and remove contaminants.
#'
#' @param proteinGroups_path Path to proteinGroups.txt
#' @examples
#' proteinGroups_path <-
#' system.file("extdata", "Conde_9508_sub.txt", package = "pumbaR")
#' pg <- load_MQ(proteinGroups_path)
#' @export
load_MQ <- function(proteinGroups_path){

  # read proteinGroups.txt
  pg <- utils::read.table(proteinGroups_path, quote="\"", row.names=NULL,
             header=TRUE, sep="\t", fill=TRUE, na.strings=c("Non Num\303\251rique"))

  # remove contaminants and keep only columns of interest
  is_contaminant <- PAFcontaminants::contaminants_MQ(pg)
  keep_columns_names <- c("Majority.protein.IDs", "Gene.names", "Fasta.headers", "Peptides", "Score", "Only.identified.by.site",
                          "Reverse", "Potential.contaminant", "id", "Peptide.IDs", "Peptide.is.razor", "Mol..weight..kDa.")
  keep_columns <- colnames(pg) %in% keep_columns_names

  keep_columns[get_intensity_columns(pg)] <- TRUE

  pg_flt <- pg[! is_contaminant, keep_columns]

  pg_flt
}

#' Get intensities
#'
#' Extract the intensity values and put slice numbers as the column names.
#'
#' @param pg ProteinGroups data.frame.
#' @examples
#' proteinGroups_path <-
#' system.file("extdata", "Conde_9508_sub.txt", package = "pumbaR")
#' pg <- load_MQ(proteinGroups_path)
#' ints <- get_intensities(pg)
#' @export
get_intensities <- function(pg){
  col_ids <- get_intensity_columns(pg)
  slice_nrs <- get_slice_numbers(pg)

  # create a data.frame with the intensities
  res <- pg[, col_ids]
  colnames(res) <- slice_nrs
  res
}

#' Get intensitie columns
#'
#' Extract the Intensity.H or Intensity. columns
#'
#' @param pg ProteinGroups data.frame.
get_intensity_columns <- function(pg){
  # we extract the Intensity.H columns
  col_ids <- grep(int_column_pattern, colnames(pg))

  # if we're not looking at Silac data
  if(length(col_ids) < 10){
    col_ids <- grep(alt_int_column_pattern, colnames(pg))
  }
  col_ids
}

#' Get the slice numbers
#'
#' Extract the slice numbers from the proteinGroups dataframe.
#'
#' @param pg ProteinGroups dataframe.
#' @examples
#' proteinGroups_path <-
#' system.file("extdata", "Conde_9508_sub.txt", package = "pumbaR")
#' pg <- load_MQ(proteinGroups_path)
#' get_slice_numbers(pg)
#' @export
get_slice_numbers <- function(pg){
  col_names <- colnames(pg)
  col_names_flt <- grep(int_column_pattern, col_names, value=TRUE)

  # in case we're not looking at Silac data
  if(length(col_names_flt) < 10){
    col_names_flt <- grep(alt_int_column_pattern, col_names, value=TRUE)
  }
  as.numeric(stringr::str_extract(col_names_flt, "[0-9]+$"))
}



