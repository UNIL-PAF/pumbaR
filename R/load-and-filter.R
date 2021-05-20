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
load_MQ <- function(proteinGroups_path, ignore_slices = NULL, sample_name = NULL){
  # read proteinGroups.txt
  pg <- utils::read.table(proteinGroups_path, quote="\"", row.names=NULL,
             header=TRUE, sep="\t", fill=TRUE, na.strings=c("Non Num\303\251rique"))

  # remove contaminants and keep only columns of interest
  is_contaminant <- PAFcontaminants::contaminants_MQ(pg)
  keep_columns_names <- c("Majority.protein.IDs", "Gene.names", "Fasta.headers", "Peptides", "Score", "Only.identified.by.site",
                          "Reverse", "Potential.contaminant", "id", "Peptide.IDs", "Peptide.is.razor", "Mol..weight..kDa.")

  keep_columns <- as.numeric(sapply(keep_columns_names, function(x){
    grep(paste0("^", x, "$"), colnames(pg))
    }))

  # we sort the intensity_columns to make sure of their order
  intensity_columns <- get_intensity_columns(pg, sample_name)
  slice_numbers <- get_slice_numbers(pg, sample_name)
  intensity_columns <- intensity_columns[order(slice_numbers)]

  # remove rows with only empty intensities
  empty_row <- apply(pg[, intensity_columns], 1, function(x){all(x == 0)})

  # remove contaminants and only keep necessary columns
  keep_columns <- c(keep_columns, intensity_columns)
  pg_flt <- pg[! (is_contaminant | empty_row), keep_columns]

  # set values to 0 in case the corresponding slice is to ignored
  if(length(ignore_slices > 0)){
    int_cols <- get_intensity_columns(pg_flt, sample_name)
    for(ignore in ignore_slices){
      pg_flt[int_cols[ignore]] <- 0
    }
  }
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
get_intensities <- function(pg, sample_name = NULL){
  col_ids <- get_intensity_columns(pg, sample_name)
  slice_nrs <- get_slice_numbers(pg, sample_name)

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
get_intensity_columns <- function(pg, sample_name = NULL){
  # we extract the Intensity.H columns
  col_ids <- grep(int_column_pattern, colnames(pg))

  # if we're not looking at Silac data
  if(length(col_ids) < 10){
    col_ids <- grep(alt_int_column_pattern, colnames(pg))
  }

  # only the provided sample_name
  if(! is.null(sample_name)){
    flt_ids <- grep(sample_name, colnames(pg)[col_ids])
    col_ids <- col_ids[flt_ids]
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
get_slice_numbers <- function(pg, sample_name = NULL){
  col_names <- colnames(pg)
  col_names_flt <- grep(int_column_pattern, col_names, value=TRUE)

  # in case we're not looking at Silac data
  if(length(col_names_flt) < 10){
    col_names_flt <- grep(alt_int_column_pattern, col_names, value=TRUE)
  }

  # only the provided sample_name
  if(! is.null(sample_name)){
    col_names_flt <- grep(sample_name, col_names_flt, value=TRUE)
  }

  as.numeric(stringr::str_extract(col_names_flt, "[0-9]+$"))
}



