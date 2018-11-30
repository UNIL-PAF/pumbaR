# normalization of experiment

#' Normalize intensities
#'
#' Normalize intensities by dividing every entry by the total sum of
#' intensities over all slices and proteins. Gives back a data-frame
#' with normalized intensities.
#'
#' @param pg ProteinGroups data.frame.
#' @examples
#' proteinGroups_path <-
#' system.file("extdata", "Conde_9508_sub.txt", package = "pumbaR")
#' pg <- load_MQ(proteinGroups_path)
#' ints <- get_intensities(pg)
#' norm_ints <- normalize_intensities(ints)
#' @export
normalize_intensities <- function(ints){
  total_sum <- sum(unlist(lapply(ints, function(x){sum(as.numeric(x))})))
  return (ints/total_sum)
}

#' Get normalize table
#'
#' Gives back the normalized intensities together with some informations
#' about the proteins (proteinAC, gene name, etc.).
#'
#' @param pg ProteinGroups data.frame.
#' @examples
#' proteinGroups_path <-
#' system.file("extdata", "Conde_9508_sub.txt", package = "pumbaR")
#' pg <- load_MQ(proteinGroups_path)
#' normalized_table <- get_normalized_table(pg)
#' @export
get_normalized_table <- function(pg){
  ints <- get_intensities(pg)
  norm_ints <- normalize_intensities(ints)

  # add information columns
  info_cols <- c("Majority.protein.IDs", "Gene.names", "Mol..weight..kDa.")
  norm_pg <- cbind(pg[,info_cols], norm_ints)

  # adapt the colnames
  colnames(norm_pg) <- c(info_cols, paste0("Intensity.Norm.", colnames(norm_ints)))
  norm_pg
}
