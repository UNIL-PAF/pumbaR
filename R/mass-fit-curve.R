# create mass fitting curve

###############################
# Typical usage:
#
## load and plot original data
# proteinGroups_path <- "/Users/admin/Work/PAF/projects/SliceSILAC/latest/data/Conde_9508/proteinGroups.txt"
# pg <- load_MQ(proteinGroups_path)
# plot_MQ(pg)

## apply the various filtering
# pg_2 <- filter_repeated_entries(pg)
# plot_MQ(pg_2)
# pg_3 <- filter_weak_intensity(pg_2)
# plot_MQ(pg_3)
# pg_4 <- filter_low_densities(pg_3)
# plot_MQ(pg_4)

## make the curve fitting and plot the result
# mass_fit <- fit_curve(pg_4)
# plot_fit(pg_4, mass_fit)
# plot_fit(pg, mass_fit)

#' Filter_and_fit
#'
#' Filter the data and apply a fit on it.
#'
#' @param pg ProteinGroups data.frame.
#' @param low_density_threshold Remove entries from regions with few data points (default is 10).
#' @examples
#' proteinGroups_path <- system.file("extdata", "Conde_9508_sub.txt", package = "pumbaR")
#' pg <- load_MQ(proteinGroups_path)
#' mass_fit <- filter_and_fit(pg, low_density_threshold = 2)
#' @export
filter_and_fit <- function(pg, low_density_threshold = 10){
  pg <- filter_repeated_entries(pg)
  pg <- filter_weak_intensity(pg)
  pg <- filter_low_densities(pg, min_nr_threshold = low_density_threshold)
  fit_curve(pg)
}


#' Fit curve
#'
#' Fit a polynomial curve to the given intensity data.
#'
#' @param pg MaxQuant ProteinGroup data.
fit_curve <- function(pg){
  # get the intensities from the proteinGroups
  ints <- get_intensities(pg)
  ints.weight <- cbind(ints, mol.weight = log10(pg$Mol..weight..kDa.))
  ints.long <- reshape2::melt(ints.weight, id="mol.weight")
  ints.flt <- ints.long[ints.long$value > 0 & ! is.na(ints.long$value),]

  ints.flt$variable <- as.numeric(ints.flt$variable)
  y.lm <- stats::lm(data = ints.flt, formula = mol.weight ~ poly(variable, 3, raw = TRUE))

  y.lm
}

#' Filter_repeated_entries
#'
#' Remove proteins which were found in too many slices. Those are very usually
#' contaminants.
#'
#' @param pg MaxQuant ProteinGroup data.
#' @param rep_threshold Max percentage of appearance. Default is 0.3.
filter_repeated_entries <- function(pg, rep_threshold = 0.3){
  ints <- get_intensities(pg)
  nr_slices <- ncol(ints)
  ints_above_zero <- apply(ints, 1, function(x) length(which(x > 0)))
  not.contaminant <- ints_above_zero/nr_slices <= rep_threshold
  pg[not.contaminant,]
}


#' Filter_low_densities
#'
#' Remove entries in regions with low density.
#'
#' @param pg MaxQuant ProteinGroup data.
#' @param min_nr_threshold Min number of entries which should be found in a
#'   slice. Default is 10.
#' @param step_nr Number of steps for the theoretical mass. Default is 500.
filter_low_densities <- function(pg, min_nr_threshold = 10, step_nr = 500){

  # get the slices
  slices <- get_slice_numbers(pg)

  # get the mol_weight range
  mol_weight <- pg$Mol..weight..kDa.
  min_mass <- min(mol_weight)
  max_mass <- max(mol_weight)
  step_size <- (max_mass - min_mass) / step_nr

  # get the intensities
  ints <- get_intensities(pg)

  # loop in slices
  for(slice in slices){
    lower_mass <- min_mass

    for(i_mass in 1:step_nr){
      higher_mass <- lower_mass + step_size

      # select points
      selected_masses <- which(mol_weight >= lower_mass & mol_weight < higher_mass)
      selected_ints <- ints[selected_masses, slice]

      # if it doesnt pass the threshold, we remove the points
      if(sum(selected_ints > 0) < min_nr_threshold){
        ints[selected_masses, slice] <- NA
      }

      lower_mass <- higher_mass
    }

  }

  col_ids <- get_intensity_columns(pg)
  pg[, col_ids] <- ints
  pg

}


#' Filter weak intensitites
#'
#' Remove entries with a weak intensity.
#'
#' @param pg MaxQuant ProteinGroup data.
#' @param int_threshold_percent Max percentage of appearance. Default is 0.3.
filter_weak_intensity <- function(pg, int_threshold_percent = 0.5){
  # work with a logarithmic scale
  ints <- get_intensities(pg)
  ints[ints == 0] <- NA
  ints <- log(ints)

  # calculate the intensity threshold
  max_int <- max(ints, na.rm = TRUE)
  min_int <- min(ints, na.rm = TRUE)
  int_span <- max_int - min_int
  int_threshold <- min_int + int_span * int_threshold_percent

  # remove too weak points from pg
  too_weak <- ints < int_threshold
  col_ids <- get_intensity_columns(pg)
  pg[, col_ids][too_weak] <- 0
  pg
}


