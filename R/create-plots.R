# Plot functions

#' plot_MQ
#'
#' plot the intensities
#'
#' @param pg ProteinGroups data.frame.
#' @examples
#' proteinGroups_path <- proteinGroups_path <- system.file("extdata", "Conde_9508_sub.txt", package = "pumbaR")
#' pg <- load_MQ(proteinGroups_path)
#' print(plot_MQ(pg))
plot_MQ <- function(pg){

  # get the intensities from the proteinGroups
  ints <- get_intensities(pg)
  ints.weight <- cbind(ints, mol.weight = log10(pg$Mol..weight..kDa.))
  ints.long <- reshape2::melt(ints.weight, id="mol.weight")
  ints.flt <- ints.long[ints.long$value > 0 & ! is.na(ints.long$value),]
  max.ints <- max(ints.flt$value)
  mn.ints <- min(ints.flt$value)

  # create a ggplot
  p <- ggplot2::ggplot(data=ints.flt, ggplot2::aes(x=variable, y=mol.weight, colour=value))
  p <- p  + ggplot2::geom_point(position="jitter", alpha=0.2)
  p <- p  + ggplot2::scale_colour_gradient2("Intensity (log)", trans="log", limits=c(min.ints, max.ints))
  p <- p  + ggplot2::xlab("slice number") + ggplot2::ylab("theoretical MW (log)")
  p <- p + ggplot2::theme_bw()

  p
}

#' plot_fit
#'
#' plot the MaxQuant data with the fitting curve
#'
#' @param pg ProteinGroups data.frame.
#' @param mass_fit Mass-slice fitting function
plot_fit <- function(pg, mass_fit){

  # get the intensities from the proteinGroups
  ints <- get_intensities(pg)
  ints.weight <- cbind(ints, mol.weight = log10(pg$Mol..weight..kDa.))
  ints.long <- reshape2::melt(ints.weight, id="mol.weight")
  ints.flt <- ints.long[ints.long$value > 0 & ! is.na(ints.long$value),]
  max.ints <- max(ints.flt$value)
  min.ints <- min(ints.flt$value)

  # the data to show the fitting curve
  plot_fit_function <- function(x){
    y <- predict(mass_fit, data.frame(variable=x))
    return(y[[1]])
  }

  # prepare the data for the fitting function
  nr_slices <- ncol(ints)
  plot_fit_data <- data.frame(slice=1:nr_slices, mass=unlist(lapply(1:nr_slices, plot_fit_function)))

  # create a ggplot
  p <- ggplot2::ggplot(data=ints.flt, ggplot2::aes(x=variable, y=mol.weight, colour=value))
  p <- p  + ggplot2::geom_point(position="jitter", alpha=0.2)
  p <- p  + ggplot2::scale_colour_gradient2("Intensity (log)", trans="log", limits=c(min.ints, max.ints))
  p <- p + ggplot2::scale_x_discrete(limits=c(1:nr_slices))
  p <- p  + ggplot2::xlab("slice number")
  p <- p  + ggplot2::ylab("theoretical MW (log)")
  p <- p + ggplot2::geom_line(data=plot_fit_data, ggplot2::aes(x=slice, y=mass), color="red", size=1.5)
  p <- p + ggplot2::theme_bw()

  p
}

