checkPackages <- function(...) {
  pkgs <- list(...)
  names(pkgs) <- pkgs
  pkgs <- vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)
  if (!all(pkgs)) {
    pkgs <- names(pkgs)[!pkgs]
    pkgs <- paste(paste0("'", pkgs, "'"), collapse = ", ")
    stop(sprintf(
      "The following packages need to be installed for this feature to work: %s",
      pkgs
    ))
  }
}

int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps^0.5]
}

# a handy function opposite to %in%
`%ni` <- Negate(`%in%`)

# Dharmesh's favorite theme
bhuvad_theme <- function(textScale = 1.1) {
  stopifnot(textScale > 0)
  ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.border = element_rect(colour = "black", fill = NA),
      panel.grid = element_blank(),
      axis.title = element_text(size = rel(textScale) * 1.1),
      axis.text = element_text(size = rel(textScale)),
      plot.title = element_text(size = rel(textScale) * 1.2),
      strip.background = element_rect(fill = NA, colour = "black"),
      strip.text = element_text(size = rel(textScale)),
      legend.text = element_text(size = rel(textScale)),
      legend.title = element_text(size = rel(textScale), face = "italic")
    )
}

orderSamples <- function(sdata, ordannots) {
  # add sample IDs
  sdata$SampleOrderID <- seq(nrow(sdata))
  
  # order samples based on provided annotations
  if (!is.null(ordannots)) {
    if (length(ordannots) == 1) {
      sdata <- sdata[order(sdata[[ordannots]]), ]
    } else {
      sdata <- sdata[do.call(order, sdata[ordannots]), ]
    }
  }
  
  return(sdata$SampleOrderID)
}
