checkPackages <- function(...) {
  pkgs = list(...)
  names(pkgs) = pkgs
  pkgs = sapply(pkgs, requireNamespace, quietly = TRUE)
  if (!all(pkgs)) {
    pkgs = names(pkgs)[!pkgs]
    pkgs = paste(paste0("'", pkgs, "'"), collapse = ', ')
    stop(sprintf(
      'The following packages need to be installed for this feature to work: %s',
      pkgs
    ))
  }
}

int_breaks <- function(x, n = 5) {
  l <- pretty(x, n)
  l[abs(l %% 1) < .Machine$double.eps ^ 0.5]
}
