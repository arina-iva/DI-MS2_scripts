#' Calculate the mass deviation between two values
#' 
#' @param refMz m/z value
#' @param measureMz m/z value
#' @returns mass deviation
massDeviation <- function(refMz, measureMz) {
  if (refMz>=200) (refMz-measureMz)/refMz*(10**6) else (refMz-measureMz)/200*(10**6)
}


#' Output matched position(s) of a certain m/z value in a peak list
#' 
#' @param mz A m/z value
#' @param pl A peak list containing m/z and intensity
#' @param dev Mass deviation
#' @param type Peak matching methods (3 options are available)
#'             "nearest": Find the closed matched peak (not recommended)
#'             "maxInt": Find a matched peak with the highest intensity
#'             "sum": Find all matched peaks.
#' @returns Position(s) of matched peak(s)
findPeak <- function(mz, pl, dev, type) {
  match <- -1
  all <- sapply(pl[, 1], function(x) {abs(massDeviation(mz, x))})
  matches <- which(all<=dev)
  if (length(matches)==1 || (type=="sum" & length(matches)>1)) {
    match <- matches
  } else if (length(matches)>1) {
    if (type=="maxInt") {
      match <- matches[which.max(pl[matches, 2])]
    } else if (type=="nearest") {
      match <- matches[which.min(matches)]
    }
  }
  match
}

