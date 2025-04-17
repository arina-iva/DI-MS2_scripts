setwd("C:/Users/aivanova/Documents/Orbitrap Data/settings comparison") # <= To change!
source("utils.R")
library(mzR)
library(msdata)
library(dplyr)


#' Get scan number of MS2 (centroided) from the mzML data
#' 
#' @param The whole mzML data read by function OpenMSfile().
#' @returns A vector of scan number.
getMS2Scan <- function(msExp) {
  Filter(function(x) {header(msExp, x)$centroided==TRUE & header(msExp, x)$msLevel==2}, 1:runInfo(msExp)$scanCount)
}


#' Filter spectra with TIC (total ion current) and outputs a list of spectra.
#' The name of each element is scan number and the element is spectrum.
#' 
#' @param msExp The whole mzML data read by function OpenMSfile().
#' @param scan_num A vector of scan numbers, which are about to be filtered and output.
#' @param thres The intensity threshold (relative to TIC). 
#'    If thres is NULL, outputs the original spectra.
#' @param int_type Intensity type of output spectrum, "unnormalized" or "normalized".
#' @returns pl A list of filtered spectra.
getSpectrumLists <- function(msExp, scan_num, thres=NULL, int_type="unnormalized") {
  pl <- list()
  for (s in scan_num) {
    if (! is.null(thres)) {
      pl[[as.character(s)]] <- TICFilter(peaks(msExp, s), header(msExp, s)$totIonCurrent, thres, int_type=int_type)
    } else {
      pl[[as.character(s)]] <- peaks(msExp, s)
    }
  }
  pl
}


#' Filter spectrum using total ion current (TIC) according to the threshold
#' 
#' @param pl Spectrum (matrix).
#' @param TIC Total ion current.
#' @param thres Intensity threshold.
#' @param int_type Intensity type of output spectrum, "unnormalized" or "normalized".
#' @returns A filtered spectrum (matrix).
TICFilter <- function(pl, TIC, thres, int_type="unnormalized") {
  pl_df <- data.frame(
    mz <- pl[, 1],
    intensity <- pl[, 2],
    norm_intensity <- pl[, 2]/TIC
  )
  colnames(pl_df) <- c("mz", "intensity", "intensity_norm")
  if (int_type == "unnormalized") {
    as.matrix(dplyr::filter(pl_df, norm_intensity >= thres)[, c("mz", "intensity")])
  } else if (int_type == "normalized") {
    as.matrix(dplyr::filter(pl_df, norm_intensity >= thres)[, c("mz", "norm_intensity")])
  }
}



#' Average spectra. Each spectrum will be first filtered by TIC filter and merged.
#' 
#' @param filtered_spectra A list of spectra, which is output by function getSpectrumListss().
#' @param scan_vector A vector of scan number, which are about to be merged and averaged
#' @param dev Mass deviation for peak matching.
#' @returns o A merged spectrum with averaged relative intensity (relative to TIC).
averageSpectra <- function(filtered_spectra, scan_vector, dev) {
  pls <- c()
  for (i in scan_vector) {
    filtered_name <- names(filtered_spectra)
    if (as.character(i) %in% filtered_name) {
      pls <- rbind(pls, filtered_spectra[[as.character(i)]])
    }
  }
  pls <- pls[order(pls[, 1]), ]
  pls_mz <- cbind(pls, -1)
  pls_int <- cbind(pls, 1:nrow(pls))
  pls_int <- pls_int[order(pls_int[, 2], decreasing = TRUE), ]
  for (i in 1:nrow(pls_int)) {
    center <- pls_int[i, 3]
    if (pls_mz[center, 3] == -1) {
      pls_mz[center, 3] <- 1
      l <- center - 1
      while (l>=1) {
        if (abs(massDeviation(pls_mz[center, 1], pls_mz[l, 1])) <= dev & pls_mz[l, 3]==-1) {
          pls_mz[center, 2] <- pls_mz[center, 2] + pls_mz[l, 2]
          pls_mz[l, 2] <- 0
          pls_mz[l, 3] <- 0
        } else {
          break
        }
        l <- l - 1
      }
      r <- center + 1
      while (r<=nrow(pls_int)) {
        if (abs(massDeviation(pls_mz[center, 1], pls_mz[r, 1])) <= dev & pls_mz[r, 3]==-1) {
          pls_mz[center, 2] <- pls_mz[center, 2] + pls_mz[r, 2]
          pls_mz[r, 2] <- 0
          pls_mz[r, 3] <- 0
        } else {
          break
        }
        r <- r + 1
      }
    }
  }
  o <- data.frame(
    mz <- pls_mz[, 1],
    intensity <- pls_mz[, 2],
    status <- pls_mz[, 3]
  )
  o <- as.matrix(dplyr::filter(o, status == 1))[, 1:2]
  o[, 2] <- o[, 2]/length(scan_vector)
  colnames(o) <- c("mz", "intensity")
  o
}


#' Output an intensity matrix. The column name is the scan number. The row name 
#' is m/z value of peaks
#'
#' @param mzs A list of m/z values.
#' @param pls A list of spectra, which is output by function getSpectrumListss(). The
#'    name (character) of each element is scan number, and the element is spectrum.
#' @param dev Mass deviation.
#' @returns o An intensity matrix. 
getIntensityMatrix <- function(mzs, pls, dev) {
  o <- c()
  scan_num <- names(pls)
  for (p in mzs) {
    profile <- rep(0, length(scan_num))
    for (i in 1:length(scan_num)) {
      spectrum <- pls[[scan_num[i]]]
      index <- findPeak(p, spectrum, dev, type="maxInt")
      if (index!=-1) {
        profile[i] <- spectrum[index, 2]
      }
    }
    o <- rbind(o, profile)
  }
  colnames(o) <- as.character(scan_num)
  rownames(o) <- as.character(mzs)
  o
}


#' Remove blank peaks from spectrum. If the spectrum and blank spectrum is empty,
#' return the empty or the original spectrum, respectively.
#' 
#' @param spec A spectrum (matrix).
#' @param blank Blank spectrum (matrix).
#' @param dev Mass deviation.
#' @returns A spectrum without blank peaks.

removeBlank_2 <- function(spec, blank, dev) {
  if (nrow(spec) > 0 & nrow(blank) > 0) {
    # sapply returns a list of row indices to delete, or NULL if no match
    toDelete <- sapply(1:nrow(spec), function(p) {
      if (findPeak(spec[p, "mz"], blank, dev, type = "maxInt") != -1) { 
        return(p)
      }
      return(NULL)
    })
    toDelete <- unlist(toDelete)
# Only attempt to delete rows if toDelete is non-empty
    if (length(toDelete) > 0) {
      spec <- spec[-toDelete, ]
    }
  }
  return(spec)
}


#' Get pearson correlation matrix
#' 
#' @param int_matrix A intensity matrix.
#' @returns A correlation matrix.
getCorrMatrix <- function(int_matrix) {
  cor(t(int_matrix), method="pearson")
}