#' Get water content based on pressure head
#'
#' @param H pressure head [cm].
#' @param WCR residual water content [cm3 cm-3].
#' @param WCS saturated water content [cm3 cm-3].
#' @param ALPHA curve shape parameter [-].
#' @param NPAR curve shape parameter [-].
#' @return return water content [cm3 cm-3].
#' @keywords internal
get_wc_MvG  <- function(H, WCR, WCS, ALPHA, NPAR) {

  # --- main part of procedure ---

  m <- 1 - 1 / NPAR
  wc <- WCR + (WCS - WCR) / ((1 + (ALPHA * H)^NPAR)^m)

  # --- return of procedure ---

  return (wc)
}

#' Get pressure head based on water content
#'
#' @param WC water content [cm3 cm-3].
#' @param WCR residual water content [cm3 cm-3].
#' @param WCS saturated water content [cm3 cm-3].
#' @param ALPHA curve shape parameter [-].
#' @param NPAR curve shape parameter [-].
#' @return return pressure head [cm].
#' @keywords internal
get_h_MvG  <- function(WC, WCR, WCS, ALPHA, NPAR) {

  # --- main part of procedure ---

  m <- 1 - 1 / NPAR
  h <- (((WCS - WCR) / (WC - WCR))^(1/m) - 1)^(1 / NPAR) / ALPHA

  # --- return of procedure ---

  return (h)
}

#' Get conductivity based on pressure head
#'
#' @param H pressure head [cm]
#' @param ALPHA curve shape parameter [-]
#' @param NPAR curve shape parameter [-]
#' @param LAMBDA exponent in hydraulic conductivity function [-]
#' @param KSAT hydraulic conductivity of saturated soil [cm d-1]
#' @return return conductivity [cm d-1]
#' @keywords internal
get_cond_MvG  <- function(H, ALPHA, NPAR, LAMBDA, KSAT) {

  # --- main part of procedure ---

  m <- 1 - 1 / NPAR
  ah <- ALPHA * H
  h1 <- (1 + ah^NPAR)^m
  h2 <- ah**(NPAR - 1)
  denom <- (1 + ah^NPAR)^(m*(LAMBDA + 2))
  cond <- KSAT * (h1 - h2)^2 / denom

  # --- return of procedure ---

  return (cond)
}

#' Create MvG-table
#'
#' @param file filename.
#' @param WCR residual water content [cm3 cm-3]
#' @param WCS saturated water content [cm3 cm-3]
#' @param ALPHA curve shape parameter [-]
#' @param NPAR curve shape parameter [-]
#' @param LAMBDA curve shape parameter [-]
#' @param KSAT hydraulic conductivity of saturated soil [cm d-1]
#' @importFrom fs path_file
#' @importFrom stringr str_c
#' @export set_table_MvG
set_table_MvG <- function(file, WCR, WCS, ALPHA, NPAR, LAMBDA, KSAT) {

  # --- main part of procedure ---

  # set sequeunce pressure head
  pF <- c(seq(from = 7, to = 4.3, by = -0.1), seq(from = 4.2, to = 0.01, by = -0.01))
  H <- c(10^pF, seq(from = 1, to = 0, by = -0.1))

  # get volumetric moisture content
  WC <- get_wc_MvG(H = H, WCR = WCR, WCS = WCS, ALPHA = ALPHA, NPAR = NPAR)

  # get hydraulic conductivity
  COND <- get_cond_MvG(H = H, ALPHA = ALPHA, NPAR = NPAR, LAMBDA = LAMBDA, KSAT = KSAT)

  # set header
  header <- c("************************************************************************",
    str_c("* Filename: ", path_file(file)),
    "* Contents: SWAP 4 - soil physical properties",
    "************************************************************************",
    "*",
    "* thetatab   volumetric moisture content [cm3/cm-3]",
    "* headtab    pressure head corresponding to theta [cm]",
    "* conductab  hydraulic conductivity corresponding to theta [cm d-1]",
    "*",
    "      thetatab         headtab      conductab")

  # set data
  data <- str_c(
    formatC(x = WC, format = "e", digits = 10, width = 15),
    formatC(x = -H, format = "e", digits = 10, width = 15),
    formatC(x = COND, format = "e", digits = 10, width = 15),
    sep = " ")

  # write table
  writeLines(text = c(header, data), con = file)
}

#' Get water content of soil profile at certain depth based on specific pressure head
#'
#' @param file_swp character string, name of swp-file.
#' @param depth character string, range of depth in [cm].
#' @param ... further arguments passed to or from other methods
#' @importFrom stringr str_c str_split
#' @importFrom tibble tibble
#' @importFrom dplyr mutate select lag
#' @importFrom rlang .data
#' @return return water content [cm3 cm-3].
#' @details
#' Pressure head can be specified by \code{h} or \code{pF}.
#' @export get_wc_h
get_wc_h <- function(file_swp, depth, ...) {

  # --- initial part of procedure ---

  # set optional arguments
  opt_param <- c("h", "pF")
  h <- pF <- NULL

  # load additional arguments
  param <- list(...)
  for (name in names(param) ) assign(name, param[[name]])

  # check unused arguments
  unused_param <- setdiff(names(param), opt_param)
  if (length(unused_param) > 0) stop(str_c("unused parameter: ", str_c(unused_param, collapse = ', ')))

  # check used arguments
  if (is.null(h) & is.null(pF)) stop("either 'h' or 'pF' should be specified")

  # set variable to extract
  variable <- c("ISOILLAY", "HSUBLAY", "ORES", "OSAT", "ALFA", "NPAR")


  # --- main part of procedure ---

  # set pressure head (in case of pF)
  if (is.null(h)) h <- -10^pF

  # set minimum and maximum depth
  depth_range <- str_split(string = depth, pattern = ":", simplify = TRUE)
  depth_min <- -as.numeric(depth_range[,1])
  depth_max <- -as.numeric(depth_range[,2])

  # extract settings SWAP
  settings <- get_records_SWAP(file = file_swp, variable = variable)

  # determine water content at given pressure head
  db_prof <- tibble(ISOILLAY = settings$ISOILLAY, HSUBLAY = settings$HSUBLAY) |>
    mutate(
      dmax = cumsum(.data$HSUBLAY),
      dmin = lag(.data$dmax, n = 1, default = 0),
      ORES = settings$ORES[.data$ISOILLAY],
      OSAT = settings$OSAT[.data$ISOILLAY],
      ALFA = settings$ALFA[.data$ISOILLAY],
      NPAR = settings$NPAR[.data$ISOILLAY],
      WC = get_wc_MvG(H = -.data$h, WCR = .data$ORES, WCS = .data$OSAT, ALPHA = .data$ALFA, NPAR = .data$NPAR)
    ) |>
    select("dmin", "dmax", "WC")

  # determine total water content over given depth
  wc <- NULL
  for (i_depth in seq_along(depth_min)) {
    wc_tmp <- 0.0
    for (i_rec in seq_along(db_prof$WC)) {
      if ((depth_min[i_depth] >= db_prof$dmin[i_rec] & depth_min[i_depth] < db_prof$dmax[i_rec]) | (depth_max[i_depth] >= db_prof$dmin[i_rec] & depth_max[i_depth] < db_prof$dmax[i_rec])) {
        wc_tmp <- wc_tmp + db_prof$WC[i_rec] * (min(depth_max[i_depth], db_prof$dmax[i_rec]) - max(depth_min[i_depth], db_prof$dmin[i_rec]))
      }
    }
    wc <- c(wc, wc_tmp)
  }

  # --- return of procedure ---

  return(wc)
}



