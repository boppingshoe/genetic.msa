
#' Sample total harvest numbers from a lognormal distribution
#'
#' @param x A vector consists of harvest mean and CV.
#' @param n Number of sample draws.
#' @param seed Optional random seed.
#'
#' @return A vector of harvest numbers drawn from a lognormal distribution with specified mean and CV.
#'
#' @examples
#' tot_harv <- harv_func(c(500, 0.05))
#'
#' @export
harv_func <- function(x, n = 5000, seed = NULL) {
  lnvar <- log(x[2]^2 + 1)
  lnmean <- log(x[1]) - lnvar / 2
  if (!is.null(seed)) set.seed(seed)
  stats::rlnorm(n, lnmean, sqrt(lnvar))
}


#' Summarize trace output for stock-specific harvest
#'
#' @param mdl_out Output of GSI model run
#' @param mdl_dat Input data for GSI model
#'
#' @return A tibble of harvest estimates for each reporting group as columns and MCMC iterations as rows.
#'
#' @examples
#' # prep input data
#' gsi_data <- prep_gsi_data(mixture_data = mix, baseline_data = baseline, pop_info = pops211)
#'
#' # run model
#' tot_harv <- harv_func(c(500, 0.05))
#' gsi_out <- gsi_mdl(gsi_data, 10, 5, 1, 1, harvest = tot_harv)
#'
#' # summarize individual assignments
#' harv_trace <- harv_summ(gsi_out, gsi_data)
#'
#' @export
harv_summ <- function(mdl_out, mdl_dat) {
  rowsum(t(mdl_out$sstc_trace), mdl_dat$group_names[mdl_dat$groups]) %>%
    t() %>%
    tibble::as_tibble()
}


