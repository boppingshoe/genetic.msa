

#' Calculate error and accuracy rates for plotting ROC
#'
#' @param scaled_like Output from self_assign_le.
#' @param threshold Cut-off for accepting individual assignment.
#' @param unit Which reporting unit to calculate.
#'
#' @return An array with error/accuracy rates for plotting ROC.
#'
#' @examples
#' \dontrun{
#' rates <- rate_calc(scaled_like, threshold = 0.5, unit = "Russia")
#' }
#'
#' @export
#'
rate_calc <- function(scaled_like, threshold, unit) {

  scaled_like %>%
    tidyr::pivot_longer(-c(indiv, collection, repunit),
                        names_to = "inferred_repunit",
                        values_to = "scaled_likelihood") %>%
    dplyr::mutate(
      assign = dplyr::case_when(repunit == inferred_repunit ~ "positive",
                                TRUE ~ "negative")
      ) %>%
    dplyr::group_by(indiv, inferred_repunit, assign) %>%
    dplyr::summarise(prob = sum(scaled_likelihood),
                     .groups = "drop") %>%
    dplyr::filter(inferred_repunit == unit) %>%
    dplyr::mutate(
      outcome = dplyr::case_when(
        assign == "positive" &
          prob >= threshold ~ "tp", # true positive
        assign == "positive" &
          prob < threshold ~ "fn", # false negative
        assign == "negative" &
          prob >= threshold ~ "fp", # false positive
        TRUE ~ "tn" # true negative
      )
    ) %>%
    dplyr::count(outcome) %>%
    tidyr::pivot_wider(names_from = outcome, values_from = n) %>%
    dplyr::mutate(tpr = tp / (tp + fn), # true positive rate (power)
                  fpr = fp / (tn + fp), # false positive rate (type 1 error)
                  acc = (tp + tn) / (tp + tn + fp + fn), # accuracy rate
                  level = threshold) # threshold level
}

#' Self-assignment using baseline data
#'
#' @param base_dat Baseline data set.
#' @param pop_info Population information for the baseline.
#'
#' @return An array of assignment probabilities for each reporting group as columns and individual fish in the baseline as rows. Each individual has a known identity.
#'
#' @examples
#' \dontrun{
#' scaled_like <- self_assign_le(baseline, pops211)
#' }
#'
#' @export
#'
self_assign_le <- function(base_dat, pop_info) {
  message("A less efficient version of rubias::self_assign()...")

  if ("SILLY_CODE" %in% names(base_dat)) base_dat <- dplyr::rename(base_dat, collection = SILLY_CODE, indiv = SillySource)

  loci <-
    dplyr::tibble(locus = names(base_dat)) %>%
    dplyr::filter(grepl("\\.1$", locus)) %>%
    dplyr::mutate(locus = substr(locus, 1, nchar(locus) - 2)) %>%
    dplyr::pull(locus)

  base <- allefreq(base_dat, base_dat, loci, collect_by = collection) %>%
    dplyr::right_join(pop_info, by = "collection", keep = FALSE) %>%
    dplyr::relocate(!dplyr::ends_with(as.character(0:9)), .after = collection) %>%
    dplyr::mutate(dplyr::across(dplyr::ends_with(as.character(0:9)), ~tidyr::replace_na(., 0)))

  mix <- allefreq(base_dat, base_dat, loci) %>%
    dplyr::right_join(dplyr::select(base_dat, collection, indiv), by = "indiv", keep = FALSE) %>%
    dplyr::relocate(dplyr::ends_with(as.character(0:9)), .after = collection)

  x <- mix %>%
    dplyr::select(dplyr::ends_with(as.character(0:9))) %>%
    dplyr::select(order(colnames(.))) %>%
    as.matrix(.)

  nalleles_tibble <- lapply(loci, function(loc) {
    dplyr::tibble(locus = loc,
                  call = base_dat %>%
                    dplyr::select(dplyr::all_of(loc), paste0(loc, ".1")) %>%
                    unlist() %>% unique() %>% .[!is.na(.)],
                  altyp = seq.int(dplyr::n_distinct(call)) %>% factor())
  }) %>% dplyr::bind_rows() %>%
    dplyr::group_by(locus) %>%
    dplyr::summarise(n_allele = max(as.numeric(altyp)), .groups = "drop")

  nalleles <- nalleles_tibble %>%
    dplyr::pull(n_allele) %>%
    stats::setNames(nalleles_tibble$locus)

  trait_fac <- factor(rep(names(nalleles), nalleles), levels = names(nalleles))

  beta <- matrix(
    rep(1 / nalleles, nalleles),
    nrow = nrow(base),
    ncol = ncol(base) - 3,
    byrow = TRUE
  )

  grp_nms <- base %>%
    dplyr::arrange(grpvec) %>%
    dplyr::pull(repunit) %>%
    unique()

  y <- base %>%
    dplyr::select(dplyr::ends_with(as.character(0:9))) %>%
    dplyr::select(order(colnames(.))) %>%
    as.matrix(.)

  t_q <- apply(y + beta, 1, function(rw) {
    unlist(tapply(rw, trait_fac, function(betty) {
      if (sum(betty)) {
        betty / sum(betty)
      } else {
        rep(1, length(betty))
      }
    }, simplify = FALSE)[names(nalleles)])
  }) # transposed (allele freq)

  # loo
  cl <- parallel::makePSOCKcluster(5)
  doParallel::registerDoParallel(cl, cores = 5)

  scaled_like <- foreach::foreach(
    m = 1:nrow(mix),
    .packages = c("magrittr", "tidyr", "dplyr", "purrr")
  ) %dorng% {
    t_q_temp <- t_q
    grp <- which(base$collection == mix$collection[m])
    loo <- y[grp, ] - x[m, ] + beta[grp, ]
    t_q_temp[, grp] <-
      unlist(tapply(loo, trait_fac, function(betty) {
        if (sum(betty)) {
          betty / sum(betty)
        } else {
          rep(1, length(betty))
        }
      }, simplify = FALSE)[names(nalleles)])

    freq <- {exp(x[m, ] %*% log(t_q_temp))} %>%
      prop.table(.)

    tapply(freq, base$grpvec, sum) %>%
      purrr::set_names(grp_nms)

  } %>% dplyr::bind_rows() %>%
    dplyr::mutate(indiv = mix$indiv,
                  collection = mix$collection) %>%
    dplyr::left_join(dplyr::select(base, collection, repunit),
                     by ="collection", keep = FALSE)

  parallel::stopCluster(cl)

  return(scaled_like)
}


#' Calculate individual assignment using likelihood
#'
#' @param gsi_dat Input data for GSI model
#'
#' @return An array of assignment probabilities for each reporting group as columns and individual fish in the mixture as rows.
#'
#' @examples
#' # prep input data
#' gsi_data <- prep_gsi_data(mixture_data = mix, baseline_data = baseline, pop_info = pops211)
#'
#' assign_like <- indiv_assign(gsi_data)
#'
#' @export
#'
indiv_assign <- function(gsi_dat) {
  message("Calculating likelihood...")

  x <- gsi_dat$x %>%
    dplyr::select(dplyr::ends_with(as.character(0:9))) %>%
    dplyr::select(order(colnames(.))) %>%
    as.matrix(.)

  y <- gsi_dat$y %>%
    dplyr::select(dplyr::ends_with(as.character(0:9))) %>%
    dplyr::select(order(colnames(.))) %>%
    as.matrix(.)

  nalleles <- gsi_dat$nalleles

  trait_fac <- factor(rep(names(nalleles), nalleles), levels = names(nalleles))

  beta <- matrix(
    rep(1 / nalleles, nalleles),
    nrow = nrow(y),
    ncol = ncol(y),
    byrow = TRUE
  )

  grp_nms <- gsi_dat$y %>%
    dplyr::arrange(grpvec) %>%
    dplyr::pull(repunit) %>%
    unique()

  t_q <- apply(y + beta, 1, function(rw) {
    unlist(tapply(rw, trait_fac, function(betty) {
      if (sum(betty)) {
        betty / sum(betty)
      } else {
        rep(1, length(betty))
      }
    }, simplify = FALSE)[names(nalleles)])
  }) # transposed (allele freq)

  freq <- exp(x %*% log(t_q))

  scaled_like <-
    apply(freq, 1, prop.table) %>%
    rowsum(., gsi_dat$y$grpvec) %>%
    t(.) %>%
    dplyr::as_tibble(., .name_repair = "minimal") %>%
    purrr::set_names(grp_nms) %>%
    dplyr::mutate(indiv = gsi_dat$x$indiv)

  return(scaled_like)
}


#' Summarize individual assignment history output
#'
#' @param mdl_out Output of GSI model run
#' @param mdl_dat Input data for GSI model
#'
#' @return A tibble of assignment probabilities for each reporting group as columns and individual fish in the mixture as rows.
#'
#' @examples
#' # prep input data
#' gsi_data <- prep_gsi_data(mixture_data = mix, baseline_data = baseline, pop_info = pops211)
#'
#' # run model
#' gsi_out <- gsi_mdl(gsi_data, 10, 5, 1, 1)
#'
#' # summarize individual assignments
#' indiv_ass <- indiv_summ(gsi_out, gsi_data)
#'
#' @export
indiv_summ <- function(mdl_out, mdl_dat) {

  pi <- apply(mdl_out$idens, 2,
              function (idens) {
                factor(idens, levels = seq(length(mdl_dat$groups))) %>%
                  table(.) %>%
                  tapply(., mdl_dat$groups, sum) %>%
                  prop.table(.)
              })

  tidyr::tibble(ID = mdl_dat$x$indiv) %>%
  # tidyr::tibble(ID = rownames(mdl_dat$x)) %>% # when mixture has no indiv id
    dplyr::bind_cols({
      t(pi) %>%
        as.data.frame() %>%
        stats::setNames(mdl_dat$group_names)
    })
}









