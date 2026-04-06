
#' Plot MCMC trace
#'
#' @param mdl_out Model output object name.
#' @param trace_obj Trace from the model output.
#' @param pop_info Population information. A tibble with columns
#'   collection (collection names), repunit (reporting unit names),
#'   and grpvec (group numbers).
#'
#' @return Trace plot in ggplot

#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # set up input data
#' gsi_data <- prep_gsi_data(mixture_data = mix, baseline_data = baseline, pop_info = pops211)
#'
#' # run model
#' gsi_out <- gsi_mdl(gsi_data, 10, 5, 1, 1)
#'
#' # trace plot
#' tr_plot(mdl_out = gsi_out, trace_obj = "trace", pop_info = gsi_out$groups)

tr_plot <- function (mdl_out, trace_obj, pop_info = NULL) {

  nburn <- as.numeric(mdl_out$specs["nburn"])
  keep_burn <- mdl_out$specs["keep_burn"] == "TRUE"
  thin <- as.numeric(mdl_out$specs["thin"])

  obj <- mdl_out[[trace_obj]] %>%
    dplyr::mutate(itr = (itr - nburn*isFALSE(keep_burn)) / thin)

  if (is.null(pop_info)) {
    name_order <- dplyr::select(obj, -c(itr, ch)) %>% colnames()
    trace <- tidyr::pivot_longer(obj, cols = -c(ch, itr),
                                 names_to = "repunit", values_to = "p")
  } else {
    name_order <- unique(pop_info$repunit)
    trace <- tidyr::pivot_longer(obj, cols = -c(ch, itr), names_to = "collection") %>%
      dplyr::left_join(pop_info, by = "collection") %>%
      dplyr::summarise(p = sum(value), .by = c(ch, itr, repunit))
  }

  trace %>%
    dplyr::mutate(repunit = factor(repunit, levels = name_order)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = itr, y = p, color = factor(ch))) +
    {if (keep_burn) {
      ggplot2::annotate("rect", fill = "red", alpha = 0.15,
                        xmin = 0, xmax = nburn/thin, ymin = -Inf, ymax = Inf)
    }} +
    ggplot2::facet_grid(repunit ~ ., scales = "free") +
    ggplot2::labs(color = "MC chain")

  # if (is.null(name_order)) {
  #   name_order <- dplyr::select(obj, -c(itr, ch)) %>% colnames()
  # }
  # tidyr::pivot_longer(obj, cols = -c(ch, itr)) %>%
  #   dplyr::mutate(name = factor(name, levels = name_order)) %>%
  #   ggplot2::ggplot() +
  #   ggplot2::geom_line(ggplot2::aes(x = itr, y = value, color = factor(ch))) +
  #   {if (nburn > 0) {
  #     ggplot2::annotate("rect", fill = "red", alpha = 0.15,
  #                       xmin = 0, xmax = nburn/thin, ymin = -Inf, ymax = Inf)
  #   }} +
  #   ggplot2::facet_grid(name ~ ., scales = "free") +
  #   ggplot2::labs(color = "MC chain")

}

utils::globalVariables(c("ch", "itr", "name", "value", "p"))




