gather_spectra_stats <- function(check_vectors, aggregated_checks){
  # check_vectors from maldipickr::check_spectra
  # aggregated_checks from Reduce(`|`, check_vectors)
  check_stats <- vapply(check_vectors, sum, FUN.VALUE = integer(1)) %>%
    tibble::as_tibble_row()
  tibble::tibble(
    "n_spectra" = length(aggregated_checks),
    "n_valid_spectra" = n_spectra - sum(aggregated_checks)
  ) %>%
    dplyr::bind_cols(check_stats) %>% 
    return()
}
