# Template for MALDI Biotyper dereplication workflow
# Charlie Pauvert
# Created: 2023-10-27

# Dereplication workflow parameters
which_plate_metadata <- here::here("data", "20230920_testRun_Sample_K0073","Report_Step3a_scdPlates_PatientID_K0073_KoelnFMT_2023.09.20_08.52.01.xlsx")
which_threshold <- 0.92

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)

# Set target options:
tar_option_set(
  packages = c("maldipickr", "tidyverse", "coop", 
               "MALDIquant","readxl","writexl"),
  # packages that your targets need to run
  format = "rds", # default storage format,
  iteration = "list"
)
options(clustermq.scheduler = "multicore")

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Workflow
list(
  tarchetypes::tar_files(
    plates,
    list.dirs(
      here::here("data","20230920_testRun_Sample_K0073/"),
      recursive = F)
  ),
  tar_target(
    spectra_raw,
    import_biotyper_spectra(plates) %>% suppressWarnings(),
    pattern = map(plates)
  ),
  tar_target(
    checks,
    check_spectra(spectra_raw, tolerance = 1),
    pattern = map(spectra_raw)
  ),
  tar_target(
    problematic_spectra,
    # Logical 'OR' combinations of checks vectors
    # src: https://stackoverflow.com/a/51140480/21085566
    Reduce(`|`, checks),
    pattern = map(checks)
  ),
  tar_target(
    spectra_stats,
    gather_spectra_stats(checks, problematic_spectra) %>%
      dplyr::mutate(maldi_plate = plates),
    pattern = map(checks, problematic_spectra, plates),
    iteration = "vector"
  ),
  tar_target( # Filter-out non empty spectra and unusual spectra
    valid_spectra,
    spectra_raw[!problematic_spectra],
    pattern = map(spectra_raw, problematic_spectra)
  ),
  tar_target(
    all_stats,
    dplyr::bind_rows(spectra_stats)
  ),
  tar_target(
    processed,
    process_spectra(valid_spectra),
    pattern = map(valid_spectra)
  ),
  tar_target(
    slow_processed, unname(c(processed))
  ),
  tar_target(
    fast_processed_path, "fast_growers/objects/fast_processed", format = "file"
  ),
  tar_target(
    fast_processed, readRDS(fast_processed_path)
  ),
  tar_target(
    fm_interpolated,
    # Named lists are problematic for dynamic branching
    # as the name are appended to the matrix rownames
    merge_processed_spectra(c(fast_processed, slow_processed))
  ),
  tar_target(
    sim_interpolated,
    coop::tcosine(fm_interpolated)
  ),
  tar_target(
    excel_metadata,
    which_plate_metadata,
    format = "file"
  ),
  tar_target(# Get metadata from excel sheet)
    metadata,
    read_excel(excel_metadata)
  ),
  tar_target(
    metadata_picking,
    metadata %>% rename(
      c("OD600" = "Well OD600_BlankCorrected_MALDI_Step2_7Tag",
        "name" = "Well SampleName",
        "picked_before" = "Well Selected_MALDI_hits")) %>%
      dplyr::mutate(
        OD600 = as.double(OD600),
        picked_before = as.logical(strtoi(picked_before)),
        well = gsub(".*_([0-9]{1,3}$)", "\\1", name) %>%
          strtoi(),
        is_edge = maldipickr::is_well_on_edge(
          well_number = well, plate_layout = 384
        )) %>%
      select(name, OD600, is_edge, picked_before)
  ),
  tar_target(
    df_interpolated,
    delineate_with_similarity(sim_interpolated, threshold = which_threshold, method = "complete")
  ),
  tar_target(
    processed_metadata,
    dplyr::bind_rows(
      lapply(c(fast_processed,slow_processed), `[[`, "metadata")
    ),
    iteration = "list"
  ),
  tar_target(
    clusters,
    set_reference_spectra(df_interpolated, processed_metadata)
  ),
  tar_target(#clean up spectra names for cluster name
    # remove trailing _B11
    clusters_clean,
    clusters %>% dplyr::mutate(
      name = gsub("_[A-Z][0-9]{1,3}$","",name)
    )
  ),
  tar_target(
    # subset metadata information
    metadata_subset,
    metadata_picking %>% semi_join(clusters_clean, by = "name")
  ),
  tar_target(
    picked,
    pick_spectra(
      cluster_df = clusters_clean,
      metadata_df = metadata_subset,
      criteria_column = "OD600",
      soft_mask_column = "is_edge",
      hard_mask_column = "picked_before")),
  tar_target(
    summary_picked,
    picked %>% filter(to_pick) %>%
      transmute(
        name = name,
        cluster_size = cluster_size,
        picked_before = picked_before,
        procedure = paste("Strejcek", which_threshold, sep = "_")
      )
  ),
  tar_target(# prep excel sheet
    prep_excel,
    metadata %>% select(-c("Well Selected_MALDI_hits")) %>% 
      left_join(
        picked %>% transmute(
          `Well SampleName` = name,
          decision = if_else(to_pick, 2, 0),
          `Well Selected_MALDI_hits` = as.integer(picked_before) + decision
        ), by = "Well SampleName"
      ) %>%
      mutate(
        `Well Selected_MALDI_hits` = replace_na(`Well Selected_MALDI_hits`, 0)
      ) %>% select(-c(decision))
  ),
  tar_target(
    excel_output,
    write_xlsx(prep_excel,
               path =
                 here::here(
                   paste0(
                     "picked_",paste("Strejcek", which_threshold, sep = "_"),
                     "_",basename(which_plate_metadata)
                   )
                 )
    )
  ),
  tar_quarto(report, "report_slow.qmd")
)