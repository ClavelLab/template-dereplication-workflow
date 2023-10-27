# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c("maldipickr", "tidyverse", "coop", 
               "MALDIquant","readxl","writexl"),
  # packages that your targets need to run
  format = "rds" # default storage format
  )
options(clustermq.scheduler = "multicore")

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Parameters
which_raw_data_directory <- "/home/cpauvert/projects/iSOMiC/MALDI/dereplication/datasets/20230915_testRun_Sample_K0073/"

list(
  tar_target(
    raw_data_dir,
    which_raw_data_directory,
    format = "file"
  ),
  tar_target(
    spectra_raw,
    import_biotyper_spectra(raw_data_dir)
  ),
  tar_target(
    checks,
    check_spectra(spectra_raw, tolerance = 1)
  ),
  tar_target(#Stats: total spectra, total empty spectra
    spectra_stats,
    tibble(
      n_spectra = length(checks$is_empty),
      n_empty = sum(checks$is_empty)
    )
  ),
  tar_target( # Filter non empty spectra and unusual spectra
    nonempty_spectra,
    spectra_raw[!checks$is_empty]
  ),
  tar_target(
    processed,
    process_spectra(nonempty_spectra)
  ),
  tar_target(
    fm_interpolated,
    merge_processed_spectra(list(processed))
  ),
  tar_target(
    sim_interpolated,
    coop::tcosine(fm_interpolated)
  )
)