
# template-dereplication-workflow

<!-- badges: start -->
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

`template-dereplication-workflow` provides a reusable template to run a reproducible workflow (using [`{targets}`](https://docs.ropensci.org/targets/)) to dereplicate and cherry-pick mass spectrometry spectra using the [`{maldipickr}`](https://cran.r-project.org/package=maldipickr), in order to reduce the redundancy of bacterial isolates.
 
## Usage

1. Clone the repository
2. Set up the dependencies using `renv::restore()`
3. Fill out the needed parameters/settings in [`setup-workflow.R`](setup-workflow.R) and run the R code to create `_targets.R`
4. Evaluate the workflow using `targets::tar_manifest()` or using `targets::tar_visnetwork()`
5. Run the workflow with `targets::tar_make()` 
