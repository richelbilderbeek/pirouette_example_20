# From https://github.com/richelbilderbeek/pirouette_article/issues/55 :
#
# Write script that shows the true and twin error for a fixed tree prior and
# reasonably small trees (say in the magnitudes of hundreds of them)
# with 10:40 taxa
suppressMessages(library(pirouette))
library(testthat)

# Constants
is_testing <- is_on_travis()
example_no <- 20
n_replicates <- 5
n_taxa <- c(10, 20, 30, 40)
crown_age <- 10
is_testing <- is_on_travis()
# Number of replicates per number of taxa
if (is_testing) {
  n_replicates <- 2
  n_taxa <- c(3, 4, 5)
}
rng_seeds <- seq(314, 314 - 1 + length(n_taxa) * n_replicates)
n_taxas <- rep(n_taxa, each = n_replicates)
expect_equal(length(rng_seeds), length(n_taxas))

# Create phylogenies
phylogenies <- list()
for (i in seq_along(rng_seeds)) {
  rng_seed <- rng_seeds[i]
  n_taxa <- n_taxas[i]
  set.seed(rng_seed)
  phylogenies[[i]] <- create_yule_tree(
    n_taxa = n_taxa,
    crown_age = crown_age
  )
}
expect_equal(length(phylogenies), length(rng_seeds))
# Create pirouette parameter sets
pir_paramses <- create_std_pir_paramses(
  n = length(phylogenies)
)
expect_equal(length(pir_paramses), length(phylogenies))
if (is_testing) {
  pir_paramses <- shorten_pir_paramses(pir_paramses)
}
# Save tree to files
for (i in seq_along(pir_paramses)) {
  expect_equal(length(pir_paramses), length(phylogenies))
  rng_seed <- pir_paramses[[i]]$alignment_params$rng_seed
  folder_name <- file.path(paste0("example_", example_no, "_", rng_seed))
  # Create folder, do not warn if it already exists
  dir.create(folder_name, showWarnings = FALSE, recursive = TRUE)
  ape::write.tree(
    phylogenies[[i]],
    file = file.path(folder_name, "true_tree.newick")
  )
}
# Do the runs
pir_outs <- pir_runs(
  phylogenies = phylogenies,
  pir_paramses = pir_paramses
)

# Save
for (i in seq_along(pir_outs)) {
  expect_equal(length(pir_paramses), length(pir_outs))
  rng_seed <- pir_paramses[[i]]$alignment_params$rng_seed
  folder_name <- file.path(paste0("example_", example_no, "_", rng_seed))
  pir_save(
    phylogeny = phylogenies[[i]],
    pir_params = pir_paramses[[i]],
    pir_out = pir_outs[[i]],
    folder_name = folder_name
  )
}

