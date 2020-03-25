# From https://github.com/richelbilderbeek/pirouette_article/issues/55 :
#
# Write script that shows the true and twin error for a fixed tree prior and
# reasonably small trees (say in the magnitudes of hundreds of them)
# with 10:40 taxa
suppressMessages(library(pirouette))
library(ggplot2)
library(testthat)

# Constants
is_testing <- is_on_ci()
example_no <- 20
n_replicates <- 5
n_taxa <- c(5, 10, 20, 30)
crown_age <- 10
folder_name <- paste0("example_", example_no)

# Number of replicates per number of taxa
if (is_testing) {
  n_replicates <- 2
  n_taxa <- c(3, 4)
}
rng_seeds <- seq(314, 314 - 1 + length(n_taxa) * n_replicates)
n_taxas <- rep(n_taxa, each = n_replicates)
expect_equal(length(rng_seeds), length(n_taxas))

# Create phylogenies
phylogenies <- list()
for (i in seq_along(rng_seeds)) {
  set.seed(rng_seeds[i])
  phylogenies[[i]] <- create_yule_tree(
    n_taxa = n_taxas[i],
    crown_age = crown_age
  )
}
expect_equal(length(phylogenies), length(rng_seeds))

# Create pirouette parameter sets
pir_paramses <- create_std_pir_paramses(
  n = length(phylogenies),
  folder_name = folder_name
)
expect_equal(length(pir_paramses), length(phylogenies))
if (is_testing) {
  pir_paramses <- shorten_pir_paramses(pir_paramses)
}

# Do the runs
pir_outs <- pir_runs(
  phylogenies = phylogenies,
  pir_paramses = pir_paramses
)

# Save summary
pir_plots(pir_outs) +
  ggtitle(paste("Number of replicates: ", n_replicates)) +
  ggsave(file.path(folder_name, "errors.png"), width = 7, height = 7)

# Save
for (i in seq_along(n_taxa)) {
  n <- n_taxa[i]
  from_index <- ((i - 1) * n_replicates) + 1
  to_index <- ((i - 1) * n_replicates) + n_replicates
  pir_plots(
    pir_outs = pir_outs[from_index:to_index]
  ) + ggtitle(paste("Number of taxa:", n)) +
    ggsave(filename = paste0("errors_", n, ".png"), width = 7, height = 7)
}
