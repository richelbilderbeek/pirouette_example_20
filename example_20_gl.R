error_vs_ntaxa <- function(n_replicates = 2, seq_length = 1e3) {
  
  # check beast
  if (beastier::is_beast2_installed() == FALSE) {
    beastier::install_beast2() 
  }
  if (mauricer::is_beast2_ns_pkg_installed() == FALSE) {
    mauricer::install_beast2_pkg("NS")
  }
  stopifnot(beastier::is_beast2_installed())
  stopifnot(mauricer::is_beast2_ns_pkg_installed())

  suppressMessages(library(pirouette))
  suppressMessages(library(ggplot2))

  # folder setting
  root_folder <- getwd()
  example_no <- 20
  example_folder <- file.path(root_folder, paste0("example_", example_no, "_", "gl"))
  dir.create(example_folder, showWarnings = FALSE, recursive = TRUE)
  setwd(example_folder)
  
  # parsetting
  l_parses <- 4
  parses <- vector("list", l_parses)
  for (i in seq_len(l_parses)) {
    parses[[i]] <- data.frame(
      lambda = 0.8,
      mu = 0,
      crown_age = 10,
      cond = 1,
      n_0 = 2,
      n_taxa = 10 * i
    )
  }
  
  # simulate trees
  sim_data <- vector("list", length(parses))
  for (i in seq_len(l_parses)) {
    pars <- parses[[i]]
    for (seed in seq_len(n_replicates)) {
      set.seed(seed)
      sim_data[[i]][[seed]] <- TESS::tess.sim.taxa.age(
        n = 1,
        nTaxa = pars$n_taxa,
        age = pars$crown_age,
        lambda = pars$lambda,
        mu = pars$mu
      )[[1]]
      ape::is.ultrametric(sim_data[[i]][[seed]])
    }
  }

  # experiments
  if (rappdirs::app_dir()$os == "win") {
    experiment <- pirouette::create_gen_experiment()
    experiment$est_evidence_mcmc$chain_length <- 1e6
    experiment$inference_model$mcmc$chain_length <- 1e6
    experiments <- list(experiment)
  } else {
    experiments <- pirouette::create_all_experiments()
  }
  
  # create pir_params
  pir_paramseses <- vector("list", length(parses))
  for (i in seq_len(l_parses)) {
    for (seed in seq_len(n_replicates)) {
      phylogeny1 <- sim_data[[i]][[1]]
      pir_paramseses[[i]][[seed]] <- pirouette::create_test_pir_params(
        alignment_params = pirouette::create_alignment_params(
          sim_tral_fun = pirouette::get_sim_tral_with_std_nsm_fun(
            mutation_rate = pirouette::create_standard_mutation_rate(
              phylogeny = phylogeny1
            ),
            site_model = beautier::create_jc69_site_model()
          ),
          root_sequence = pirouette::create_blocked_dna(length = seq_length),
        ),
        twinning_params = pirouette::create_twinning_params(
          rng_seed_twin_tree = seed,
          rng_seed_twin_alignment = seed,
          sim_twin_tree_fun = pirouette::get_sim_bd_twin_tree_fun(),
          sim_twal_fun = pirouette::get_sim_twal_with_std_nsm_fun(
            mutation_rate = pirouette::create_standard_mutation_rate(
              phylogeny = phylogeny1
            ),
            site_model = beautier::create_jc69_site_model()
          )
        ), 
        experiments = experiments,
        error_measure_params = pirouette::create_error_measure_params(
          error_fun = pirouette::get_nltt_error_fun()
        )
      )
    }
  }
  
  # pir run!
  pir_outs <- vector("list", length(parses))
  for (i in seq_along(parses)){
    phylogenies <- sim_data[[i]]
    pir_paramses <- pir_paramseses[[i]]
    pir_outs[[i]] <- pirouette::pir_runs(
      phylogenies = phylogenies,
      pir_paramses = pir_paramses
    )
  }

  # return and save out
  save(
    pir_outs,
    file = file.path(
      example_folder,
      paste0("error_vs_ntaxa_", n_replicates, ".RData")
    )
  )
  pir_outs
}
error_vs_ntaxa(
  n_replicates = 2,
  seq_length = 1e2 * (rappdirs::app_dir()$os == "win") +
    1e3 * (rappdirs::app_dir()$os != "win")
)
