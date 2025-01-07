## Purpose: Run second-step models
## Project: FireBioOptima

run_2step <- function(model_path, ## where to find models (these models are too large for a standard github repo)
                      out_path, ## where to save models
                      taxa,
                      non_fdis = F) ## T when running models other than fdis-focused?
{
  library(tidyverse)
  library(brms)
  library(spOccupancy)

  
  ## read in msom models
  ml <- read_rds(model_path)
  
  ## pull out needed components
  d <- ml$d
  m <- ml$m
  
  ## set priors
  bprior <- prior(normal(0, 10), class = "b")
  
  ## If not plants do some more work
  if(taxa != "plants") {
    ## pull out covariate matrix to predict back to
    x.0 <- m$X
    
    ## posterior predictions back to observed sites
    post <- predict(m, x.0, ignore.RE = T)
    
    rich <- apply(post$psi.0.samples, c(1, 3), sum)
    
    rich.mn <- apply(rich, 2, mean)
    rich.sd <- apply(rich, 2, sd)
    rich.pi <- apply(rich, 2, quantile, probs = c(.05, .95))
    
    ## put it back together
    dat <- bind_cols(d$occ.covs, 
                     data.frame(rich_mn = rich.mn,
                                rich_sd = rich.sd,
                                rich_lower = rich.pi[1,],
                                rich_upper = rich.pi[2,])) %>% 
      ## lets standardize fdis, sev_mn, fri_mn, pat_mn, and elev; make new columns with 's' suffix
      mutate(fdis_s = scale(fdis),
             fdis_ls = scale(fdis_l),
             sev_mn_s = scale(sev_mn),
             fri_mn_s = scale(fri_mn),
             pat_mn_s = scale(pat_mn),
             elev_s = scale(elev),
             ## create a four class factor based on sev_mn_s breaks: breaks <- c(0, .1, 1.25, 2.25, 3)
             sev4 = cut(sev_mn, breaks = c(0, .1, 1.25, 2.25, 3), 
                        labels = c("1", "2", "3", "4"),
                        include.lowest = T)
             )
  }

    #### Bird models ####
  if(taxa == "birds") {
    
    bm.fdis.q <- brm(rich_mn | se(rich_sd) ~
                      fdis_s + I(fdis_s^2) +
                      elev_s + I(elev_s^2) +
                      (1 | tr + pt),
                    prior = bprior,
                    control = list(adapt_delta = 0.95),
                    chains = 3, data = dat)
    
    ## do we run the models below?
    if(non_fdis) {
      ## models below only need to be run for the 'best' MSOM predictions
      bm.sev.q <- brm(rich_mn | se(rich_sd) ~ 
                      sev_mn_s + I(sev_mn_s^2) +
                      fri_mn_s + I(fri_mn_s^2) +
                      elev_s + I(elev_s^2) +
                      (1 | tr + pt),
                    prior = bprior,
                    control = list(adapt_delta = 0.95),
                    chains = 3, data = dat)

      bm.fri.q <- brm(rich_mn | se(rich_sd) ~ 
                      fri_mn_s + I(fri_mn_s^2) +
                      elev_s + I(elev_s^2) +
                      (1 | tr + pt),
                    prior = bprior,
                    control = list(adapt_delta = 0.95),
                    chains = 3, data = dat)
      
      ## bundle models and save
      msom_models <- list(fdis.q = bm.fdis.q,
                          sev.q = bm.sev.q,
                          fri.q = bm.fri.q,
                          dat = dat)
      write_rds(msom_models, file = out_path, compress = "gz")
    } else {
      msom_models <- list(fdis = bm.fdis.q, dat = dat)
      write_rds(msom_models, file = out_path, compress = "gz")
    }
  }
  
  #### Bat models ####
  if(taxa == "bats") {
    bm.fdis.q <- brm(rich_mn | se(rich_sd) ~ 
                     fdis_s + I(fdis_s^2) +
                     elev_s + I(elev_s^2) +
                     (1 | pt),
                   prior = bprior,
                   control = list(adapt_delta = 0.95),
                   chains = 3, data = dat)
    
    ## do we run the models below?
    if(non_fdis) {

      bm.sev.q <- brm(rich_mn | se(rich_sd) ~ 
                        sev_mn_s + I(sev_mn_s^2) +
                        fri_mn_s + I(fri_mn_s^2) +
                        elev_s + I(elev_s^2) +
                        (1 | pt),
                      prior = bprior,
                      control = list(adapt_delta = 0.95),
                      chains = 3, data = dat)

      bm.fri.q <- brm(rich_mn | se(rich_sd) ~ 
                        fri_mn_s + I(fri_mn_s^2) +
                        elev_s + I(elev_s^2) +
                        (1 | pt),
                      prior = bprior,
                      control = list(adapt_delta = 0.95),
                      chains = 3, data = dat)
      
      ## bundle models and save
      msom_models <- list(fdis.q = bm.fdis.q,
                          sev.q = bm.sev.q,
                          fri.q = bm.fri.q,
                          dat = dat)
      write_rds(msom_models, file = out_path, compress = "gz")
    } else {
      msom_models <- list(fdis = bm.fdis.q,  
                          dat = dat)
      write_rds(msom_models, file = out_path, compress = "gz")
    }
  }
  
  

  
  #### Plant models ####

  if(taxa == "plants") {

    bm.fdis.q <- brm(n ~ 
                     fdis_s + I(fdis_s^2) +
                     elev_s + I(elev_s^2) +
                     (1 | soil),
                   prior = bprior,
                   family = 'poisson', 
                   control = list(adapt_delta = 0.95),
                   chains = 3, data = d)
    
    ## bundle models and save
    plant_models <- list(fdis = bm.fdis.q, 
                         dat = d)
    write_rds(plant_models, file = out_path, compress = "gz")

    ## do we run the models below?
    if(non_fdis) {

      bm.sev.q <- brm(n ~ 
                      sev_s + I(sev_s^2) +
                      fri_s + I(fri_s^2) +
                      elev_s + I(elev_s^2) +
                      (1 | soil),
                    prior = bprior,
                    family = 'poisson',
                    control = list(adapt_delta = 0.95),
                    chains = 3, data = d)

      bm.fri.q <- brm(n ~ 
                      fri_s + I(fri_s^2) +
                      elev_s + I(elev_s^2) +
                      (1 | soil),
                    prior = bprior,
                    family = 'poisson',
                    control = list(adapt_delta = 0.95),
                    chains = 3, data = d)
      
      ## bundle models and save
      plant_models <- list(fdis.q = bm.fdis.q,
                           sev.q = bm.sev.q,
                           fri.q = bm.fri.q,
                           dat = d)
      write_rds(plant_models, file = out_path, compress = "gz")

    }
    }
}

  
