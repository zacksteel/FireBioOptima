## Purpose: Run MSOM models with different data inputs; Full model
## Project: FireBioOptima

msom_full <- function(decay_rate,
                     taxa,
                     buff_width,
                     tws, # trait weights in percent, must be in order of severity, fri, patch size
                     n.fold = NULL, # number of folds used in k-fold CV
                     n.core = 1, # number of cores used for k-fold CV
                     save_model = T, # logical, do we save the full model or just the summary output
                     save_waic = T,
                     out_base = NULL, # base name for output files, no file suffix
                     return = F # return the model object?
) 
{
  library(tidyverse)
  library(spOccupancy)
  library(tidybayes)
  
  ## Slightly different code for birds and bats
  #### Birds ####
  if(taxa == 'birds')
  {
    ## read in y and covariate data
    data.msom <- read_rds('Data/bird_y_data.rds')
    covs <- read_rds(paste0('Data/bird_covs', buff_width, '.rds'))
    
    ## Pull out severity data
    sev <- filter(covs$sev, buffer == buff_width, 
                  decay == decay_rate) %>% 
      pivot_wider(id_cols = site, names_from = stat, values_from = cbi) %>% 
      mutate(point = str_sub(site, 1L, -6L)) %>% 
      filter(site %in% data.msom$sites)
    
    ## pull average frequency
    fri <- filter(covs$freq, buffer == buff_width, 
                  decay == decay_rate, stat == "mean") %>% 
      filter(site %in% data.msom$sites) %>% 
      select(site, fri_mn = value)
    
    ## average patch size
    pat <- filter(covs$pat, buffer == buff_width, 
                  decay == decay_rate, stat == "mean") %>% 
      filter(site %in% data.msom$sites) %>% 
      select(site, pat_mn = log_ha)
    
    ## pull out pyrodiversity
    pd <- filter(covs$pd, radius == buff_width, decay == decay_rate,
                 sev_w == tws[1], fri_w == tws[2], pat_w == tws[3]) %>% #, sea_w == tws[4]/100) %>% 
      filter(site %in% data.msom$sites) %>% 
      dplyr::select(FDis, site) %>% 
      arrange(site)
    
    ## get unique pts
    pts <- dplyr::select(sev, point, site) %>% 
      unique()
    
    ## add in some veg detection covariates
    veg <- read_csv('Data/veg_simple.csv') %>% 
      filter(point %in% pts$point) %>%
      mutate(can_cov = ifelse(is.na(can_cov), mean(can_cov, na.rm = T), can_cov),
             shr_cov = ifelse(is.na(shr_cov), mean(shr_cov, na.rm = T), shr_cov)) %>%
      ## expand to account for point x year stacking
      expand_grid(year = c(2021, 2022)) %>%
      mutate(site = paste0(point, ".", year)) %>% 
      unique() %>% 
      filter(site %in% data.msom$sites) %>%
      arrange(site)
    veg <- veg[with(veg, order(site)),]
    
    ## Check to make sure we align; although shouldn't matter with joining below
    if(!identical(sev$site, dimnames(data.msom$y)[[2]])) stop('points dont match')
    if(!identical(fri$site, dimnames(data.msom$y)[[2]])) stop('points dont match')
    if(!identical(pd$site, dimnames(data.msom$y)[[2]])) stop('points dont match')
    if(!identical(pat$site, dimnames(data.msom$y)[[2]])) stop('points dont match')
    
    
    ## put it all together
    occ.covs <- full_join(pd, covs$elev, by = "site") %>% 
      full_join(sev, by = c("site", "point")) %>% 
      full_join(fri, by = "site") %>% 
      full_join(pat, by = "site") %>%
      mutate(fdis_l = log(FDis + 0.01),
             year = as.integer(year), 
             # ## random effects must be integers
             pt = as.integer(as.factor(point)),
             tran = str_extract(point, "^\\w+"),
             tr = as.integer(as.factor(tran))
             ) %>%
      dplyr::select(site, point, pt, tran, tr, year, 
                    fdis = FDis, fdis_l, sev_mn = mean, sev_sd = sd, fri_mn, pat_mn,
                    elev = elevation)
    
    if(!identical(occ.covs$site, dimnames(data.msom$y)[[2]])) stop('points dont match')
    
    ## add to y and detection data
    data.msom$occ.covs <- occ.covs
    
    ## add veg to detection
    if(!identical(veg$site, dimnames(data.msom$det.covs$day)[[1]])) stop('points dont match')
    data.msom$det.covs$can_cov <- veg$can_cov
    data.msom$det.covs$shr_cov <- veg$shr_cov
    
    print(str(data.msom))
    
    ## Specifying the formula first
    occ.ms.formula <- ~ 
      scale(fdis) + I(scale(fdis)^2) + 
      scale(sev_mn) + I(scale(sev_mn)^2) +
      scale(fri_mn) + I(scale(fri_mn)^2) +
      scale(pat_mn) + I(scale(pat_mn)^2) +
      scale(elev) + I(scale(elev)^2) + 
      (1 | tr) + (1 | pt) 
    
    det.ms.formula <- ~ scale(day) + scale(tod) + I(scale(day)^2) +
      scale(can_cov) + scale(shr_cov)
    
    ## Setting up initial values
    N <- dim(data.msom$y)[1]
    ms.inits <- list(alpha.comm = 0, 
                     beta.comm = 0, 
                     beta = 0, 
                     alpha = 0,
                     tau.sq.beta = 1, 
                     tau.sq.alpha = 1, 
                     z = apply(data.msom$y, c(1, 2), max, na.rm = TRUE))
    
    ## default vaguley informative priors
    ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                      alpha.comm.normal = list(mean = 0, var = 2.72), 
                      tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                      tau.sq.alpha.ig = list(a = 0.1, b = 0.1),                      
                      ## the defaults are super flat, not helpful with sparse replication in groups
                      sigma.sq.psi.ig = list(a = c(10), b = c(1, 3)))
    
    # Approx. run time:  ~10 min w/ 3000 samples
    if(is.null(n.fold)) {
      out.ms <- msPGOcc(occ.formula = occ.ms.formula, 
                        det.formula = det.ms.formula, 
                        data = data.msom, 
                        inits = ms.inits, 
                        n.samples = 30000, 
                        priors = ms.priors, 
                        n.omp.threads = 1, 
                        verbose = TRUE, 
                        n.report = 2000,
                        # n.burn = 10000,
                        n.chains = 3,
                        n.thin = 20)
    } else {
      out.ms <- msPGOcc(occ.formula = occ.ms.formula, 
                        det.formula = det.ms.formula, 
                        data = data.msom, 
                        inits = ms.inits, 
                        n.samples = 20000, 
                        priors = ms.priors, 
                        n.omp.threads = 1, 
                        verbose = TRUE, 
                        n.report = 2000,
                        # n.burn = 10000,
                        n.thin = 20, 
                        n.chains = 3,
                        k.fold = n.fold, 
                        k.fold.threads = n.core,
                        k.fold.seed = 1234,
                        k.fold.only = F)
    }
    summary(out.ms, level = 'community')

    ## Calculate waic
    m.waic <- waicOcc(out.ms) %>% 
      as.data.frame()
    
    ## sum up species-level deviance 
    if("k.fold.deviance" %in% names(out.ms)){
      sum_deviance = sum(out.ms$k.fold.deviance)
      mean_deviance = mean(out.ms$k.fold.deviance)
    } else {
      sum_deviance = NULL
      mean_deviance = NULL
    }
    
    ## bundle data with model
    out <- list(d = data.msom,
                m = out.ms,
                waic = m.waic,
                taxa = taxa,
                decay_rate = decay_rate,
                buff_width = buff_width,
                tws = tws)
    
    if(save_model) {

      
      if(is.null(out_base)) {
        fn <- paste0("models/bird_msom_d", decay_rate, "_w", buff_width, 
                     '_s', tws[1], 'f', tws[2], 'p', tws[3], #'ss', tws[4], 
                     "_full.rds")
      } else {
        fn = paste0("models/", out_base, ".rds")
      }
      
      write_rds(out, fn,
                compress = "bz") 
    }
    
    ## Save waic score
    if(save_waic) {
      out.waic <- data.frame(taxa = taxa, decay_rate = decay_rate, buff_width = buff_width, 
                             tws = paste(tws, collapse = " ")) %>% 
        bind_cols(t(as.data.frame(m.waic))) %>% 
        mutate(sum_deviance = sum_deviance,
               mean_deviance = mean_deviance)
      
      ## append to existing dataframe
      ap <- ifelse(file.exists('Data/bird_full_waic.csv'), T, F)
      write_csv(out.waic, file = "Data/bird_full_waic.csv", append = ap)
    }
    
    if(is.null(out_base)) {
      capture.output(summary(out.ms), file = paste0("models/summaries/bird_msom_d", decay_rate, "_w", buff_width, 
                                                    '_s', tws[1], 'f', tws[2], 'p', tws[3], 
                                                    "_full_summary.txt"))
    } else
    {
      capture.output(summary(out.ms), file = paste0("models/summaries/", out_base, "_summary.txt"))
    }
    
    if(return) {return(out)}

  }

  #### Bats ####
  if(taxa == "bats")
  {
    ## read in y and covariate data
    data.msom <- read_rds('Data/bat_y_data.rds')
    covs <- read_rds(paste0('Data/ucb_covs', buff_width, '.rds'))
    
    sev <- filter(covs$sev, decay == decay_rate, win == 35) %>% 
      pivot_wider(id_cols = point, names_from = stat, values_from = cbi) %>% 
      ## expand to account for point x year stacking
      expand_grid(year = c(2021, 2022)) %>%
      ## make sure all sampled single digits have zeros before
      mutate(point = case_when(point == "ynp_2" ~ "YNP_02",
                               point == "ynp_1" ~ "YNP_01",
                               point == "ynp_4" ~ "YNP_04",
                               TRUE ~ toupper(point)),
             site = paste0(point, "_", year)) %>%  
      filter(site %in% data.msom$sites) 
    ## sort and arrange don't work exactly the same...
    sev <- sev[with(sev, order(site)),]
    
    ## Check to make sure we align
    if(!identical(sev$site, dimnames(data.msom$y)[[2]])) stop('points dont match')
    
    ## pull average TSLF
    fri <- filter(covs$freq, buffer == buff_width, 
                  decay == decay_rate, stat == "mean", win == 35) %>% 
      select(point, fri_mn = value) %>% 
      ## expand to account for point x year stacking
      expand_grid(year = c(2021, 2022)) %>%
      ## make sure all sampled single digits have zeros before
      mutate(point = case_when(point == "ynp_2" ~ "YNP_02",
                               point == "ynp_1" ~ "YNP_01",
                               point == "ynp_4" ~ "YNP_04",
                               TRUE ~ toupper(point)),
             site = paste0(point, "_", year)) %>%    
      filter(site %in% data.msom$sites) 
    fri <- fri[with(fri, order(site)),]
    
    if(!identical(fri$site, dimnames(data.msom$y)[[2]])) stop('points dont match')
    
    ## average patch size
    pat <- filter(covs$pat, buffer == buff_width, 
                  decay == decay_rate, stat == "mean", win == 35) %>% 
      select(point, pat_mn = log_ha) %>% 
      ## expand to account for point x year stacking
      expand_grid(year = c(2021, 2022)) %>%
      ## make sure all sampled single digits have zeros before
      mutate(point = case_when(point == "ynp_2" ~ "YNP_02",
                               point == "ynp_1" ~ "YNP_01",
                               point == "ynp_4" ~ "YNP_04",
                               TRUE ~ toupper(point)),
             site = paste0(point, "_", year)) %>%    
      filter(site %in% data.msom$sites) 
    pat <- pat[with(pat, order(site)),]
    
    if(!identical(pat$site, dimnames(data.msom$y)[[2]])) stop('points dont match')
    
    ## pull out pyrodiversity
    pd <- filter(covs$pd, radius == buff_width, decay == decay_rate, wind == 35,
                 sev_w == tws[1], fri_w == tws[2], pat_w == tws[3]) %>% 
      dplyr::select(FDis, point = name, year) %>% 
      mutate(point = case_when(point == "ynp_2" ~ "YNP_02",
                               point == "ynp_1" ~ "YNP_01",
                               point == "ynp_4" ~ "YNP_04",
                               TRUE ~ toupper(point)),
             site = paste0(point, "_", year)) %>%    
      filter(site %in% data.msom$sites) 
    pd <- pd[with(pd, order(site)),]
    
    if(!identical(pd$site, dimnames(data.msom$y)[[2]])) stop('points dont match')
    
    ## get unique pts
    pts <- dplyr::select(sev, point, site)
    
    ## add in some veg detection covariates
    veg <- read_csv('Data/Intermediate/model_data/veg_simple.csv') %>% 
      select(point = ucb_name, can_cov, shr_cov) %>% 
      ## expand to account for point x year stacking
      expand_grid(year = c(2021, 2022)) %>%
      # mutate(point = toupper(point)) %>% 
      mutate(point = case_when(point == "ynp_2" ~ "YNP_02",
                               point == "ynp_1" ~ "YNP_01",
                               point == "ynp_4" ~ "YNP_04",
                               TRUE ~ toupper(point)),
             site = paste0(point, "_", year)) %>% 
      unique() %>% 
      filter(site %in% data.msom$sites) %>% 
      mutate(can_cov = ifelse(is.na(can_cov), mean(can_cov, na.rm = T), can_cov),
             shr_cov = ifelse(is.na(shr_cov), mean(shr_cov, na.rm = T), shr_cov))
    veg <- veg[with(veg, order(site)),]
    
    ## add veg to detection
    if(!identical(veg$site, dimnames(data.msom$det.covs$day)[[1]])) stop('points dont match')
    data.msom$det.covs$can_cov <- veg$can_cov
    data.msom$det.covs$shr_cov <- veg$shr_cov
    
    elev <- covs$elev %>% 
      ## expand to account for point x year stacking
      expand_grid(year = c(2021, 2022)) %>%
      ## make sure all sampled single digits have zeros before
      mutate(point = case_when(point == "ynp_2" ~ "YNP_02",
                               point == "ynp_1" ~ "YNP_01",
                               point == "ynp_4" ~ "YNP_04",
                               TRUE ~ toupper(point)),
             site = paste0(point, "_", year)) %>%  
      filter(site %in% data.msom$sites) %>% 
      arrange(site)
    elev <- elev[with(elev, order(site)),]
    
    ## Check to make sure we align
    if(!identical(elev$site, dimnames(data.msom$y)[[2]])) stop('points dont match')
    
    ## put it all together
    occ.covs <- full_join(pd, elev, by = c("site", 'point', 'year')) %>% 
      full_join(sev, by = c("site", "point", 'year')) %>% 
      full_join(fri, by = c("site", "point", 'year')) %>% 
      ## same for pat
      full_join(pat, by = c("site", "point", 'year')) %>%
      mutate(fdis_l = log(FDis + 0.01),
             year = as.integer(year), 
             # ## random effects must be integers
             pt = as.integer(as.factor(point)),) %>%
      dplyr::select(site, point, pt, year, 
                    fdis = FDis, fdis_l, sev_mn = mean, sev_sd = sd, fri_mn, pat_mn,
                    elev = elevation)
    
    ## add to y and detection data
    data.msom$occ.covs <- occ.covs

        ## Specifying the formula first
    ## different formulas based on 'form' specification
    occ.ms.formula <- ~ 
      scale(fdis) + I(scale(fdis)^2) +
      scale(sev_mn) + I(scale(sev_mn)^2) +
      scale(fri_mn) + I(scale(fri_mn)^2) +
      scale(pat_mn) + I(scale(pat_mn)^2) +
      scale(elev) + I(scale(elev)^2) +
      (1 | pt)
    
    det.ms.formula <- ~ scale(day) + I(scale(day)^2) +
      scale(temp_avg) +
      scale(noise) + scale(illum) +
      scale(can_cov) + scale(shr_cov) 
    
    
    ## Setting up initial values
    N <- dim(data.msom$y)[1]
    ms.inits <- list(alpha.comm = 0, 
                     beta.comm = 0, 
                     beta = 0, 
                     alpha = 0,
                     tau.sq.beta = 1, 
                     tau.sq.alpha = 1, 
                     z = apply(data.msom$y, c(1, 2), max, na.rm = TRUE))
    
    ## default vaguely informative priors (quite flat)
    ms.priors <- list(beta.comm.normal = list(mean = 0, var = 2.72),
                      alpha.comm.normal = list(mean = 0, var = 2.72), 
                      tau.sq.beta.ig = list(a = 0.1, b = 0.1), 
                      tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                      ## the defaults are super flat, not helpful with sparse replication in groups
                      sigma.sq.psi = list(a = 1, b = 1),
                      sigma.sq.p = list(a = 1, b = 1),
                      sigma.sq.psi.ig = list(a = c(10), b = c(10)))
    
    # Approx. run time:  ~10 min w/ 3000 samples
    if(is.null(n.fold)) {
      
      out.ms <- msPGOcc(occ.formula = occ.ms.formula, 
                        det.formula = det.ms.formula, 
                        data = data.msom, 
                        inits = ms.inits, 
                        n.samples = 30000, # usually 30000
                        priors = ms.priors, 
                        n.omp.threads = 1, 
                        verbose = TRUE, 
                        n.report = 3000,
                        # n.burn = 10000,
                        n.thin = 20, 
                        n.chains = 3) #usually 3
    } else {

      ## Setting up initial values
      N <- dim(data.msom$y)[1]
      ms.inits <- list(alpha.comm = 0, 
                       beta.comm = 0, 
                       beta = 0, 
                       alpha = 0,
                       tau.sq.beta = 1, 
                       tau.sq.alpha = 1, 
                       z = apply(data.msom$y, c(1, 2), max, na.rm = TRUE))
      
      out.ms <- msPGOcc(occ.formula = occ.ms.formula, 
                        det.formula = det.ms.formula, 
                        data = data.msom, 
                        inits = ms.inits, 
                        n.samples = 30000, 
                        priors = ms.priors, 
                        n.omp.threads = 1, 
                        verbose = TRUE, 
                        n.report = 3000,
                        # n.burn = 10000,
                        n.thin = 20, 
                        n.chains = 3,
                        k.fold = n.fold, 
                        k.fold.threads = n.core,
                        k.fold.seed = 1234,
                        k.fold.only = F)
    }

    
    summary(out.ms, level = 'community')
    
    ## sum up species-level deviance 
    if("k.fold.deviance" %in% names(out.ms)){
      sum_deviance = sum(out.ms$k.fold.deviance)
      mean_deviance = mean(out.ms$k.fold.deviance, na.rm = T)
    } else {
      sum_deviance = NULL
      mean_deviance = NULL
    }
    
    ## Calculate waic
    m.waic <- waicOcc(out.ms) %>% 
      as.data.frame()
    
    ## bundle data with model
    out <- list(d = data.msom,
                m = out.ms,
                waic = m.waic,
                taxa = taxa,
                decay_rate = decay_rate,
                buff_width = buff_width,
                tws = tws)
    
    if(save_model) {
      
      ## file name dependent on out_base
      if(is.null(out_base)) {
        fn <- paste0("models/bat_msom_d", decay_rate, "_w", buff_width, 
                     '_s', tws[1], 'f', tws[2], 'p', tws[3], 
                     "_full.rds")
      } else {
        fn = paste0("models/", out_base, ".rds")
      }
      write_rds(out, fn, 
                compress = "bz") 
    }
    
    ## Save waic score
    if(save_waic) {
      out.waic <- data.frame(taxa = taxa, decay_rate = decay_rate, buff_width = buff_width, 
                             tws = paste(tws, collapse = " ")) %>% 
        bind_cols(t(as.data.frame(m.waic))) %>% 
        mutate(sum_deviance = sum_deviance,
               mean_deviance = mean_deviance)
      
      ## append to existing dataframe
      ap <- ifelse(file.exists('Data/Results/bat_full_waic.csv'), T, F)
      write_csv(out.waic, file = "Data/Results/bat_full_waic.csv", append = ap)
    }

    ## save summary
    if(is.null(out_base)) {
      capture.output(summary(out.ms), file = paste0("models/summaries/bat_msom_d", decay_rate, "_w", buff_width, 
                                                    '_s', tws[1], 'f', tws[2], 'p', tws[3],  
                                                    "_full_summary.txt"))
    } else
    {
      capture.output(summary(out.ms), file = paste0("models/summaries/", out_base, "_summary.txt"))
    }
    
    if(return) {return(out)}

  }
 
}
