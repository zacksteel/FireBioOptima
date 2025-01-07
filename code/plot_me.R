## Purpose: plot marginal effects on the community
## Project: FireBioOptima

plot_me <- function() 
{
  library(tidyverse)
  library(brms)
  library(colorspace)
  library(tidybayes)
  library(ggdist)
  library(patchwork)
  library(modelr)
  
  #### Birds Fdis ####
  ## First let's plot the marginal effects of fdis on richness
  ## Read in all the bird models
  bi.bal <- read_rds('models/step2/bird_2step_s2f1p0_fdis.rds')[c('fdis', 'dat')]

  ## get posterior epreds 
  bi.bal.pred.q <- bi.bal$dat %>% 
    data_grid(fdis_s = seq_range(fdis_s, n = 51),
              rich_sd = mean(rich_sd),
              elev_s = c(-1, 0, 1),
              tr = "new",
              pt = "new") %>%
    add_epred_draws(bi.bal$fdis, re_formula = NA, allow_new_levels = T) %>% 
    group_by(fdis_s, rich_sd, .draw) %>%
    summarize(.epred = mean(.epred)) %>%
    mutate(taxa = 'Birds',
           model = "Multi-trait",
           fit = "Quad") %>%
    mutate(## reverse scaling of fdis
      fdis = fdis_s * sd(bi.bal$dat$fdis) + mean(bi.bal$dat$fdis)
    ) 
  
  
  #### Bats Fdis ####
  ## Now for bats
  ba.fdis <- read_rds('models/step2/bat_2step_s0f1p2_fdis.rds')[c('fdis', 'dat')]
    
  ba.fdis.pred.q <- ba.fdis$dat %>%
    data_grid(fdis_s = seq_range(fdis_s, n = 51),
              rich_sd = mean(rich_sd),
              elev_s = c(-1, 0, 1),
              pt = "new") %>%
    add_epred_draws(ba.fdis$fdis, re_formula = NULL, allow_new_levels = T) %>% 
    group_by(fdis_s, rich_sd, .draw) %>%
    summarize(.epred = mean(.epred)) %>%
    mutate(taxa = 'Bats',
           model = "Patch Size",
           fit = "Quadratic") %>%
    mutate(## reverse scaling of fdis
      fdis = fdis_s * sd(ba.fdis$dat$fdis) + mean(ba.fdis$dat$fdis)
    )
  
  
  ## Now for plants
  pl.bal <- read_rds('models/step2/plant_2step_s2f1p0_fdis.rds')[c('fdis', 'dat')]
  
  ## get posterior epreds 
  pl.bal.pred.q <- pl.bal$dat %>% 
    data_grid(fdis_s = seq_range(fdis_s, n = 51),
              elev_s = c(-1, 0, 1),
              soil = unique(soil)) %>%
    add_epred_draws(pl.bal$fdis, re_formula = NULL) %>% 
    group_by(fdis_s, .draw) %>%
    summarize(.epred = mean(.epred)) %>%
    mutate(taxa = 'Plants',
           model = "Multi-trait",
           fit = "Quadratic") %>% 
    mutate(## reverse scaling of fdis
      fdis = fdis_s * sd(pl.bal$dat$fdis) + mean(pl.bal$dat$fdis)
    ) 
  
   
  dat2 <- bind_rows(
    mutate(bi.bal$dat, taxa = "Birds"),
    mutate(ba.fdis$dat, taxa = "Bats"),
    mutate(pl.bal$dat, taxa = "Plants",
           year = as.numeric(year),
           rich_mn = n)
    ) %>% 
    ## for each taxa, scale fdis by it's max
    group_by(taxa) %>%
    mutate(fdis = fdis / max(fdis)) %>%
    ungroup() 
  
  ## combine the best bird, bat, and plant models into a single panel
  pp.1 <- bind_rows(bi.bal.pred.q, ba.fdis.pred.q, pl.bal.pred.q) %>% 
    ## for each taxa, scale fdis by it's max
    group_by(taxa) %>%
    mutate(fdis = fdis / max(fdis)) %>%
    ungroup() %>%
    ggplot(aes(x = fdis, y = rich_mn, fill = taxa, color = taxa)) +
    geom_linerange(data = dat2[dat2$taxa == "Birds",], aes(ymin = rich_lower, ymax = rich_upper), alpha = 0.15) + 
    geom_linerange(data = dat2[dat2$taxa == "Bats",], aes(ymin = rich_lower, ymax = rich_upper), alpha = 0.25) + 
    geom_point(data = dat2, aes(color = taxa), size = 0.8, alpha = 0.7) +
    stat_lineribbon(aes(y = .epred),
                    .width = .95, alpha = 0.8) + #, color = 'grey30') +
    scale_fill_manual(values = c("#9C76A8", "#9CD3B3", "#594409"),
                      breaks = c("Birds", "Plants", "Bats"),
                      name = 'Taxa') +
    scale_color_manual(values = c("#9C76A8", "#9CD3B3", "#594409"),
                      breaks = c("Birds", "Plants", "Bats"),
                      name = 'Taxa') +
    theme_bw() + ylab('Species Richness') + xlab('Pyrodiversity (scaled)') +
    coord_cartesian(ylim = c(5, 40))
  
  
  ## Now let's plot the marginal effects of severity and frequency on richness
  ## for each taxa and model, pull out the posterior epreds and bind
  ## these into a single data frame
  #### Birds all ####
  bi.all <- read_rds('models/step2/bird_2step_s1f0p2_all.rds')[c('sev.q', 'fri.q', 'pat.q', 
                                                                 'dat')]
  ## Severity
  bi.sev.q <- bi.all$sev$data %>% 
    data_grid(sev_mn_s = seq_range(sev_mn_s, n = 51),
              rich_sd = mean(rich_sd),
              fri_mn_s = c(-1, 0, 1), #0,
              elev_s = c(-1, 0, 1), 
              pt = "new",
              tr = "new") %>%  #median(elev_s)) %>%
    add_epred_draws(bi.all$sev.q, re_formula = NULL, allow_new_levels = T) %>% 
    ## average over fri
    group_by(sev_mn_s, rich_sd, .draw) %>%
    summarize(.epred = mean(.epred)) %>% 
    mutate(taxa = 'Birds',
           model = "Mean Severity",
           fit = "Quadratic",
           ## reverse scaling of fdis
           sev_mn = sev_mn_s * sd(bi.all$dat$sev_mn) + mean(bi.all$dat$sev_mn)) 
  
  ## Frequency
  bi.fri.q <- bi.all$dat %>% 
    data_grid(fri_mn_s = seq_range(fri_mn_s, n = 51),
              rich_sd = mean(rich_sd),
              elev_s = c(-1, 0, 1), 
              pt = "new",
              tr = "new") %>%
    add_epred_draws(bi.all$fri.q, re_formula = NULL, allow_new_levels = T) %>% 
    ## average over fri
    group_by(fri_mn_s, rich_sd, .draw) %>%
    summarize(.epred = mean(.epred))   %>%
    mutate(taxa = 'Birds',
           model = "Mean Frequency",
           fit = "Quadratic",
           ## reverse scaling of fdis
           fri_mn = fri_mn_s * sd(bi.all$dat$fri_mn) + mean(bi.all$dat$fri_mn))

  #### Bats - All ####
  ## NOw the same for bat models
  ba.all <- read_rds('models/step2/bat_2step_s0f2p1_all.rds')
  
  ## Severity
  ba.sev.q <- ba.all$dat %>% 
    data_grid(sev_mn_s = seq_range(sev_mn_s, n = 51),
              rich_sd = mean(rich_sd),
              fri_mn_s = c(-1, 0, 1),
              elev_s = c(-1, 0, 1),
              pt = "new") %>%
    add_epred_draws(ba.all$sev.q, re_formula = NULL, allow_new_levels = T) %>% 
    group_by(sev_mn_s, rich_sd, .draw) %>%
    summarize(.epred = mean(.epred)) %>%
    mutate(taxa = 'Bats',
           model = "Mean Severity",
           fit = "Quadratic",
           ## reverse scaling of fdis
           sev_mn = sev_mn_s * sd(ba.all$dat$sev_mn) + mean(ba.all$dat$sev_mn))
  
  ## Frequency
  ba.fri.q <- ba.all$dat %>% 
    data_grid(fri_mn_s = seq_range(fri_mn_s, n = 51),
              rich_sd = mean(rich_sd),
              elev_s = c(-1, 0, 1)) %>%
    add_epred_draws(ba.all$fri.q, re_formula = NULL, allow_new_levels = T) %>%
    group_by(fri_mn_s, rich_sd, .draw) %>%
    summarize(.epred = mean(.epred)) %>%
    mutate(taxa = 'Bats',
           model = "Mean Frequency",
           fit = "Quadratic",
           ## reverse scaling of fdis
           fri_mn = fri_mn_s * sd(ba.all$dat$fri_mn) + mean(ba.all$dat$fri_mn))
  
    #### Plants - All ####
  ## repeat for plant models
  pl.all <- read_rds('models/step2/plant_2step_s0f0p3_all.rds')
  
  ## Severity
  pl.sev.q <- pl.all$dat %>% 
    data_grid(sev_s = seq_range(sev_s, n = 51),
              fri_s = c(-1, 0, 1),
              elev_s = c(-1, 0, 1),
              soil = soil) %>%
    add_epred_draws(pl.all$sev.q, re_formula = NULL) %>% 
    ## average over fri
    group_by(sev_s, .draw) %>%
    summarize(.epred = mean(.epred)) %>% 
    mutate(taxa = 'Plants',
           model = "Mean Severity",
           fit = "Quadratic",
           ## reverse scaling of fdis
           sev_mn = sev_s * sd(pl.all$dat$sev_mn) + mean(pl.all$dat$sev_mn))
  
  
  ## Frequency
  pl.fri.q <- pl.all$dat %>% 
    data_grid(fri_s = seq_range(fri_s, n = 51),
              elev_s = c(-1, 0, 1),
              soil = soil) %>%
    add_epred_draws(pl.all$fri.q, re_formula = NULL) %>%
    ## average over fri
    group_by(fri_s, .draw) %>%
    summarize(.epred = mean(.epred)) %>% 
    mutate(taxa = 'Plants',
           model = "Mean Frequency",
           fit = "Quadratic",
           ## reverse scaling of fdis
           fri_mn = fri_s * sd(pl.all$dat$fri) + mean(pl.all$dat$fri))
  
  ## combine into a single dataframe
  dat3 <- bind_rows(
    mutate(bi.all$dat, taxa = 'Birds'),
    mutate(ba.all$dat, taxa = 'Bats'),
    mutate(pl.all$dat, taxa = 'Plants',
           year = as.numeric(year),
           rich_mn = n,
           fri_mn = fri))
  
  ## combine the best bird, bat, and plant models into a single panel
  p.sev <- bind_rows(bi.sev.q, ba.sev.q, pl.sev.q) %>% 
    ggplot(aes(x = sev_mn, y = rich_mn, fill = taxa, color = taxa)) +
    geom_linerange(data = dat3[dat3$taxa == "Birds",], aes(ymin = rich_lower, ymax = rich_upper), alpha = 0.1) + 
    geom_linerange(data = dat3[dat3$taxa == "Bats",], aes(ymin = rich_lower, ymax = rich_upper), alpha = 0.25) + 
    geom_point(data = dat3, aes(color = taxa), size = 0.8, alpha = 0.7) +
    stat_lineribbon(aes(y = .epred),
                    .width = .95, alpha = 0.8) + 
    scale_fill_manual(
      values = c("#9C76A8", "#9CD3B3", "#594409"),
      breaks = c("Birds", "Plants", "Bats"),
      name = 'Taxa') +
    scale_color_manual(
      values = c("#9C76A8", "#9CD3B3", "#594409"),
      breaks = c("Birds", "Plants", "Bats"),
      name = 'Taxa') +
    theme_bw() + ylab('Species Richness') + xlab('Severity (CBI)') +
    coord_cartesian(ylim = c(5, 40), xlim = c(0, 3))  +
    ## remove legends
    theme(legend.position="none") +
    ## annotate taxa directly instead of legend
    annotate(geom = "text", x = 2.5, y = 28.5, label = "Birds") +
    annotate(geom = "text", x = 2.7, y = 21.5, label = "Plants") +
    annotate(geom = "text", x = 2.3, y = 14, label = "Bats")
  
  p.fri <- bind_rows(bi.fri.q, ba.fri.q, pl.fri.q) %>%
    ## order taxa so that plants are on top
    mutate(taxa = factor(taxa, levels = c("Plants", "Birds", "Bats"))) %>% 
    ggplot(aes(x = fri_mn, y = rich_mn, fill = taxa, color = taxa)) +
    ## add lineranges
    geom_linerange(data = dat3[dat3$taxa == "Birds",], aes(ymin = rich_lower, ymax = rich_upper), alpha = 0.1) +
    geom_linerange(data = dat3[dat3$taxa == "Bats",], aes(ymin = rich_lower, ymax = rich_upper), alpha = 0.25) +
    geom_point(data = dat3, aes(color = taxa), size = 0.8, alpha = 0.7) +
    stat_lineribbon(aes(y = .epred),
                    .width = .95, alpha = 0.8) + 
    scale_fill_manual(
      values = c("#9C76A8", "#9CD3B3", "#594409"),
      breaks = c("Birds", "Plants", "Bats"),
      name = 'Taxa') +
    scale_color_manual(
      values = c("#9C76A8", "#9CD3B3", "#594409"),
      breaks = c("Birds", "Plants", "Bats"),
      name = 'Taxa') +
    theme_bw() + ylab('Species Richness') + xlab('Frequency (years)') +
    coord_cartesian(ylim = c(5, 40), xlim = c(0, 40))  +
    ## remove legends
    theme(legend.position="none") 
  
  ## combine p.sev, p.fri, p.pat with patchwork
  pp2 <- 
    p.sev + 
    p.fri + 
    pp.1 +
    ## remove legends
    theme(legend.position="none") +
    # p.pat + 
    plot_annotation(tag_levels = 'A', tag_suffix = ')') +
    plot_layout(guides = "collect", axes = "collect") 
  
  ggsave("figures/me_all.png", pp2, width = 7, height = 5)
  
  ## save as tiff for publication quality
  ggsave("figures/Fig3.tiff", pp2, width = 7, height = 5)


}

