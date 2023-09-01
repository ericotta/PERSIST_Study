######################################################################
#                                                                    #
#      Analysis of anti-S IgG MSD data for PERSIST 1 paper           #
#                                                                    #
######################################################################





# Download packages ----



library(here)
library(ggplot2)
library(dplyr)
library(magrittr)
library(patchwork)
library(tidyr)
library(lubridate)
library(readxl)
library(ggdist)
library(forcats)
library(scales)
library(mgcv)
library(ggside)
library(confintr)
library(rsample)
library(purrr)
library(viridis)
library(coxed)
library(beepr)
library(ungeviz)
library(ggpubr)
library(rcompanion)
library(writexl)
library(grid)
library(beepr)



# User defined functions ----



# Not in
`%!in%` = Negate(`%in%`)

# Table 2
table2 <- function(x){
  vec <- length(x)
  names(vec) <- "Total"
  c(table(x, useNA = "always"), vec)
}

# Custom theme
theme_fausto <- function(){
  theme_bw() + 
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(size = 18, face = "bold"),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 15),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14, face = "bold"),
          legend.margin = margin(unit(0, "cm")),
          strip.text.x = element_text(margin = margin(0.5, 0, 0.5, 0, "mm")),
          strip.background = element_blank())
}

# Plotting GAMs in ggplot2 with random effects
# Adapted from https://gist.github.com/richardbeare/b679c38dcb644ec50ea34ac061200504
# The original function assumes a Gaussian model so it's not general
# I kept the top part of the github function that defined newdata based on presence of randomid
# I also kept the middle part to use the mgcv predict command rather than stats::predict
# Importantly, the bottom part of the github function assumed an identity link, so it 
# wouldn't work with binary data. Here I c/p the original code from ggplot2:::predictdf.glm
# Thus, it works normally when there's no randomid and incorporates randomid when there is one
# The bottom environment part redirects predictdf.glm, a
# built-in function within ggplot2, when it's called to use predictdf.gam

predictdf.gam <- function(model, xseq, se, level) {
  
  # Define newdata based on presents of a random id variable
  olddata <- model.frame(model)
  if (is.null(olddata$randomid)) {
    newdata = tibble(x = xseq)
  } else {
    newdata = tibble(x = xseq, randomid = olddata$randomid[1])
  }
  
  # Use mgcv to make a prediction
  pred <- predict(model, newdata = newdata,
                  se.fit = se, type = "link", exclude = "s(randomid)")
  
  # Original code from ggplot2:::predictdf.glm
  if (se) {
    std <- stats::qnorm(level/2 + 0.5)
    base::data.frame(x = xseq, y = model$family$linkinv(as.vector(pred$fit)), 
                     ymin = model$family$linkinv(as.vector(pred$fit - 
                                                             std * pred$se.fit)), ymax = model$family$linkinv(as.vector(pred$fit + 
                                                                                                                          std * pred$se.fit)), se = as.vector(pred$se.fit))
  }
  else {
    base::data.frame(x = xseq, y = model$family$linkinv(as.vector(pred)))
  }
  
}
environment(predictdf.gam) <- environment(ggplot2:::predictdf.glm)

# Custom ggplot panel borders from http://egret.psychol.cam.ac.uk/statistics/R/extensions/rnc_ggplot2_border_themes_2013_01.r
.pt <- 1 / 0.352777778 / 2

len0_null <- function(x) {
  if (length(x) == 0)  NULL
  else                 x
}

theme_border <- function(
  type = c("left", "right", "bottom", "top", "none"),
  colour = "black", size = 1, linetype = 1) {
  type <- match.arg(type, several.ok=TRUE)
  structure(
    list(type = type, colour = colour, size = size, linetype = linetype),
    class = c("theme_border", "element_blank", "element")
  )
}

element_grob.theme_border <- function(
  element, x = 0, y = 0, width = 1, height = 1,
  type = NULL,
  colour = NULL, size = NULL, linetype = NULL,
  ...) {
  if (is.null(type)) type = element$type
  xlist <- c()
  ylist <- c()
  idlist <- c()
  if ("bottom" %in% type) { # bottom
    xlist <- append(xlist, c(x, x+width))
    ylist <- append(ylist, c(y, y))
    idlist <- append(idlist, c(1,1))
  }
  if ("top" %in% type) { # top
    xlist <- append(xlist, c(x, x+width))
    ylist <- append(ylist, c(y+height, y+height))
    idlist <- append(idlist, c(2,2))
  }
  if ("left" %in% type) { # left
    xlist <- append(xlist, c(x, x))
    ylist <- append(ylist, c(y, y+height))
    idlist <- append(idlist, c(3,3))
  }
  if ("right" %in% type) { # right
    xlist <- append(xlist, c(x+width, x+width))
    ylist <- append(ylist, c(y, y+height))
    idlist <- append(idlist, c(4,4))
  }
  if (length(type)==0 || "none" %in% type) { # blank; cannot pass absence of coordinates, so pass a single point and use an invisible line
    xlist <- c(x,x)
    ylist <- c(y,y)
    idlist <- c(5,5)
    linetype <- "blank"
  }
  gp <- gpar(lwd = len0_null(size * .pt), col = colour, lty = linetype)
  element_gp <- gpar(lwd = len0_null(element$size * .pt), col = element$colour, lty = element$linetype)
  polylineGrob(
    x = xlist, y = ylist, id = idlist, ..., default.units = "npc",
    gp = modifyList(element_gp, gp),
  )
}



# Load the data ----



load("C:/Users/bustosfa/Rotation 2/Data/Final_MSD_IgG_data.RData")
load("C:/Users/bustosfa/Rotation 2/Data/Final_MSD_inhib_data.RData")

# Set LOD variables
LOD_IgG <- 21.05 * 0.00901
LOD_inhib <- 5



#### Figure S5 - IgG concentration over time by subgroup ----



# Make the labels
fig_s5_labels <- c(expression("10"^0), expression("10"^2), 
                   expression("10"^4), expression("10"^6))

# Make the plot
all_IgG2 %>%
  filter(pos_before == 0 & sample_bad == 0) %>%
  ggplot() +
  aes(x = antigen, y = binding, color = antigen, fill = antigen) +
  
  # Add in grid lines
  geom_hline(yintercept = c(1/10, 10, 1000, 100000), 
             linewidth = rel(0.5), color = "grey92") +

  stat_boxplot(aes(x = antigen, y = binding, color = antigen), 
             outlier.shape = NA, show.legend = FALSE, 
             geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = antigen, y = binding, color = antigen),
              width = 0.2, alpha = 0.4, size = 2.5) + 
  stat_summary(aes(x = antigen, y = binding, color = antigen),
               fun = "median", linewidth = 0.5, geom = "crossbar", show.legend = FALSE) + 
  
  facet_wrap(~time) +
  
  #ggtitle("Anti-spike IgG across antigens by time for all participants") +
  ylab("Anti-S IgG (BAU/ml)") + 
  xlab("SARS-CoV-2 variants") +
  theme_fausto() + 
  scale_color_manual(name = "SARS-CoV-2 variants",
                     values = c("Ancestral" = "#7209b7B2",
                                "Alpha" = "#048B90B2",
                                "Beta" = "#000000",
                                "Gamma" = "#DB9303B2",
                                "Delta" = "#0001AAB2",
                                "Omicron BA.1" = "#e00065")) + 
  scale_fill_manual(name = "SARS-CoV-2 variants",
                    values = c("Ancestral" = "#7209b7B2",
                               "Alpha" = "#048B90B2",
                               "Beta" = "#000000",
                               "Gamma" = "#DB9303B2",
                               "Delta" = "#0001AAB2",
                               "Omicron BA.1" = "#e00065")) + 
  scale_y_continuous(trans = "log10",
                     breaks = c(1, 100, 10000, 1000000),
                     labels = fig_s5_labels, 
                     limits = c(NA, 470000)) + 
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) + 
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 13.5),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) + 
  # good reference: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  geom_hline(yintercept = LOD_IgG, linetype = 2) + 
  stat_compare_means(ref.group = "Omicron BA.1", 
                     hide.ns = T, label = "p.signif",
                     size = 7, label.y.npc = 0.965) -> binding_pic_time

binding_pic_time



#### Figure S6 - ACE2 Inhibition over time by subgroup ----



all_inhib2 %>%
  filter(pos_before == 0 & sample_bad == 0) %>%
  ggplot() +
  aes(x = antigen, y = binding, color = antigen, fill = antigen) +

  stat_boxplot(aes(x = antigen, y = binding, color = antigen), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = antigen, y = binding, color = antigen),
              width = 0.2, alpha = 0.4, size = 2.5) + 
  stat_summary(aes(x = antigen, y = binding, color = antigen),
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) + 
  
  facet_wrap(~time) +
  
  ylab("ACE2 inhibition (%)") + 
  xlab("SARS-CoV-2 variants") +
  theme_fausto() + 
  scale_color_manual(name = "SARS-CoV-2 variants",
                     values = c("Ancestral" = "#7209b7B2",
                                "Alpha" = "#048B90B2",
                                "Beta" = "#000000",
                                "Gamma" = "#DB9303B2",
                                "Delta" = "#0001AAB2",
                                "Omicron BA.1" = "#e00065")) + 
  scale_fill_manual(name = "SARS-CoV-2 variants",
                    values = c("Ancestral" = "#7209b7B2",
                               "Alpha" = "#048B90B2",
                               "Beta" = "#000000",
                               "Gamma" = "#DB9303B2",
                               "Delta" = "#0001AAB2",
                               "Omicron BA.1" = "#e00065")) + 
  scale_y_continuous(breaks = seq(0, 100, 25),
                     labels = seq(0, 100, 25),
                     limits = c(0, 104)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) + 
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 13.5),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) + 
  # good reference: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  geom_hline(yintercept = LOD_inhib, linetype = 2) + 
  stat_compare_means(ref.group = "Omicron BA.1", 
                     hide.ns = T, label = "p.signif",
                     size = 7, label.y.npc = 0.965) -> inhib_pic_time

inhib_pic_time



#### Compare IgG concentration and ACE2 inhibition at pd2 and 6mpd2 ---- 



# Decrease from pd2 to 6mpds: IgG concentration
all_IgG2 %>%
  group_by(time) %>%
  summarise(titers = median(binding, na.rm = T)) %>% 
  filter(time %in% c("Post-dose 2", "6 months post-dose 2")) %>%
  pivot_wider(names_from = time, 
              values_from = titers) %>%
  mutate(decrease = paste0(round(100 * (1 - `6 months post-dose 2`/`Post-dose 2`), 1), "%"),
         fold = `Post-dose 2` / `6 months post-dose 2`)

# Decrease from pd2 to 6mpds: ACE2 inhibition
all_inhib2 %>%
  group_by(time) %>%
  summarise(titers = median(binding)) %>% 
  filter(time %in% c("Post-dose 2", "6 months post-dose 2")) %>%
  pivot_wider(names_from = time, 
              values_from = titers) %>%
  mutate(decrease = paste0(round(100 * (1 - `6 months post-dose 2`/`Post-dose 2`), 1), "%"),
         fold = `Post-dose 2` / `6 months post-dose 2`)



#### Are more Omicron than Wuhan samples before the LOD? ----



df_LOD <- data.frame("Above LOD" = c(nrow(all_inhib2 %>%
                                            filter(pos_before == 0 & sample_bad == 0) %>%
                                            filter(antigen == "Ancestral")) - 
                                       nrow(all_inhib2 %>%
                                              filter(pos_before == 0 & sample_bad == 0) %>%
                                              filter(antigen == "Ancestral") %>%
                                              filter(binding < LOD_inhib)),
                                     (nrow(all_inhib2 %>%
                                             filter(pos_before == 0 & sample_bad == 0) %>%
                                             filter(antigen == "Omicron BA.1"))) - 
                                       nrow(all_inhib2 %>%
                                              filter(pos_before == 0 & sample_bad == 0) %>%
                                              filter(antigen == "Omicron BA.1") %>%
                                              filter(binding < LOD_inhib))),
                     "Below LOD" = c(nrow(all_inhib2 %>%
                                            filter(pos_before == 0 & sample_bad == 0) %>%
                                            filter(antigen == "Ancestral") %>%
                                            filter(binding < LOD_inhib)),
                                     nrow(all_inhib2 %>%
                                            filter(pos_before == 0 & sample_bad == 0) %>%
                                            filter(antigen == "Omicron BA.1") %>%
                                            filter(binding < LOD_inhib))), 
                     row.names = c("Ancestral", "Omicron BA.1"))
df_LOD
fisher.test(df_LOD) # p-value = 0.5881...if you remove the 2 bad conditions pvalue=0.1529
rm(df_LOD, fig_s5_labels)

# Make correlation plot
corr_pre_IgG <- all_IgG2 %>%
  filter(pos_before == 0 & sample_bad == 0) %>%
  filter(time != "Post-dose 4") %>%
  arrange(id, time, antigen) %>%
  as.data.frame()

corr_pre_inhib <- all_inhib2 %>%
  filter(pos_before == 0 & sample_bad == 0) %>%
  filter(time != "Post-dose 4")%>%
  arrange(id, time, antigen) %>%
  rename(inhib = binding) %>%
  as.data.frame()

all.equal(corr_pre_IgG$id, corr_pre_inhib$id)
all.equal(corr_pre_IgG$time, corr_pre_inhib$time)
all.equal(corr_pre_IgG$antigen, corr_pre_inhib$antigen)

corr_df <- cbind(corr_pre_IgG, inhib = corr_pre_inhib$inhib) # 984 rows
corr_df <- corr_df[corr_df$binding < 4e+06 * 0.00901,] # 981, to drop outliers
corr_df <- na.omit(corr_df) # 975, can't do correlation with missing values



#### Overall correlation ----



# Note: This section is commented out because the code takes a long time to run. Some results and given in the comments. Uncomment the code to run the analysis.



# # Assuming independence at the row level
# ci_cor(corr_df$binding, corr_df$inhib, seed = 1989,
#        method = "kendall", type = "bootstrap", boot_type = "bca") # 0.53, (0.49, 0.56)
# # 
# # 	Two-sided 95% bootstrap confidence interval for the true Kendall correlation
# # 	coefficient based on 9999 bootstrap replications and the bca method
# # 
# # Sample estimate: 0.5269546
# # Confidence interval:
# #      2.5%      97.5%
# # 0.4920645    0.5585753
# beep(sound = 3)

# # Doing a cluster bootstrap to figure out the correlation
# # Step 1, make the appropriate bootstrap datasets, all contained in bs
# # From: https://www.r-bloggers.com/2018/08/bootstrapping-clustered-data/
# set.seed(666)
# corr_df %>%
#   mutate(id2 = factor(id)) %>%
#   dplyr::select(id2, binding, inhib) %>%
#   nest(data = -id2) %>%
#   bootstraps(times = 9999) -> bs
# 
# # For step 2, calculate the correlation over the bootstrap datasets
# # From: https://www.r-bloggers.com/2018/08/bootstrapping-clustered-data/
# bs_cor <- map(bs$splits, ~as_tibble(.) %>% unnest(cols = c(id2, data)) %>%
#                 dplyr::summarize(cor = cor(.$binding, .$inhib, method = "kendall"))) %>%
#   bind_rows(.id = 'boots')
# 
# # For step 3, calculate BCa confidence intervals
# # From: https://stackoverflow.com/questions/55401615/r-calculate-bca-from-vector-of-bootstrapped-results
# round(bca(bs_cor$cor), 2) #0.46, 0.59 -> Final answer: 0.53 (0.46, 0.59)
# beep(sound = 3)


# # Now repeat for the Ancestral variant
# ci_cor(corr_df[corr_df$antigen == "Ancestral",]$binding,
#        corr_df[corr_df$antigen == "Ancestral",]$inhib, seed = 1989,
#        method = "kendall", type = "bootstrap", boot_type = "bca") # 0.60
# beep(sound = 3)
# 
# set.seed(666)
# corr_df %>%
#   filter(antigen == "Ancestral") %>%
#   mutate(id2 = factor(id)) %>%
#   dplyr::select(id2, binding, inhib) %>%
#   nest(data = -id2) %>%
#   bootstraps(times = 9999) -> bs_wuhan
# 
# bs_cor_wuhan <- map(bs_wuhan$splits, ~as_tibble(.) %>% unnest(cols = c(id2, data)) %>%
#                 dplyr::summarize(cor = cor(.$binding, .$inhib, method = "kendall"))) %>%
#   bind_rows(.id = 'boots')
# 
# round(bca(bs_cor_wuhan$cor), 2) #0.53, 0.67 -> Final answer: 0.60 (0.53, 0.67)
# wuhan_int <- bca(bs_cor_wuhan$cor)
# beep(sound = 3)



# # Now repeat for the Omicron variant
# ci_cor(corr_df[corr_df$antigen == "Omicron BA.1",]$binding,
#        corr_df[corr_df$antigen == "Omicron BA.1",]$inhib, seed = 1989,
#        method = "kendall", type = "bootstrap", boot_type = "bca") # 0.25
# 
# set.seed(666)
# corr_df %>%
#   filter(antigen == "Omicron BA.1") %>%
#   mutate(id2 = factor(id)) %>%
#   dplyr::select(id2, binding, inhib) %>%
#   nest(data = -id2) %>%
#   bootstraps(times = 9999) -> bs_omi
# 
# bs_cor_omi <- map(bs_omi$splits, ~as_tibble(.) %>% unnest(cols = c(id2, data)) %>%
#                 dplyr::summarize(cor = cor(.$binding, .$inhib, method = "kendall"))) %>%
#   bind_rows(.id = 'boots')
# 
# round(bca(bs_cor_omi$cor), 2) # 0.14, 0.34 -> Final answer: 0.25 (0.14, 0.34)
# beep(sound = 3)
# 
# # Calculate empirical p-value for Omicron vs Wuhan
# omi_int <- bca(bs_cor_omi$cor)
# sum(between(bs_cor_omi$cor, wuhan_int[1], wuhan_int[2])) / length(bs_cor_omi$cor) #0
# 1 - sum(bs_cor_wuhan$cor > omi_int[2]) / length(bs_cor_wuhan$cor) #0



# # Now repeat for the Beta variant
# ci_cor(corr_df[corr_df$antigen == "Beta",]$binding,
#        corr_df[corr_df$antigen == "Beta",]$inhib, seed = 1989,
#        method = "kendall", type = "bootstrap", boot_type = "bca") # 0.52
# 
# set.seed(666)
# corr_df %>%
#   filter(antigen == "Beta") %>%
#   mutate(id2 = factor(id)) %>%
#   dplyr::select(id2, binding, inhib) %>%
#   nest(data = -id2) %>%
#   bootstraps(times = 9999) -> bs_beta
# 
# bs_cor_beta <- map(bs_beta$splits, ~as_tibble(.) %>% unnest(cols = c(id2, data)) %>%
#                 dplyr::summarize(cor = cor(.$binding, .$inhib, method = "kendall"))) %>%
#   bind_rows(.id = 'boots')
# 
# round(bca(bs_cor_beta$cor), 2) #0.44  0.59 -> Final answer: 0.52 (0.44, 0.59)
# beep(sound = 3)
# 
# # Calculate empirical p-value for Beta vs Omicron
# sum(between(bs_cor_beta$cor, omi_int[1], omi_int[2])) / length(bs_cor_beta$cor) #0
# 1 - sum(bs_cor_beta$cor > omi_int[2]) / length(bs_cor_beta$cor) #0



# # Now repeat for the Gamma variant
# ci_cor(corr_df[corr_df$antigen == "Gamma",]$binding,
#        corr_df[corr_df$antigen == "Gamma",]$inhib, seed = 1989,
#        method = "kendall", type = "bootstrap", boot_type = "bca") # 0.53
# 
# set.seed(666)
# corr_df %>%
#   filter(antigen == "Gamma") %>%
#   mutate(id2 = factor(id)) %>%
#   dplyr::select(id2, binding, inhib) %>%
#   nest(data = -id2) %>%
#   bootstraps(times = 9999) -> bs_gamma
# 
# bs_cor_gamma <- map(bs_gamma$splits, ~as_tibble(.) %>% unnest(cols = c(id2, data)) %>%
#                 dplyr::summarize(cor = cor(.$binding, .$inhib, method = "kendall"))) %>%
#   bind_rows(.id = 'boots')
# 
# round(bca(bs_cor_gamma$cor), 2) #0.45  0.61 -> Final answer: 0.53 (0.45, 0.61)
# beep(sound = 3)
# 
# # Calculate empirical p-value for Gamma vs Omicron
# sum(between(bs_cor_gamma$cor, omi_int[1], omi_int[2])) / length(bs_cor_gamma$cor) #0
# 1 - sum(bs_cor_gamma$cor > omi_int[2]) / length(bs_cor_gamma$cor) #0



#### Table S8 ----



# # Make it a for loop across the different variants of concern!
# antigens <- levels(corr_df$antigen)
# 
# corrs_df <- data.frame(Antigen = antigens, Corr = NA, LL = NA, UL = NA)
# 
# for(i in 1:length(antigens)){
#   corrs_df$Corr[i] <- round(ci_cor(corr_df[corr_df$antigen == antigens[i],]$binding,
#                                    corr_df[corr_df$antigen == antigens[i],]$inhib,
#                                    seed = 1989, method = "kendall",
#                                    type = "bootstrap", boot_type = "bca")$estimate, 2)
# 
#   set.seed(666)
#   corr_df %>%
#   filter(antigen == antigens[i]) %>%
#   mutate(id2 = factor(id)) %>%
#   dplyr::select(id2, binding, inhib) %>%
#   nest(data = -id2) %>%
#   bootstraps(times = 9999) -> bs_antigen
# 
#   bs_cor_antigen <- map(bs_antigen$splits, ~as_tibble(.) %>% unnest(cols = c(id2, data)) %>%
#                 dplyr::summarize(cor = cor(.$binding, .$inhib, method = "kendall"))) %>%
#   bind_rows(.id = 'boots')
# 
#   answer <- bca(bs_cor_antigen$cor)
# 
#   corrs_df$LL[i] <- round(answer[1], 2)
#   corrs_df$UL[i] <- round(answer[2], 2)
# }
# beep(sound = 3)
# 
# save(corrs_df, file = "correlations_variant.RData")
# 
# Calculate the width of the CI
# corrs_df$width = corrs_df$UL - corrs_df$LL



#### Table S9 ----



# # Now do this across immune groups within variants of concern
# immunegroups <- levels(corr_df$new_group)
# immunegroups <- immunegroups[immunegroups != "Combined immunodeficiency"]
# 
# corrs_immune_df <- expand.grid(immunegroups, antigens)
# names(corrs_immune_df) <- c("Immune group", "Antigen")
# corrs_immune_df$UL <- corrs_immune_df$LL <- corrs_immune_df$Corr <- NA
# 
# for(i in 1:length(antigens)){
#   for(j in 1:length(immunegroups)){
#     n = length(immunegroups) # n = number of immune groups
#   corrs_immune_df$Corr[(i*n + j) - n] <-
#       round(ci_cor(corr_df[corr_df$antigen == antigens[i] &
#                            corr_df$new_group == immunegroups[j],]$binding,
#                    corr_df[corr_df$antigen == antigens[i] &
#                            corr_df$new_group == immunegroups[j],]$inhib,
#                    seed = 1989, method = "kendall",
#                    type = "bootstrap", boot_type = "bca")$estimate, 2)
# 
#   set.seed(666)
#   corr_df %>%
#   filter(antigen == antigens[i] & new_group == immunegroups[j]) %>%
#   mutate(id2 = factor(id)) %>%
#   dplyr::select(id2, binding, inhib) %>%
#   nest(data = -id2) %>%
#   bootstraps(times = 9999) -> bs_antigen
# 
#   bs_cor_antigen <- map(bs_antigen$splits, ~as_tibble(.) %>% unnest(cols = c(id2, data)) %>%
#                 dplyr::summarize(cor = cor(.$binding, .$inhib, method = "kendall"))) %>%
#   bind_rows(.id = 'boots')
# 
#   answer <- bca(bs_cor_antigen$cor)
# 
#   corrs_immune_df$LL[(i*n + j) - n] <- round(answer[1], 2)
#   corrs_immune_df$UL[(i*n + j) - n] <- round(answer[2], 2)
# 
#   print(paste0((i*n + j) - n, " of ", length(antigens) * length(immunegroups)))
# 
#   }
# }
# 
# beep(sound = 3)
# save(corrs_immune_df, file = "correlations_variant_immunegroup.RData")
# 
# # for loop for the combined immunity group
# comb_group_df <- corrs_immune_df[1:length(antigens),]
# comb_group_df$Antigen <- antigens
# comb_group_df$Corr <- comb_group_df$LL <- comb_group_df$UL <- NA
# comb_group_df$`Immune group` <- "Combined immunodeficiency"
# 
# for(i in 1:length(antigens)){
#       answer <- ci_cor(corr_df[corr_df$antigen == antigens[i] &
#                                corr_df$new_group == "Combined immunodeficiency",]$binding,
#                        corr_df[corr_df$antigen == antigens[i] &
#                                corr_df$new_group == "Combined immunodeficiency",]$inhib,
#                        seed = 1989, method = "kendall",
#                        type = "bootstrap", boot_type = "bca")
# 
#   comb_group_df$Corr[i] <- round(answer$estimate, 2)
#   comb_group_df$LL[i] <- round(answer$interval[1], 2)
#   comb_group_df$UL[i] <- round(answer$interval[2], 2)
# 
#   print(i)
#   }
# 
# # Put it all together
# together <- rbind(corrs_immune_df, comb_group_df)
# together$`Immune group` <- factor(together$`Immune group`, 
#                                   levels = c("Healthy volunteer", "Antibody deficiency",
#                                              "PIRD", "Combined immunodeficiency",
#                                              "Other IEI", "Other immune disorder"))
# together$Antigen <- factor(together$Antigen, 
#                            levels = c("Ancestral", "Alpha", "Beta",
#                                       "Gamma", "Delta", "Omicron BA.1"))
# 
# together %>%
#   arrange(`Immune group`, Antigen) -> together



#### Table S6 ----



# # Now do this across immune groups 
# immunegroups <- levels(corr_df$new_group)
# 
# coors_immunegroups_df <- data.frame("Immune group" = immunegroups)
# coors_immunegroups_df$UL <- coors_immunegroups_df$LL <- coors_immunegroups_df$Corr <- NA
# 
# for(i in 1:length(immunegroups)){
#   coors_immunegroups_df$Corr[i] <- round(ci_cor(corr_df[corr_df$new_group ==
#                                                           immunegroups[i],]$binding,
#                                    corr_df[corr_df$new_group == immunegroups[i],]$inhib,
#                                    seed = 1989, method = "kendall",
#                                    type = "bootstrap", boot_type = "bca")$estimate, 2)
# 
#   set.seed(666)
#   corr_df %>%
#   filter(new_group == immunegroups[i]) %>%
#   mutate(id2 = factor(id)) %>%
#   dplyr::select(id2, binding, inhib) %>%
#   nest(data = -id2) %>%
#   bootstraps(times = 9999) -> bs_antigen
# 
#   bs_cor_antigen <- map(bs_antigen$splits, ~as_tibble(.) %>% unnest(cols = c(id2, data)) %>%
#                 dplyr::summarize(cor = cor(.$binding, .$inhib, method = "kendall"))) %>%
#   bind_rows(.id = 'boots')
# 
#   answer <- bca(bs_cor_antigen$cor)
# 
#   coors_immunegroups_df$LL[i] <- round(answer[1], 2)
#   coors_immunegroups_df$UL[i] <- round(answer[2], 2)
# 
#   print(paste0(i, " of ", length(immunegroups)))
# }
# beep(sound = 3)
# 
# save(coors_immunegroups_df, file = "correlations_immunegroups.RData")

# Comparing correlations for HVs and combined immunodeficiencies
## This was useful with data6, but with data7, all CIs overlap substantially
# # Calculate corr for the healthy volunteers
#   corr_hv <- ci_cor(corr_df[corr_df$new_group == "Healthy volunteer",]$binding,
#                                    corr_df[corr_df$new_group == "Healthy volunteer",]$inhib,
#                                    seed = 1989, method = "kendall",
#                                    type = "bootstrap", boot_type = "bca")$estimate
# 
# # Calculate CIs for immune dysreg
#   set.seed(666)
#   corr_df %>%
#   filter(new_group == "Immune dysregulation") %>%
#   mutate(id2 = factor(id)) %>%
#   dplyr::select(id2, binding, inhib) %>%
#   nest(data = -id2) %>%
#   bootstraps(times = 9999) -> bs_ID
# 
#   bs_cor_ID <- map(bs_ID$splits, ~as_tibble(.) %>% unnest(cols = c(id2, data)) %>%
#                 dplyr::summarize(cor = cor(.$binding, .$inhib, method = "kendall"))) %>%
#   bind_rows(.id = 'boots')
# 
# # Now calculate an empirical p-value
# sum(bs_cor_ID$cor <= corr_hv) / length(bs_cor_ID$cor)
# 
#   # Calculate CIs for combined immunodeficiency
#   set.seed(666)
#   corr_df %>%
#   filter(new_group == "Combined immunodeficiency") %>%
#   mutate(id2 = factor(id)) %>%
#   dplyr::select(id2, binding, inhib) %>%
#   nest(data = -id2) %>%
#   bootstraps(times = 9999) -> bs_CI
# 
#   bs_cor_CI <- map(bs_CI$splits, ~as_tibble(.) %>% unnest(cols = c(id2, data)) %>%
#                 dplyr::summarize(cor = cor(.$binding, .$inhib, method = "kendall"))) %>%
#   bind_rows(.id = 'boots')
#   
# # Now calculate an empirical p-value
# sum(bs_cor_CI$cor >= corr_hv) / length(bs_cor_CI$cor)
# 



#### Table S7 ----



# # Now do this across timepoints
# timess <- levels(corr_df$time)
# 
# coors_time_df <- data.frame("time" = timess)
# coors_time_df$UL <- coors_time_df$LL <- coors_time_df$Corr <- NA
# 
# for(i in 1:length(timess)){
#   coors_time_df$Corr[i] <-         round(ci_cor(corr_df[corr_df$time == timess[i],]$binding,
#                                    corr_df[corr_df$time == timess[i],]$inhib,
#                                    seed = 1989, method = "kendall",
#                                    type = "bootstrap", boot_type = "bca")$estimate, 2)
# 
#   set.seed(666)
#   corr_df %>%
#   filter(time == timess[i]) %>%
#   mutate(id2 = factor(id)) %>%
#   dplyr::select(id2, binding, inhib) %>%
#   nest(data = -id2) %>%
#   bootstraps(times = 9999) -> bs_time
# 
#   bs_cor_time <- map(bs_time$splits, ~as_tibble(.) %>% unnest(cols = c(id2, data)) %>%
#                 dplyr::summarize(cor = cor(.$binding, .$inhib, method = "kendall"))) %>%
#   bind_rows(.id = 'boots')
# 
#   answer <- bca(bs_cor_time$cor)
# 
#   coors_time_df$LL[i] <- round(answer[1], 2)
#   coors_time_df$UL[i] <- round(answer[2], 2)
# 
#   print(paste0("Done with ", i, " out of ", length(timess)))
# 
# }
# 
# beep(sound = 4)
# coors_time_df$width <- coors_time_df$UL - coors_time_df$LL
# save(coors_time_df, file = "correlations_time.RData")


# # Now do right after vs 6 months after
# coors_time2_df <- data.frame("time2" = c("Right after", "6 months after"))
# coors_time2_df$UL <- coors_time2_df$LL <- coors_time2_df$Corr <- NA
# 
# # Right away
# coors_time2_df$Corr[1] <- round(ci_cor(corr_df[corr_df$time %in%
#                                                 c("Post-dose 2", "Post-dose 3"),]$binding,
#                                        corr_df[corr_df$time  %in%
#                                                 c("Post-dose 2", "Post-dose 3"),]$inhib,
#                                        seed = 1989, method = "kendall",
#                                        type = "bootstrap", boot_type = "bca")$estimate, 2)
# 
#   set.seed(666)
#   corr_df %>%
#   filter(time %in% c("Post-dose 2", "Post-dose 3")) %>%
#   mutate(id2 = factor(id)) %>%
#   dplyr::select(id2, binding, inhib) %>%
#   nest(data = -id2) %>%
#   bootstraps(times = 9999) -> bs_time2
# 
#   bs_cor_time2 <- map(bs_time2$splits, ~as_tibble(.) %>% unnest(cols = c(id2, data)) %>%
#                 dplyr::summarize(cor = cor(.$binding, .$inhib, method = "kendall"))) %>%
#   bind_rows(.id = 'boots')
# 
#   answer_rightaway <- bca(bs_cor_time2$cor)
# 
#   coors_time2_df$LL[1] <- round(answer_rightaway[1], 2)
#   coors_time2_df$UL[1] <- round(answer_rightaway[2], 2)
# 
# 
# # 6 months later
#   coors_time2_df$Corr[2] <- round(ci_cor(corr_df[corr_df$time %in%
#                                                 c("6 months post-dose 2",
#                                                   "6 months post-dose 3"),]$binding,
#                                        corr_df[corr_df$time  %in%
#                                                 c("6 months post-dose 2",
#                                                   "6 months post-dose 3"),]$inhib,
#                                        seed = 1989, method = "kendall",
#                                        type = "bootstrap", boot_type = "bca")$estimate, 2)
# 
#   set.seed(666)
#   corr_df %>%
#   filter(time %in% c("6 months post-dose 2", "6 months post-dose 3")) %>%
#   mutate(id2 = factor(id)) %>%
#   dplyr::select(id2, binding, inhib) %>%
#   nest(data = -id2) %>%
#   bootstraps(times = 9999) -> bs_time23
# 
#   bs_cor_time23 <- map(bs_time23$splits, ~as_tibble(.) %>% unnest(cols = c(id2, data)) %>%
#                 dplyr::summarize(cor = cor(.$binding, .$inhib, method = "kendall"))) %>%
#   bind_rows(.id = 'boots')
# 
#   answer_later <- bca(bs_cor_time23$cor)
# 
#   coors_time2_df$LL[2] <- round(answer_later[1], 2)
#   coors_time2_df$UL[2] <- round(answer_later[2], 2)
# 
# Now calculate an empirical p-value
# sum(bs_cor_time23$cor >= coors_time2_df$Corr[1]) / length(bs_cor_time23$cor)

# Remove correlation data objects
rm(bs, bs_antigen, bs_beta, bs_cor, bs_cor_antigen, bs_cor_beta, bs_cor_gamma,
   bs_cor_omi, bs_cor_time, bs_cor_time2, bs_cor_time23, bs_cor_wuhan, bs_gamma,
   bs_omi, bs_time, bs_time2, bs_time23, bs_wuhan)
rm(comb_group_df, coors_immunegroups_df, coors_time_df, coors_time2_df,
   corr_pre_IgG, corr_pre_inhib, corrs_df, corrs_immune_df, together)
rm(answer, answer_later, answer_rightaway, antigens, i, j, immunegroups,
   million_labels, n, omi_int, timess, wuhan_int)


# Now let's do some plotting
corr_df$inhib2 <- corr_df$inhib/100

million_labels <- c(expression("0x10"^6), expression("1x10"^6), 
                    expression("2x10"^6), expression("3x10"^6))



#### Figure S7 ----



overall_corr_re_voc2 <- ggplot(corr_df, 
                               aes(x = binding, y = inhib2)) +
  geom_hline(yintercept = LOD_inhib/100, linetype = 2) +
  geom_vline(xintercept = LOD_IgG, linetype = 2) +
  geom_point(alpha = 0.25, aes(color = factor(antigen)), size = 2.5) + 
  theme_fausto() + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     labels = seq(0, 1, 0.25)*100,
                     limits = c(0,1),
                     expand = c(0.01, 0)) +
  scale_x_continuous(breaks = seq(0, 25000, 5000),
                     labels = c("0", "5,000", "10,000", "15,000", "20,000", "25,000"),
                     limits = c(0, 28000),
                     expand = c(0.015, 0)) +
  xlab("Anti-S IgG (BAU/ml)") +
  ylab("ACE2 inhibition (%)") + 
  stat_smooth(method = "gam", 
              se = F, aes(group = factor(antigen), color = factor(antigen),
                          linetype = factor(antigen), size = factor(antigen),
                          randomid = factor(id)),
              formula = y ~ s(x, bs = "tp") + s(randomid, bs = "re"),
              method.args = list(family = "quasibinomial", method = "REML")) + 
  scale_color_manual(name = "SARS-CoV-2 variants",
                     values = c("Ancestral" = "#7209b7B2",
                                "Alpha" = "#048B90B2",
                                "Beta" = "#000000",
                                "Gamma" = "#DB9303B2",
                                "Delta" = "#0001AAB2",
                                "Omicron BA.1" = "#e00065")) + 
  scale_size_manual(values = c("Ancestral" = 2,
                               "Alpha" = 2,
                               "Beta" = 2,
                               "Gamma" = 2,
                               "Delta" = 2,
                               "Omicron BA.1" = 2),
                    guide = "none") + 
  scale_linetype_manual(values = c("Ancestral" = 1,
                                   "Alpha" = 1,
                                   "Beta" = 1,
                                   "Gamma" = 1,
                                   "Delta" = 1,
                                   "Omicron BA.1" = 1),
                        guide = "none") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4), 
                               nrow = 1, byrow = T)) + 
  theme(legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, -10, -5, -10),
        legend.key.width = unit(1,"cm"), 
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) 

overall_corr_re_voc2



#### Compare the median binding and inhibition values of Omicron to each variant ----



wilcox.test(corr_df[corr_df$antigen == "Ancestral",]$binding, # 8.446e-09
            corr_df[corr_df$antigen == "Omicron BA.1",]$binding)
wilcox.test(corr_df[corr_df$antigen == "Ancestral",]$inhib2, # 2.78e-09
            corr_df[corr_df$antigen == "Omicron BA.1",]$inhib2)

wilcox.test(corr_df[corr_df$antigen == "Alpha",]$binding, # 2.833e-06
            corr_df[corr_df$antigen == "Omicron BA.1",]$binding)
wilcox.test(corr_df[corr_df$antigen == "Alpha",]$inhib2, # 1.193e-08
            corr_df[corr_df$antigen == "Omicron BA.1",]$inhib2)

wilcox.test(corr_df[corr_df$antigen == "Beta",]$binding, # 3.498e-05
            corr_df[corr_df$antigen == "Omicron BA.1",]$binding)
wilcox.test(corr_df[corr_df$antigen == "Beta",]$inhib2, # 0.006416
            corr_df[corr_df$antigen == "Omicron BA.1",]$inhib2)

wilcox.test(corr_df[corr_df$antigen == "Gamma",]$binding, # 0.0006184
            corr_df[corr_df$antigen == "Omicron BA.1",]$binding)
wilcox.test(corr_df[corr_df$antigen == "Gamma",]$inhib2, # 0.007937
            corr_df[corr_df$antigen == "Omicron BA.1",]$inhib2)

wilcox.test(corr_df[corr_df$antigen == "Delta",]$binding, # 9.853e-06
            corr_df[corr_df$antigen == "Omicron BA.1",]$binding)
wilcox.test(corr_df[corr_df$antigen == "Delta",]$inhib2, # 3.188e-07
            corr_df[corr_df$antigen == "Omicron BA.1",]$inhib2)

# Make some GAM-based predictions
corr_df_factor <- corr_df
corr_df_factor$id <- factor(corr_df_factor$id)

pred.dat <- data.frame(binding = c(73000, 425000, 600000, 975000,
                                   1200000, 1600000, 1700000, 1750000, 
                                   1800000, 2000000, 2500000, 1664817, 1609323),
                       id = factor(1))
pred.dat$binding <- pred.dat$binding * 0.00901

gam_omi <- gam(inhib2 ~ s(binding, bs = "tp") + s(id, bs = "re"),
                   data = corr_df_factor[corr_df_factor$antigen == "Omicron BA.1",],
                   family = "quasibinomial", method = "REML")
preds_omi <- predict(gam_omi,
                 newdata = pred.dat,
                 exclude = "s(id)",
                 se = T, type = 'response')
preds_omi_results <- round(data.frame(at_binding = pred.dat$binding,
                            estimate_inhib = preds_omi$fit,
                            LL = preds_omi$fit - qnorm(0.975)*preds_omi$se.fit,
                            UL = preds_omi$fit + qnorm(0.975)*preds_omi$se.fit), 2)
preds_omi_results$variant <- "Omicron BA.1"
conf_omi <- confidence_band(gam_omi, pred.dat, unconditional = TRUE, exclude = "s(id)")
conf_omi$variant <- "Omicron BA.1"



#### Figure S9 ----



# Binding efficiency and ACE2 inhibition by time
dose.labs <- c("Post-dose 2", "6 months post-dose 2",
               "Post-dose 3", "6 months post-dose 3")
names(dose.labs) <- c("Post-dose 2", "6 months post-dose 2", 
                      "Post-dose 3", "6 months post-dose 3")

# Binding efficiency and ACE2 inhibition by time, with random effects
time_corr_re <- ggplot(corr_df, aes(x = binding, y = inhib2)) +
  geom_hline(yintercept = LOD_inhib/100, linetype = 2) +
  geom_vline(xintercept = LOD_IgG, linetype = 2) +
  geom_point(alpha = 0.4, size = 2.5) + 
  theme_fausto() + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     labels = seq(0, 1, 0.25)*100,
                     limits = c(0,1),
                     expand = c(0.001, 0)) +
  scale_x_continuous(breaks = seq(0, 30000, 5000),
                     labels = c("0", "5,000", "10,000", "15,000", "20,000", "25,000", "30,000"),
                     limits = c(0, 30000,
                                expand = c(0.0001, 0))) +
  coord_cartesian(xlim = c(50000 * 0.00901, NA)) +
  xlab("Anti-S IgG (BAU/ml)") +
  ylab("ACE2 inhibition (%)") + 
  stat_smooth(method = "gam", color = "#6a17de", fill = "#6a17de",
              se = T, aes(group = 1, randomid = factor(id)),
              formula = y ~ s(x, bs = "tp") + s(randomid, bs = "re"),
              method.args = list(family = "quasibinomial", method = "REML")) + 
  facet_wrap(~ time, labeller = labeller(time = dose.labs)) +
  theme(axis.text = element_text(size = 13),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left")))

time_corr_re



#### Compare correlations at 75% at different time points ----



df_time <- data.frame("Above 75%" = c(nrow(corr_df %>%
                                             filter(time == "Post-dose 2") %>%
                                             filter(inhib2 >= 0.75)),
                                      (nrow(corr_df %>%
                                              filter(time == "Post-dose 3") %>%
                                              filter(inhib2 >= 0.75)))),
                      "Below 75%" = c(nrow(corr_df %>%
                                             filter(time == "Post-dose 2") %>%
                                             filter(inhib2 < 0.75)),
                                      (nrow(corr_df %>%
                                              filter(time == "Post-dose 3") %>%
                                              filter(inhib2 < 0.75)))), 
                      row.names = c("Post-dose 2", "Post-dose 3"))
df_time
fisher.test(df_time) # p-value = 3.741e-10

df_time2 <- data.frame("Above 75%" = c(nrow(corr_df %>%
                                              filter(time == "6 months post-dose 2") %>%
                                              filter(inhib2 >= 0.75)),
                                       (nrow(corr_df %>%
                                               filter(time == "6 months post-dose 3") %>%
                                               filter(inhib2 >= 0.75)))),
                       "Below 75%" = c(nrow(corr_df %>%
                                              filter(time == "6 months post-dose 2") %>%
                                              filter(inhib2 < 0.75)),
                                       (nrow(corr_df %>%
                                               filter(time == "6 months post-dose 3") %>%
                                               filter(inhib2 < 0.75)))), 
                       row.names = c("6 months post-dose 2", "6 months post-dose 3"))
df_time2
fisher.test(df_time2) # p-value = 0.0007479

df_time3 <- data.frame("Above 75%" = c(nrow(corr_df %>%
                                              filter(time == "Post-dose 3") %>%
                                              filter(inhib2 >= 0.75)),
                                       (nrow(corr_df %>%
                                               filter(time == "6 months post-dose 3") %>%
                                               filter(inhib2 >= 0.75)))),
                       "Below 75%" = c(nrow(corr_df %>%
                                              filter(time == "Post-dose 3") %>%
                                              filter(inhib2 < 0.75)),
                                       (nrow(corr_df %>%
                                               filter(time == "6 months post-dose 3") %>%
                                               filter(inhib2 < 0.75)))), 
                       row.names = c("Post-dose 3", "6 months post-dose 3"))
df_time3
fisher.test(df_time3) # p-value = 9.034e-13



#### Figure S10 ----



group_corr_re <- ggplot(corr_df, aes(x = binding, y = inhib2)) +
  geom_hline(yintercept = LOD_inhib/100, linetype = 2) +
  geom_vline(xintercept = LOD_IgG, linetype = 2) +
  geom_point(alpha = 0.4, size = 2.5) + 
  theme_fausto() + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     labels = seq(0, 1, 0.25)*100,
                     limits = c(0, 1),
                     expand = c(0.01, 0)) +
  scale_x_continuous(breaks = seq(0, 30000, 5000),
                     labels = c("0", "5,000", "10,000", "15,000", "20,000", "25,000", "30,000"),
                     limits = c(0, 27000,
                                expand = c(0.0001, 0))) +
  coord_cartesian(xlim = c(50000 * 0.00901, 27000)) +
  xlab("Anti-S IgG (BAU/ml)") +
  ylab("ACE2 inhibition (%)") + 
  stat_smooth(method = "gam", color = "#6a17de", fill = "#6a17de",
              se = T, aes(group = 1, randomid = factor(id)),
              formula = y ~ s(x, bs = "tp") + s(randomid, bs = "re"),
              method.args = list(family = "quasibinomial", method = "REML")) + 
  facet_wrap(~ new_group) + 
  theme(axis.text = element_text(size = 13),
        panel.grid.major.x = element_blank(),
        #axis.text.x = element_text(angle = 45/2),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left")))


group_corr_re



#### Figure S8 ----



corr_df_mega <- rbind(corr_df, corr_df)
corr_df_mega$antigen <- c(rep("Overall", nrow(corr_df)), as.character(corr_df$antigen))
corr_df_mega$antigen <- factor(corr_df_mega$antigen,
                               levels = c("Overall",
                                          "Ancestral", "Alpha", "Beta", "Gamma", "Delta",
                                          "Omicron BA.1"))

# Add in a density plot to the variant-specific graphs of binding vs neutralization
# Had to log the x-axis to get reasonable density visualizations...
# need a reasonably regular grid in order to get sensible contour plots
# Did log10 transformation so that 5 (log plot) = 100,000 (MSD reg value) ~= 1 ug/ml on ELISA
corr_df_mega %>%
  filter(antigen != "Overall") %>%
  group_by(antigen) %>%
  summarise(Medianx = median(log10(binding)),
            Mediany = median(inhib2)) -> medians

rcomp_data <- corr_df_mega %>%
  filter(antigen != "Overall") %>%
  mutate(logx = log10(binding))

medians_x <- groupwiseMedian(logx ~ antigen,
                             data       = rcomp_data,
                             conf       = 0.95,
                             R          = 10000,
                             bca        = TRUE,
                             digits     = 7)

medians_y <- groupwiseMedian(inhib2 ~ antigen,
                             data       = rcomp_data,
                             conf       = 0.95,
                             R          = 10000,
                             bca        = TRUE,
                             digits     = 7)

medians_xy <- medians_x %>%
  dplyr::select(antigen, Bca.lower, Bca.upper)
medians_xy <- left_join(medians_xy, medians, by = "antigen")

medians_y %>%
  dplyr::select(antigen, lowery = Bca.lower, uppery = Bca.upper) -> medians_y2
medians_xy <- left_join(medians_xy, medians_y2, by = "antigen")

# Panel plot with densities
logistic_panel_plot <- corr_df_mega %>%
  filter(antigen != "Overall") %>%
  ggplot(aes(x = log10(binding), y = inhib2, group = factor(antigen))) +
  geom_hline(yintercept = LOD_inhib/100, linetype = 2) +
  geom_vline(xintercept = log10(LOD_IgG), linetype = 2) +
  stat_density_2d(geom = "polygon",
                  contour = TRUE,
                  aes(alpha = ..level.., fill = after_stat(level)),
                  bins = 5) + 
  scale_fill_gradient(low = "#e5e5e5", high = "#7f7f7f") +# gray90 to gray50
  geom_point(aes(color = factor(antigen)), alpha = 0.45, size = 2.5) + 
  theme_fausto() + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     labels = seq(0, 1, 0.25)*100,
                     limits = c(0, 1),
                     expand = c(0.02, 0)) +
  scale_x_continuous(breaks = seq(-1, 5, 1),
                     labels = expression("10"^-1, "10"^0, "10"^1, "10"^2, "10"^3,
                                         "10"^4, "10"^5),
                     limits = c(-0.9, 4.5)) +
  xlab("Anti-S IgG (BAU/ml)") +
  ylab("ACE2 inhibition (%)") + 
  facet_wrap(~factor(antigen)) + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(size = 13),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) +
  stat_smooth(method = "gam", size = 2.5,
              se = F, aes(group = factor(antigen), color = factor(antigen),
                          randomid = factor(id)),
              formula = y ~ s(x, bs = "tp") + s(randomid, bs = "re"),
              method.args = list(family = "quasibinomial", method = "REML")) +
  geom_point(data = medians, aes(x = Medianx, y = Mediany), 
             fill = "#C7510F", color = c(rep("black", 2), "white", rep("black", 3)), 
             size = 4, shape = 22, alpha = 1) + 
  scale_color_manual(values = c("Ancestral" = "#7209b7B2",
                                "Alpha" = "#048B90B2",
                                "Beta" = "#000000",
                                "Gamma" = "#DB9303B2",
                                "Delta" = "#0001AAB2",
                                "Omicron BA.1" = "#e00065"))

logistic_panel_plot



#### Figure 2C ----



rainbow_lines_plot_ms <- corr_df_mega %>%
  filter(antigen != "Overall") %>%
  ggplot(aes(x = log10(binding), y = inhib2, group = factor(antigen))) +
  geom_hline(yintercept = LOD_inhib / 100, linetype = 2) + 
  geom_vline(xintercept = log10(LOD_IgG), linetype = 2) +
  geom_point(aes(color = factor(antigen)), alpha = 0.30, size = 1.25, show.legend = F) + 
  theme_fausto() + 
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     labels = seq(0, 1, 0.25)*100,
                     limits = c(0, 1),
                     expand = c(0.02, 0)) +
  scale_x_continuous(breaks = seq(-1, 5, 1),
                     labels = expression("10"^-1, "10"^0, "10"^1, "10"^2, "10"^3,
                                         "10"^4, "10"^5),
                     limits = c(-0.9, 4.9)) +
  xlab("Anti-S IgG (BAU/ml)") +
  ylab("ACE2 inhibition (%)") + 
  stat_smooth(method = "gam", size = 2, alpha = 0.9,
              se = F, aes(group = factor(antigen), color = factor(antigen),
                          randomid = factor(id)),
              formula = y ~ s(x, bs = "tp") + s(randomid, bs = "re"),
              method.args = list(family = "quasibinomial", method = "REML")) +
  geom_point(data = medians, aes(x = Medianx, y = Mediany, fill = factor(antigen)), 
             color = "black", size = 4, shape = 22, alpha = 1) + 
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) + 
  scale_color_manual(name = "SARS-CoV-2 variants",
                     values = c("Ancestral" = "#7209b7B2",
                                "Alpha" = "#048B90B2",
                                "Beta" = "#000000B2",
                                "Gamma" = "#DB9303B2",
                                "Delta" = "#0001AAB2",
                                "Omicron BA.1" = "#e00065")) + 
  scale_fill_manual(name = "SARS-CoV-2 variants",
                    values = c("Ancestral" = "#7209b7B2",
                               "Alpha" = "#048B90B2",
                               "Beta" = "#000000B2",
                               "Gamma" = "#DB9303B2",
                               "Delta" = "#0001AAB2",
                               "Omicron BA.1" = "#e00065")) + 
  theme(axis.text = element_text(size = 14),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left")))

rainbow_lines_plot_ms



#### Figure 2A ---- 



powers_of_ten_labels <- c(expression("10"^2), expression("10"^3), 
                          expression("10"^4), expression("10"^5),
                          expression("10"^6))

powers_of_ten_labels2 <- c(expression("10"^-1), expression("10"^0), 
                           expression("10"^1), expression("10"^2), 
                           expression("10"^3), expression("10"^4))

corr_df %>%
  filter(antigen == "Ancestral") %>% 
  ggplot() + 
  aes(x = new_group, y = binding, shape = time,
      color = new_group, fill = new_group, dodge = time) + 
  
  geom_hline(yintercept = LOD_IgG, linetype = 2) + 

  stat_boxplot(aes(x = new_group, y = binding, color = new_group), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = new_group, y = binding, color = new_group),
              alpha = 0.4, size = 2.5, position = position_dodge(width = 0.75)) + 
  stat_summary(aes(x = new_group, y = binding, color = new_group),
               position = position_dodge(width = 0.75), width = 0.65,
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) + 
  
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3),
                               nrow = 1,
                               direction = "horizontal",
                               title.position = "top")) + #Default alpha for the legend
  guides(shape = guide_legend(direction = "horizontal",
                              title.position = "top")) + 
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000),
                labels = powers_of_ten_labels2) +
  coord_cartesian(ylim = c(17, 3000000) * 0.00901) +
  
  ylab("Anti-S IgG (BAU/ml)") + 
  theme_fausto() +
  theme(legend.position = "bottom", 
        legend.box = "vertical",
        panel.grid.major.x = element_blank()) + 
  
  scale_shape_manual(name = "Time point",
                     breaks = c("Post-dose 2",
                                "6 months post-dose 2",
                                "Post-dose 3",
                                "6 months post-dose 3"),
                     values = c("Post-dose 2" = 25,
                                "6 months post-dose 2" = 21,
                                "Post-dose 3" = 24,
                                "6 months post-dose 3" = 22),
                     labels = c("Post-dose 2" = "Post-dose 2",
                                "6 months post-dose 2" = "6 months\npost-dose 2",
                                "Post-dose 3" = "Post-dose 3",
                                "6 months post-dose 3" = "6 months\npost-dose 3")) +
  scale_color_manual(name = "Subgroup",
                     breaks = c("Healthy volunteer",
                                "Antibody deficiency",
                                "PIRD",
                                "Combined immunodeficiency",
                                "Other IEI",
                                "Other immune disorder"),
                     values = c("Healthy volunteer" = "black",
                                "Antibody deficiency" = "#7209b7",
                                "PIRD" = "#048B90",
                                "Combined immunodeficiency" = "#CF9A03",
                                "Other IEI" = "#D10106",
                                "Other immune disorder" = "#0006F2"),
                     labels = c("Healthy volunteer" = "Healthy\nvolunteer",
                                "Antibody deficiency" = "Antibody\ndeficiency",
                                "PIRD" = "PIRD",
                                "Combined immunodeficiency" = "Combined\nimmunodeficiency",
                                "Other IEI" = "Other\nIEI",
                                "Other immune disorder" = "Other immune\ndisorder")) +
  scale_fill_manual(name = "Immunological subgroup",
                    breaks = c("Healthy volunteer",
                               "Antibody deficiency",
                               "PIRD",
                               "Combined immunodeficiency",
                               "Other IEI",
                               "Other immune disorder"),
                    values = c("Healthy volunteer" = "black",
                               "Antibody deficiency" = "#7209b7",
                               "PIRD" = "#048B90",
                               "Combined immunodeficiency" = "#CF9A03",
                               "Other IEI" = "#D10106",
                               "Other immune disorder" = "#0006F2")) +
  
  guides(fill = "none") +  # Remove redundant color/fill legends
  theme(axis.text.x = element_blank(), #remove labels &  ticks
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) -> temporal_MESO_IgG_wuhan_ms

temporal_MESO_IgG_wuhan_ms



#### Figure S4A ---- 



corr_df %>%
  filter(antigen == "Omicron BA.1") %>% 
  ggplot() + 
  aes(x = new_group, y = binding, shape = time,
      color = new_group, fill = new_group, dodge = time) + 
  
  geom_hline(yintercept = LOD_IgG, linetype = 2) + 

  stat_boxplot(aes(x = new_group, y = binding, color = new_group), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = new_group, y = binding, color = new_group),
              alpha = 0.4, size = 2.5, position = position_dodge(width = 0.75)) + 
  stat_summary(aes(x = new_group, y = binding, color = new_group),
               position = position_dodge(width = 0.75), width = 0.65,
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) + 
  
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3),
                               nrow = 1,
                               direction = "horizontal",
                               title.position = "top")) + #Default alpha for the legend
  guides(shape = guide_legend(direction = "horizontal",
                              title.position = "top")) + 
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000),
                labels = powers_of_ten_labels2) +
  coord_cartesian(ylim = c(17, 3000000) * 0.00901) +
  
  ylab("Anti-S IgG (BAU/ml)") + 
  theme_fausto() +
  theme(legend.position = "bottom", 
        legend.box = "vertical",
        panel.grid.major.x = element_blank()) + 
  
  scale_shape_manual(name = "Time point",
                     breaks = c("Post-dose 2",
                                "6 months post-dose 2",
                                "Post-dose 3",
                                "6 months post-dose 3"),
                     values = c("Post-dose 2" = 25,
                                "6 months post-dose 2" = 21,
                                "Post-dose 3" = 24,
                                "6 months post-dose 3" = 22),
                     labels = c("Post-dose 2" = "Post-dose 2",
                                "6 months post-dose 2" = "6 months\npost-dose 2",
                                "Post-dose 3" = "Post-dose 3",
                                "6 months post-dose 3" = "6 months\npost-dose 3")) +
  scale_color_manual(name = "Subgroup",
                     breaks = c("Healthy volunteer",
                                "Antibody deficiency",
                                "PIRD",
                                "Combined immunodeficiency",
                                "Other IEI",
                                "Other immune disorder"),
                     values = c("Healthy volunteer" = "black",
                                "Antibody deficiency" = "#7209b7",
                                "PIRD" = "#048B90",
                                "Combined immunodeficiency" = "#CF9A03",
                                "Other IEI" = "#D10106",
                                "Other immune disorder" = "#0006F2"),
                     labels = c("Healthy volunteer" = "Healthy\nvolunteer",
                                "Antibody deficiency" = "Antibody\ndeficiency",
                                "PIRD" = "PIRD",
                                "Combined immunodeficiency" = "Combined\nimmunodeficiency",
                                "Other IEI" = "Other\nIEI",
                                "Other immune disorder" = "Other immune\ndisorder")) +
  scale_fill_manual(name = "Immunological subgroup",
                    breaks = c("Healthy volunteer",
                               "Antibody deficiency",
                               "PIRD",
                               "Combined immunodeficiency",
                               "Other IEI",
                               "Other immune disorder"),
                    values = c("Healthy volunteer" = "black",
                               "Antibody deficiency" = "#7209b7",
                               "PIRD" = "#048B90",
                               "Combined immunodeficiency" = "#CF9A03",
                               "Other IEI" = "#D10106",
                               "Other immune disorder" = "#0006F2")) +
  
  guides(fill = "none") +  # Remove redundant color/fill legends
  theme(axis.text.x = element_blank(), #remove labels &  ticks
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) -> temporal_MESO_IgG_omicron_ms

temporal_MESO_IgG_omicron_ms



#### Figure 2B ----


corr_df %>%
  filter(antigen == "Ancestral") %>% 
  ggplot() + 
  aes(x = new_group, y = inhib, shape = time,
      color = new_group, fill = new_group, dodge = time) + 
  
  geom_hline(yintercept = LOD_inhib, linetype = 2) + 

  stat_boxplot(aes(x = new_group, y = inhib, color = new_group), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = new_group, y = inhib, color = new_group),
              alpha = 0.4, size = 2.5, position = position_dodge(width = 0.75)) + 
  stat_summary(aes(x = new_group, y = inhib, color = new_group),
               position = position_dodge(width = 0.75), width = 0.65,
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) +
  
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3),
                               nrow = 1,
                               direction = "horizontal",
                               title.position = "top")) + #Default alpha for the legend
  guides(shape = guide_legend(direction = "horizontal",
                              title.position = "top")) + 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100),
                     labels = c("0", "20", "40", "60", "80", "100"),
                     expand = c(0.02, 0)) +
  coord_cartesian(ylim = c(0, 100)) +
  
  ylab("ACE2 inhibition (%)") + 
  theme_fausto() +
  theme(legend.position = 'bottom', 
        legend.box = "vertical",
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left")),
        axis.text.x = element_blank(), #remove labels &  ticks
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) + 
  
  scale_shape_manual(name = "Time point",
                     breaks = c("Post-dose 2",
                                "6 months post-dose 2",
                                "Post-dose 3",
                                "6 months post-dose 3"),
                     values = c("Post-dose 2" = 25,
                                "6 months post-dose 2" = 21,
                                "Post-dose 3" = 24,
                                "6 months post-dose 3" = 22),
                     labels = c("Post-dose 2" = "Post-dose 2",
                                "6 months post-dose 2" = "6 months\npost-dose 2",
                                "Post-dose 3" = "Post-dose 3",
                                "6 months post-dose 3" = "6 months\npost-dose 3")) +
  scale_color_manual(name = "Subgroup",
                     breaks = c("Healthy volunteer",
                                "Antibody deficiency",
                                "PIRD",
                                "Combined immunodeficiency",
                                "Other IEI",
                                "Other immune disorder"),
                     values = c("Healthy volunteer" = "black",
                                "Antibody deficiency" = "#7209b7",
                                "PIRD" = "#048B90",
                                "Combined immunodeficiency" = "#CF9A03",
                                "Other IEI" = "#D10106",
                                "Other immune disorder" = "#0006F2"),
                     labels = c("Healthy volunteer" = "Healthy\nvolunteer",
                                "Antibody deficiency" = "Antibody\ndeficiency",
                                "PIRD" = "PIRD",
                                "Combined immunodeficiency" = "Combined\nimmunodeficiency",
                                "Other IEI" = "Other\nIEI",
                                "Other immune disorder" = "Other immune\ndisorder")) +
  scale_fill_manual(name = "Immunological subgroup",
                    breaks = c("Healthy volunteer",
                               "Antibody deficiency",
                               "PIRD",
                               "Combined immunodeficiency",
                               "Other IEI",
                               "Other immune disorder"),
                    values = c("Healthy volunteer" = "black",
                               "Antibody deficiency" = "#7209b7",
                               "PIRD" = "#048B90",
                               "Combined immunodeficiency" = "#CF9A03",
                               "Other IEI" = "#D10106",
                               "Other immune disorder" = "#0006F2")) +
  
  guides(fill = "none") -> temporal_MESO_inhib_wuhan_ms 

temporal_MESO_inhib_wuhan_ms



#### Figure 2 ----



layout <- "
AAAAAAAA
#BBBBBB#
"
x <- ((temporal_MESO_IgG_wuhan_ms + temporal_MESO_inhib_wuhan_ms) & 
        guides(shape = guide_legend(override.aes = list(fill = "black", size = 2.5, alpha = 1),
                                    direction = "horizontal", title.position = "top")) &
        scale_x_discrete(expand = c(-0.7, 0)) &
        theme(legend.position = 'bottom', 
              legend.box = "horizontal",
              legend.spacing.y = unit(0.05, "cm"),
              legend.spacing.x = unit(0.15, "cm"),
              legend.margin = margin(0, 0, 0, 0),
              legend.box.margin = margin(-10, -10, -5, -10),
              legend.text = element_text(size = 12.5),
              legend.title = element_text(size = 12, face = "bold"),
              axis.text = element_text(size = 12.5),
              axis.title = element_text(size = 12.5, face = "bold"))) +
  plot_layout(guides = 'collect') + 
  theme(legend.margin = unit(x = c(0, 25, 0, 0), units = "cm"))

y <- (rainbow_lines_plot_ms & 
        theme(legend.position = 'right', 
              legend.box = "horizontal",
              legend.key.width = unit(2, "line"),
              panel.grid.major.x = element_blank(),
              plot.margin = margin(0.4, 0, -0.75, 0, "cm"),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12.5, face = "bold"),
              axis.text = element_text(size = 12.5),
              axis.title = element_text(size = 12.5, face = "bold")) & 
        guides(color = guide_legend(ncol = 1, byrow = TRUE))) 

fig2 <- (x / y) + 
  plot_layout(design = layout) + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(face = "bold", hjust = 0, vjust = 0.75))

fig2
save_plot("fig2.pdf", fig2, base_height = 6.5, base_width = 13.25, dpi = 300) # Excellent



#### Figure S4B ----



corr_df %>%
  filter(antigen == "Omicron BA.1") %>% 
  ggplot() + 
  aes(x = new_group, y = inhib, shape = time,
      color = new_group, fill = new_group, dodge = time) + 
  
  geom_hline(yintercept = LOD_inhib, linetype = 2) + 

  stat_boxplot(aes(x = new_group, y = inhib, color = new_group), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = new_group, y = inhib, color = new_group),
              alpha = 0.4, size = 2.5, position = position_dodge(width = 0.75)) + 
  stat_summary(aes(x = new_group, y = inhib, color = new_group),
               position = position_dodge(width = 0.75), width = 0.65,
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) +
  
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3),
                               nrow = 1,
                               direction = "horizontal",
                               title.position = "top")) + #Default alpha for the legend
  guides(shape = guide_legend(direction = "horizontal",
                              title.position = "top")) + 
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100),
                     labels = c("0", "20", "40", "60", "80", "100"),
                     expand = c(0.02, 0)) +
  coord_cartesian(ylim = c(0, 100)) +
  
  ylab("ACE2 inhibition (%)") + 
  theme_fausto() +
  theme(legend.position = 'bottom', 
        legend.box = "vertical",
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left")),
        axis.text.x = element_blank(), #remove labels &  ticks
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) + 
  
  scale_shape_manual(name = "Time point",
                     breaks = c("Post-dose 2",
                                "6 months post-dose 2",
                                "Post-dose 3",
                                "6 months post-dose 3"),
                     values = c("Post-dose 2" = 25,
                                "6 months post-dose 2" = 21,
                                "Post-dose 3" = 24,
                                "6 months post-dose 3" = 22),
                     labels = c("Post-dose 2" = "Post-dose 2",
                                "6 months post-dose 2" = "6 months\npost-dose 2",
                                "Post-dose 3" = "Post-dose 3",
                                "6 months post-dose 3" = "6 months\npost-dose 3")) +
  scale_color_manual(name = "Subgroup",
                     breaks = c("Healthy volunteer",
                                "Antibody deficiency",
                                "PIRD",
                                "Combined immunodeficiency",
                                "Other IEI",
                                "Other immune disorder"),
                     values = c("Healthy volunteer" = "black",
                                "Antibody deficiency" = "#7209b7",
                                "PIRD" = "#048B90",
                                "Combined immunodeficiency" = "#CF9A03",
                                "Other IEI" = "#D10106",
                                "Other immune disorder" = "#0006F2"),
                     labels = c("Healthy volunteer" = "Healthy\nvolunteer",
                                "Antibody deficiency" = "Antibody\ndeficiency",
                                "PIRD" = "PIRD",
                                "Combined immunodeficiency" = "Combined\nimmunodeficiency",
                                "Other IEI" = "Other\nIEI",
                                "Other immune disorder" = "Other immune\ndisorder")) +
  scale_fill_manual(name = "Subgroup",
                    breaks = c("Healthy volunteer",
                               "Antibody deficiency",
                               "PIRD",
                               "Combined immunodeficiency",
                               "Other IEI",
                               "Other immune disorder"),
                    values = c("Healthy volunteer" = "black",
                               "Antibody deficiency" = "#7209b7",
                               "PIRD" = "#048B90",
                               "Combined immunodeficiency" = "#CF9A03",
                               "Other IEI" = "#D10106",
                               "Other immune disorder" = "#0006F2")) +
  
  guides(fill = "none") -> temporal_MESO_inhib_omicron_ms 

temporal_MESO_inhib_omicron_ms



#### Fig S4 ----



fig_s4 <- ((temporal_MESO_IgG_omicron_ms + temporal_MESO_inhib_omicron_ms) & 
             guides(shape = guide_legend(override.aes = list(fill = "black", size = 2, alpha = 1),
                                         direction = "horizontal", title.position = "top")) &
             scale_x_discrete(expand = c(-0.7, 0)) &
             theme(legend.position = 'bottom', 
                   legend.box = "horizontal",
                   legend.spacing.y = unit(0.05, "cm"),
                   legend.spacing.x = unit(0.15, "cm"),
                   legend.margin = margin(0, 0, 0, 0),
                   legend.box.margin = margin(-10, -10, -5, -10),
                   legend.text = element_text(size = 12.5),
                   legend.title = element_text(size = 12, face = "bold"),
                   axis.text = element_text(size = 12.5),
                   axis.title = element_text(size = 12.5, face = "bold"))) +
  plot_layout(guides = 'collect') 

fig_s4 + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position = c(0, 1), 
        plot.tag = element_text(face = "bold", hjust = 0, vjust = 0.75)) 



#### Fig S18 ----



# Load the ELISA data
load("C:/Users/bustosfa/Rotation 2/Data/Final_ELISA_data.RData")

# Comparing ELISA IgG to MSD
powers_of_ten_labels3 <- c(expression("10"^-3), expression("10"^-2), 
                           expression("10"^-1), expression("10"^0),
                           expression("10"^1), expression("10"^2), 
                           expression("10"^3), expression("10"^4), 
                           expression("10"^5), expression("10"^6))

elisa <- df %>%
  filter(analyte == "Spike IgG" & timepoint %in% as.character(unique(all_IgG2$time)))  %>%
  dplyr::select(id = `Subject ID`, group = new_group, value = conc.ug_ml, time = timepoint) %>%
  mutate(source = "ELISA") %>%
  mutate(id2 = paste0(id, "_", time))

msd <- all_IgG2
msd <- msd %>%
  mutate(id2 = paste0(id, "_", time)) %>%
  filter(id2 %in% elisa$id2 & antigen == "Ancestral") %>%
  mutate(source = "MSD") %>%
  dplyr::select(id, group = new_group, value = binding, time, source, id2)

elisa <- elisa %>%
  filter(id2 %in% msd$id2) %>%
  dplyr::select(-id2)

msd <- msd %>%
  dplyr::select(-id2)

elisa$time <- factor(elisa$time, levels = c("Post-dose 2", "6 months post-dose 2",
                                            "Post-dose 3", "6 months post-dose 3"))

msd$time <- factor(msd$time, levels = c("Post-dose 2", "6 months post-dose 2",
                                        "Post-dose 3", "6 months post-dose 3"))

together <- rbind(elisa, msd) %>%
  arrange(id, source, time)

elisa_plot <- together %>% filter(group == "Healthy volunteer" & source == "ELISA") %>%
  ggplot() + 
  aes(y = value, shape = time, dodge = time) + 
  
  stat_boxplot(aes(x = time, y = value), color = "#EC449B", fill = "#EC449BB2",
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = time, y = value), color = "#EC449B", fill = "#EC449BB2",
              alpha = 0.4, size = 3.5, position = position_dodge(width = 0.75)) + 
  stat_summary(aes(x = time, y = value), color = "#EC449B", fill = "#EC449BB2",
               position = position_dodge(width = 0.75), width = 0.65,
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) + 
  
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3),
                               nrow = 1,
                               direction = "horizontal",
                               title.position = "top")) + #Default alpha for the legend
  guides(shape = guide_legend(direction = "horizontal",
                              title.position = "top")) + 
  
  scale_y_log10(breaks = c(0.0001, 0.001, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000),
                labels = powers_of_ten_labels3) +
  
  ylab("Anti-S IgG (ug/ml)") + 
  theme_fausto() +
  theme(legend.position = "bottom", 
        legend.box = "vertical",
        panel.grid.major.x = element_blank()) + 
  
  scale_shape_manual(name = "Time point",
                     breaks = c("Post-dose 2",
                                "6 months post-dose 2",
                                "Post-dose 3",
                                "6 months post-dose 3"),
                     values = c("Post-dose 2" = 25,
                                "6 months post-dose 2" = 21,
                                "Post-dose 3" = 24,
                                "6 months post-dose 3" = 22),
                     labels = c("Post-dose 2" = "Post-dose 2",
                                "6 months post-dose 2" = "6 months\npost-dose 2",
                                "Post-dose 3" = "Post-dose 3",
                                "6 months post-dose 3" = "6 months\npost-dose 3")) +
  
  guides(fill = "none") +  # Remove redundant color/fill legends
  theme(axis.text.x = element_blank(), #remove labels &  ticks
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) + 
  ggtitle("ELISA")


msd_plot <- together %>% filter(group == "Healthy volunteer" & source == "MSD") %>%
  ggplot() + 
  aes(y = value, shape = time, dodge = time) + 
  
  stat_boxplot(aes(x = time, y = value), color = "#00008B", fill = "#00008BB2",
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = time, y = value), color = "#00008B", fill = "#00008BB2",
              alpha = 0.4, size = 3.5, position = position_dodge(width = 0.75)) + 
  stat_summary(aes(x = time, y = value), color = "#00008B", fill = "#00008BB2",
               position = position_dodge(width = 0.75), width = 0.65,
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) + 
  
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3),
                               nrow = 1,
                               direction = "horizontal",
                               title.position = "top")) + #Default alpha for the legend
  guides(shape = guide_legend(direction = "horizontal",
                              title.position = "top")) + 
  
  scale_y_log10(breaks = c(0.0001, 0.001, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000),
                labels = powers_of_ten_labels3) +
  
  ylab("Anti-S IgG (BAU/ml)") + 
  theme_fausto() +
  theme(legend.position = "bottom", 
        legend.box = "vertical",
        panel.grid.major.x = element_blank()) + 
  
  scale_shape_manual(name = "Time point",
                     breaks = c("Post-dose 2",
                                "6 months post-dose 2",
                                "Post-dose 3",
                                "6 months post-dose 3"),
                     values = c("Post-dose 2" = 25,
                                "6 months post-dose 2" = 21,
                                "Post-dose 3" = 24,
                                "6 months post-dose 3" = 22),
                     labels = c("Post-dose 2" = "Post-dose 2",
                                "6 months post-dose 2" = "6 months\npost-dose 2",
                                "Post-dose 3" = "Post-dose 3",
                                "6 months post-dose 3" = "6 months\npost-dose 3")) +
  
  guides(fill = "none") +  # Remove redundant color/fill legends
  theme(axis.text.x = element_blank(), #remove labels &  ticks
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) + 
  ggtitle("Electrochemiluminescence")

together_plot <- ((elisa_plot + msd_plot) & 
                    guides(shape = guide_legend(override.aes = list(fill = "black", size = 4, 
                                                                    alpha = 1, color = "black"),
                                                direction = "horizontal", title.position = "top")) &
                    scale_x_discrete(expand = c(-0.7, 0)) &
                    theme(legend.position = 'bottom', 
                          legend.box = "horizontal",
                          legend.spacing.y = unit(0.05, "cm"),
                          legend.spacing.x = unit(0.15, "cm"),
                          legend.margin = margin(0, 0, 0, 0),
                          legend.box.margin = margin(-10, -10, -5, -10),
                          legend.text = element_text(size = 12.5),
                          legend.title = element_text(size = 12, face = "bold"),
                          axis.text = element_text(size = 13.5),
                          axis.title = element_text(size = 12.5, face = "bold"))) +
  plot_layout(guides = 'collect') 

together_plot
