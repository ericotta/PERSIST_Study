######################################################################
#                                                                    #
#      Analysis of anti-S IgG ELISA data for PERSIST 1 paper         #
#                                                                    #
######################################################################



#### Data cleaning ----



# Load packages
library(here)
library(dplyr)
library(tidyverse)
library(tableone)
library(kableExtra)
library(mice) 
library(lubridate)
library(RColorBrewer)
library(ggpubr)
library(asht)
library(scales)
library(logistf)
library(RColorBrewer)
library(lubridate)
library(ggthemes)
library(patchwork)
library(grid)
library(cowplot)

# Load the data...might need to change this as needed
load("C:/Users/bustosfa/Rotation 2/Data/PERSIST_data_v7_2022-12-09.RData")
variants <- read.csv(here("C:/Users/bustosfa/Rotation 2/Data", "variants.csv"))

# Custom FB functions
`%!in%` <- Negate(`%in%`)

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

# Start cleaning the data
data7 %<>%
  group_by(ks_patid) %>%
  mutate(N_pos_before_sample = case_when(date_first_N_pos < Collection.date ~ 1,
                                         TRUE ~ 0)) %>%   
  mutate(N_pos = case_when(nucleocapsidOD >= 0.9 ~ 1,
                           TRUE ~ 0),
         conc_cutoff = case_when(conc.ug_ml >= 0.0009426984 & ks_cutoff == 1 ~ "Positive",
                                 conc.ug_ml >= 0.0009426984 & ks_cutoff == 0 ~ "Indeterminate",
                                 TRUE ~ "Negative"))

# FB: Combine 6 months post-dose 2 and pre-dose 3
check <- data7 %>%
  dplyr::select(`Subject ID`, ks_patid, timepoint) %>%
  filter(timepoint %in% c("6 months post-dose 2", "Pre-dose 3")) %>%
  distinct() %>%
  arrange(`Subject ID`) %>%
  group_by(`Subject ID`) %>%
  filter(n() > 1)
check <- as.vector(unique(check$`Subject ID`)) #ID with both timepoints

data7$timepoint2 <- data7$timepoint
for(i in 1:nrow(data7)){
  if(data7$`Subject ID`[i] %!in% check & data7$timepoint[i] == "Pre-dose 3"){
    data7$timepoint2[i] <- "6 months post-dose 2"
  }
  if(data7$`Subject ID`[i] %in% check & data7$timepoint[i] == "Pre-dose 3"){
    data7$timepoint2[i] <- "Remove"
  }
}
check <- data7 %>%
  dplyr::select(`Subject ID`, timepoint, timepoint2)

data7 <- data7[data7$timepoint2 != "Remove",]

data7$timepoint <- NULL
data7$timepoint <- data7$timepoint2
data7$timepoint2 <- NULL
rm(check)

# Continue with MZ's data cleaning
df <- data7 %>% 
  mutate(
    control = case_when(
      new_group == "Healthy volunteer" ~ "Healthy volunteer",
      TRUE ~ "PID"), 
    nucleocapsid_ks_cutoff = as.factor(
      case_when(
        nucleocapsidOD >= 0.9 ~ "Positive",
        TRUE ~ "Negative")), 
    timepoint = case_when(
      timepoint == "6 months post-dose 2" ~ "6 months post-dose 2",
      timepoint == "Post-dose 1" ~ "Post-dose 1",
      timepoint == "Post-dose 2" ~ "Post-dose 2",
      timepoint == "Post-dose 3" ~ "Post-dose 3",
      timepoint == "Post-dose 4" ~ "Post-dose 4",
      timepoint == "Post-dose 5" ~ "Post-dose 5",
      timepoint == "Baseline" ~ "Baseline",
      timepoint == "6 months post-dose 3" ~ "6 months post-dose 3",
      timepoint == "6 months post-dose 4" ~ "6 months post-dose 4", 
      timepoint == "Pre-dose 4" ~"Pre-dose 4"
    )
  ) %>%
  subset(analyte == "Spike IgG") %>% 
  filter(timepoint %in% c("Baseline", "Post-dose 1", "Post-dose 2", "Post-dose 3",
                          "6 months post-dose 2", "6 months post-dose 3"))

df$new_group <- factor(df$new_group, levels = c("Healthy volunteer",
                                                "Antibody deficiency",
                                                "PIRD",
                                                "Combined immunodeficiency",
                                                "Other IEI",
                                                "Other immune disorder"))

df$control <- factor(df$control, levels = c("Healthy volunteer", "PID"))

df$timepoint <- factor(df$timepoint, levels = c("Baseline", "Post-dose 1", 
                                                "Post-dose 2", "6 months post-dose 2", 
                                                "Post-dose 3", "6 months post-dose 3"))

# FB: Change the data's lower limit from 0 to LOD/sqrt(2) and ungroup the data
table(df$conc.ug_ml >= 0.0009426984)
df[df$conc.ug_ml < 0.0009426984,]$conc.ug_ml <- 0.0009426984 / sqrt(2)

df <- df %>%
  ungroup()



#### Figure 1a - IgG titers for HV and all IDP ----



# Make labels for the figure
time.labs <- c("Baseline", "Post-dose 1",
               "Post-dose 2", "6 months\npost-dose 2",
               "Post-dose 3", "6 months\npost-dose 3")
names(time.labs) <- c("Baseline", "Post-dose 1",
                      "Post-dose 2", "6 months post-dose 2", 
                      "Post-dose 3", "6 months post-dose 3")

df %>%
  ggplot(aes(x = control, y = conc.ug_ml)) +
  geom_hline(yintercept = 0.0009426984, linetype = 2) +
  stat_boxplot(aes(x = control, y = conc.ug_ml, color = control), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = control, y = conc.ug_ml, color = control),
              width = 0.2, alpha = 0.4, size = 2.5) + 
  stat_summary(aes(x = control, y = conc.ug_ml, color = control),
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = seq(0, 40, 10)) +
  stat_compare_means(method = "wilcox.test", # To mirror option,s for fig 1b
                     label = "p.signif",
                     label.y.npc = 0.99,
                     label.x.npc = 0.5,
                     hide.ns = T, 
                     size = 5.5) +
  coord_cartesian(ylim = c(0.001, 45)) + 
  facet_wrap(~timepoint, ncol = 6) + 
  ylab("Anti-S IgG (ug/ml)") +
  scale_color_manual(name = "Subgroup",
                     values = c("Healthy volunteer" = "black",
                                "PID" = "#d60270"),
                     labels = c("Healthy volunteer" = "Healthy volunteer",
                                "PID" = "All IDP")) + 
  theme_fausto() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4.5), nrow = 1)) -> fig_1a

fig_1a



#### Figure 1b - IgG titers for HV and IDP subgroups ----



df %>%
  ggplot(aes(x = new_group, y = conc.ug_ml)) +
  geom_hline(yintercept = 0.0009426984, linetype = 2) +
  stat_boxplot(aes(x = new_group, y = conc.ug_ml, color = new_group), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = new_group, y = conc.ug_ml, color = new_group),
              width = 0.2, alpha = 0.4, size = 2.5) + 
  stat_summary(aes(x = new_group, y = conc.ug_ml, color = new_group),
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = seq(0, 40, 10)) +
  coord_cartesian(ylim = c(0, 45)) + 
  stat_compare_means(ref.group = "Healthy volunteer", 
                     method = "wilcox.test",
                     label = "p.signif",
                     label.y.npc = 0.99,
                     hide.ns = T, 
                     size = 5.5) +
  facet_wrap(~timepoint, ncol = 6) +
  ylab("Anti-S IgG (ug/ml)") +
  scale_color_manual(name = "Subgroup",
                     values = c("Healthy volunteer" = "black",
                                "Antibody deficiency" = "#7209b7",
                                "PIRD" = "#048B90",
                                "Combined immunodeficiency" = "#CF9A03",
                                "Other IEI" = "#D10106",
                                "Other immune disorder" = "#0006F2")) + 
  theme_fausto() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4.5), nrow = 1)) -> fig_1b

fig_1b



#### Figure 1 - IgG titers ----



(fig_1a / fig_1b) & 
  theme(legend.spacing.y = unit(0.05, "cm"),
        legend.spacing.x = unit(0.15, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, -10, -5, -10),
        legend.text = element_text(size = 12.5),
        legend.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 13, face = "bold"),
        plot.tag.position = c(0, 1),
        plot.tag = element_text(face = "bold", hjust = 0, vjust = 0.75),
        strip.text.x = element_text(size = 14)) & 
  plot_annotation(tag_levels = 'A') -> fig_1

fig_1 

save_plot("fig1.pdf", fig_1, base_height = 6.5, base_width = 13.25, dpi = 300) # Excellent



#### Figure 3a - Sensitivity analysis: IgRT ----



df1 <- df %>% 
  mutate(igrt = case_when(IgRT == "Yes" ~ "Yes",
                          TRUE ~ "No"))

df1$igrt <- factor(df1$igrt, levels = c("Yes","No"))


df1 %>% 
  filter(control == "PID") %>% 
  ggplot(aes(x = igrt, y = conc.ug_ml)) +
  geom_hline(yintercept = 0.0009426984, linetype = 2) +
  stat_boxplot(aes(x = igrt, y = conc.ug_ml, color = igrt), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = igrt, y = conc.ug_ml, color = igrt),
              width = 0.2, alpha = 0.4, size = 3.5) + 
  stat_summary(aes(x = igrt, y = conc.ug_ml, color = igrt),
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) + 
  scale_color_manual(values = c("Yes" = "royal blue", 
                                "No" = "salmon")) + 
  facet_wrap(~timepoint, ncol = 6) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     label.y.npc = 1,
                     label.x.npc = 0.5,
                     hide.ns = T, size = 4.75) +
  ylab("Anti-S IgG (ug/ml)") +
  theme_fausto() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4.5))) -> fig_3a

fig_3a



#### Figure 3b - Sensitivity analysis: Immunosuppressants ----



df1 <- df1 %>% 
  mutate(immunosuppressants = replace_na(immunosuppressants, "No"))

df1$immunosuppressants = factor(df1$immunosuppressants, levels = c("Yes","No"))

df1 %>%
  filter(control == "PID") %>%
  ggplot(aes(x = immunosuppressants, y = conc.ug_ml)) +
  geom_hline(yintercept = 0.0009426984, linetype = 2) +
  stat_boxplot(aes(x = immunosuppressants, y = conc.ug_ml, color = immunosuppressants), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = immunosuppressants, y = conc.ug_ml, color = immunosuppressants),
              width = 0.2, alpha = 0.4, size = 3.5) + 
  stat_summary(aes(x = immunosuppressants, y = conc.ug_ml, color = immunosuppressants),
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) + 
  scale_color_manual(values = c("Yes" = "royal blue", 
                                "No" = "salmon")) + 
  facet_wrap(~timepoint, ncol = 6) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  ylab("Anti-S IgG (ug/ml)") +
  theme_fausto() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4.5)))  -> fig_3b

fig_3b


#### Figure 3c - Sensitivity analysis: Evusheld use ----



df1 <- df1 %>%
  mutate(evu_plot = case_when(`Subject ID` %in% c("E0026", "E0079",
                                                  "E0101", "E0062",
                                                  "E0181", "E0201",
                                                  "E0090") ~ "Yes",
                              TRUE ~ "No"))

df1$evu_plot = factor(df1$evu_plot, levels = c("Yes","No"))


df1 %>% 
  dplyr::select(`Subject ID`, Any_evusheld, date_evusheld, Collection.date, 
                evusheld_before_sample, timepoint, conc.ug_ml, new_group, evu_plot) %>%
  subset(new_group != "Healthy volunteer") %>%
  subset(date_evusheld < ymd("2022-04-29") | is.na(date_evusheld) == TRUE) -> df3c 


df3c %>%                      
  ggplot(aes(x = evu_plot, y = conc.ug_ml)) +
  geom_hline(yintercept = 0.0009426984, linetype = 2) +
  stat_boxplot(aes(x = evu_plot, y = conc.ug_ml, color = evu_plot), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = evu_plot, y = conc.ug_ml, color = evu_plot),
              width = 0.2, alpha = 0.4, size = 3.5) + 
  stat_summary(aes(x = evu_plot, y = conc.ug_ml, color = evu_plot),
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) +   
  scale_color_manual(values = c("Yes" = "royal blue", 
                                "No" = "salmon")) + 
  facet_wrap(~timepoint, ncol = 6) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  ylab("Anti-S IgG (ug/ml)") +
  theme_fausto() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4.5))) -> fig_3c

fig_3c



#### Figure 3 - Sensitivity analysis ---



text1 <- data.frame(igrt = 
                      factor("No", levels = c("Yes", "No")), 
                    conc.ug_ml = 47, lab = "Text",
                    timepoint = factor("Baseline", 
                                       levels = c("Baseline", "Post-dose 1",
                                                  "Post-dose 2", "6 months post-dose 2",
                                                  "Post-dose 3", "6 months post-dose 3")))

text2 <- data.frame(immunosuppressants = 
                      factor("No", levels = c("Yes", "No")), 
                    conc.ug_ml = 36, lab = "Text",
                    timepoint = factor("Baseline", 
                                       levels = c("Baseline", "Post-dose 1",
                                                  "Post-dose 2", "6 months post-dose 2",
                                                  "Post-dose 3", "6 months post-dose 3")))

text3 <- data.frame(evu_plot = 
                      factor("No", levels = c("Yes", "No")), 
                    conc.ug_ml = 47, lab = "Text",
                    timepoint = factor("Baseline", 
                                       levels = c("Baseline", "Post-dose 1",
                                                  "Post-dose 2", "6 months post-dose 2",
                                                  "Post-dose 3", "6 months post-dose 3")))

((fig_3a + scale_color_manual(values = c("Yes" = "royal blue", "No" = "salmon")) + 
    theme(legend.position = "none",
          strip.text.x = element_text(size = 14)) + 
    geom_text(data = text1, label = "IgRT use",
              hjust = 1.6, fontface = "bold", color = "#e00065", size = 5)) / 
    (fig_3b + scale_color_manual(values = c("Yes" = "royal blue", "No" = "salmon")) + 
       theme(legend.position = "none",
             strip.background = element_blank(),
             strip.text.x = element_blank()) + 
       geom_text(data = text2, label = "Immunosup- \npressant use",
                 hjust = 1.055, fontface = "bold", color = "#e00065", size = 5)) / 
    (fig_3c + scale_color_manual(values = c("Yes" = "royal blue", "No" = "salmon")) + 
       theme(strip.background = element_blank(),
             strip.text.x = element_blank()) + 
       geom_text(data = text3, label = "AZD7442 use",
                 hjust = 1.0525, fontface = "bold", color = "#e00065", size = 5))) & 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = seq(0, 40, 10),
                     labels = seq(0, 40, 10)) &
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     label.y.npc = 0.94,
                     label.x.npc = 0.5,
                     hide.ns = T, 
                     size = 7) &
  theme(legend.spacing.y = unit(0.05, "cm"),
        legend.spacing.x = unit(0.15, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, -10, -5, -10),
        legend.text = element_text(size = 12.5),
        legend.title = element_blank(),
        axis.text = element_text(size = 12.5),
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.minor.y = element_blank(),
        plot.tag.position = c(0, 1),
        plot.tag = element_text(face = "bold", hjust = 0, vjust = 0.5)) & 
  coord_cartesian(ylim = c(0, 50)) &
  plot_annotation(tag_levels = 'A') -> fig_3

fig_3

save_plot("fig3.pdf", fig_3, base_height = 6.5, base_width = 13.25, dpi = 300) # Excellent



#### Figure S3 - Days between doses ----



# Make the dataset for this figure
dfs3 <- data7 %>%
  ungroup() %>%
  dplyr::select(`Subject ID`, vaccine_1_date, vaccine_2_date, vaccine_3_date,
                IDorHV) %>%
  distinct() 

dfs3$vaccine_3_date2 <- dfs3$vaccine_3_date
dfs3[dfs3$vaccine_3_date2 >= "2022-04-28" & 
       !is.na(dfs3$vaccine_3_date2),]$vaccine_3_date2 <- NA

dfs3 %>%
  mutate(diff1 = as.integer(difftime(vaccine_2_date, vaccine_1_date, "days")),
         diff2 = as.integer(difftime(vaccine_3_date2, vaccine_2_date, "days"))) %>%
  dplyr::select(`Subject ID`, diff1, diff2, IDorHV) -> dfs3

# Make the figures
fig_s3a <- ggplot(dfs3, aes(x = diff1)) + 
  geom_histogram(binwidth = 7) + 
  theme_fausto() + 
  scale_x_continuous(breaks = seq(0, 300, 50),
                     limits = c(0, 300)) + 
  scale_y_continuous(limits = c(0, 150)) +
  ylab("Frequency") + 
  xlab("Difference in days\nbetween doses 1 and 2")

fig_s3b <- ggplot(dfs3, aes(x = diff2)) + 
  geom_histogram(binwidth = 7) + 
  theme_fausto() + 
  scale_x_continuous(breaks = seq(0, 450, 75),
                     limits = c(0, 450)) + 
  ylab("Frequency") + 
  xlab("Difference in days\nbetween doses 2 and 3")

# Calculate means and medians for the figure caption
mean(dfs3$diff1, na.rm = T); median(dfs3$diff1, na.rm = T)
mean(dfs3$diff2, na.rm = T); median(dfs3$diff2, na.rm = T)

# Are the times to doses different between IDP and HV?
median(dfs3[dfs3$IDorHV == "HV",]$diff2, na.rm = T)
median(dfs3[dfs3$IDorHV == "PID",]$diff2, na.rm = T)

wilcox.test(dfs3[dfs3$IDorHV == "HV",]$diff2, dfs3[dfs3$IDorHV == "PID",]$diff2)

# Make the final figure
fig_s3 <- (fig_s3a + fig_s3b) &
  theme(panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left")),
        plot.tag.position = c(0, 1),
        plot.tag = element_text(face = "bold", hjust = 0, vjust = 0.5)) & 
  plot_annotation(tag_levels = 'A')

fig_s3



#### Figure S11 - Moderna vs Pfizer ----



# Make the labels for facet_grid()
time.labs2 <- c("Baseline", "Post-dose 1",
                "Post-dose 2", "6 months post-dose 2",
                "Post-dose 3", "6 months post-dose 3")
names(time.labs2) <- c("Baseline", "Post-dose 1",
                       "Post-dose 2", "6 months post-dose 2", 
                       "Post-dose 3", "6 months post-dose 3")

group.labs <- c("Healthy\nvolunteer", "Antibody\ndeficiency",
                "PIRD", "Combined\nimmuno-\ndeficiency",
                "Other\nIEI", "Other\nimmune\ndisorder")
names(group.labs) <- c("Healthy volunteer", "Antibody deficiency",
                       "PIRD", "Combined immunodeficiency", 
                       "Other IEI", "Other immune disorder")

# Get the data that you need for this analysis 
df_s11 <- df1 %>%
  filter(vaccine_1 != "Janssen" & 
           (vaccine_2 != "Janssen" | is.na(vaccine_2)) & 
           (vaccine_3 != "Janssen" | is.na(vaccine_3))) %>%
  filter((switch_1_2 != "Switch" | is.na(switch_1_2)) & 
           (switch_2_3 != "Switch" | is.na(switch_2_3))) %>%
  filter(timepoint != "Baseline") %>%
  dplyr::select(`Subject ID`, timepoint, conc.ug_ml, vaccine_1, new_group)

# Make the plot  
df_s11 %>%                      
  ggplot(aes(x = as.factor(vaccine_1), y = conc.ug_ml)) +
  geom_hline(yintercept = 0.0009426984, linetype = 2) +
  stat_boxplot(aes(x = as.factor(vaccine_1), y = conc.ug_ml, color = as.factor(vaccine_1)), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = as.factor(vaccine_1), y = conc.ug_ml, color = as.factor(vaccine_1)),
              width = 0.2, alpha = 0.4, size = 3.5) + 
  stat_summary(aes(x = as.factor(vaccine_1), y = conc.ug_ml, color = as.factor(vaccine_1)),
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) + 
  scale_color_manual(values = c("Moderna" = "royal blue", 
                                "Pfizer" = "salmon"),
                     labels = c("Moderna" = "mRNA-1273 (Moderna)",
                                "Pfizer" = "BNT162b2 (Pfizer-BioNTech)")) + 
  facet_grid(new_group ~ timepoint, 
             labeller = labeller(timepoint = time.labs2, new_group = group.labs)) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = seq(0, 60, 30)) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     label.y.npc = 0.85, label.x.npc = 0.5,
                     hide.ns = T, size = 7) +
  ylab("Anti-S IgG (ug/ml)") +
  theme_fausto() +
  coord_cartesian(ylim = c(0, 50)) + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.spacing.y = unit(0.05, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, -10, -5, -10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12.5),
        panel.grid.minor.y = element_line(linewidth = rel(0.5), color = "grey92"),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size = 13),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4.5))) -> fig_s11

fig_s11



#### Figure S12 - Vaccine switch ----



# Get the data that you need for this analysis
# This dataset excludes Janssen 
df_s12 <- df1 %>%
  filter(timepoint %in% c("Post-dose 3", "6 months post-dose 3")) %>%
  filter((vaccine_1 != "Janssen" & !is.na(vaccine_1)) & 
           (vaccine_2 != "Janssen" & !is.na(vaccine_2)) & 
           (vaccine_3 != "Janssen" & !is.na(vaccine_3))) %>% 
  mutate(switch = ifelse (vaccine_1 == vaccine_2 & vaccine_2 == vaccine_3, 
                          "Same vaccine brand for doses 1-3", 
                          "Switched vaccine brand across doses 1-3")) %>%
  dplyr::select(`Subject ID`, timepoint, conc.ug_ml, new_group, switch,
                vaccine_1, vaccine_2, vaccine_3)
table(df_s12$switch, useNA = "ifany") # 182 Same, 12 Switch...6% Yes, P(switch) = 0.039

# This dataset allows Janssen, same results either way
df_s12 <- df1 %>%
  filter(timepoint %in% c("Post-dose 3", "6 months post-dose 3")) %>%
  filter((!is.na(vaccine_1)) & 
           (!is.na(vaccine_2)) & 
           (!is.na(vaccine_3))) %>% 
  mutate(switch = ifelse (vaccine_1 == vaccine_2 & vaccine_2 == vaccine_3, 
                          "Same vaccine brand for doses 1-3", 
                          "Switched vaccine brand across doses 1-3")) %>%
  dplyr::select(`Subject ID`, timepoint, conc.ug_ml, new_group, switch,
                vaccine_1, vaccine_2, vaccine_3)
table(df_s12$switch, useNA = "ifany") # 182 Same, 19 Switch...9.5% Yes, P(switch) = 0.039


# Make the plot  
df_s12 %>%                      
  ggplot(aes(x = as.factor(switch), y = conc.ug_ml)) +
  geom_hline(yintercept = 0.0009426984, linetype = 2) +
  stat_boxplot(aes(x = as.factor(switch), y = conc.ug_ml, color = as.factor(switch)), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = as.factor(switch), y = conc.ug_ml, color = as.factor(switch)),
              width = 0.2, alpha = 0.4, size = 3.5) + 
  stat_summary(aes(x = as.factor(switch), y = conc.ug_ml, color = as.factor(switch)),
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) + 
  scale_color_manual(values = c("Same vaccine brand for doses 1-3" = "royal blue", 
                                "Switched vaccine brand across doses 1-3" = "salmon")) + 
  facet_grid(new_group ~ timepoint, 
             labeller = labeller(new_group = group.labs)) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = seq(0, 60, 30)) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     label.y.npc = 0.85, label.x.npc = 0.5,
                     hide.ns = T, size = 7) +
  ylab("Anti-S IgG (ug/ml)") +
  theme_fausto() +
  coord_cartesian(ylim = c(0, 50)) + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.spacing.y = unit(0.05, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, -10, -5, -10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12.5),
        panel.grid.minor.y = element_line(linewidth = rel(0.5), color = "grey92"),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size = 13),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) +  
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4,5))) -> fig_s12

fig_s12



#### Figure S13 - Figure 1b but only for those not on IgRT ----



df1 %>%
  subset(igrt == "No") %>% 
  ggplot(aes(x = new_group, y = conc.ug_ml)) +
  geom_hline(yintercept = 0.0009426984, linetype = 2) +
  stat_boxplot(aes(x = new_group, y = conc.ug_ml, color = new_group), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = new_group, y = conc.ug_ml, color = new_group),
              width = 0.2, alpha = 0.4, size = 2.5) + 
  stat_summary(aes(x = new_group, y = conc.ug_ml, color = new_group),
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = seq(0, 40, 10)) +
  coord_cartesian(ylim = c(0, 35)) + 
  stat_compare_means(ref.group = "Healthy volunteer", 
                     method = "wilcox.test",
                     label = "p.signif",
                     label.y.npc = 1,
                     label.x.npc = 0.5,
                     hide.ns = T, 
                     size = 6.5)+
  facet_wrap(~timepoint, ncol = 6) +
  ylab("Anti-S IgG (ug/ml)") +
  scale_color_manual(name = "Subgroup",
                     values = c("Healthy volunteer" = "black",
                                "Antibody deficiency" = "#7209b7",
                                "PIRD" = "#048B90",
                                "Combined immunodeficiency" = "#CF9A03",
                                "Other IEI" = "#D10106",
                                "Other immune disorder" = "#0006F2")) + 
  theme_fausto() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.spacing.x = unit(0.15, "cm"),
        legend.text = element_text(size = 13.5),
        legend.title = element_text(size = 13.5, face = "bold"),
        strip.text = element_text(size = 14),
        plot.tag.position = c(0, 1),
        plot.tag = element_text(face = "bold", hjust = 0, vjust = 0.75),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4.5), nrow = 1)) -> fig_s13

fig_s13



#### Figure S14 - Rituximab ----



# Get the data for this plot
df6 <- df %>%
  subset(Rituximab == "Yes") %>%
  subset(timepoint != "Post-dose 4" & !is.na(timepoint)) %>%
  mutate(time.diff = ymd(vaccine_1_date) - ymd(last_ritux_before_vaccination),
         rituximab.group = 
           case_when(
             time.diff <= 180 ~ "<6 months",
             time.diff > 180 & time.diff <= 365 ~ "6-12 months",
             time.diff > 365 ~ ">12 months")) 

df6$rituximab.group <- factor(df6$rituximab.group, levels = c("<6 months",
                                                              "6-12 months",
                                                              ">12 months"))

df6 <- df6 %>%
  select(c("Subject ID", timepoint, rituximab.group, conc.ug_ml, Rituximab)) %>% 
  subset(Rituximab == "Yes"  & !is.na(timepoint))
  

# Make the plot
set.seed(1)
df6 %>%
  ggplot(aes(x = rituximab.group, y = conc.ug_ml)) +
  geom_hline(yintercept = 0.0009426984, linetype = 2) +
  stat_boxplot(aes(x = rituximab.group, y = conc.ug_ml, color = rituximab.group), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = rituximab.group, y = conc.ug_ml, color = rituximab.group),
              width = 0.2, alpha = 0.4, size = 3.5) + 
  stat_summary(aes(x = rituximab.group, y = conc.ug_ml, color = rituximab.group),
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  facet_wrap(~timepoint, ncol = 3) + 
  ylab("Anti-S IgG (ug/ml)") +
  scale_color_manual(name = "Time since last receipt of rituximab relative to date of dose 1",
                     values = c("<6 months" = "#d60270",
                                "6-12 months" = "#7209b7",
                                ">12 months" = "#048B90")) + 
  theme_fausto() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.tag.position = c(0, 1),
        plot.tag = element_text(face = "bold", hjust = 0, vjust = 0.75),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4.5), nrow = 1)) -> fig_s14

fig_s14



#### Figure S15 - N-positivity ----



df1 %>% 
  ggplot(aes(x = factor(N_pos), y = conc.ug_ml)) +
  geom_hline(yintercept = 0.0009426984, linetype = 2) +
  stat_boxplot(aes(x = factor(N_pos), y = conc.ug_ml, color = factor(N_pos)), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = factor(N_pos), y = conc.ug_ml, color = factor(N_pos)),
              width = 0.2, alpha = 0.4, size = 3.5) + 
  stat_summary(aes(x = factor(N_pos), y = conc.ug_ml, color = factor(N_pos)),
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) + 
  scale_color_manual(values = c("0" = "royal blue", "1" = "salmon"),
                     labels = c("0" = "Anti-N IgG negative", "1" = "Anti-N IgG positive")) + 
  facet_wrap(~timepoint, ncol = 6) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     label.y.npc = 1,
                     label.x.npc = 0.5,
                     hide.ns = T, size = 7.5)+
  ylab("Anti-S IgG (ug/ml)") +
  theme_fausto() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 14),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4.5))) -> fig_s15

fig_s15



#### Figure S16a - Adverse Events: Local ----



# Get the data sand format it
fp_survey <- read.csv(here("C:/Users/bustosfa/Rotation 2/Data", "fp_survey.csv"))

fp_survey$PATIENT.TYPE[fp_survey$PATIENT.TYPE == "Patient Type: Healthy Volunteer"] <- "HV"
fp_survey$PATIENT.TYPE[fp_survey$PATIENT.TYPE == "Patient Type: Immune Deficient"] <- "IDP"

fp_survey$PATIENT.TYPE <- factor(fp_survey$PATIENT.TYPE, levels = c("HV", "IDP"))

# Calculate the numbers you need for the graph
fp_survey %>%
  mutate(any_local_AE = factor(any_local_AE, levels=c("0", "1"), labels = c("None", "Any Local AE"))) %>%
  filter(any_local_AE != "") %>%
  filter(TREATMENT.GROUP != "") %>%
  filter(X3..Which.dose.is.this.for != "",
         X3..Which.dose.is.this.for != "Dose 4",
         X3..Which.dose.is.this.for != "Dose 5",
         X3..Which.dose.is.this.for != "NA") %>%
  dplyr::select(id = "Participant.ID", control = PATIENT.TYPE, any_local_AE, dose = X3..Which.dose.is.this.for) %>%
  mutate(local_ae = ifelse(any_local_AE == "None", 0, 1)) %>%
  dplyr::select(-any_local_AE) %>%
  group_by(dose, control) %>%
  summarize(n = n(), x = sum(local_ae == 1)) %>%
  group_by(dose, control) %>%
  mutate(mean = round(binom.test(x, n)$estimate, 3)) %>%
  mutate(LL = round(binom.test(x, n)$conf.int[1], 3)) %>%
  mutate(UL = round(binom.test(x, n)$conf.int[2], 3)) %>%
  ungroup() -> local_ae_df

# Make the plot
ggplot(local_ae_df, aes(x = as.factor(control), y = mean, 
                        ymin = LL, ymax = UL,
                        color = as.factor(control))) + 
  geom_linerange(linewidth = 1.25) + 
  geom_point(size = 4) +
  theme_fausto() + 
  facet_wrap(~dose) + 
  scale_y_continuous("Frequency (%)",
                     breaks = seq(0, 1, 0.2),
                     labels = paste0(seq(0, 1, 0.2)*100, "%"),
                     limits = c(0, 1)) + 
  scale_color_manual(values = c("HV" = "royal blue", "IDP" = "salmon")) + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) -> fig_s16a

fig_s16a



#### Figure S16b - Adverse Events: Systemic ----



# Get the data that I need
fp_survey %>%
  mutate(any_systemic_AE = factor(any_systemic_AE, levels=c("0", "1"), labels = c("None", "Any Systemic AE"))) %>%
  filter(any_systemic_AE != "") %>%
  filter(TREATMENT.GROUP != "") %>%
  filter(X3..Which.dose.is.this.for != "",
         X3..Which.dose.is.this.for != "Dose 4",
         X3..Which.dose.is.this.for != "Dose 5",
         X3..Which.dose.is.this.for != "NA") %>%
  dplyr::select(id = "Participant.ID", control = PATIENT.TYPE, any_systemic_AE, dose = X3..Which.dose.is.this.for) %>%
  mutate(systemic_ae = ifelse(any_systemic_AE == "None", 0, 1)) %>%
  dplyr::select(-any_systemic_AE) %>%
  group_by(dose, control) %>%
  summarize(n = n(), x = sum(systemic_ae == 1)) %>%
  group_by(dose, control) %>%
  mutate(mean = round(binom.test(x, n)$estimate, 3)) %>%
  mutate(LL = round(binom.test(x, n)$conf.int[1], 3)) %>%
  mutate(UL = round(binom.test(x, n)$conf.int[2], 3)) %>%
  ungroup() -> systemic_ae_df

# Make the plot
ggplot(systemic_ae_df, aes(x = as.factor(control), y = mean, 
                           ymin = LL, ymax = UL,
                           color = as.factor(control))) + 
  geom_linerange(linewidth = 1.25) + 
  geom_point(size = 4) +
  theme_fausto() + 
  facet_wrap(~dose) + 
  scale_y_continuous("Frequency (%)",
                     breaks = seq(0, 1, 0.2),
                     labels = paste0(seq(0, 1, 0.2)*100, "%"),
                     limits = c(0, 1)) + 
  scale_color_manual(values = c("HV" = "royal blue", "IDP" = "salmon")) + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) -> fig_s16b

fig_s16b



#### Figure S16 - Adverse Events ----




local_ae_df$type <- "Local\nadverse events"
systemic_ae_df$type <- "Systemic\nadverse events"

ae_df <- rbind(local_ae_df, systemic_ae_df)

ggplot(ae_df, aes(x = as.factor(control), y = mean, 
                  ymin = LL, ymax = UL,
                  color = as.factor(control))) + 
  geom_linerange(linewidth = 1.25) + 
  geom_point(size = 4) +
  theme_fausto() + 
  facet_grid(type ~ dose) + 
  scale_y_continuous("Frequency (%)",
                     breaks = seq(0, 1, 0.2),
                     labels = paste0(seq(0, 1, 0.2)*100, "%"),
                     limits = c(0, 1)) + 
  scale_color_manual(values = c("HV" = "black", 
                                "IDP" = "#d60270"),
                     labels = c("HV" = "Healthy volunteers",
                                "IDP" = "Immune-deficient/disordered people")) + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) -> fig_s16

fig_s16



### Seroconversion rates ----



# Seroconversion rates by HV status and time point
df %>%
  group_by(control, timepoint) %>%
  summarize(detected = mean(conc_cutoff == "Positive"))
  


### Describe titers over time, Tables S4-S5 ----



# Median titers by HV status and timepoint
df %>%
  filter(timepoint != "Baseline") %>%
  group_by(control, timepoint) %>%
  summarize(medians = median(conc.ug_ml)) -> seroconversion

seroconversion %>%
  pivot_wider(id_cols = timepoint, names_from = control, values_from = medians) -> sero2

sero2 %>%
  ungroup() %>%
  group_by(timepoint) %>%
  mutate(across(everything(), ~ . / `Healthy volunteer`)*100) -> sero3

# Look at titers over time (Table S4)
df %>%
  group_by(control, timepoint) %>%
  summarize(titers = round(median(conc.ug_ml), 3))

# Look at titers over time, raw numbers (Table S4)
df %>%
  group_by(new_group, timepoint) %>%
  filter(timepoint != "Baseline") %>%
  summarize(titers = round(median(conc.ug_ml), 3)) %>%
  pivot_wider(id_cols = timepoint, names_from = new_group, values_from = titers) -> check

# Look at titers over time, relative numbers (Table S3)
check %>%
  ungroup() %>%
  filter(timepoint != "6 months post-dose 3") %>%
  group_by(timepoint) %>%
  mutate(across(everything(), ~ . / `Healthy volunteer`)*100) -> check2 

# Calculate decrease in titers from pd2 to 6mpds
df %>%
  group_by(control, timepoint) %>%
  summarize(titers = median(conc.ug_ml)) %>% 
  filter(timepoint %in% c("Post-dose 2", "6 months post-dose 2")) %>%
  pivot_wider(id_cols = control, 
              names_from = timepoint, 
              values_from = titers) %>%
  mutate(decrease = paste0(round(100 * (1 - `6 months post-dose 2`/`Post-dose 2`), 1), "%"),
         fold = `Post-dose 2` / `6 months post-dose 2`)
  
rm(seroconversion, sero2, sero3, check, check2)
rm(text1, text2, text3, group.labs, i, time.labs, time.labs2)
rm(check_wide, check_wide2, df_s11, df_s12, df3c, df6)



### Compare IDP titers at pd3 and 6mpd3 ----



# Identify the breakthrough infections
df$breakthru = case_when(ymd(df$date_first_COVID_pos) > ymd(df$vaccine_1_date) & 
                           ymd(df$date_first_COVID_pos) <= ymd("2022-04-28") ~ "Breakthru",
                         ymd(df$date_first_COVID_pos) < ymd(df$vaccine_1_date) ~ "COVID pre-dose 1",
                         TRUE ~ "No breakthru")

# Check if the titer levels among IDP differ between post-dose 3 and 6 months later
check <- df %>%
  ungroup() %>% 
  dplyr::select(`Subject ID`, timepoint, conc.ug_ml, control, breakthru) %>%
  filter(control == "PID") %>%
  group_by(`Subject ID`) %>%
  filter(timepoint %in% c("Post-dose 3", "6 months post-dose 3")) %>%
  filter(all(c("Post-dose 3", "6 months post-dose 3") %in% timepoint)) %>%
  arrange(`Subject ID`) %>%
  ungroup() %>%
  dplyr::select(-control)

check$timepoint <- droplevels(check$timepoint)

check_wide <- check %>%
  pivot_wider(id_cols = c(`Subject ID`, breakthru), 
              names_from = timepoint, values_from = conc.ug_ml)

wilcox.test(check_wide$`Post-dose 3`, check_wide$`6 months post-dose 3`, paired = TRUE)
# p-value = 0.5988

# Now remove the info for the breakthrough infections and check again
check_wide2 <- check_wide %>%
  filter(breakthru == "No breakthru")

wilcox.test(check_wide2$`Post-dose 3`, check_wide2$`6 months post-dose 3`, paired = TRUE)
# p-value = 0.5466



### Detectability of titers for rituximab users ----


df6 <- df %>%
  subset(Rituximab == "Yes") %>%
  subset(timepoint != "Post-dose 4" & !is.na(timepoint)) %>%
  mutate(time.diff = ymd(vaccine_1_date) - ymd(last_ritux_before_vaccination),
         rituximab.group = 
           case_when(
             time.diff <= 180 ~ "<6 months",
             time.diff > 180 & time.diff <= 365 ~ "6-12 months",
             time.diff > 365 ~ ">12 months")) 

df6$rituximab.group <- factor(df6$rituximab.group, levels = c("<6 months",
                                                              "6-12 months",
                                                              ">12 months"))

df6 <- df6 %>%
  select(c("Subject ID", timepoint, rituximab.group, conc.ug_ml, conc_cutoff, Rituximab)) %>% 
  subset(Rituximab == "Yes"  & !is.na(timepoint)) %>%
  mutate(Pos = ifelse(conc_cutoff == "Positive", 1, 0)) %>%
  group_by(`Subject ID`) %>%
  mutate(Pos_at_least_once = max(Pos)) %>%
  arrange(`Subject ID`)

check <- df6 %>%
  filter(timepoint != "Baseline") %>%
  slice(1)
table(check$Pos_at_least_once)
rm(check, df6, df1)



#### Figure 4: Variant wave & Table S10 ----



# Process the variants data
names(variants)[1] <- "region"
variants2 <- variants %>%
  filter(region == "USA") %>%
  dplyr::select(week_ending, variant, share)

variants2$week_ending <- format(as.Date(as.POSIXct(variants2$week_ending, format = '%m/%d/%Y %H:%M'),
                                        format='%m/%d/%Y'))

variants2[variants2 == "B.1.1.529"] <- "Omicron BA.1"
variants2[variants2 == "BA.1.1"] <- "Omicron BA.1"
variants2[variants2 == "BA.2.12.1"] <- "Omicron BA.2"
variants2[variants2 == "BA.2"] <- "Omicron BA.2"
variants2[variants2 == "BA.4"] <- "Omicron BA.4/5"
variants2[variants2 == "BA.4.6"] <- "Omicron BA.4/5"
variants2[variants2 == "BA.5"] <- "Omicron BA.4/5"
variants2[variants2 == "B.1.617.2"] <- "Delta"
variants2[variants2 == "AY.1"] <- "Delta"
variants2[variants2 == "AY.2"] <- "Delta"
variants2[variants2 == "AY.3"] <- "Delta"
variants2[variants2 == "B.1.1.7"] <- "Alpha"
variants2[variants2 == "P.1"] <- "Gamma"
variants2[variants2 == "B.1.621"] <- "Mu"
variants2[variants2 == "B.1.526"] <- "Iota"
variants2[variants2 == "B.1.628"] <- "Other"
variants2[variants2 == "B.1.621.1"] <- "Mu"
variants2[variants2 == "A.2.5"] <- "Other"
variants2[variants2 == "B.1.351"] <- "Beta"
variants2[variants2 == "B.1.626"] <- "Other"
variants2[variants2 == "B.1.617.1"] <- "Kappa"
variants2[variants2 == "B.1.525"] <- "Eta"
variants2[variants2 == "B.1.617.3"] <- "Other"
variants2[variants2 == "B.1.526.2"] <- "Iota"
variants2[variants2 == "B.1.526.1"] <- "Iota"
variants2[variants2 == "B.1"] <- "Other"
variants2[variants2 == "B.1.1.519"] <- "Other"
variants2[variants2 == "B.1.429"] <- "Epsilon"
variants2[variants2 == "B.1.427"] <- "Epsilon"
variants2[variants2 == "B.1.2"] <- "Other"
variants2[variants2 == "B.1.596"] <- "Other"
variants2[variants2 == "P.2"] <- "Zeta"
variants2[variants2 == "B.1.617"] <- "Other"
variants2[variants2 == "B.1.637"] <- "Other"
variants2[variants2 == "B.1.1"] <- "Other"
variants2[variants2 == "B.1.1.194"] <- "Other"
variants2[variants2 == "B.1.243"] <- "Other"
variants2[variants2 == "B.1.234"] <- "Other"
variants2[variants2 == "R.1"] <- "Other"

variants3 <- variants2 %>%
  group_by(week_ending, variant) %>%
  summarize(percent = sum(share)) %>%
  ungroup()

min(variants3$week_ending)

unique(variants3$variant)
variants3$variant <- factor(variants3$variant,
                            levels = c("Alpha", "Beta", "Gamma", "Delta", "Epsilon",
                                       "Zeta", "Eta", "Iota", "Kappa", "Mu",
                                       "Omicron BA.1", "Omicron BA.2", 
                                       "Omicron BA.4/5", "Other"))

variants3 <- arrange(variants3, week_ending, variant)

# Deal with the x-axis labels
variants3$label <- format(as.Date(variants3$week_ending),"%b %d")
variants3$label2 <- format(as.Date(variants3$week_ending),"%b %d %y")
variants3$label3 <- format(as.Date(variants3$week_ending),"%b %d %y")

every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}

labels <- unique(variants3$label3)
labels <- labels[c(T, F, F, F)]

labels3 <- variants3 %>%
  filter(label3 %in% labels) %>%
  distinct(label, .keep_all = T)
labels_vec <- labels3$label3
labels_vec <- gsub("(\\d{2}) (\\d+)","\\1, \\2", labels_vec)

# Get the dates of all saturdays
saturdays <-  seq(as.Date(min(variants3$week_ending)), 
                  as.Date(max(variants3$week_ending)), 
                  "7 days")

## Get the data set up for the line version of the final graph

# Get the dates of the breakthru cases
bt <- df %>%
  filter(evusheld_before_sample == 0) %>%
  filter(breakthru != "COVID pre-dose 1") %>%
  distinct(`Subject ID`, .keep_all = T)
bt$pos_before <- ifelse(bt$COVID_pos_before_sample == 1 |
                          bt$N_pos_before_sample == 1, 1, 0)
bt$Healthy_binary <- ifelse(bt$new_group == "Healthy volunteer", 1, 0)
bt$bt_binary <- ifelse(bt$breakthru == "Breakthru", 1, 0)

newvaxdf <- bt %>%
  filter(bt_binary==1) %>%
  dplyr::select(`Subject ID`, vaccine_2_date, vaccine_3_date, date_first_COVID_pos) %>%
  mutate(between_vax23 = ifelse(date_first_COVID_pos >= vaccine_2_date & 
                                  date_first_COVID_pos <= vaccine_3_date, 1, 0)) %>%
  mutate(after_vax3 = ifelse(date_first_COVID_pos > vaccine_3_date, 1, 0)) %>%
  mutate(after_vax2_no_vax3 = ifelse(date_first_COVID_pos >= vaccine_2_date & 
                                       is.na(vaccine_3_date), 1, 0))

# Make sure dates are formatted
newvaxdf$vaccine_2_date <- format(as.Date(newvaxdf$vaccine_2_date, 
                                          format = '%m/%d/%Y'))
newvaxdf$vaccine_3_date <- format(as.Date(newvaxdf$vaccine_3_date, 
                                          format = '%m/%d/%Y'))
newvaxdf$date_first_COVID_pos <- format(as.Date(newvaxdf$date_first_COVID_pos, 
                                                format = '%m/%d/%Y'))

# Get the saturday for the second vaccine date
newvaxdf$saturday_vax2 <- NA
for(i in 1:nrow(newvaxdf)){
  for(j in 1:length(saturdays)){
    if(newvaxdf$vaccine_2_date[i] > saturdays[j]){
      newvaxdf$saturday_vax2[i] <- format(as.Date(saturdays[j + 1], format = '%m/%d/%Y'))
    }
  }
}

# Get the saturday for the date of infection
newvaxdf$saturday_covid <- NA
for(i in 1:nrow(newvaxdf)){
  for(j in 1:length(saturdays)){
    if(newvaxdf$date_first_COVID_pos[i] > saturdays[j]){
      newvaxdf$saturday_covid[i] <- format(as.Date(saturdays[j + 1], format = '%m/%d/%Y'))
    }
  }
}

# Get the saturday for the third vaccine date
newvaxdf$saturday_vax3 <- NA
for(i in 1:nrow(newvaxdf)){
  for(j in 1:length(saturdays)){
    if(!is.na(newvaxdf$vaccine_3_date[i]) & newvaxdf$vaccine_3_date[i] > saturdays[j]){
      newvaxdf$saturday_vax3[i] <- format(as.Date(saturdays[j + 1], format = '%m/%d/%Y'))
    }
  }
}

# Simplify the dataset
newvaxdf2 <- newvaxdf %>%
  dplyr::select(id = `Subject ID`, saturday_vax2, saturday_vax3, saturday_covid)

# Format the dataset appropriately
newvaxdf2$variant = "Kappa"
newvaxdf2 <- newvaxdf2 %>%
  arrange(saturday_covid)
newvaxdf2$id <- factor(newvaxdf2$id, levels = newvaxdf2$id)
newvaxdf2$percent = seq(1, 0, 
                        length.out = nrow(newvaxdf) + 2)[c(F, rep(T, nrow(newvaxdf)), F)]

# Get the most common variant for each breakthrough infection
likely_var <- newvaxdf2 %>%
  dplyr::select(id, saturday_covid)

variants4 <- variants3 %>%
  group_by(week_ending) %>%
  mutate(most_common = variant[which.max(percent)]) %>%
  ungroup() %>%
  dplyr::select(week_ending, most_common) %>%
  distinct()

likely_var <- left_join(likely_var, variants4, by = c("saturday_covid" = "week_ending"))
likely_var$most_common <- as.character(likely_var$most_common)

likely_var <- likely_var %>%
  mutate(confirmed = case_when(id == "E0012" ~ "Omicron BA.2",
                               id == "E0013" ~ "Omicron BA.1",
                               id == "E0014" ~ "Omicron BA.1",
                               id == "E0055" ~ "Delta",
                               id == "E0061" ~ "Omicron BA.2",
                               
                               id == "E0082" ~ "Omicron BA.2",
                               id == "E0087" ~ "Omicron BA.1",
                               id == "E0107" ~ "Omicron BA.1",
                               id == "E0118" ~ "Omicron BA.2",
                               id == "E0140" ~ "Omicron BA.2",
                               
                               id == "E0169" ~ "Omicron BA.1",
                               id == "E0209" ~ "Omicron BA.2",
                               id == "H0015" ~ "Omicron BA.2",
                               id == "H0022" ~ "Omicron BA.1",
                               id == "H0024" ~ "Omicron BA.1",
                               
                               TRUE ~ "None"))

likely_var$confirmed <- ifelse(likely_var$confirmed == "None", NA, likely_var$confirmed)

likely_var$final <- ifelse(!is.na(likely_var$confirmed), 
                           likely_var$confirmed, likely_var$most_common)

table(likely_var[!is.na(likely_var$confirmed),]$final) # confirmed numbers
table(likely_var[is.na(likely_var$confirmed),]$final) # assigned numbers

likely_var$confirmed_yn <- ifelse(!is.na(likely_var$confirmed), 1, 0)

integrate <- likely_var %>%
  dplyr::select(id, most_common, final, confirmed_yn)

# Integrate back into newvaxdf2
newvaxdf2 <- left_join(newvaxdf2, integrate, by = "id")

# For the 2 people who have a 3rd shot after "2022-06-04", set that to NA
newvaxdf2[newvaxdf2$saturday_vax3 >= "2022-06-04" & 
            !is.na(newvaxdf2$saturday_vax3),]$saturday_vax3 <- NA

# Add in variable for HV or not
controlvar <- df %>%
  dplyr::select(id = `Subject ID`, control) %>%
  filter(id %in% newvaxdf2$id) %>%
  distinct()

newvaxdf2 <- left_join(newvaxdf2, controlvar, by = "id")
newvaxdf2 <- as.data.frame(newvaxdf2)

# Make final a factor variable
newvaxdf2$final <- factor(newvaxdf2$final, levels = levels(variants3$variant))

# Make a line dataset
newvaxdf2_line <- newvaxdf2 %>%
  dplyr::select(-saturday_vax3)

newvaxdf2_line <- gather(newvaxdf2_line, saturday, week_ending, 
                         saturday_vax2:saturday_covid, factor_key = TRUE)

newvaxdf2_line$id <- factor(newvaxdf2_line$id, levels = newvaxdf2$id)

newvaxdf2_line <- newvaxdf2_line %>%
  arrange(id, week_ending)

# FOR THE ZOOM WINDOW
# Now make the final version of the variant_plot
variant_plot_swimmer_plot <- ggplot(variants3[variants3$week_ending <= "2022-06-04",],  
                                    aes(fill = variant, y = percent, x = week_ending)) + 
  geom_bar(position = "fill", stat = "identity", width = 1) + 
  ylab("National percentage of variants in USA") + 
  xlab("Week ending") +
  theme_fausto() +
  coord_cartesian(ylim = c(0, 1), expand = F) + 
  scale_y_continuous(labels = scales::percent_format(scale = 100), 
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  theme(legend.position = "right") + 
  scale_fill_manual(name = "SARS-CoV-2 variant", #B2 adds alpha of 0.7, see alpha("red", .3)
                    values = c("Alpha" = "#A732F6B2",
                               "Beta" = "#590004B2",
                               "Epsilon" = "#E3A803B2",
                               "Eta" = "#D10106B2",
                               "Gamma" = "#2103D7B2",
                               "Iota" = "#00BCD0B2",
                               "Zeta" = "#F17427B2",
                               "Delta" = "#048B90B2",
                               "Kappa" = "#FF9C85B2",
                               "Mu" = "#29A91BB2",
                               "Omicron BA.1" = "#FFCDFFB2",
                               "Omicron BA.2" = "#FF2FFFB2",
                               "Omicron BA.4/5" = "#8E008EB2",
                               "Other" = "#808080B2")) + 
scale_color_manual(name = "SARS-CoV-2 variant", #B2 adds alpha of 0.7, see alpha("red", .3)
                   values = c("Alpha" = "#A732F6B2",
                              "Beta" = "#590004B2",
                              "Epsilon" = "#E3A803B2",
                              "Eta" = "#D10106B2",
                              "Gamma" = "#2103D7B2",
                              "Iota" = "#00BCD0B2",
                              "Zeta" = "#F17427B2",
                              "Delta" = "#048B9B2",
                              "Kappa" = "#FF9C85B2",
                              "Mu" = "#29A91BB2",
                              "Omicron BA.1" = "#FFCDFFB2",
                              "Omicron BA.2" = "#FF2FFFB2",
                              "Omicron BA.4/5" = "#8E008EB2",
                              "Other" = "#808080B2")) +
scale_x_discrete(breaks = every_nth(n = 4),
                 labels = labels_vec) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.y.right = element_text(angle = 90)) + 
  geom_vline(xintercept = "2021-05-01", color = "black", 
             linewidth = 1.75, linetype = "dashed") + 
  guides(color = guide_legend(override.aes = list(linetype = rep(1, 14), color = "black"))) + 
  geom_line(inherit.aes = F, data = newvaxdf2_line, # Add swimmer lines
            aes(group = id, y = percent, x = week_ending),
            linewidth = 0.3, color = "black") + 
  geom_point(inherit.aes = F, data = newvaxdf2, # Add dose 2
             aes(group = id, y = percent, x = saturday_vax2),
             size = 3.5, color = "black", fill = "gray28", shape = 21) + 
  geom_point(inherit.aes = F, data = newvaxdf2, # Add dose 3
             aes(group = id, y = percent, x = saturday_vax3),
             size = 3.5, color = "black", fill = "gray28", shape = 22) + 
  geom_point(inherit.aes = F, # Add infection for unconfirmed PID
             data = newvaxdf2[newvaxdf2$confirmed_yn == 0 & newvaxdf2$control == "PID",],
             aes(group = id, y = percent, x = saturday_covid),
             size = 3.5, color = "black", fill = "gray28", shape = 24) + 
  geom_point(inherit.aes = F, # Add infection for unconfirmed HIV
             data = newvaxdf2[newvaxdf2$confirmed_yn == 0 & newvaxdf2$control != "PID",],
             aes(group = id, y = percent, x = saturday_covid),
             size = 3.5, color = "black", fill = "gray28", shape = 25) +
  geom_point(inherit.aes = F, # Add infection for confirmed PID
             data = newvaxdf2[newvaxdf2$confirmed_yn == 1 & newvaxdf2$control == "PID",],
             aes(group = id, y = percent, x = saturday_covid, fill = final),
             size = 3.5, color = "black", shape = 24, show.legend = FALSE, alpha = 1) + 
  geom_point(inherit.aes = F, # Add infection for confirmed HV
             data = newvaxdf2[newvaxdf2$confirmed_yn == 1 & newvaxdf2$control != "PID",],
             aes(group = id, y = percent, x = saturday_covid, fill = final),
             size = 3.5, color = "black", shape = 25, show.legend = FALSE, alpha = 1)


### FOR THE PDF
# Now make the final version of the variant_plot
variant_plot_swimmer_plot <- ggplot(variants3[variants3$week_ending <= "2022-06-04",],  
                                    aes(fill = variant, y = percent, x = week_ending)) + 
  geom_bar(position = "fill", stat = "identity", width = 0.92, color = NA, size = 0) + # This combination works when you crop thingssss!!!!
  ylab("National percentage of variants in USA") + 
  xlab("Week ending") +
  theme_fausto() +
  coord_cartesian(ylim = c(0, 1), expand = F) + 
  scale_y_continuous(labels = scales::percent_format(scale = 100), 
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
  theme(legend.position = "right") + 
  scale_fill_manual(name = "SARS-CoV-2 variant", #B2 adds alpha of 0.7, see alpha("red", .3)
                    values = c("Alpha" = "#A732F6B2",
                               "Beta" = "#590004B2",
                               "Epsilon" = "#E3A803B2",
                               "Eta" = "#D10106B2",
                               "Gamma" = "#2103D7B2",
                               "Iota" = "#00BCD0B2",
                               "Zeta" = "#F17427B2",
                               "Delta" = "#048B90B2",
                               "Kappa" = "#FF9C85B2",
                               "Mu" = "#29A91BB2",
                               "Omicron BA.1" = "#FFCDFFB2",
                               "Omicron BA.2" = "#FF2FFFB2",
                               "Omicron BA.4/5" = "#8E008EB2",
                               "Other" = "#808080B2")) + 
  scale_x_discrete(breaks = every_nth(n = 4),
                   labels = labels_vec) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.y.right = element_text(angle = 90)) + 
  geom_vline(xintercept = "2021-05-01", color = "black", 
             linewidth = 1.75, linetype = "dashed") + 
  guides(color = guide_legend(override.aes = list(linetype = rep(1, 14), color = "black"))) + 
  geom_line(inherit.aes = F, data = newvaxdf2_line, # Add swimmer lines
            aes(group = id, y = percent, x = week_ending),
            linewidth = 0.3, color = "black") + 
  geom_point(inherit.aes = F, data = newvaxdf2, # Add dose 2
             aes(group = id, y = percent, x = saturday_vax2),
             size = 3.5, color = "black", fill = "gray28", shape = 21) + 
  geom_point(inherit.aes = F, data = newvaxdf2, # Add dose 3
             aes(group = id, y = percent, x = saturday_vax3),
             size = 3.5, color = "black", fill = "gray28", shape = 22) + 
  geom_point(inherit.aes = F, # Add infection for unconfirmed PID
             data = newvaxdf2[newvaxdf2$confirmed_yn == 0 & newvaxdf2$control == "PID",],
             aes(group = id, y = percent, x = saturday_covid),
             size = 3.5, color = "black", fill = "gray28", shape = 24) + 
  geom_point(inherit.aes = F, # Add infection for unconfirmed HIV
             data = newvaxdf2[newvaxdf2$confirmed_yn == 0 & newvaxdf2$control != "PID",],
             aes(group = id, y = percent, x = saturday_covid),
             size = 3.5, color = "black", fill = "gray28", shape = 25) +
  geom_point(inherit.aes = F, # Add infection for confirmed PID
             data = newvaxdf2[newvaxdf2$confirmed_yn == 1 & newvaxdf2$control == "PID",],
             aes(group = id, y = percent, x = saturday_covid, fill = final),
             size = 3.5, color = "black", shape = 24, show.legend = FALSE, alpha = 1) + 
  geom_point(inherit.aes = F, # Add infection for confirmed HV
             data = newvaxdf2[newvaxdf2$confirmed_yn == 1 & newvaxdf2$control != "PID",],
             aes(group = id, y = percent, x = saturday_covid, fill = final),
             size = 3.5, color = "black", shape = 25, show.legend = FALSE, alpha = 1)


variant_plot_swimmer_plot

save_plot("fig4.pdf", variant_plot_swimmer_plot, base_height = 6.5, base_width = 13.25, dpi = 300) # Excellent

rm(bt, controlvar, integrate, labels3, likely_var, newvaxdf2_line,
   variants, variants2, variants3, variants4, by, i, j, labels,
   labels_vec, saturdays)



#### Median time bw most recent dose and breakthrough infection ----



# Assign NA to dose 3 date for those after infection or realm of consideration
newvaxdf[newvaxdf$vaccine_3_date >= "2022-06-04" & 
            !is.na(newvaxdf$vaccine_3_date),]$vaccine_3_date <- NA

newvaxdf[newvaxdf$vaccine_3_date >= newvaxdf$date_first_COVID_pos & 
           !is.na(newvaxdf$vaccine_3_date),]$vaccine_3_date <- NA

# Calculate date of doses to breakthrough infection
newvaxdf$date_dose2_to_inf <- as.numeric(difftime(newvaxdf$date_first_COVID_pos, 
                                                  newvaxdf$vaccine_2_date, 
                                                  units = "days"))

newvaxdf$date_dose3_to_inf <- as.numeric(difftime(newvaxdf$date_first_COVID_pos, 
                                                  newvaxdf$vaccine_3_date, 
                                                  units = "days"))

# Calculate days from covid infection to most recent vaccination
newvaxdf %>%
  group_by(`Subject ID`) %>%
  mutate(most_recent_dose_to_inf = min(c(date_dose2_to_inf, date_dose3_to_inf),
                                       na.rm = T)) %>%
  ungroup() -> check

check <- left_join(check, newvaxdf2 %>% dplyr::select(id, control), 
                   by = c(`Subject ID` = "id"))

# Get the medians for IDP and HV
median(check[check$control == "PID",]$most_recent_dose_to_inf) / 30
median(check[check$control != "PID",]$most_recent_dose_to_inf) / 30

wilcox.test(check[check$control == "PID",]$most_recent_dose_to_inf, 
            check[check$control != "PID",]$most_recent_dose_to_inf)

# By the time of infection, how many had 2 doses? 3 doses?
sum(!is.na(check$date_dose2_to_inf))
sum(!is.na(check$date_dose3_to_inf))

rm(check)
                   