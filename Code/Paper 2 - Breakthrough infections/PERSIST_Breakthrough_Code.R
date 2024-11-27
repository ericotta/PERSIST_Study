##################################################################################
#                                                                                #
#           Analysis of breakthrough infections in the PERSIST cohort            #
#                                                                                #
##################################################################################



# **** PRELIMINARIES **** ----



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
library(ggfortify)
library(beepr)
library(ggpubr)
library(rcompanion)
library(writexl)
library(grid)
library(beepr)
library(stringr)
library(varhandle)
library(coin)
library(waffle)
library(survival)
library(survminer)
remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(tableone)
library(kableExtra)
library(MatchIt)



# User-defined functions ----



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
          axis.title = element_text(size = 17, face = "bold"),
          axis.text = element_text(size = 16),
          strip.text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15, face = "bold"),
          legend.margin = margin(unit(0, "cm")),
          strip.text.x = element_text(margin = margin(0.5, 0, 0.5, 0, "mm")),
          strip.background = element_blank())
}


# LOD for the ELISA
LOD_elisa <- 0.0009426984


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
    xlist <- append(xlist, c(x, x + width))
    ylist <- append(ylist, c(y, y))
    idlist <- append(idlist, c(1, 1))
  }
  if ("top" %in% type) { # top
    xlist <- append(xlist, c(x, x + width))
    ylist <- append(ylist, c(y + height, y + height))
    idlist <- append(idlist, c(2, 2))
  }
  if ("left" %in% type) { # left
    xlist <- append(xlist, c(x, x))
    ylist <- append(ylist, c(y, y + height))
    idlist <- append(idlist, c(3, 3))
  }
  if ("right" %in% type) { # right
    xlist <- append(xlist, c(x + width, x + width))
    ylist <- append(ylist, c(y, y + height))
    idlist <- append(idlist, c(4, 4))
  }
  if (length(type) == 0 || "none" %in% type) { # blank; cannot pass absence of coordinates, so pass a single point and use an invisible line
    xlist <- c(x, x)
    ylist <- c(y, y)
    idlist <- c(5, 5)
    linetype <- "blank"
  }
  gp <- gpar(lwd = len0_null(size * .pt), col = colour, lty = linetype)
  element_gp <- gpar(lwd = len0_null(element$size * .pt), col = element$colour, lty = element$linetype)
  polylineGrob(
    x = xlist, y = ylist, id = idlist, ..., default.units = "npc",
    gp = modifyList(element_gp, gp),
  )
}



# Load publicly available datasets ----



# Load the data from the PERSIST study
behavior_raw_public <- read.csv("C:/Users/bustosfa/Box/PIDvax Study/Queries/Analysis Datasets/Final data and code for breakthrough paper/Publicly released data and code/behavior_raw_public.csv")
inf_public <- read.csv("C:/Users/bustosfa/Box/PIDvax Study/Queries/Analysis Datasets/Final data and code for breakthrough paper/Publicly released data and code/inf_public.csv")
new_elisa_public <- read.csv("C:/Users/bustosfa/Box/PIDvax Study/Queries/Analysis Datasets/Final data and code for breakthrough paper/Publicly released data and code/new_elisa_public.csv")
outcomes_public <- read.csv("C:/Users/bustosfa/Box/PIDvax Study/Queries/Analysis Datasets/Final data and code for breakthrough paper/Publicly released data and code/outcomes_public.csv")
pp_raw_public <- read.csv("C:/Users/bustosfa/Box/PIDvax Study/Queries/Analysis Datasets/Final data and code for breakthrough paper/Publicly released data and code/pp_raw_public.csv")
reinfection_raw_public <- read.csv("C:/Users/bustosfa/Box/PIDvax Study/Queries/Analysis Datasets/Final data and code for breakthrough paper/Publicly released data and code/reinfection_raw_public.csv")
sdat_public <- read.csv("C:/Users/bustosfa/Box/PIDvax Study/Queries/Analysis Datasets/Final data and code for breakthrough paper/Publicly released data and code/sdat_public.csv")
visit_public <- read.csv("C:/Users/bustosfa/Box/PIDvax Study/Queries/Analysis Datasets/Final data and code for breakthrough paper/Publicly released data and code/visit_public.csv")


# Load the variant data from the CDC
variants <- read.csv("C:/Users/bustosfa/Box/PIDvax Study/Queries/Analysis Datasets/Final data and code for breakthrough paper/SARS-CoV-2_Variant_Proportions 2023-09-28.csv")


# Rename the public files to correspond with the name of the true datasets
behavior_raw <- behavior_raw_public
inf <- inf_public
new_elisa <- new_elisa_public
outcomes <- outcomes_public
pp_raw <- pp_raw_public
reinfection_raw <- reinfection_raw_public
sdat <- sdat_public
visit <- visit_public


# Remove unnecessary data objects
rm(behavior_raw_public, inf_public, new_elisa_public, outcomes_public,
   pp_raw_public, reinfection_raw_public, sdat_public, visit_public)



# **** PARAGRAPH 1: PARTICIPANT CHARACTERISTICS **** ----



# Sankey plot ----



# Put everything you need in this dataset
sankey1 <- new_elisa %>%
    filter(!is.na(vaccine_1_date)) %>%
    filter(vaccine_1 != "Janssen" | is.na(vaccine_1)) %>% 
    filter(vaccine_2 != "Janssen" | is.na(vaccine_2)) %>% 
    filter(vaccine_3 != "Janssen" | is.na(vaccine_3)) %>% 
    filter(vaccine_4 != "Janssen" | is.na(vaccine_4)) %>% 
    filter(vaccine_5 != "Janssen" | is.na(vaccine_5)) %>% 
    dplyr::select(fb_init_bt_date, vaccine_1_date:vaccine_6_date, CovidPos_1) %>%
    mutate_all(list(~ as.Date(., "%Y-%m-%d")))

check <- new_elisa %>%
    filter(!is.na(vaccine_1_date)) %>%
    filter(vaccine_1 != "Janssen" | is.na(vaccine_1)) %>% 
    filter(vaccine_2 != "Janssen" | is.na(vaccine_2)) %>% 
    filter(vaccine_3 != "Janssen" | is.na(vaccine_3)) %>% 
    filter(vaccine_4 != "Janssen" | is.na(vaccine_4)) %>% 
    filter(vaccine_5 != "Janssen" | is.na(vaccine_5)) %>% 
    dplyr::select(Subject.ID, vaccine_1:vaccine_6)

sankey1 <- cbind(check, sankey1)

sankey1 <- sankey1 %>%
    group_by(Subject.ID) %>%
    slice(1) %>%
    ungroup()


# Make new variables for the columns of the Sankey diagram
sankey1$Dose6 <- sankey1$Dose5 <- sankey1$Dose4 <- 
    sankey1$Dose3 <- sankey1$Dose2 <- sankey1$Dose1 <- sankey1$baseline <- NA
sankey1$Last <- NA


# Extract whether the vaccines were monovalent or bivalent
table2(sankey1$vaccine_1) # 0 NA
table2(sankey1$vaccine_2) # 1 NA
table2(sankey1$vaccine_3) # 44 NA
table2(sankey1$vaccine_4) # 120 NA
table2(sankey1$vaccine_5) # 192 NA
table2(sankey1$vaccine_6) # 243 NA
table2(sankey1$Last) # 269 NA


# Redo vaccine_N to focus on bivalent or monovalent
check <- sankey1 %>%
    dplyr::select(Subject.ID:vaccine_6)

check$vaccine_6v2 <- check$vaccine_5v2 <- check$vaccine_4v2 <- 
    check$vaccine_3v2 <- check$vaccine_2v2 <- check$vaccine_1v2 <- NA

names <- c("vaccine_1", "vaccine_2", "vaccine_3", "vaccine_4", "vaccine_5", "vaccine_6")

for(i in 1:nrow(check)){
    for(j in names){
        
        if (is.na(check[[j]][i])) 
        {check[[paste0(j, "v2")]][i] <- NA}
        
        else if(check[[j]][i] %in% c("Moderna (monovalent)", "Pfizer (monovalent)"))
        {check[[paste0(j, "v2")]][i] <- "Monovalent"}
        
        else if(check[[j]][i] %in% c("Moderna (bivalent)", "Pfizer (bivalent)"))
        {check[[paste0(j, "v2")]][i] <- "Bivalent"}
        
    }}

check$vaccine_1 <- check$vaccine_2 <- check$vaccine_3 <- check$vaccine_4 <- 
    check$vaccine_5 <- check$vaccine_6 <- NULL

names(check) <- c("id", names)

sankey1$vaccine_1 <- sankey1$vaccine_2 <- sankey1$vaccine_3 <- sankey1$vaccine_4 <- 
    sankey1$vaccine_5 <- sankey1$vaccine_6 <- NULL

names(sankey1)[1] <- "id"

sankey1 <- left_join(sankey1, check, by = "id")

rm(check)


# Define the status of each Dose variable
for(i in 1:nrow(sankey1)){
    
    # Baseline possibilities
    if(!is.na(sankey1$CovidPos_1[i]) & 
       !is.na(sankey1$vaccine_1_date[i]) & 
       (sankey1$CovidPos_1[i] < sankey1$vaccine_1_date[i]))
    {sankey1$baseline[i] <- "Infected\npre-dose 1"}
    
    else(sankey1$baseline[i] <- "Uninfected\npre-dose 1")
    
    # Dose 1 possibilities
    if(!is.na(sankey1$vaccine_1[i])) {sankey1$Dose1[i] <- sankey1$vaccine_1[i]}
    
    else(sankey1$Dose1[i] <- "None")
    
    
    # Dose 2 possibilities
    if(!is.na(sankey1$fb_init_bt_date[i]) & 
       !is.na(sankey1$vaccine_2_date[i]) &
       (sankey1$fb_init_bt_date[i] < sankey1$vaccine_2_date[i]))
    {sankey1$Dose2[i] <- sankey1$Dose3[i] <- sankey1$Dose4[i] <- 
        sankey1$Dose5[i] <- sankey1$Dose6[i] <- "Breakthrough"}
    
    else if(!is.na(sankey1$fb_init_bt_date[i]) &
            !is.na(sankey1$vaccine_1_date[i]) &
            is.na(sankey1$vaccine_2_date[i]) & 
            (sankey1$fb_init_bt_date[i] > sankey1$vaccine_1_date[i]))
    {sankey1$Dose2[i] <- sankey1$Dose3[i] <- sankey1$Dose4[i] <- 
        sankey1$Dose5[i] <- sankey1$Dose6[i] <- "Breakthrough"}  
    
    else if(!is.na(sankey1$vaccine_2[i])) {sankey1$Dose2[i] <- sankey1$vaccine_2[i]}
    
    else if(is.na(sankey1$Dose2[i])) {sankey1$Dose2[i] <- "None"}
    
    
    # Dose 3 possibilities
    if(!is.na(sankey1$fb_init_bt_date[i]) & 
       !is.na(sankey1$vaccine_3_date[i]) &
       (sankey1$fb_init_bt_date[i] < sankey1$vaccine_3_date[i]))
    {sankey1$Dose3[i] <- sankey1$Dose4[i] <- 
        sankey1$Dose5[i] <- sankey1$Dose6[i] <- "Breakthrough"}
    
    else if(!is.na(sankey1$fb_init_bt_date[i]) &
            !is.na(sankey1$vaccine_2_date[i]) &
            is.na(sankey1$vaccine_3_date[i]) & 
            (sankey1$fb_init_bt_date[i] > sankey1$vaccine_2_date[i]))
    {sankey1$Dose3[i] <- sankey1$Dose4[i] <- 
        sankey1$Dose5[i] <- sankey1$Dose6[i] <- "Breakthrough"}  
    
    else if(!is.na(sankey1$vaccine_3[i])) {sankey1$Dose3[i] <- sankey1$vaccine_3[i]}
    
    else if(is.na(sankey1$Dose3[i])) {sankey1$Dose3[i] <- "None"}
    
    
    # Dose 4 possibilities
    if(!is.na(sankey1$fb_init_bt_date[i]) &
       !is.na(sankey1$vaccine_4_date[i]) &
       (sankey1$fb_init_bt_date[i] < sankey1$vaccine_4_date[i]))
    {sankey1$Dose4[i] <- sankey1$Dose5[i] <- sankey1$Dose6[i] <- "Breakthrough"}
    
    else if(!is.na(sankey1$fb_init_bt_date[i]) &
            !is.na(sankey1$vaccine_3_date[i]) &
            is.na(sankey1$vaccine_4_date[i]) & 
            (sankey1$fb_init_bt_date[i] > sankey1$vaccine_3_date[i]))
    {sankey1$Dose4[i] <- sankey1$Dose5[i] <- sankey1$Dose6[i] <- "Breakthrough"}
    
    else if(!is.na(sankey1$vaccine_4[i])) {sankey1$Dose4[i] <- sankey1$vaccine_4[i]}
    
    else if(is.na(sankey1$Dose4[i])) {sankey1$Dose4[i] <- "None"}
    
    
    # Dose 5 possibilities
    if(!is.na(sankey1$fb_init_bt_date[i]) &
       !is.na(sankey1$vaccine_5_date[i]) &
       (sankey1$fb_init_bt_date[i] < sankey1$vaccine_5_date[i]))
    {sankey1$Dose5[i] <- sankey1$Dose6[i] <- "Breakthrough"}
    
    else if(!is.na(sankey1$fb_init_bt_date[i]) &
            !is.na(sankey1$vaccine_4_date[i]) &
            is.na(sankey1$vaccine_5_date[i]) & 
            (sankey1$fb_init_bt_date[i] > sankey1$vaccine_4_date[i]))
    {sankey1$Dose5[i] <- sankey1$Dose6[i] <- "Breakthrough"}
    
    else if(!is.na(sankey1$vaccine_5[i])) {sankey1$Dose5[i] <- sankey1$vaccine_5[i]}
    
    else if(is.na(sankey1$Dose5[i])) {sankey1$Dose5[i] <- "None"}
    
    
    # Dose 6 possibilities
    if(!is.na(sankey1$fb_init_bt_date[i]) &
       !is.na(sankey1$vaccine_6_date[i]) &
       (sankey1$fb_init_bt_date[i] < sankey1$vaccine_6_date[i]))
    {sankey1$Dose6[i] <- "Breakthrough"}
    
    else if(!is.na(sankey1$fb_init_bt_date[i]) &
            !is.na(sankey1$vaccine_5_date[i]) &
            is.na(sankey1$vaccine_6_date[i]) & 
            (sankey1$fb_init_bt_date[i] > sankey1$vaccine_5_date[i]))
    {sankey1$Dose6[i] <- "Breakthrough"}
    
    else if(!is.na(sankey1$vaccine_6[i])) {sankey1$Dose6[i] <- sankey1$vaccine_6[i]}
    
    else if(is.na(sankey1$Dose6[i])) {sankey1$Dose6[i] <- "None"}
}


# Deal with the End of study variable
for(i in 1:nrow(sankey1)){
    if(sankey1$Dose6[i] == "Breakthrough") 
    {sankey1$Last[i] <- "Breakthrough"}
    
    else if(!is.na(sankey1$fb_init_bt_date[i]) &
            !is.na(sankey1$vaccine_6_date[i]) &
            (sankey1$fb_init_bt_date[i] > sankey1$vaccine_6_date[i]))
    {sankey1$Last[i] <- "Breakthrough"}
    
    else(sankey1$Last[i] <- "None")
}


# Turning this on keeps the pre-dose 1 infected independent from everything else. 
# Turning this off just keeps them independent before dose 1 since they all got dose 1
# # Deal with the Infected pre-dose 1 participants
# for(i in 1:nrow(sankey1)){
#   if(sankey1$baseline[i] == "Infected\npre-dose 1")
#    {sankey1$Dose1[i] <- sankey1$Dose2[i] <- sankey1$Dose3[i] <-
#     sankey1$Dose4[i] <- sankey1$Dose5[i] <- sankey1$Dose6[i] <-
#     sankey1$Last[i] <- "Infected\npre-dose 1"}
# }


# Get it in the right format for the Sankey plot
sankey2 <- sankey1 %>%
    dplyr::select(Baseline = baseline, `Dose 1` = Dose1, `Dose 2` = Dose2, 
                  `Dose3` = Dose3, `Dose 4` = Dose4,
                  `Dose 5` = Dose5, `Dose 6` = Dose6, `End of study` = Last) %>%
    mutate_all(str_replace_all, "None", "Uninfected, no\nmore doses") %>%
    make_long(Baseline:`End of study`) %>%
    rename(`Vaccine and Breakthrough Infection Status` = x)


# Factor the columns
sankey2$node <- factor(sankey2$node, 
                       levels = rev(c("Uninfected\npre-dose 1", 
                                      "Uninfected, no\nmore doses", "Monovalent", 
                                      "Bivalent", "Breakthrough", 
                                      "Infected\npre-dose 1")))

sankey2$next_node <- factor(sankey2$next_node,
                            levels = rev(c("Uninfected\npre-dose 1", 
                                           "Uninfected, no\nmore doses", "Monovalent", 
                                           "Bivalent", "Breakthrough", 
                                           "Infected\npre-dose 1")))


# Make the full plot
ggplot(sankey2, aes(x = `Vaccine and Breakthrough Infection Status`, 
                    next_x = next_x, 
                    node = node,
                    next_node = next_node,
                    fill = node,
                    label = node)) +
    geom_sankey(flow.alpha = 0.5, node.color = 1) +
    scale_fill_manual(values = c("Monovalent" = "#d60270",
                                 "Bivalent" = "#7209b7",
                                 "Breakthrough" = "#000000",
                                 "Uninfected, no\nmore doses" = "#048B90",
                                 "Infected\npre-dose 1" = "#067FD0",
                                 "Uninfected\npre-dose 1" = "#ffc719")) +
    geom_sankey_label(size = 3.5, color = 1, fill = "white") +
    theme_sankey(base_size = 16) + 
    theme(legend.position = "none")



# Basic numbers ----



# Overall N
length(unique(new_elisa$Subject.ID)) # 282
length(unique(new_elisa[new_elisa$IDorHV == "PID",]$Subject.ID)) # 219
length(unique(new_elisa[new_elisa$IDorHV == "HV",]$Subject.ID)) # 63

# Age distribution
new_elisa %>%
    group_by(Subject.ID) %>%
    slice(1) %>%
    ungroup() %>%
    summarise(mean = mean(AGE))

range(new_elisa$AGE)

new_elisa %>%
    group_by(Subject.ID) %>%
    slice(1) %>%
    ungroup() %>%
    pull(AGE) %>%
    quantile(., probs = c(0, 0.25, 0.5, 0.75, 1))


# Sex distribution
new_elisa %>%
    group_by(Subject.ID) %>%
    slice(1) %>%
    ungroup() %>%
    summarise(count = n(),
              .by = GENDER)


# N participants with pre-vaccine infections
length(unique(new_elisa[new_elisa$fb_predose1_inf == 1,]$Subject.ID)) # 10


# N initial breakthrough infections before and after the Omicron era
check <- new_elisa %>%
    dplyr::select(Subject.ID, IDorHV, fb_init_bt_date) %>%
    group_by(Subject.ID) %>%
    slice(1) %>%
    ungroup

sum(check$fb_init_bt_date <= ymd("2021-12-15"), na.rm = T) # 5
table2(check[check$fb_init_bt_date <= ymd("2021-12-15"),]$IDorHV) # all IDP

sum(check$fb_init_bt_date > ymd("2021-12-15"), na.rm = T) # 111


# N pre-dose 1 infections among 116 breakthrough infections
length(unique(new_elisa[new_elisa$fb_bt == 1 & new_elisa$fb_predose1_inf == 1,]$Subject.ID))


# N secondary and tertiary breakthrough infections
inf_min <- inf %>%
    select(SUBJECT.ID., covid_date, infection_timing, vaccine_1_date)

inf_min <- inf_min %>%
    filter(infection_timing == "Breakthrough")

inf_min <- inf_min %>%
    group_by(SUBJECT.ID.) %>%
    arrange(covid_date) %>%
    mutate(bt_count = row_number()) %>%
    ungroup()

table(inf_min$bt_count)


# Remove useless data objects
rm(sankey1, sankey2, i, j, names, check, inf_min)



# Table 1 ----



# Filter to only the breakthroughs 
all_bt <- inf %>%
    filter(infection_timing == "Breakthrough")

all_bt <- all_bt %>%
    select(SUBJECT.ID., covid_date, infection_number) %>%
    group_by(SUBJECT.ID.) %>%
    arrange(SUBJECT.ID., covid_date) %>% # Select only the first of the breakthrough infections
    slice_head()

all_bt <- all_bt %>%
    select(!infection_number) %>%
    rename(date_first_bt_infection = covid_date)

BT_IDs <- unique(all_bt$SUBJECT.ID.)


# Create more variables for the table 
## Get the antibody titers after dose 2 for everyone
df4 <- new_elisa %>% 
    subset(time == "Post-dose 2") %>% 
    dplyr::select(Subject.ID, Collection.date, conc.ug_ml) %>% 
    rename(
        "postdose2.collection.date" = Collection.date,
        "postdose2.titer" = conc.ug_ml
    ) %>% 
    unique()

df4b <- new_elisa %>% 
    subset(time == "Post-dose 3") %>% 
    dplyr::select(Subject.ID, Collection.date, conc.ug_ml) %>% 
    rename(
        "postdose3.collection.date" = Collection.date,
        "postdose3.titer" = conc.ug_ml
    ) %>% 
    unique()

df4c <- new_elisa %>% 
    subset(time == "Post-dose 4") %>% 
    dplyr::select(Subject.ID, Collection.date, conc.ug_ml) %>% 
    rename(
        "postdose4.collection.date" = Collection.date,
        "postdose4.titer" = conc.ug_ml
    ) %>% 
    unique()


## Get the number of vaccinations for each person. 
df5a <- new_elisa %>% 
    dplyr::select(Subject.ID, vaccine_1_date:vaccine_6_date) %>% 
    unique() %>% 
    pivot_longer(vaccine_1_date:vaccine_6_date,
                 names_to = "vaccine_num", 
                 values_to = "date") %>% 
    drop_na(date) %>% 
    group_by(Subject.ID) %>% 
    nest() %>% 
    mutate(
        total.number.vaccines = map(data, 
                                    ~ length(.$vaccine_num))[[1]]
    ) %>% 
    dplyr::select(Subject.ID, total.number.vaccines)


# Join the information back into the elisa dataset
ELISA <- new_elisa %>%
    left_join(df4, by = "Subject.ID") %>% # antibody titer for postdose 2
    left_join(df4b, by = "Subject.ID") %>% # antibody titer for postdose 3
    left_join(df4c, by = "Subject.ID") %>% # antibody titer for postdose 4
    left_join(df5a, by = "Subject.ID")


# Make a participant-level dataset
part <- ELISA %>%
    group_by(Subject.ID) %>%
    slice_head()


# Filter ELISA data to only breakthrough infections using the selected IDs
part_b <- part %>%
    filter(Subject.ID %in% BT_IDs)

table(part_b$IDorHV) #116 participants with a breakthrough with ELISA data: 28 HV, 88 IDP


# Get the variables of interest
myVars <- c("AGE", 
            "GENDER", 
            "race", 
            "ethnicity", 
            "new_group",
            "IgRT",
            "immunosuppressants",
            "vaccine_2",
            "any_bivalent",
            "total.number.vaccines",
            "CovidPos1_variant",
            "number_doses_before_covid",
            "days_between_vax_covid",
            "bivalent_before_covid",
            "evusheld_before_covid",
            "postdose2.titer",  
            "postdose3.titer", 
            "postdose4.titer")

catVars <- c("GENDER", 
             "race", 
             "ethnicity",
             "new_group",
             "vaccine_2",
             "number_doses_before_covid",
             "any_bivalent",
             "bivalent_before_covid",
             "IgRT",
             "evusheld_before_covid",
             "immunosuppressants",
             "CovidPos1_variant",
             "total.number.vaccines")


# Generate the table of interest
tab1 <- CreateTableOne(vars = myVars,
                       strata = "IDorHV",
                       data = part_b,
                       factorVars = catVars,
                       test = T)

print(tab1, showAllLevels = TRUE, catDigits = 1, contDigits = 1) %>%
    kableone(longtable = T) %>%
    kable_styling(latex_options = c("hold_position", "repeat_header")) %>% 
    column_spec(1:6, width = c("4cm", "3cm", "2cm","2cm","1.25cm", "0.5cm"))

tab1_excel <- print(tab1, showAllLevels = TRUE, catDigits = 1, contDigits = 1, printToggle = FALSE)


# Remove useless data objects
rm(tab1, BT_IDs, catVars, myVars, part, part_b, ELISA, all_bt,
   df4, df4b, df4c, df5a, tab1_excel)



# Swimmer plot ----



# Process the variants data
names(variants)[1] <- "region"

variants2 <- variants %>%
    filter(region == "USA") 

variants2$week_ending <- format(as.Date(as.POSIXct(variants2$week_ending, 
                                                   format = '%m/%d/%Y %H:%M'), format = '%m/%d/%Y'))

variants2$published_date <- format(as.Date(as.POSIXct(variants2$published_date, 
                                                      format = '%m/%d/%Y %H:%M'), format = '%m/%d/%Y'))


# Take the most recent published date
variants2$published_date <- gsub("\\ .*", "", variants2$published_date)
variants2$published_date <- ymd(variants2$published_date)
variants2$published_date <- as.Date(format(ymd(variants2$published_date), format = "%Y-%m-%d"))

variants2$week_ending <- gsub("\\ .*", "", variants2$week_ending)
variants2$week_ending <- ymd(variants2$week_ending)
variants2$week_ending <- as.Date(format(ymd(variants2$week_ending), format = "%Y-%m-%d"))


# Look at the dates
table(variants2$published_date, useNA = "always")

variants2 <- variants2 %>%
    group_by(week_ending) %>%
    slice_max(published_date, n = 1) %>%
    ungroup()


# Now take only variants of interest
variants2 <- variants2 %>%
    dplyr::select(week_ending, variant, share)

table(variants2$variant)


# Look at frequencies
sort <- variants2 %>%
    group_by(variant) %>%
    tally()


# Change the format of the Date variable; otherwise an error occurs below
variants2$week_ending <- as.character(variants2$week_ending)


# Simplify the variants in the CDC data
variants2[variants2 == "XBB"] <- "Omicron XBB"
variants2[variants2 == "XBB.1.16"] <- "Omicron XBB"
variants2[variants2 == "FD.2"] <- "Omicron XBB"
variants2[variants2 == "XBB.1.5"] <- "Omicron XBB"
variants2[variants2 == "XBB.1.9.1"] <- "Omicron XBB"
variants2[variants2 == "XBB.1.9.2"] <- "Omicron XBB"
variants2[variants2 == "XBB.1.16.1"] <- "Omicron XBB"
variants2[variants2 == "XBB.1.16.11"] <- "Omicron XBB"
variants2[variants2 == "XBB.1.16.6"] <- "Omicron XBB"
variants2[variants2 == "XBB.1.5.10"] <- "Omicron XBB"
variants2[variants2 == "XBB.1.5.59"] <- "Omicron XBB"
variants2[variants2 == "XBB.1.16.6"] <- "Omicron XBB"
variants2[variants2 == "XBB.1.5.68"] <- "Omicron XBB"
variants2[variants2 == "XBB.1.5.70"] <- "Omicron XBB"
variants2[variants2 == "XBB.1.5.72"] <- "Omicron XBB"
variants2[variants2 == "XBB.2.3"] <- "Omicron XBB"
variants2[variants2 == "XBB.2.3"] <- "Omicron XBB"
variants2[variants2 == "XBB.2.3.8"] <- "Omicron XBB"
variants2[variants2 == "XBB.1.42.2"] <- "Omicron XBB"
variants2[variants2 == "XBB.1.5.1"] <- "Omicron XBB"

variants2[variants2 == "EG.5"] <- "Omicron EG.5/6"
variants2[variants2 == "EG.6.1"] <- "Omicron EG.5/6"

variants2[variants2 == "B.1.1.529"] <- "Omicron BA.1"
variants2[variants2 == "BA.1.1"] <- "Omicron BA.1"

variants2[variants2 == "BA.2.12.1"] <- "Omicron BA.2"
variants2[variants2 == "BA.2"] <- "Omicron BA.2"
variants2[variants2 == "BA.2.75"] <- "Omicron BA.2"
variants2[variants2 == "BA.2.75.2"] <- "Omicron BA.2"
variants2[variants2 == "BN.1"] <- "Omicron BA.2"
variants2[variants2 == "CH.1.1"] <- "Omicron BA.2"

variants2[variants2 == "BA.4"] <- "Omicron BA.4/5"
variants2[variants2 == "BA.4.6"] <- "Omicron BA.4/5"
variants2[variants2 == "BA.5"] <- "Omicron BA.4/5"
variants2[variants2 == "BF.7"] <- "Omicron BA.4/5"
variants2[variants2 == "BQ.1"] <- "Omicron BA.4/5"
variants2[variants2 == "BQ.1.1"] <- "Omicron BA.4/5"
variants2[variants2 == "BA.5.2.6"] <- "Omicron BA.4/5"
variants2[variants2 == "BF.11"] <- "Omicron BA.4/5"

variants2[variants2 == "B.1.617"] <- "Delta"
variants2[variants2 == "B.1.617.2"] <- "Delta"
variants2[variants2 == "AY.1"] <- "Delta"
variants2[variants2 == "AY.2"] <- "Delta"
variants2[variants2 == "AY.3"] <- "Delta"

variants2[variants2 == "B.1.1.7"] <- "Alpha"
variants2[variants2 == "P.1"] <- "Other"
variants2[variants2 == "B.1.621"] <- "Other"
variants2[variants2 == "B.1.526"] <- "Other"
variants2[variants2 == "B.1.621.1"] <- "Other"
variants2[variants2 == "B.1.351"] <- "Other"
variants2[variants2 == "B.1.617.1"] <- "Other"
variants2[variants2 == "B.1.525"] <- "Other"
variants2[variants2 == "B.1.526.2"] <- "Other"
variants2[variants2 == "B.1.526.1"] <- "Other"
variants2[variants2 == "B.1.429"] <- "Other"
variants2[variants2 == "B.1.427"] <- "Other"
variants2[variants2 == "P.2"] <- "Other"

variants2[variants2 == "B.1.637"] <- "Other"
variants2[variants2 == "B.1.1"] <- "Other"
variants2[variants2 == "B.1.1.194"] <- "Other"
variants2[variants2 == "B.1.243"] <- "Other"
variants2[variants2 == "B.1.234"] <- "Other"
variants2[variants2 == "R.1"] <- "Other"
variants2[variants2 == "B.1.628"] <- "Other"
variants2[variants2 == "A.2.5"] <- "Other"
variants2[variants2 == "B.1.626"] <- "Other"
variants2[variants2 == "B.1.617.3"] <- "Other"
variants2[variants2 == "B.1"] <- "Other"
variants2[variants2 == "B.1.1.519"] <- "Other"
variants2[variants2 == "B.1.2"] <- "Other"
variants2[variants2 == "B.1.596"] <- "Other"

variants2[variants2 == "EU.1.1"] <- "Other"
variants2[variants2 == "FD.1.1"] <- "Other"
variants2[variants2 == "FE.1.1"] <- "Other"
variants2[variants2 == "FL.1.5.1"] <- "Other"
variants2[variants2 == "GE.1"] <- "Other"
variants2[variants2 == "HV.1"] <- "Other"


# look at frequencies again
sort2 <- variants2 %>%
    group_by(variant) %>%
    tally()


# New dataset with percentages
variants3 <- variants2 %>%
    group_by(week_ending, variant) %>%
    summarise(percent = sum(share)) %>%
    ungroup()

min(variants3$week_ending)

unique(variants3$variant)

variants3$variant <- factor(variants3$variant,
                            levels = c("Alpha", "Delta", "Omicron BA.1", "Omicron BA.2", 
                                       "Omicron BA.4/5", "Omicron XBB", "Omicron EG.5/6", "Other"))

variants3 <- arrange(variants3, week_ending, variant)


# Create a new row with values
new_row <- data.frame(
    week_ending = "2021-01-09",
    variant = "Other",
    percent = 100.00
)

new_row_2 <- data.frame(
    week_ending = "2021-01-23",
    variant = "Other",
    percent = 100.00
)


# Add the new row to the dataset
variants3 <- rbind(variants3, new_row, new_row_2)


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
names(new_elisa)
table(new_elisa$breakthru)


# Get the dates of the breakthru cases
bt <- new_elisa %>%
    filter(breakthru == 1) %>%
    # filter(evusheld_before_sample == 0) %>%
    #  filter(breakthru != "COVID pre-dose 1") %>%
    distinct(Subject.ID, .keep_all = T)

bt$pos_before <- ifelse(bt$COVID_pos_before_sample == 1 |
                            bt$N_pos_before_sample == 1, 1, 0)

bt$Healthy_binary <- ifelse(bt$new_group == "Healthy volunteer", 1, 0)


# Dataset for the breakthroughs
newvaxdf <- bt %>%
    filter(breakthru == 1) %>%
    dplyr::select(Subject.ID, vaccine_2_date, vaccine_3_date, vaccine_4_date, 
                  vaccine_5_date, vaccine_6_date, CovidPos_1)


# Make sure dates are formatted
newvaxdf$vaccine_2_date <- as.Date(newvaxdf$vaccine_2_date)
newvaxdf$vaccine_3_date <- as.Date(newvaxdf$vaccine_3_date)
newvaxdf$vaccine_4_date <- as.Date(newvaxdf$vaccine_4_date)
newvaxdf$vaccine_5_date <- as.Date(newvaxdf$vaccine_5_date)
newvaxdf$vaccine_6_date <- as.Date(newvaxdf$vaccine_6_date)
newvaxdf$CovidPos_1 <- as.Date(newvaxdf$CovidPos_1)


# Drop the one person who only has gotten one vaccination
newvaxdf <- newvaxdf %>%
    filter(Subject.ID != "E0080")


# Drop the people who did not finish primary series or got a covidpos 1 beforehand
newvaxdf <- newvaxdf %>%
    filter(!CovidPos_1 < vaccine_2_date)


# Identify the date of last vaccine dose before brekthrough
newvaxdf <- newvaxdf %>%
    group_by(Subject.ID) %>%
    mutate(closest_date = case_when(
        vaccine_2_date < CovidPos_1 & vaccine_3_date > CovidPos_1 ~ vaccine_2_date,
        vaccine_2_date < CovidPos_1 & is.na(vaccine_3_date) ~ vaccine_2_date,
        vaccine_3_date < CovidPos_1 & vaccine_4_date > CovidPos_1 ~ vaccine_3_date,
        vaccine_3_date < CovidPos_1 & is.na(vaccine_4_date) ~ vaccine_3_date,
        vaccine_4_date < CovidPos_1 & vaccine_5_date > CovidPos_1 ~ vaccine_4_date,
        vaccine_4_date < CovidPos_1 & is.na(vaccine_5_date) ~ vaccine_4_date,
        vaccine_5_date < CovidPos_1 & vaccine_6_date > CovidPos_1 ~ vaccine_5_date,
        vaccine_5_date < CovidPos_1 & is.na(vaccine_6_date) ~ vaccine_5_date,
        vaccine_6_date < CovidPos_1 ~ vaccine_6_date,
        TRUE ~ NA)
    )


# Get the Saturday for the second vaccine date
newvaxdf$saturday_vax2 <- NA
for(i in 1:nrow(newvaxdf)){
    for(j in 1:length(saturdays)){
        if(newvaxdf$vaccine_2_date[i] > saturdays[j]){
            newvaxdf$saturday_vax2[i] <- format(as.Date(saturdays[j + 1], format = '%m/%d/%Y'))
        }
    }
}


# Get the Saturday for the date of infection
newvaxdf$saturday_covid <- NA
for(i in 1:nrow(newvaxdf)){
    for(j in 1:length(saturdays)){
        if(newvaxdf$CovidPos_1[i] > saturdays[j]){
            newvaxdf$saturday_covid[i] <- format(as.Date(saturdays[j + 1], format = '%m/%d/%Y'))
        }
    }
}

# Get the Saturday for the most recent vaccine date
newvaxdf$saturday_closest_date <- NA
for(i in 1:nrow(newvaxdf)){
    for(j in 1:length(saturdays)){
        if(newvaxdf$closest_date[i] > saturdays[j]){
            newvaxdf$saturday_closest_date[i] <- format(as.Date(saturdays[j + 1], format = '%m/%d/%Y'))
        }
    }
}
# a few need to be manually edited because of how CDC didn't always post data weekly, sometimes it was bi-weekly
# this happens in May 2023 then it goes bi-weekly
newvaxdf <- newvaxdf %>%
    mutate(saturday_covid = case_when(Subject.ID == "H0038" ~ "2023-07-22",
                                      TRUE ~ saturday_covid))


# Simplify the dataset
newvaxdf2 <- newvaxdf %>%
    dplyr::select(id = Subject.ID, saturday_vax2, saturday_closest_date, saturday_covid)


# Format the dataset appropriately
newvaxdf2 <- newvaxdf2 %>%
    arrange(saturday_covid)

newvaxdf2$id <- factor(newvaxdf2$id, levels = newvaxdf2$id)

newvaxdf2$percent = seq(1, 0, length.out = nrow(newvaxdf) + 2)[c(F, rep(T, nrow(newvaxdf)), F)]


# Add in variable for HV or not
controlvar <- new_elisa %>%
    filter(Subject.ID != "E0080") %>%
    dplyr::select(id = Subject.ID, IDorHV) %>%
    filter(id %in% newvaxdf2$id) %>%
    distinct()

newvaxdf2 <- left_join(newvaxdf2, controlvar, by = "id")
newvaxdf2 <- as.data.frame(newvaxdf2)


# Make a line dataset
newvaxdf2_line <- newvaxdf2

newvaxdf2_line <- gather(newvaxdf2_line, saturday, week_ending, 
                         saturday_vax2, saturday_closest_date, saturday_covid, factor_key = TRUE)

newvaxdf2_line$id <- factor(newvaxdf2_line$id, levels = newvaxdf2$id)

newvaxdf2_line <- newvaxdf2_line %>%
    arrange(id, week_ending)


# Arrange data appropriately
var_highest <- variants2 %>%
    arrange(variant, week_ending) %>%
    filter(share > 50)

highest_date <- variants3 %>%
    group_by(week_ending) %>%
    mutate(most_comm_variant = variant[which.max(percent)])


# Now make the final version of the variant_plot
# mrz version
swimmer_simple <- ggplot(variants3[variants3$week_ending <= "2023-08-08",], 
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
    scale_fill_manual(name = "SARS-CoV-2 \nvariant", #B2 adds alpha of 0.7, see alpha("red", .3)
                      values = c("Alpha" = "#dadaeb",
                                 "Delta" = "#bcbddc",
                                 "Omicron BA.1" = "#fff7f3",
                                 "Omicron BA.2" = "#fde0dd",
                                 "Omicron BA.4/5" = "#fcc5c0",
                                 "Omicron XBB" = "#fa9fb5",
                                 "Omicron EG.5/6" = "#f768a1",
                                 "Other" = "#f0f0f0")) + 
    scale_x_discrete(breaks = every_nth(n = 4),
                     labels = labels_vec) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.title.y.right = element_text(angle = 90)) + 
    geom_vline(xintercept = "2021-05-01", color = "black", 
               linewidth = 2.25, linetype = "dotted") + 
    guides(color = guide_legend(override.aes = list(linetype = rep(1, 14), color = "black"))) + 
    geom_line(inherit.aes = F, data = newvaxdf2_line, # Add swimmer lines
              aes(group = id, y = percent, x = week_ending),
              linewidth = 0.3, color = "black") + 
    geom_point(inherit.aes = F, data = newvaxdf2, # Add dose 2
               aes(group = id, y = percent, x = saturday_vax2),
               size = 2.5, color = "black", shape = 1) + 
    geom_point(inherit.aes = F, data = newvaxdf2, # Add closest date
               aes(group = id, y = percent, x = saturday_closest_date),
               size = 2.5, color = "black", shape = 10) + 
    geom_point(inherit.aes = F, # Add infection for confirmed HV
               data = newvaxdf2[newvaxdf2$IDorHV == "HV",],
               aes(group = id, y = percent, x = saturday_covid),
               size = 3, color = "#41b6c4", fill = "#41b6c4", shape = 17,
               show.legend = FALSE, alpha = 1) +
    geom_point(inherit.aes = F, # Add infection for confirmed PID
               data = newvaxdf2[newvaxdf2$IDorHV == "PID",],
               aes(group = id, y = percent, x = saturday_covid),
               size = 3, color = "#7a0177", fill = "#7a0177", shape = 17, 
               show.legend = FALSE, alpha = 1)

swimmer_simple


# Save it
ggsave(filename = "swimmer_simple_2024-01-19.png", width = 2.14, height = 1.51, dpi = 300)


# Remove useless data objects
rm(variants, variants2, sort2, controlvar, highest_date, labels3, new_row, new_row_2,
   newvaxdf, newvaxdf2, newvaxdf2_line, sort, swimmer_simple, var_highest, variants3, 
   i, j, labels, labels_vec, saturdays, bt)



# **** PARAGRAPH 2: INCIDENCE OF INITIAL BREAKTHROUGH INFECTIONS **** ----



# N breakthroughs among IDP and HV ----



# Count them up
length(unique(new_elisa[new_elisa$IDorHV == "PID" & new_elisa$fb_bt == 1,]$Subject.ID)) # 88
length(unique(new_elisa[new_elisa$IDorHV == "HV" & new_elisa$fb_bt == 1,]$Subject.ID)) # 28



# Incidence rate: IDP vs HV, overall ----



# Simplify the dataset to what I need
bt <- new_elisa %>%
  dplyr::select(Subject.ID, IDorHV, new_group, fb_bt, end,
                vaccine_2_date, vaccine_1_date, fb_init_bt_date, vaccine_1) %>%
  group_by(Subject.ID) %>%
  slice(1) %>%
  ungroup


# Calculate person-time at risk
bt %<>%
  mutate(start_date = case_when(vaccine_1 == "Janssen" ~ ymd(vaccine_1_date),
                                TRUE ~ ymd(vaccine_2_date))) %>%
  mutate(end_date = case_when(!is.na(fb_init_bt_date) ~ ymd(fb_init_bt_date),
                              is.na(fb_init_bt_date) ~ ymd(end))) %>%
  
  drop_na(start_date, end_date) %>% # Drop H0060 with only 1 vaccine
  mutate(persondays = as.numeric(ymd(end_date) - ymd(start_date))) %>%
  mutate(personyears = as.numeric(persondays / 365))


# Calculate IRs
IR <- bt %>%
  dplyr::select(Subject.ID, IDorHV, fb_bt, persondays, personyears) %>%
  summarise(cases = sum(fb_bt),
            PT_d = sum(persondays),
            PT_y = sum(personyears),
            IR_y = round(cases / PT_y, 2),
            IR_y100 = round(100 * cases / PT_y, 2),
            .by = IDorHV)
IR


# Add in the CIs, numbers from p. 434 of Epidemiology Beyond The Basics, 3rd ed.
IR$`95% UL` <- IR$`95% LL` <- NA

IR[IR$IDorHV == "PID",]$`95% LL` <- round(100 * IR[IR$IDorHV == "PID",]$cases * 
                                  (0.789 + 0.8 * (0.809 - 0.789)) / IR[IR$IDorHV == "PID",]$PT_y, 2)
IR[IR$IDorHV == "PID",]$`95% UL` <- round(100 * IR[IR$IDorHV == "PID",]$cases * 
                                  (1.250 - 0.2 * (1.250 - 1.240)) / IR[IR$IDorHV == "PID",]$PT_y, 2)

IR[IR$IDorHV == "HV",]$`95% LL` <- round(100 * IR[IR$IDorHV == "HV",]$cases * 
                                           0.665 / IR[IR$IDorHV == "HV",]$PT_y, 2)
IR[IR$IDorHV == "HV",]$`95% UL` <- round(100 * IR[IR$IDorHV == "HV",]$cases * 
                                           1.450 / IR[IR$IDorHV == "HV",]$PT_y, 2)


# Calculate IRR
IRR <- round((IR[IR$IDorHV == "PID",]$cases / IR[IR$IDorHV == "PID",]$PT_y) / 
               (IR[IR$IDorHV == "HV",]$cases / IR[IR$IDorHV == "HV",]$PT_y), 2)


# Calculate IRR 95% CIs
Phat <- IR[IR$IDorHV == "PID",]$cases / sum(IR$cases)

PL <- Phat - (qnorm(0.975) * sqrt(Phat * (1 - Phat) / sum(IR$cases)))
PU <- Phat + (qnorm(0.975) * sqrt(Phat * (1 - Phat) / sum(IR$cases)))

RRL <- round((PL / (1 - PL)) * (IR[IR$IDorHV == "HV",]$PT_y / IR[IR$IDorHV == "PID",]$PT_y), 2)
RRU <- round((PU / (1 - PU)) * (IR[IR$IDorHV == "HV",]$PT_y / IR[IR$IDorHV == "PID",]$PT_y), 2)

IRR; RRL; RRU


# Get pvalue
E1 <- sum(IR$cases) * IR[IR$IDorHV == "PID",]$PT_y / sum(IR$PT_y)
E2 <- sum(IR$cases) * IR[IR$IDorHV == "HV",]$PT_y / sum(IR$PT_y)

ts <- ((IR[IR$IDorHV == "PID",]$cases - E1)^2 / E1) + ((IR[IR$IDorHV == "HV",]$cases - E2)^2 / E2)

pvalue <- pchisq(q = ts, df = 1, lower.tail = F)
pvalue


# Put everything together
IR2 <- IR
IR2$pvalue <- IR2$RRU <- IR2$RRL <- IR2$IRR <- NA
IR2$IRR[1] <- IRR
IR2$RRL[1] <- RRL
IR2$RRU[1] <- RRU
IR2$pvalue[1] <- round(pvalue, 3)

IR2

rm(Phat, PL, PU, ts, E1, E2, IRR, RRL, RRU, pvalue, IR, IR2)



# Incidence rate: pre- vs post-Omicron, overall ----



# Calculate cases and person-time at risk in the pre-Omicron and post-Omicron era
sum(bt$fb_init_bt_date <= ymd("2021-12-15"), na.rm = T) # Only 5 pre-Omicron era breakthroughs
table(bt[bt$fb_init_bt_date <= ymd("2021-12-15"),]$IDorHV) # all in IDP

bt_pre <- bt %>%
    mutate(cases_pre = case_when(fb_bt == 1 & 
                                     !is.na(fb_init_bt_date) & 
                                     bt$fb_init_bt_date <= ymd("2021-12-15") ~ 1,
                                 TRUE ~ 0)) %>% 
    mutate(start_date_pre = case_when(vaccine_1 == "Janssen" ~ ymd(vaccine_1_date),
                                      TRUE ~ ymd(vaccine_2_date))) %>%
    mutate(end_date_pre = 
               case_when(!is.na(fb_init_bt_date) & 
                             (ymd(fb_init_bt_date) <= ymd("2021-12-15")) ~ ymd(fb_init_bt_date),
                         !is.na(fb_init_bt_date) & 
                             (ymd(fb_init_bt_date) > ymd("2021-12-15")) ~ ymd("2021-12-15")))


# Deal with end date without invoking min(end, ymd("2021-12-15")) because that will get the end from the entire variable, not the minimum within the person. The rest is just formatting to make sure the dates stay formatted as date throughout
bt_pre$end_date_pre <- ifelse(is.na(bt_pre$fb_init_bt_date) & bt_pre$end <= ymd("2021-12-15"), 
                              as.Date(bt_pre$end, format = "%Y-%m-%d"), ymd(bt_pre$end_date_pre))
bt_pre$end_date_pre <- as.Date(bt_pre$end_date_pre, format = "%Y-%m-%d")
table2(bt_pre$end_date_pre)

bt_pre$end_date_pre <- ifelse(is.na(bt_pre$fb_init_bt_date) & bt_pre$end > ymd("2021-12-15"), 
                              ymd("2021-12-15"), bt_pre$end_date_pre)
bt_pre$end_date_pre <- as.Date(bt_pre$end_date_pre, format = "%Y-%m-%d")
table2(bt_pre$end_date_pre)


# Continue as before
bt_pre <- bt_pre %>%   
    drop_na(start_date_pre, end_date_pre) %>% # Remove any NAs
    filter(Subject.ID != "H0060") %>% # Drop H0060 with only 1 vaccine
    mutate(persondays_pre = as.numeric(ymd(end_date_pre) - ymd(start_date_pre))) %>%
    mutate(personyears_pre = as.numeric(persondays_pre / 365)) %>%
    filter(persondays_pre > 0) # Remove people who didn't conclude vax sequence till after Omicron era
    
bt_post <- bt %>%
    mutate(cases_pre = case_when(fb_bt == 1 & 
                                     !is.na(fb_init_bt_date) & 
                                     bt$fb_init_bt_date <= ymd("2021-12-15") ~ 1,
                                 TRUE ~ 0)) %>% 
    mutate(cases_post = case_when(fb_bt == 1 & 
                                      !is.na(fb_init_bt_date) & 
                                      bt$fb_init_bt_date > ymd("2021-12-15") ~ 1,
                                  TRUE ~ 0)) %>% 
    mutate(start_date_post = 
               case_when(vaccine_1 == "Janssen" & 
                             cases_pre == 0 &
                             ymd(vaccine_1_date) <= ymd("2021-12-15") ~ ymd("2021-12-15"),
                         vaccine_1 == "Janssen" & 
                             cases_pre == 0 &
                             ymd(vaccine_1_date) > ymd("2021-12-15") ~ ymd(vaccine_1_date),
                         vaccine_1 != "Janssen" & 
                             cases_pre == 0 &
                             ymd(vaccine_2_date) <= ymd("2021-12-15") ~ ymd("2021-12-15"),
                         vaccine_1 != "Janssen" & 
                             cases_pre == 0 &
                             ymd(vaccine_2_date) > ymd("2021-12-15") ~ ymd(vaccine_2_date),
                         TRUE ~ NA)) %>%
    mutate(end_date_post = 
               case_when(!is.na(fb_init_bt_date) & 
                             (ymd(fb_init_bt_date) <= ymd("2021-12-15")) ~ NA,
                         !is.na(fb_init_bt_date) & 
                             (ymd(fb_init_bt_date) > ymd("2021-12-15")) ~ ymd(fb_init_bt_date),
                         is.na(fb_init_bt_date) ~ ymd(end))) %>%
    
    drop_na(start_date_post, end_date_post) %>% # Drop ineligibles from post-Omicron consideration
    filter(Subject.ID != "H0060") %>% # Drop H0060 with only 1 vaccine
    mutate(persondays_post = as.numeric(ymd(end_date_post) - ymd(start_date_post))) %>%
    mutate(personyears_post = as.numeric(persondays_post / 365)) 
    

# Calculate the IRs
IR_pre <- bt_pre %>%
    dplyr::select(Subject.ID, IDorHV, cases_pre, persondays_pre, personyears_pre) %>%
    summarise(cases = sum(cases_pre),
              PT_d = sum(persondays_pre),
              PT_y = sum(personyears_pre),
              IR_y = round(cases / PT_y, 2),
              IR_y100 = round(100 * cases / PT_y, 2)) %>%
    mutate(era = "pre")
IR_pre

IR_post <- bt_post %>%
    dplyr::select(Subject.ID, IDorHV, cases_post, persondays_post, personyears_post) %>%
    summarise(cases = sum(cases_post),
              PT_d = sum(persondays_post),
              PT_y = sum(personyears_post),
              IR_y = round(cases / PT_y, 2),
              IR_y100 = round(100 * cases / PT_y, 2)) %>%
    mutate(era = "post")
IR_post


# Combine the incidence information across the two eras
IR_era <- rbind(IR_pre, IR_post)
IR_era


# Add in the CIs, numbers from p. 434 of Epidemiology Beyond The Basics, 3rd ed.
IR_era$`95% UL` <- IR_era$`95% LL` <- NA

IR_era[IR_era$era == "post",]$`95% LL` <- 
    round(100 * IR_era[IR_era$era == "post",]$cases *  (0.818 + (11/20) * (0.833 - 0.818)) /
              IR_era[IR_era$era == "post",]$PT_y, 2)

IR_era[IR_era$era == "post",]$`95% UL` <- 
    round(100 * IR_era[IR_era$era == "post",]$cases * (1.22 - (9/20) * (1.22 - 1.2)) /
              IR_era[IR_era$era == "post",]$PT_y, 2)

IR_era[IR_era$era == "pre",]$`95% LL` <- round(100 * IR_era[IR_era$era == "pre",]$cases * 
                                             0.324 / IR_era[IR_era$era == "pre",]$PT_y, 2)
IR_era[IR_era$era == "pre",]$`95% UL` <- round(100 * IR_era[IR_era$era == "pre",]$cases * 
                                             2.33 / IR_era[IR_era$era == "pre",]$PT_y, 2)


# Calculate IRR (PID -> post, HV -> pre)
IRR <- round((IR_era[IR_era$era == "post",]$cases / IR_era[IR_era$era == "post",]$PT_y) / 
                 (IR_era[IR_era$era == "pre",]$cases / IR_era[IR_era$era == "pre",]$PT_y), 2)


# Calculate IRR 95% CIs
Phat <- IR_era[IR_era$era == "post",]$cases / sum(IR_era$cases)

PL <- Phat - (qnorm(0.975) * sqrt(Phat * (1 - Phat) / sum(IR_era$cases)))
PU <- Phat + (qnorm(0.975) * sqrt(Phat * (1 - Phat) / sum(IR_era$cases)))

RRL <- round((PL / (1 - PL)) * (IR_era[IR_era$era == "pre",]$PT_y / 
                                    IR_era[IR_era$era == "post",]$PT_y), 2)
RRU <- round((PU / (1 - PU)) * (IR_era[IR_era$era == "pre",]$PT_y / 
                                    IR_era[IR_era$era == "post",]$PT_y), 2)

IRR; RRL; RRU


# Get pvalue
E1 <- sum(IR_era$cases) * IR_era[IR_era$era == "post",]$PT_y / sum(IR_era$PT_y)
E2 <- sum(IR_era$cases) * IR_era[IR_era$era == "pre",]$PT_y / sum(IR_era$PT_y)

ts <- ((IR_era[IR_era$era == "post",]$cases - E1)^2 / E1) + 
      ((IR_era[IR_era$era == "pre",]$cases - E2)^2 / E2)

pvalue <- pchisq(q = ts, df = 1, lower.tail = F)
pvalue


# Put everything together
IR2 <- IR_era
IR2$pvalue <- IR2$RRU <- IR2$RRL <- IR2$IRR <- NA
IR2$IRR[1] <- IRR
IR2$RRL[1] <- RRL
IR2$RRU[1] <- RRU
IR2$pvalue[1] <- round(pvalue, 3)

IR2

rm(Phat, PL, PU, ts, E1, E2, IRR, RRL, RRU, pvalue, IR_era, IR2, bt_pre)



# Incidence rate: IDP vs HV, during Omicron ----



# Get the basic information that we need
IR_post2 <- bt_post %>%
    dplyr::select(Subject.ID, IDorHV, cases_post, persondays_post, personyears_post) %>%
    summarise(cases = sum(cases_post),
              PT_d = sum(persondays_post),
              PT_y = sum(personyears_post),
              IR_y = round(cases / PT_y, 2),
              IR_y100 = round(100 * cases / PT_y, 2),
              .by = IDorHV)
IR_post2

# Add in the CIs, numbers from p. 434 of Epidemiology Beyond The Basics, 3rd ed.
IR_post2$`95% UL` <- IR_post2$`95% LL` <- NA

IR_post2[IR_post2$IDorHV == "PID",]$`95% LL` <- 
    round(100 * IR_post2[IR_post2$IDorHV == "PID",]$cases * (0.798 + 0.3 * (0.809 - 0.798)) / 
              IR_post2[IR_post2$IDorHV == "PID",]$PT_y, 2)
IR_post2[IR_post2$IDorHV == "PID",]$`95% UL` <- 
    round(100 * IR_post2[IR_post2$IDorHV == "PID",]$cases * (1.24 - 0.7 * (1.25 - 1.24)) / 
              IR_post2[IR_post2$IDorHV == "PID",]$PT_y, 2)

IR_post2[IR_post2$IDorHV == "HV",]$`95% LL` <- 
    round(100 * IR_post2[IR_post2$IDorHV == "HV",]$cases * 0.665 / 
              IR_post2[IR_post2$IDorHV == "HV",]$PT_y, 2)
IR_post2[IR_post2$IDorHV == "HV",]$`95% UL` <- 
    round(100 * IR_post2[IR_post2$IDorHV == "HV",]$cases * 1.450 / 
              IR_post2[IR_post2$IDorHV == "HV",]$PT_y, 2)

# Calculate IRR
IRR <- 
    round((IR_post2[IR_post2$IDorHV == "PID",]$cases / IR_post2[IR_post2$IDorHV == "PID",]$PT_y) / 
          (IR_post2[IR_post2$IDorHV == "HV",]$cases / IR_post2[IR_post2$IDorHV == "HV",]$PT_y), 2)

# Calculate IRR 95% CIs
Phat <- IR_post2[IR_post2$IDorHV == "PID",]$cases / sum(IR_post2$cases)

PL <- Phat - (qnorm(0.975) * sqrt(Phat * (1 - Phat) / sum(IR_post2$cases)))
PU <- Phat + (qnorm(0.975) * sqrt(Phat * (1 - Phat) / sum(IR_post2$cases)))

RRL <- round((PL / (1 - PL)) * (IR_post2[IR_post2$IDorHV == "HV",]$PT_y / 
                                    IR_post2[IR_post2$IDorHV == "PID",]$PT_y), 2)
RRU <- round((PU / (1 - PU)) * (IR_post2[IR_post2$IDorHV == "HV",]$PT_y /
                                    IR_post2[IR_post2$IDorHV == "PID",]$PT_y), 2)

IRR; RRL; RRU

# Get pvalue
E1 <- sum(IR_post2$cases) * IR_post2[IR_post2$IDorHV == "PID",]$PT_y / sum(IR_post2$PT_y)
E2 <- sum(IR_post2$cases) * IR_post2[IR_post2$IDorHV == "HV",]$PT_y / sum(IR_post2$PT_y)

ts <- ((IR_post2[IR_post2$IDorHV == "PID",]$cases - E1)^2 / E1) + 
      ((IR_post2[IR_post2$IDorHV == "HV",]$cases - E2)^2 / E2)

pvalue <- pchisq(q = ts, df = 1, lower.tail = F)
pvalue


# Put everything together
IR_post3 <- IR_post2
IR_post3$pvalue <- IR_post3$RRU <- IR_post3$RRL <- IR_post3$IRR <- NA
IR_post3$IRR[1] <- IRR
IR_post3$RRL[1] <- RRL
IR_post3$RRU[1] <- RRU
IR_post3$pvalue[1] <- round(pvalue, 3)

IR_post3



# Incidence rate: pre- vs post-Omicron, just IDP ----



# Get the data that I need, but filter to IDP
bt_pre <- bt %>%
    filter(IDorHV == "PID") %>%
    mutate(cases_pre = case_when(fb_bt == 1 & 
                                     !is.na(fb_init_bt_date) & 
                                     bt[bt$IDorHV == "PID",]$fb_init_bt_date <= ymd("2021-12-15") ~ 1,
                                 TRUE ~ 0)) %>% 
    mutate(start_date_pre = case_when(vaccine_1 == "Janssen" ~ ymd(vaccine_1_date),
                                      TRUE ~ ymd(vaccine_2_date))) %>%
    mutate(end_date_pre = 
               case_when(!is.na(fb_init_bt_date) & 
                             (ymd(fb_init_bt_date) <= ymd("2021-12-15")) ~ ymd(fb_init_bt_date),
                         !is.na(fb_init_bt_date) & 
                             (ymd(fb_init_bt_date) > ymd("2021-12-15")) ~ ymd("2021-12-15")))


# Deal with end date without invoking min(end, ymd("2021-12-15")) because that will get the end from the entire variable, not the minimum within the person. The rest is just formatting to make sure the dates stay formatted as date throughout
bt_pre$end_date_pre <- ifelse(is.na(bt_pre$fb_init_bt_date) & bt_pre$end <= ymd("2021-12-15"), 
                              as.Date(bt_pre$end, format = "%Y-%m-%d"), ymd(bt_pre$end_date_pre))
bt_pre$end_date_pre <- as.Date(bt_pre$end_date_pre, format = "%Y-%m-%d")
table2(bt_pre$end_date_pre)

bt_pre$end_date_pre <- ifelse(is.na(bt_pre$fb_init_bt_date) & bt_pre$end > ymd("2021-12-15"), 
                              ymd("2021-12-15"), bt_pre$end_date_pre)
bt_pre$end_date_pre <- as.Date(bt_pre$end_date_pre, format = "%Y-%m-%d")
table2(bt_pre$end_date_pre)


# Continue as before
bt_pre <- bt_pre %>%   
    drop_na(start_date_pre, end_date_pre) %>% # Remove any NAs
    filter(Subject.ID != "H0060") %>% # Drop H0060 with only 1 vaccine
    mutate(persondays_pre = as.numeric(ymd(end_date_pre) - ymd(start_date_pre))) %>%
    mutate(personyears_pre = as.numeric(persondays_pre / 365)) %>%
    filter(persondays_pre > 0) # Remove people who didn't conclude vax sequence till after Omicron era

bt_post <- bt %>%
    filter(IDorHV == "PID") %>%
    mutate(cases_pre = case_when(fb_bt == 1 & 
                                     !is.na(fb_init_bt_date) & 
                                     bt[bt$IDorHV == "PID",]$fb_init_bt_date <= ymd("2021-12-15") ~ 1,
                                 TRUE ~ 0)) %>% 
    mutate(cases_post = case_when(fb_bt == 1 & 
                                      !is.na(fb_init_bt_date) & 
                                      bt[bt$IDorHV == "PID",]$fb_init_bt_date > ymd("2021-12-15") ~ 1,
                                  TRUE ~ 0)) %>% 
    mutate(start_date_post = 
               case_when(vaccine_1 == "Janssen" & 
                             cases_pre == 0 &
                             ymd(vaccine_1_date) <= ymd("2021-12-15") ~ ymd("2021-12-15"),
                         vaccine_1 == "Janssen" & 
                             cases_pre == 0 &
                             ymd(vaccine_1_date) > ymd("2021-12-15") ~ ymd(vaccine_1_date),
                         vaccine_1 != "Janssen" & 
                             cases_pre == 0 &
                             ymd(vaccine_2_date) <= ymd("2021-12-15") ~ ymd("2021-12-15"),
                         vaccine_1 != "Janssen" & 
                             cases_pre == 0 &
                             ymd(vaccine_2_date) > ymd("2021-12-15") ~ ymd(vaccine_2_date),
                         TRUE ~ NA)) %>%
    mutate(end_date_post = 
               case_when(!is.na(fb_init_bt_date) & 
                             (ymd(fb_init_bt_date) <= ymd("2021-12-15")) ~ NA,
                         !is.na(fb_init_bt_date) & 
                             (ymd(fb_init_bt_date) > ymd("2021-12-15")) ~ ymd(fb_init_bt_date),
                         is.na(fb_init_bt_date) ~ ymd(end))) %>%
    
    drop_na(start_date_post, end_date_post) %>% # Drop ineligibles from post-Omicron consideration
    filter(Subject.ID != "H0060") %>% # Drop H0060 with only 1 vaccine
    mutate(persondays_post = as.numeric(ymd(end_date_post) - ymd(start_date_post))) %>%
    mutate(personyears_post = as.numeric(persondays_post / 365)) 


# Calculate the IRs
IR_pre <- bt_pre %>%
    dplyr::select(Subject.ID, IDorHV, cases_pre, persondays_pre, personyears_pre) %>%
    summarise(cases = sum(cases_pre),
              PT_d = sum(persondays_pre),
              PT_y = sum(personyears_pre),
              IR_y = round(cases / PT_y, 2),
              IR_y100 = round(100 * cases / PT_y, 2)) %>%
    mutate(era = "pre")
IR_pre

IR_post <- bt_post %>%
    dplyr::select(Subject.ID, IDorHV, cases_post, persondays_post, personyears_post) %>%
    summarise(cases = sum(cases_post),
              PT_d = sum(persondays_post),
              PT_y = sum(personyears_post),
              IR_y = round(cases / PT_y, 2),
              IR_y100 = round(100 * cases / PT_y, 2)) %>%
    mutate(era = "post")
IR_post


# Combine the incidence information across the two eras
IR_era <- rbind(IR_pre, IR_post)
IR_era


# Add in the CIs, numbers from p. 434 of Epidemiology Beyond The Basics, 3rd ed.
IR_era$`95% UL` <- IR_era$`95% LL` <- NA

IR_era[IR_era$era == "post",]$`95% LL` <- 
    round(100 * IR_era[IR_era$era == "post",]$cases *  (0.798 + (3/10) * (0.809 - 0.798)) /
              IR_era[IR_era$era == "post",]$PT_y, 2)

IR_era[IR_era$era == "post",]$`95% UL` <- 
    round(100 * IR_era[IR_era$era == "post",]$cases * (1.24 - (7/10) * (1.25 - 1.24)) /
              IR_era[IR_era$era == "post",]$PT_y, 2)

IR_era[IR_era$era == "pre",]$`95% LL` <- round(100 * IR_era[IR_era$era == "pre",]$cases * 
                                                   0.324 / IR_era[IR_era$era == "pre",]$PT_y, 2)
IR_era[IR_era$era == "pre",]$`95% UL` <- round(100 * IR_era[IR_era$era == "pre",]$cases * 
                                                   2.33 / IR_era[IR_era$era == "pre",]$PT_y, 2)


# Calculate IRR (PID -> post, HV -> pre)
IRR <- round((IR_era[IR_era$era == "post",]$cases / IR_era[IR_era$era == "post",]$PT_y) / 
                 (IR_era[IR_era$era == "pre",]$cases / IR_era[IR_era$era == "pre",]$PT_y), 2)


# Calculate IRR 95% CIs
Phat <- IR_era[IR_era$era == "post",]$cases / sum(IR_era$cases)

PL <- Phat - (qnorm(0.975) * sqrt(Phat * (1 - Phat) / sum(IR_era$cases)))
PU <- Phat + (qnorm(0.975) * sqrt(Phat * (1 - Phat) / sum(IR_era$cases)))

RRL <- round((PL / (1 - PL)) * (IR_era[IR_era$era == "pre",]$PT_y / 
                                    IR_era[IR_era$era == "post",]$PT_y), 2)
RRU <- round((PU / (1 - PU)) * (IR_era[IR_era$era == "pre",]$PT_y / 
                                    IR_era[IR_era$era == "post",]$PT_y), 2)

IRR; RRL; RRU


# Get pvalue
E1 <- sum(IR_era$cases) * IR_era[IR_era$era == "post",]$PT_y / sum(IR_era$PT_y)
E2 <- sum(IR_era$cases) * IR_era[IR_era$era == "pre",]$PT_y / sum(IR_era$PT_y)

ts <- ((IR_era[IR_era$era == "post",]$cases - E1)^2 / E1) + 
    ((IR_era[IR_era$era == "pre",]$cases - E2)^2 / E2)

pvalue <- pchisq(q = ts, df = 1, lower.tail = F)
pvalue


# Put everything together
IR2 <- IR_era
IR2$pvalue <- IR2$RRU <- IR2$RRL <- IR2$IRR <- NA
IR2$IRR[1] <- IRR
IR2$RRL[1] <- RRL
IR2$RRU[1] <- RRU
IR2$pvalue[1] <- round(pvalue, 3)

IR2

rm(Phat, PL, PU, ts, E1, E2, IRR, RRL, RRU, pvalue, IR_era, IR2, bt_pre)
rm(IR_post, IR_post2, IR_post3, IR_pre, bt_post)



# Kaplan-Meier curves: HV vs IDP ----



# Generate survival curves for HV and IDP
survival_IDPHV <- survfit(Surv(time = persondays, event = fb_bt) ~ IDorHV, data = bt)


# KM plot for HV and IDP
# Get pvalue by switching pval argument on and off
KM_1 <- ggsurvplot(survival_IDPHV,
                   data = bt,
                   conf.int = F,
                   #pval = T,
                   #pval.method = T,
                   size = 1,
                   censor.size = 7,
                   legend.title = "",
                   legend.labs = c("Healthy volunteer", "Immune-deficient person"),
                   legend = c(0.15, 0.20),
                   censor.shape = "l",
                   fontsize = 6,
                   
                   # Add risk table
                   risk.table = "nrisk_cumcensor",
                   tables.height = 0.2,
                   tables.theme = theme_cleantable(),
                   tables.y.text = FALSE,
                   palette = c("black", "#8856a7"),
                   xlab = "Days from primary series completion to first breakthrough",
                   ylab = "Survival probability (%)") 


# Fix the y-axis and the legend
y <- which(KM_1$plot$scales$find("y"))
KM_1$plot$scales$scales[[y]] <- scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                                                   labels = c(0, 25, 50, 75, 100))

KM_1$plot <- KM_1$plot + 
    theme(legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18)) + 
    annotate("text", x = 800, y = 0.9, label = "p = 0.49", size = 14/.pt) 
KM_1$plot <- KM_1$plot + 
    guides(color = guide_legend(override.aes = list(shape = NA)))

KM_1



# Kaplan-Meier curves: HV vs IDP subgroups ----



# Generate survival curves by immune group
bt$new_group <- factor(bt$new_group, 
                       levels = c("Healthy volunteer", "Antibody deficiency", 
                                  "PIRD", "Combined immunodeficiency",
                                  "Other IEI", "Other immune disorder"))

survival_newgroup <- survfit(Surv(time = persondays, event = fb_bt) ~ new_group, data = bt)


# KM plot by immune subgroup
# Get pvalue by switching pval argument on and off
KM_2 <- ggsurvplot(survival_newgroup,
                   data = bt,
                   conf.int = F,
                   #pval = T,
                   #pval.method = T,
                   size = 1,
                   censor.size = 7,
                   legend.title = "",
                   legend.labs = c("Healthy volunteer", "Antibody deficiency", 
                                   "PIRD", "Combined immunodeficiency",
                                   "Other IEI", "Other immune disorder"),
                   legend = c(0.15, 0.20),
                   censor.shape = "l",
                   fontsize = 6,
                   
                   # Add risk table
                   risk.table = "nrisk_cumcensor",
                   tables.height = 0.2,
                   tables.theme = theme_cleantable(),
                   tables.y.text = FALSE,
                   palette = c("black", "#785EF0", "#048B90", "#CF9A03", "#D10106", "#648fff"),
                   xlab = "Days from primary series completion to first breakthrough",
                   ylab = "Survival probability (%)")


# Fix the y-axis and the legend
y2 <- which(KM_2$plot$scales$find("y"))
KM_2$plot$scales$scales[[y2]] <- scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                                                    labels = c(0, 25, 50, 75, 100))

KM_2$plot <- KM_2$plot + 
    theme(legend.title = element_text(size = 18),
          legend.text = element_text(size = 18),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18)) + 
    annotate("text", x = 800, y = 0.9, label = "p = 0.19", size = 14/.pt) 
KM_2$plot <- KM_2$plot + 
    guides(color = guide_legend(override.aes = list(shape = NA)))

KM_2


# Remove unless data objects
rm(KM_1, KM_2, bt, y, y2, survival_IDPHV, survival_newgroup)



# **** PARAGRAPH 3: DISEASE OCCURRENCE **** ----



# Comparing COVID symptoms ----


# NOTE: Several symptoms and medical conditions we examined (i.e., rash, blood clots, stroke, low blood oxygen, disorientation, kidney injury, shock, myocarditis, overnight hospital stays, supplemental oxygen, ECMO, and death) were quite rare in the cohort. To maintain participants' privacy, we have omitted these variables from the publicly available dataset. To be transparent about our analsyis, we have left the symptom analysis as-is. We have commented the code out to prevent the occurrence of errors upon running the full code with a variety of missing symptom data. 


# # Formatting of the infection dataset
# inf$symptom_onset_date <- as.Date(inf$symptom_onset_date)
# inf$recovery_date <- as.Date(inf$recovery_date)
# 
# inf <- inf %>%
#   group_by(SUBJECT.ID.) %>%
#   mutate(symptom_days = recovery_date - symptom_onset_date)
# 
# 
# # Edit some of the symptom data
# inf <- inf %>%
#   mutate(asymp_infection = case_when(felt_feverish != "Yes" & 
#                                      chills != "Yes" & 
#                                      cough != "Yes" & 
#                                      sorethroat != "Yes" & 
#                                      runnynose  != "Yes" &                                    
#                                      difficult_breathing!= "Yes" & 
#                                      musclepain!= "Yes" & 
#                                      chestpain!= "Yes" & 
#                                      abdominalpain != "Yes" & 
#                                      nausea_vomiting!= "Yes" & 
#                                      diarrhea!= "Yes" & 
#                                      headache!= "Yes" & 
#                                      fatigue!= "Yes" & 
#                                      lossofsmelltaste!= "Yes" & 
#                                      rash!= "Yes" & 
#                                      bloodclots!= "Yes" & 
#                                      stroke!= "Yes" & 
#                                      lowbloodoxlevel != "Yes" & 
#                                      disorientation != "Yes" & 
#                                      kidneyinjury != "Yes" & 
#                                      shock!= "Yes" &
#                                      myocarditis != "Yes"  ~ "Yes",
#                                      TRUE ~ "No"))
# 
# columns_with_idk <- c("felt_feverish", "chills", "cough", 
#                       "sorethroat", "runnynose", "difficult_breathing", "musclepain", 
#                       "chestpain", "abdominalpain", "nausea_vomiting", "diarrhea", 
#                       "headache", "fatigue", "lossofsmelltaste", "rash", "bloodclots", 
#                       "stroke", "lowbloodoxlevel", "disorientation", "kidneyinjury", 
#                       "shock", "myocarditis", "other_viral", "other_bacterial", 
#                       "doctor_ER", "influenza_positive", 
#                       "hosp_overnight", "supplemental_oxygen", "ECMO", "death")
# 
# 
# # Replace idks with NA
# inf <- inf %>%
#   mutate(across(all_of(columns_with_idk), ~str_replace_all(., "Don't know", NA_character_)))
# 
# 
# # Filter dataset just to breakthrough infections
# inf_bk <- inf %>%
#   filter(infection_timing == "Breakthrough")
# 
# 
# # How many of these resulted in hospitalization?
# table2(inf_bk$hosp_admission_date)
# table2(inf_bk$hosp_days)
# table2(inf_bk$hosp_discharge_date)
# 
# check <- inf_bk[!is.na(inf_bk$hosp_admission_date) | 
#                     !is.na(inf_bk$hosp_days) | 
#                     (!is.na(inf_bk$hosp_overnight) & inf_bk$hosp_overnight == "Yes"),]
# table2(check$IDorHV) # 5 
# 
# 
# # How many deaths?
# table2(inf_bk$death) # 0 among the breakthroughs
# table2(inf$death) # 0 overall
# 
# 
# # Make a list of the IDs and breakthrough dates for the initial breakthroughs
# list <- new_elisa %>%
#   filter(fb_bt == 1) %>%
#   ungroup() %>%
#   dplyr::select(Subject.ID, fb_init_bt_date) %>%
#   distinct() %>%
#   mutate(newvar = paste0(Subject.ID, "__", fb_init_bt_date))
# 
# 
# # How many people of the initial breakthroughs completed their symptoms survey?
# check <- inf_bk %>%
#     mutate(newvar = paste0(SUBJECT.ID., "__", covid_date)) %>%
#     ungroup() %>%
#     filter(newvar %in% list$newvar) %>%
#     filter(!is.na(symptom_onset_date))
# 
# nrow(check) / nrow(list)
# 
# 
# # Make a visual, will require reformatting the data
# bk_viz <- inf_bk %>%
#   mutate(newvar = paste0(SUBJECT.ID., "__", covid_date)) %>%
#   ungroup() %>%
#   filter(newvar %in% list$newvar) %>%
#   filter(!is.na(symptom_onset_date)) %>% # Filter to people (62/116) who completed symptoms survey
#   dplyr::select("SUBJECT.ID.", "IDorHV", "infection_number", "felt_feverish", "chills", "cough", 
#                 "sorethroat", "runnynose", "difficult_breathing", "musclepain", 
#                 "chestpain", "abdominalpain", "nausea_vomiting", "diarrhea", 
#                 "headache", "fatigue", "lossofsmelltaste", "rash", "bloodclots", 
#                 "stroke", "lowbloodoxlevel", "disorientation", "kidneyinjury", 
#                 "shock", "myocarditis", "other_viral", "other_bacterial", 
#                 "doctor_ER", 
#                 "hosp_overnight", "supplemental_oxygen", "ECMO", 
#                 "death", "asymp_infection") %>%
#   distinct(SUBJECT.ID., infection_number, .keep_all = T)
# 
# 
# # Restrict to the useful variables
# bk_viz <- bk_viz %>%
#   pivot_longer(cols = 4:33, names_to = "symptom_type", values_to = "value")
# 
# 
# # Ensure that there are rows for 0 Yes answers among HVs
# check <- bk_viz %>%
#   expand(IDorHV, symptom_type, value) 
# 
# check2 <- check %>% 
#   mutate(count = 0) %>%
#   mutate(newvar = paste0(symptom_type, "__", value, "__", IDorHV))
# 
# bk_viz_1fb_pre <- bk_viz %>%
#   group_by(IDorHV, symptom_type, value) %>%
#   summarise(count = n()) %>%
#   mutate(newvar = paste0(symptom_type, "__", value, "__", IDorHV))
# 
# check2 <- check2 %>%
#   filter(newvar %!in% bk_viz_1fb_pre$newvar)
# 
# bk_viz_1fb_post <- rbind(bk_viz_1fb_pre, check2)
# 
# 
# # Remove some symptoms and add count variable
# bk_viz_1fb_post2 <- bk_viz_1fb_post %>%
#   filter(symptom_type %!in% c("asymp_infection", "doctor_ER", "ECMO", "hosp_overnight",
#                               "other_bacterial", "other_viral", "supplemental_oxygen")) %>%
#   group_by(symptom_type, IDorHV) %>%
#   mutate(sum = sum(count)) %>%
#   arrange(symptom_type, IDorHV) %>%
#   filter(value != "No") %>%
#   group_by(symptom_type)
# 
# 
# # Add pvalues and proportion variables
# bk_viz_1fb_post2$pval <- NA
# bk_viz_1fb_post2$prop <- NA
# for(i in 1:nrow(bk_viz_1fb_post2)){
#   if(i %% 2 == 1){bk_viz_1fb_post2$pval[i] <- 
#     round(prop.test(x = c(bk_viz_1fb_post2$count[i], bk_viz_1fb_post2$count[i + 1]),
#                     n = c(bk_viz_1fb_post2$sum[i], bk_viz_1fb_post2$sum[i + 1]))$p.value, 3)}
#   if(i %% 2 != 1){bk_viz_1fb_post2$pval[i] <- bk_viz_1fb_post2$pval[i - 1]}
#   bk_viz_1fb_post2$prop[i] <- round(bk_viz_1fb_post2$count[i] / bk_viz_1fb_post2$sum[i], 3)
#   
# } # Warnings ok
# 
# bk_viz_1fb_post2 %>%
#   filter(pval < 0.05)
# 
# 
# # Make the heatmap
# heatmap <- bk_viz_1fb_post2 %>%
#   filter(symptom_type %!in% c("death", "stroke", "myocarditis", "shock"))
# 
# heatmap$symptom_type <- case_when(heatmap$symptom_type == "abdominalpain" ~ "Abdominal pain",
#                                   heatmap$symptom_type == "bloodclots" ~ "Blood clot",
#                                   heatmap$symptom_type == "chestpain" ~ "Chest pain",
#                                   heatmap$symptom_type == "chills" ~ "Chills",
#                                   heatmap$symptom_type == "cough" ~ "Cough",
#                                   heatmap$symptom_type == "diarrhea" ~ "Diarrhea",
#                                   heatmap$symptom_type == "difficult_breathing" ~ "Difficulty breathing",
#                                   heatmap$symptom_type == "disorientation" ~ "Disorientation",
#                                   heatmap$symptom_type == "fatigue" ~ "Fatigue",
#                                   heatmap$symptom_type == "felt_feverish" ~ "Fever",
#                                   heatmap$symptom_type == "headache" ~ "Headache",
#                                   heatmap$symptom_type == "kidneyinjury" ~ "Kidney injury",
#                                   heatmap$symptom_type == "lossofsmelltaste" ~ "Loss of smell/taste",
#                                   heatmap$symptom_type == "lowbloodoxlevel" ~ "Low blood oxygen",
#                                   heatmap$symptom_type == "musclepain" ~ "Muscle pain",
#                                   heatmap$symptom_type == "nausea_vomiting" ~ "Nausea/vomiting",
#                                   heatmap$symptom_type == "rash" ~ "Rash",
#                                   heatmap$symptom_type == "runnynose" ~ "Runny nose",
#                                   heatmap$symptom_type == "sorethroat" ~ "Sore throat",
#                                   TRUE ~ heatmap$symptom_type)
# 
# ggplot(heatmap, aes(y = IDorHV, x = symptom_type, fill = prop)) + 
#   theme_fausto() + 
#   geom_tile() + 
#   scale_y_discrete(breaks = c("HV", "PID"),
#                    labels = c("HV" = "HV",
#                               "PID" = "IDP"),
#                    expand = c(0, 0)) +
#   scale_x_discrete(expand = c(0, 0)) + 
#   scale_fill_viridis(option = "A",
#                      begin = 0, end = 1,
#                      name = "Prevalence (%) ",
#                      limits = c(0, 1),
#                      breaks = c(0, 0.25, 0.5, 0.75, 1),
#                      labels = c("0", "25", "50", "75", "100")) + 
#   guides(fill = guide_colourbar(barheight = 1,
#                                 barwidth = 40)) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", 
#                                    size = 25),
#         axis.text.y = element_text(angle = 90, hjust = 0.5, face = "bold", 
#                                    color = "black", size = 25),
#         legend.text = element_text(angle = 90, hjust = 1, size = 25),
#         legend.title = element_text(size = 25),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank()) # Export as PDF w/ length 19in length, 6.5in width landscape
# 
# 
# # Look at most and least common symptoms among IDP
# View(heatmap)
# 
# 
# # Remove useless data objects
# rm(heatmap, inf_bk, list, columns_with_idk, i, check, check2,
#    bk_viz, bk_viz_1fb_pre, bk_viz_1fb_post, bk_viz_1fb_post2)



# **** PARAGRAPH 4: VACCINATION HISTORIES AND IMPACT **** ----



# N vaccines prior to initial breakthrough ----



# Put all breakthrough infections in one dataset and check for accuracy
fb_dataset <- new_elisa %>%
  group_by(Subject.ID) %>%
  slice(1) %>%
  filter(fb_bt == 1) %>%
  dplyr::select(Subject.ID, IDorHV, fb_predose1_inf, CovidPos_1:CovidPos_3, 
                fb_init_bt_date, vaccine_1_date:vaccine_6_date,
                vaccine_1:vaccine_6) %>%
  ungroup()


# Since monovalent and bivalent vaccines are different and most infections were Omicron, focus on the number of monovalent and bivalent vaccines separately
# # Number of vaccines before initial breakthrough infection by immune status
# vaxx <- fb_dataset %>% 
#   rowwise() %>%
#   mutate(N_vax = sum(across(vaccine_1_date:vaccine_6_date, ~ !is.na(.) & . < fb_init_bt_date))) %>%
#   ungroup()
# 
# vaxx %>%
#   summarise(mean = round(mean(N_vax), 1),
#             median = round(median(N_vax), 1), 
#             .by = IDorHV) 
# 
# 
# # Test the difference of the mean
# wilcox.test(vaxx[vaxx$IDorHV == "PID",]$N_vax, vaxx[vaxx$IDorHV == "HV",]$N_vax)
# 
# 
# # Remove useless data objects
# rm(vaxx)


# Set post-breakthrough vaccines to NA
fb_dataset2 <- fb_dataset

for(i in 1:nrow(fb_dataset2)){
    dat <- fb_dataset2[i, c("vaccine_1_date", "vaccine_2_date", "vaccine_3_date", 
                            "vaccine_4_date", "vaccine_5_date", "vaccine_6_date",
                            "fb_init_bt_date")]
    index <- which(dat[1:6] >= dat$fb_init_bt_date)
    fb_dataset2[i, 13 + index] <- NA
}
    

# Count the number of mono and bivalent vaccines before breakthrough infection by immune group
vaxx2 <- fb_dataset2 %>% 
     rowwise() %>%
     mutate(N_mono_moderna = sum(across(vaccine_1:vaccine_6, 
                                        ~ !is.na(.) &  . == "Moderna (monovalent)")), 
            N_mono_pfizer = sum(across(vaccine_1:vaccine_6, 
                                       ~ !is.na(.) &  . == "Pfizer (monovalent)")),
            N_mono = sum(N_mono_moderna, N_mono_pfizer)) %>%
    
     mutate(N_biva_moderna = sum(across(vaccine_1:vaccine_6, 
                                       ~ !is.na(.) &  . == "Moderna (bivalent)")), 
           N_biva_pfizer = sum(across(vaccine_1:vaccine_6, 
                                      ~ !is.na(.) &  . == "Pfizer (bivalent)")),
           N_biva = sum(N_biva_moderna, N_biva_pfizer)) %>%
     ungroup()


# Summarize the data
vaxx2 %>%
    summarise(mean_mono = round(mean(N_mono), 1),
              mean_biva = round(mean(N_biva), 1),
              .by = IDorHV) 


# Check whether the numbers are different across immune groups
wilcox.test(vaxx2[vaxx2$IDorHV == "PID",]$N_mono, vaxx2[vaxx2$IDorHV == "HV",]$N_mono)



# N 1st gen mRNA vaccines before Omicron ----



# Generate the dataset for this analysis
vaxx_omi <- new_elisa %>%
    group_by(Subject.ID) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(vaccine_1_date = 
               case_when(vaccine_1 %in% c("Moderna (monovalent)", "Pfizer (monovalent)") ~ 
                             vaccine_1_date, TRUE ~ NA)) %>%
    mutate(vaccine_2_date = 
               case_when(vaccine_2 %in% c("Moderna (monovalent)", "Pfizer (monovalent)") ~ 
                             vaccine_2_date, TRUE ~ NA)) %>%
    mutate(vaccine_3_date = 
               case_when(vaccine_3 %in% c("Moderna (monovalent)", "Pfizer (monovalent)") ~ 
                             vaccine_3_date, TRUE ~ NA)) %>%
    mutate(vaccine_4_date = 
               case_when(vaccine_4 %in% c("Moderna (monovalent)", "Pfizer (monovalent)") ~ 
                             vaccine_4_date, TRUE ~ NA)) %>%
    mutate(vaccine_5_date = 
               case_when(vaccine_5 %in% c("Moderna (monovalent)", "Pfizer (monovalent)") ~ 
                             vaccine_5_date, TRUE ~ NA)) %>%
    mutate(vaccine_6_date = 
               case_when(vaccine_6 %in% c("Moderna (monovalent)", "Pfizer (monovalent)") ~ 
                             vaccine_6_date, TRUE ~ NA)) %>%
    rowwise() %>%
    mutate(N_vax = sum(across(vaccine_1_date:vaccine_6_date, ~ !is.na(.) & . < ymd("2021-12-15")))) %>%
    ungroup()


# Summarize the data
check <- vaxx_omi %>%
    summarise(atleast3 = sum(N_vax >= 3),
              n = n(),
              percent = round(100*sum(N_vax >= 3) / n, 1),
              .by = IDorHV)


# Test the proportions for equality
prop.test(x = check$atleast3, n = check$n)



# Percent uptake of the bivalent vaccine ----



# Overall, how many received a bivalent vaccine?
new_elisa %>%
  dplyr::select(Subject.ID, IDorHV, vaccine_1:vaccine_6) %>%
  group_by(Subject.ID) %>%
  slice(1) %>%
  mutate(bivalent = case_when(vaccine_1 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") ~ 1,
                              vaccine_2 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") ~ 1,
                              vaccine_3 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") ~ 1,
                              vaccine_4 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") ~ 1,
                              vaccine_5 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") ~ 1,
                              vaccine_6 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") ~ 1,
                              TRUE ~ 0)) %>%
  ungroup() %>%
  summarise(mean = round(100*mean(bivalent), 2),
            count = sum(bivalent),
            .by = IDorHV)

# Among those who had a breakthrough infection, how many received a bivalent vaccine?
new_elisa %>%
  filter(fb_bt == 1) %>%
  dplyr::select(Subject.ID, IDorHV, vaccine_1:vaccine_6) %>%
  group_by(Subject.ID) %>%
  slice(1) %>%
  mutate(bivalent = case_when(vaccine_1 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") ~ 1,
                              vaccine_2 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") ~ 1,
                              vaccine_3 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") ~ 1,
                              vaccine_4 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") ~ 1,
                              vaccine_5 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") ~ 1,
                              vaccine_6 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") ~ 1,
                              TRUE ~ 0)) %>%
  ungroup() %>%
  summarise(mean = round(100*mean(bivalent), 2),
            count = sum(bivalent),
            .by = IDorHV)



# N received bivalent show within 2 months of release ----



# Bivalent vaccines were approved/recommended on Sept 1 2022
# https://www.hhs.gov/sites/default/files/secretarial-directive-covid-19-bivalent-vaccine-boosters.pdf#:~:text=On%20September%201,%202022,%20the%20Advisory%20Committee%20on%20Immunization%20Practices

# Make a dataset
ELISA_ind <- new_elisa %>%
    group_by(Subject.ID) %>%
    slice_head() 


# Make a plot
ELISA_ind %>%
    filter(!is.na(date_first_bivalent)) %>%
    ggplot(aes(x = date_first_bivalent)) + geom_bar()

table(ELISA_ind$date_first_bivalent, useNA = "always")
152 / 282


# Define who received a bivalent vaccine at any point
# NOTE: Because noise was introduced into date variables, results here disagree slightly with what is reported in the paper, the latter of which is based on the full, true data.
bivalent <- ELISA_ind %>%
    filter(!is.na(date_first_bivalent)) %>%
    select(Subject.ID, date_first_bivalent, IDorHV) %>%
    mutate(bivalent_2mo = case_when(date_first_bivalent < ymd("2022-11-01") ~ 1,
                                    TRUE ~ 0))


# Get the statistics
table(bivalent$IDorHV)
table(bivalent$bivalent_2mo)
table(bivalent$bivalent_2mo, bivalent$IDorHV)

70/130
12/29
58/101

prop.test(x = c(12, 58), n = c(29, 101))

# 53.8% of all participants
# 57% of IDP received it in 2 months
# 41% of HV received it in 2 months


# Remove useless data objects
rm(bivalent, ELISA_ind)



# N first breakthrough after bivalent vaccine ----



# Define the necessary data needed to answer this question
bt_bi <- fb_dataset %>%
  mutate(inf_after_bivalent = 
           case_when(vaccine_1 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") & 
                       ymd(fb_init_bt_date) >= ymd(vaccine_1_date) ~ 1,
                     vaccine_2 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") & 
                       ymd(fb_init_bt_date) >= ymd(vaccine_2_date) ~ 1,
                     vaccine_3 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") & 
                       ymd(fb_init_bt_date) >= ymd(vaccine_3_date) ~ 1,
                     vaccine_4 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") & 
                       ymd(fb_init_bt_date) >= ymd(vaccine_4_date) ~ 1,
                     vaccine_5 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") & 
                       ymd(fb_init_bt_date) >= ymd(vaccine_5_date) ~ 1,
                     vaccine_6 %in% c("Pfizer (bivalent)", "Moderna (bivalent)") & 
                       ymd(fb_init_bt_date) >= ymd(vaccine_6_date) ~ 1,
                     TRUE ~ 0)) %>%
  ungroup()


# Identify when the first bivalent vaccine was received among the breakthroughs
bt_bi$first_bivalent_date <- NA
vars <- names(bt_bi)[c(7, 15:19)]
for(i in 1:nrow(bt_bi)){
  for(j in 2:6){
    if((bt_bi[i, vars[j]] %in% c("Pfizer (bivalent)", "Moderna (bivalent)")) & 
       (bt_bi[i, vars[j - 1]] %!in% c("Pfizer (bivalent)", "Moderna (bivalent)"))){
      bt_bi$first_bivalent_date[i] <- bt_bi[i, 8:13][j]
    }
  }
}
bt_bi$first_bivalent_date <- ymd(bt_bi$first_bivalent_date) # Warning OK


# Define those with a pre-bivalent breakthrough as a check on the prior code
bt_bi$pre_bivalent_bt <- 0
for(i in 1:nrow(bt_bi)){
  if(!is.na(bt_bi$first_bivalent_date)[i] & 
     !is.na(bt_bi$fb_init_bt_date)[i] & 
     bt_bi$first_bivalent_date[i] >= ymd(bt_bi$fb_init_bt_date)[i]){
    bt_bi$pre_bivalent_bt[i] <- 1
  }
}


# So how many got initially infected after a bivalent vaccine?
nrow(bt_bi[bt_bi$inf_after_bivalent == 1,])
sum(bt_bi[bt_bi$inf_after_bivalent == 1,]$pre_bivalent_bt) # Check their pre-valent status


# Summarize the data
bt_bi %>%
  summarise(n = sum(inf_after_bivalent),
            inf_after_bivalent = round(100*mean(inf_after_bivalent), 2),
            .by = IDorHV)


# Test for a difference
prop.test(x = c(sum(bt_bi[bt_bi$IDorHV == "PID",]$inf_after_bivalent), 
                sum(bt_bi[bt_bi$IDorHV == "HV",]$inf_after_bivalent)),
          n = c(sum(bt_bi$IDorHV == "PID"), 
                sum(bt_bi$IDorHV == "HV")))


# Check what happened with the pre-dose 1 infections
check <- new_elisa %>%
  filter(fb_predose1_inf == 1) %>%
  dplyr::select(names(new_elisa)[names(new_elisa) %in% names(fb_dataset)]) %>%
  group_by(Subject.ID) %>%
  slice(1) %>%
  ungroup


# Remove useless data objects
rm(check, dat, fb_dataset2, vaxx_omi, vaxx2, i, j, bt_bi, fb_dataset)



# Incidence before and after the bivalent booster ----



bt <- new_elisa %>%
    dplyr::select(Subject.ID, IDorHV, new_group, fb_bt, end,
                  vaccine_1_date:vaccine_6_date, fb_init_bt_date, vaccine_1:vaccine_6) %>%
    group_by(Subject.ID) %>%
    slice(1) %>%
    ungroup

bt$first_bivalent_date <- "2024-01-01"
for(i in 1:nrow(bt)){
    for(j in 2:6){
        if((bt[i, vars[j]] %in% c("Pfizer (bivalent)", "Moderna (bivalent)")) & 
           (bt[i, vars[j - 1]] %!in% c("Pfizer (bivalent)", "Moderna (bivalent)"))){
            bt$first_bivalent_date[i] <- bt[i, 6:11][j]
        }
    }
}

bt <- as.data.frame(bt)

bt_pre <- bt %>%
    mutate(cases_pre = case_when(fb_bt == 1 & 
                                     !is.na(fb_init_bt_date) & 
                                     bt$fb_init_bt_date <= ymd(first_bivalent_date) ~ 1,
                                 TRUE ~ 0)) %>% 
    mutate(start_date_pre = case_when(vaccine_1 == "Janssen" ~ ymd(vaccine_1_date),
                                      TRUE ~ ymd(vaccine_2_date))) %>%
    mutate(end_date_pre = 
               case_when(!is.na(fb_init_bt_date) & 
                             (ymd(fb_init_bt_date) <= ymd(first_bivalent_date)) ~ ymd(fb_init_bt_date),
                         !is.na(fb_init_bt_date) & 
                             (ymd(fb_init_bt_date) > ymd(first_bivalent_date)) ~ ymd(first_bivalent_date)))


# Deal with end date without invoking min(end, ymd(first_bivalent_date)) because that will get the end from the entire variable, not the minimum within the person. The rest is just formatting to make sure the dates stay formatted as date throughout
bt_pre$end_date_pre <- ifelse(is.na(bt_pre$fb_init_bt_date) & bt_pre$end <= bt_pre$first_bivalent_date,
                              as.Date(bt_pre$end, format = "%Y-%m-%d"), 
                              ymd(bt_pre$end_date_pre))
bt_pre$end_date_pre <- as.Date(bt_pre$end_date_pre, format = "%Y-%m-%d")
table2(bt_pre$end_date_pre)

bt_pre$first_bivalent_date <- as.character(bt_pre$first_bivalent_date)
bt_pre$end_date_pre <- ifelse(is.na(bt_pre$fb_init_bt_date) & bt_pre$end > bt_pre$first_bivalent_date, 
                              as.Date(bt_pre$first_bivalent_date, format = "%Y-%m-%d"), 
                              ymd(bt_pre$end_date_pre))
bt_pre$end_date_pre <- as.Date(bt_pre$end_date_pre, format = "%Y-%m-%d")
table2(bt_pre$end_date_pre)

    
# Continue as before
bt_pre <- bt_pre %>%
    drop_na(start_date_pre, end_date_pre) %>% # Remove any NAs
    filter(Subject.ID != "H0060") %>% # Drop H0060 with only 1 vaccine
    mutate(persondays_pre = as.numeric(ymd(end_date_pre) - ymd(start_date_pre))) %>%
    mutate(personyears_pre = as.numeric(persondays_pre / 365))

table2(bt_pre$cases_pre)

# Now do the calculation for the post-bivalent period
bt_post <- bt %>%
    mutate(cases_pre = case_when(fb_bt == 1 & 
                                     !is.na(fb_init_bt_date) & 
                                     bt$fb_init_bt_date <= ymd(first_bivalent_date) ~ 1,
                                 TRUE ~ 0)) %>% 
    mutate(cases_post = case_when(fb_bt == 1 & 
                                      !is.na(fb_init_bt_date) & 
                                      bt$fb_init_bt_date > ymd(first_bivalent_date) ~ 1,
                                  TRUE ~ 0)) %>% 
    mutate(start_date_post = 
               case_when(vaccine_1 == "Janssen" & 
                             cases_pre == 0 &
                             ymd(vaccine_1_date) <= ymd(first_bivalent_date) ~ 
                             ymd(first_bivalent_date),
                         vaccine_1 == "Janssen" & 
                             cases_pre == 0 &
                             ymd(vaccine_1_date) > ymd(first_bivalent_date) ~ NA,
                         vaccine_1 != "Janssen" & 
                             cases_pre == 0 &
                             ymd(vaccine_2_date) <= ymd(first_bivalent_date) ~ 
                             ymd(first_bivalent_date),
                         vaccine_1 != "Janssen" & 
                             cases_pre == 0 &
                             ymd(vaccine_2_date) > ymd(first_bivalent_date) ~ NA,
                         TRUE ~ NA)) %>%
    mutate(end_date_post = 
               case_when(!is.na(fb_init_bt_date) & 
                             (ymd(fb_init_bt_date) <= ymd(first_bivalent_date)) ~ NA,
                         !is.na(fb_init_bt_date) & 
                             (ymd(fb_init_bt_date) > ymd(first_bivalent_date)) ~
                             ymd(fb_init_bt_date),
                         is.na(fb_init_bt_date) ~ ymd(end))) %>%
    
    drop_na(start_date_post, end_date_post) %>% # Drop ineligibles from post-Omicron consideration
    filter(Subject.ID != "H0060") %>% # Drop H0060 with only 1 vaccine
    mutate(persondays_post = as.numeric(ymd(end_date_post) - ymd(start_date_post))) %>%
    mutate(personyears_post = as.numeric(persondays_post / 365)) %>% 
    # n=177, number of people without a case before their bivalent vaccine
    filter(persondays_post > 0) # Remove people who didn't get bivalent booster
    # n=87 after removing those without a bivalent vaccine


# Calculate the IRs
IR_pre <- bt_pre %>%
    dplyr::select(Subject.ID, IDorHV, cases_pre, persondays_pre, personyears_pre) %>%
    summarise(cases = sum(cases_pre),
              PT_d = sum(persondays_pre),
              PT_y = sum(personyears_pre),
              IR_y = round(cases / PT_y, 2),
              IR_y100 = round(100 * cases / PT_y, 2)) %>%
    mutate(era = "pre")
IR_pre

IR_post <- bt_post %>%
    dplyr::select(Subject.ID, IDorHV, cases_post, persondays_post, personyears_post) %>%
    summarise(cases = sum(cases_post),
              PT_d = sum(persondays_post),
              PT_y = sum(personyears_post),
              IR_y = round(cases / PT_y, 2),
              IR_y100 = round(100 * cases / PT_y, 2)) %>%
    mutate(era = "post")
IR_post



# Combine the incidence information across the two eras
IR_era <- rbind(IR_pre, IR_post)
IR_era


# Add in the CIs, numbers from p. 434 of Epidemiology Beyond The Basics, 3rd ed.
IR_era$`95% UL` <- IR_era$`95% LL` <- NA

IR_era[IR_era$era == "pre",]$`95% LL` <- 
    round(100 * IR_era[IR_era$era == "pre",]$cases *  (0.818 + (14/20) * (0.833 - 0.818)) /
              IR_era[IR_era$era == "pre",]$PT_y, 2)

IR_era[IR_era$era == "pre",]$`95% UL` <- 
    round(100 * IR_era[IR_era$era == "pre",]$cases * (1.22 - (6/20) * (1.22 - 1.2)) /
              IR_era[IR_era$era == "pre",]$PT_y, 2)

IR_era[IR_era$era == "post",]$`95% LL` <- round(100 * IR_era[IR_era$era == "post",]$cases * 
                                                   0.517 / IR_era[IR_era$era == "post",]$PT_y, 2)
IR_era[IR_era$era == "post",]$`95% UL` <- round(100 * IR_era[IR_era$era == "post",]$cases * 
                                                   1.75 / IR_era[IR_era$era == "post",]$PT_y, 2)


# Calculate IRR (PID -> post, HV -> pre)
IRR <- round((IR_era[IR_era$era == "post",]$cases / IR_era[IR_era$era == "post",]$PT_y) / 
                 (IR_era[IR_era$era == "pre",]$cases / IR_era[IR_era$era == "pre",]$PT_y), 2)


# Calculate IRR 95% CIs
Phat <- IR_era[IR_era$era == "post",]$cases / sum(IR_era$cases)

PL <- Phat - (qnorm(0.975) * sqrt(Phat * (1 - Phat) / sum(IR_era$cases)))
PU <- Phat + (qnorm(0.975) * sqrt(Phat * (1 - Phat) / sum(IR_era$cases)))

RRL <- round((PL / (1 - PL)) * (IR_era[IR_era$era == "pre",]$PT_y / 
                                    IR_era[IR_era$era == "post",]$PT_y), 2)
RRU <- round((PU / (1 - PU)) * (IR_era[IR_era$era == "pre",]$PT_y / 
                                    IR_era[IR_era$era == "post",]$PT_y), 2)

IRR; RRL; RRU


# Get pvalue
E1 <- sum(IR_era$cases) * IR_era[IR_era$era == "post",]$PT_y / sum(IR_era$PT_y)
E2 <- sum(IR_era$cases) * IR_era[IR_era$era == "pre",]$PT_y / sum(IR_era$PT_y)

ts <- ((IR_era[IR_era$era == "post",]$cases - E1)^2 / E1) + 
    ((IR_era[IR_era$era == "pre",]$cases - E2)^2 / E2)

pvalue <- pchisq(q = ts, df = 1, lower.tail = F)
pvalue


# Put everything together
IR2 <- IR_era
IR2$pvalue <- IR2$RRU <- IR2$RRL <- IR2$IRR <- NA
IR2$IRR[1] <- IRR
IR2$RRL[1] <- RRL
IR2$RRU[1] <- RRU
IR2$pvalue[1] <- round(pvalue, 3)

IR2


# Remove useless data objects
rm(bt, bt_post, bt_pre, IR_era, IR_post, IR_pre, IRR, E1, E2, i, index, IR2, j, Phat, 
   PL, PU, pvalue, RRL, RRU, ts, vars)



# **** PARAGRAPH 5: ANTI-SPIKE IGG TITER AND BREAKTHROUGH RISK **** ----



# Traditional analysis: IgG titer one-month after the most recent vaccination pre-infection----



# Adapting Viviane's code to get new variables into the ELISA dataset
Viv <- new_elisa %>% 
  dplyr::select(Subject.ID, fb_bt, fb_init_bt_date, vaccine_1_date:vaccine_6_date) %>% 
  unique() %>% 
  group_by(Subject.ID) %>% 
  arrange(Subject.ID, fb_init_bt_date) 


# Get date of 2nd-5th doses
pd2_date <- new_elisa %>%
  subset(time == "Post-dose 2") %>%
  dplyr::select(Subject.ID, Collection.date) %>%
  rename("postdose2.collection.date" = Collection.date) %>%
  unique()

pd3_date <- new_elisa %>%
  subset(time == "Post-dose 3") %>%
  dplyr::select(Subject.ID, Collection.date) %>%
  rename("postdose3.collection.date" = Collection.date) %>%
  unique()

pd4_date <- new_elisa %>%
  subset(time == "Post-dose 4") %>%
  dplyr::select(Subject.ID, Collection.date) %>%
  rename("postdose4.collection.date" = Collection.date) %>%
  unique()

pd5_date <- new_elisa %>%
  subset(time == "Post-dose 5") %>%
  dplyr::select(Subject.ID, Collection.date) %>%
  rename("postdose5.collection.date" = Collection.date) %>%
  unique()


# Combine everything
Viv <- plyr::join_all(list(Viv, pd2_date, pd3_date, pd4_date, pd5_date),
                      by = "Subject.ID", type = 'left')


# Identify the people whose outcome data will be compared at pd2 time point
# Not looking pd1 so don't need to consider the 5 bts after a Janssen dose 1
Viv$pd2_analysis <- if_else(!is.na(Viv$postdose2.collection.date) & # Must have pd2 outcome data 
                              ((Viv$fb_bt == 0) | # Either never had breakthrough OR
                                 (!is.na(Viv$fb_init_bt_date) & # They had a breakthrough
                                    (!is.na(Viv$vaccine_2_date) & ymd(Viv$fb_init_bt_date) > ymd(Viv$vaccine_2_date)) & # Post vaxx 2
                                    # And either didn't have vaccine 3 or they did and the breakthrough was before vaccine 3
                                    (is.na(Viv$vaccine_3_date) | (ymd(Viv$fb_init_bt_date) < ymd(Viv$vaccine_3_date))))), 1, 0)


# Identify the people whose MSD data will be compared at pd3 time point
Viv$pd3_analysis <- if_else(!is.na(Viv$postdose3.collection.date) & # Must have pd3 outcome data 
                              ((Viv$fb_bt == 0) | # Either never had breakthrough OR
                                 (!is.na(Viv$fb_init_bt_date) & # They had a breakthrough
                                    (!is.na(Viv$vaccine_3_date) & ymd(Viv$fb_init_bt_date) > ymd(Viv$vaccine_3_date)) & # Post vaxx 3
                                    # And either didn't have vaccine 4 or they did and the breakthrough was before vaccine 4
                                    (is.na(Viv$vaccine_4_date) | (ymd(Viv$fb_init_bt_date) < ymd(Viv$vaccine_4_date))))), 1, 0)


# Identify the people whose MSD data will be compared at pd4 time point
Viv$pd4_analysis <- if_else(!is.na(Viv$postdose4.collection.date) & # Must have pd2 outcome data 
                              ((Viv$fb_bt == 0) | # Either never had breakthrough OR
                                 (!is.na(Viv$fb_init_bt_date) & # They had a breakthrough
                                    (!is.na(Viv$vaccine_4_date) & ymd(Viv$fb_init_bt_date) > ymd(Viv$vaccine_4_date)) & # Post vaxx 4
                                    # And either didn't have vaccine 5 or they did and the breakthrough was before vaccine 5
                                    (is.na(Viv$vaccine_5_date) | (ymd(Viv$fb_init_bt_date) < ymd(Viv$vaccine_5_date))))), 1, 0)


# Identify the people whose MSD data will be compared at pd5 time point
Viv$pd5_analysis <- if_else(!is.na(Viv$postdose2.collection.date) & # Must have pd2 outcome data 
                              ((Viv$fb_bt == 0) | # Either never had breakthrough OR
                                 (!is.na(Viv$fb_init_bt_date) & # They had a breakthrough
                                    (!is.na(Viv$vaccine_5_date) & ymd(Viv$fb_init_bt_date) > ymd(Viv$vaccine_5_date)) & # Post vaxx 5
                                    # And either didn't have vaccine 6 or they did and the breakthrough was before vaccine 6
                                    (is.na(Viv$vaccine_6_date) | (ymd(Viv$fb_init_bt_date) < ymd(Viv$vaccine_6_date))))), 1, 0)


# Reduce the dataset to the needed variables
Viv <- Viv %>%
  dplyr::select(id = Subject.ID, fb_bt, pd2_analysis:pd5_analysis) 


# Get elisa outcome data
check <- new_elisa %>%
  dplyr::select(id = Subject.ID, time, conc.ug_ml, InforN) %>%
  filter(time %in% c("Post-dose 2", "Post-dose 3", "Post-dose 4", "Post-dose 5"))


# Add in ELISA outcome data
Viv2_elisa <- left_join(Viv, check, by = c("id"))


# Remove unnecessary variables and categories from the MSD dataset
Viv2_elisa <- Viv2_elisa %>%
  ungroup() %>%
  pivot_longer(!c(id, fb_bt, time, conc.ug_ml, InforN),
               names_to = "analysis", values_to = "in_analysis") %>%
  filter(in_analysis == 1) %>%
  dplyr::select(-in_analysis) %>%
  mutate(flag = if_else((analysis == "pd2_analysis" & time == "Post-dose 2") | 
                          (analysis == "pd3_analysis" & time == "Post-dose 3") | 
                          (analysis == "pd4_analysis" & time == "Post-dose 4") |
                          (analysis == "pd5_analysis" & time == "Post-dose 5"), 1, 0)) %>%
  filter(flag == 1) %>%
  dplyr::select(-flag)


# Check that due to ^^ time var now doubles as an indicator of them being in the right analysis
table(Viv2_elisa$time, Viv2_elisa$analysis)

Viv2_elisa <- Viv2_elisa %>%
  dplyr::select(-analysis)


# Add labels for the traditional analysis of the ELISA data
time.labs <- c("Post-dose 2", "Post-dose 3", 
               "Post-dose 4", "Post-dose 5")
names(time.labs) <- c("Post-dose 2", "Post-dose 3", 
                      "Post-dose 4", "Post-dose 5")


# Make a new datset that removes those individuals where n<5 in the relevant panel
Viv2_elisa_new <- Viv2_elisa

Viv2_elisa_new[Viv2_elisa_new$time == "Post-dose 5",]$conc.ug_ml <- NA


# Old way of plotting the ELISA data
# Note that pvalues fail to generate ugh, do it separately instead
set.seed(123)
Viv2_elisa_new %>%
  filter(InforN == 0) %>%
  filter(time != "Post-dose 5") %>%
  mutate(breakthru = factor(fb_bt)) %>%
  ggplot(aes(x = fb_bt, y = conc.ug_ml)) +
  
  geom_hline(yintercept = LOD_elisa, linetype = 2) +
  stat_boxplot(aes(x = breakthru, y = conc.ug_ml, color = breakthru), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = breakthru, y = conc.ug_ml, color = breakthru),
              width = 0.2, alpha = 0.4, size = 2.5) + 
  stat_summary(aes(x = breakthru, y = conc.ug_ml, color = breakthru),
               fun = "median", linewidth = 0.5, geom = "crossbar", show.legend = FALSE) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = seq(0, 65, 10)) +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",
                     label.y.npc = 0.95,
                     label.x.npc = 0.45,
                     hide.ns = F, 
                     size = 5.5) +
  coord_cartesian(ylim = c(0.001, 31)) + 
  facet_wrap(~time, nrow = 3) + 
  ylab("Anti-S IgG (\U003BCg/ml)") +
  scale_color_manual(name = "Status between this dose and the next",
                     values = c("1" = "#d60270",
                                "0" = "black"),
                     labels = c("1" = "Breakthrough",
                                "0" = "No breakthrough")) + 
  theme_fausto() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4.5), nrow = 1))


# Now make three panels instead and combine everything
# Plot the ELISA data
set.seed(123)
plot1_viv <- Viv2_elisa_new %>%
  filter(InforN == 0) %>%
  filter(time == "Post-dose 2") %>%
  mutate(breakthru = factor(fb_bt)) %>%
  ggplot(aes(x = breakthru, y = conc.ug_ml)) +
  
  geom_hline(yintercept = LOD_elisa, linetype = 2) +
  stat_boxplot(aes(x = breakthru, y = conc.ug_ml, color = breakthru), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = breakthru, y = conc.ug_ml, color = breakthru),
              width = 0.2, alpha = 0.4, size = 2.5) + 
  stat_summary(aes(x = breakthru, y = conc.ug_ml, color = breakthru), width = 0.65,
               fun = "median", linewidth = 0.5, geom = "crossbar", show.legend = FALSE) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = seq(0, 65, 10)) +
  scale_x_discrete(breaks = c(0, 1),
                   labels = c("0" = "No breakthrough",
                              "1" = "Breakthrough")) +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",
                     label.y.npc = 0.98,
                     label.x.npc = 0.38,
                     hide.ns = F, 
                     size = 5.5) +
  coord_cartesian(ylim = c(0.001, 31)) + 
  facet_wrap(~time, nrow = 3) + 
  ylab("Anti-S IgG (\U003BCg/ml)") +
  xlab("Between doses 2 and 3") +
  scale_color_manual(name = "Between doses 2 and 3",
                     values = c("1" = "#d60270",
                                "0" = "black"),
                     labels = c("1" = "Breakthrough",
                                "0" = "No breakthrough")) + 
  theme_fausto() +
  theme(legend.position = "none",
        #axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size = 16, face = "plain"),
        axis.text = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 16, face = "bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4.5), nrow = 1))
plot1_viv

set.seed(123)
plot2_viv <- Viv2_elisa_new %>%
  filter(InforN == 0) %>%
  filter(time == "Post-dose 3") %>%
  mutate(breakthru = factor(fb_bt)) %>%
  ggplot(aes(x = breakthru, y = conc.ug_ml)) +
  
  geom_hline(yintercept = LOD_elisa, linetype = 2) +
  stat_boxplot(aes(x = breakthru, y = conc.ug_ml, color = breakthru), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = breakthru, y = conc.ug_ml, color = breakthru),
              width = 0.2, alpha = 0.4, size = 2.5) + 
  stat_summary(aes(x = breakthru, y = conc.ug_ml, color = breakthru), width = 0.65,
               fun = "median", linewidth = 0.5, geom = "crossbar", show.legend = FALSE) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = seq(0, 65, 10)) +
  scale_x_discrete(breaks = c(0, 1),
                   labels = c("0" = "No breakthrough",
                              "1" = "Breakthrough")) +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",
                     label.y.npc = 0.75,
                     label.x.npc = 0.38,
                     hide.ns = F, 
                     size = 5.5) +
  coord_cartesian(ylim = c(0.001, 31)) + 
  facet_wrap(~time, nrow = 3) + 
  ylab("Anti-S IgG (\U003BCg/ml)") +
  xlab("Between doses 3 and 4") +
  scale_color_manual(name = "Between doses 3 and 4",
                     values = c("1" = "#d60270",
                                "0" = "black"),
                     labels = c("1" = "Breakthrough",
                                "0" = "No breakthrough")) + 
  theme_fausto() +
  theme(legend.position = "none",
        #axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size = 16, face = "plain"),
        axis.text = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 16, face = "bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4.5), nrow = 1))
plot2_viv

set.seed(123)
plot3_viv <- Viv2_elisa_new %>%
  filter(InforN == 0) %>%
  filter(time == "Post-dose 4") %>%
  mutate(breakthru = factor(fb_bt)) %>%
  ggplot(aes(x = breakthru, y = conc.ug_ml)) +
  
  geom_hline(yintercept = LOD_elisa, linetype = 2) +
  stat_boxplot(aes(x = breakthru, y = conc.ug_ml, color = breakthru), 
               outlier.shape = NA, show.legend = FALSE, 
               geom ='boxplot', width = 0, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(x = breakthru, y = conc.ug_ml, color = breakthru),
              width = 0.2, alpha = 0.4, size = 2.5) + 
  stat_summary(aes(x = breakthru, y = conc.ug_ml, color = breakthru), width = 0.65,
               fun = "median", linewidth = 0.5, geom = "crossbar", show.legend = FALSE) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = seq(0, 65, 10)) +
  scale_x_discrete(breaks = c(0, 1),
                   labels = c("0" = "No breakthrough",
                              "1" = "Breakthrough")) +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.format",
                     label.y.npc = 0.82,
                     label.x.npc = 0.38,
                     hide.ns = F, 
                     size = 5.5) +
  coord_cartesian(ylim = c(0.001, 31)) + 
  facet_wrap(~time, nrow = 3) + 
  ylab("Anti-S IgG (\U003BCg/ml)") +
  xlab("Between doses 4 and 5") +
  scale_color_manual(name = "Between doses 4 and 5",
                     values = c("1" = "#d60270",
                                "0" = "black"),
                     labels = c("1" = "Breakthrough",
                                "0" = "No breakthrough")) + 
  theme_fausto() +
  theme(legend.position = "none",
        #axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title = element_text(size = 16, face = "plain"),
        axis.text = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 16, face = "bold"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing.x = unit(0, "line"),
        panel.border = theme_border(type = c("bottom", "left"))) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 4.5), nrow = 1))
plot3_viv


# Final plot
ggarrange(
  plot1_viv, plot2_viv, plot3_viv,
  nrow = 3
)


# Remove useless data objects
rm(pd2_date, pd3_date, pd4_date, pd5_date, Viv, Viv2_elisa, time.labs,
   check, plot1_viv, plot2_viv, plot3_viv, Viv2_elisa_new, LOD_elisa)



# Emulated trial: IgG titer one-month after the most recent vaccination pre-infection----


# NOTE: Because noise was introduced into date variables, results here disagree slightly with what is reported in the paper, the latter of which is based on the full, true data.


# Fix the Collected.date variable in the sdat dataset
sdat <- left_join(sdat %>% 
                    select(-c("Collection.date")) %>% 
                    mutate(idt = paste(Subject.ID, time, sep = "-")), 
                  visit %>% 
                    #mutate(idt = paste(`Subject ID`, time, sep = "-")) %>% 
                    select(c("idt", "Collection.date" = "Date")),
                  by = "idt")


# Filter out anyone who had COVID for the first time before BA1
sdat %<>% 
    filter(CovidPos_1 > as.Date("2021-12-15") | 
               is.na(CovidPos_1)) %>% 
    filter(time %in% c("Post-dose 1", "Post-dose 2", "Post-dose 3",
                       "Post-dose 4", "Post-dose 5"))


# Data preparation: Create start and end dates for variant waves, # days from vax to start of window, 
# categorical variables for KM plots
sdat %<>% 
    mutate(tstart = 0,
           start_BA1 = as.Date("2021-12-15"),
           end_BA1 = as.Date("2022-03-25"),
           start_BA2 = as.Date("2022-03-26"),
           end_BA2 = as.Date("2022-07-01"),
           start_BA4.5 = as.Date("2022-07-02"),
           end_BA4.5 = as.Date("2023-01-27"),
           t_vax_to_BA1 = start_BA1 - as.Date(date_vax_before_BA1_begins),
           t_vax_to_BA2 = start_BA2 - as.Date(date_vax_before_BA2_begins),
           t_vax_to_BA4.5 = start_BA4.5 - as.Date(date_vax_before_BA45_begins),
           ab_cat = factor(case_when(conc.ug_ml < 1 ~ "Negative/Low",
                                     conc.ug_ml >= 1 & 
                                         conc.ug_ml < 10 ~ "Medium",
                                     TRUE ~ "High"),
                           labels = c("Ab negative/low", "Ab positive - medium", "Ab positive - high")),
           IDorHV_cat = factor(IDorHV,
                               levels = c("HV", "PID"), 
                               labels = c("Healthy", "IDP")),
           age_cat = factor(case_when(AGE < 16 ~ "<16",
                                      AGE >= 16 &
                                          AGE < 65 ~ "16-64",
                                      TRUE ~ "65+"),
                            labels = c("<16", "16-64", "65+")),
           tv_cat_BA1 = factor(case_when(t_vax_to_BA1 < 184 ~ "1 to <6 months",
                                         t_vax_to_BA1 >= 184 &
                                             t_vax_to_BA1 < 365 ~ "6 months to <12 months",
                                         TRUE ~ "12 or more months"),
                               labels = c("1 to <6 months", "6 months to <12 months", "12 or more months")),
           tv_cat_BA2 = factor(case_when(t_vax_to_BA2 < 184 ~ "1 to <6 months",
                                         t_vax_to_BA2 >= 184 &
                                             t_vax_to_BA2 < 365 ~ "6 months to <12 months",
                                         TRUE ~ "12 or more months"),
                               labels = c("1 to <6 months", "6 months to <12 months", "12 or more months")),
           tv_cat_BA4.5 = factor(case_when(t_vax_to_BA4.5 < 184 ~ "1 to <6 months",
                                           t_vax_to_BA4.5 >= 184 &
                                               t_vax_to_BA4.5 < 365 ~ "6 months to <12 months",
                                           TRUE ~ "12 or more months"),
                                 labels = c("1 to <6 months", "6 months to <12 months", "12 or more months")))


sdat %<>% 
    mutate(CovidPos_1 = as.Date(CovidPos_1),
           BA1_event = case_when(CovidPos_1 > start_BA1 &
                                     CovidPos_1 <= end_BA1 ~ 1,
                                 TRUE ~ 0),
           BA1_event_time = case_when(BA1_event == 1 ~ as.numeric(CovidPos_1 - start_BA1),
                                      TRUE ~ as.numeric(end_BA1 - start_BA1)),
           BA1_t_vax_countup = t_vax_to_BA1 + BA1_event_time,
           BA2_event = case_when(CovidPos_1 > start_BA2 &
                                     CovidPos_1 <= end_BA2 ~ 1,
                                 BA1_event == 1 ~ NA_real_,
                                 TRUE ~ 0),
           BA2_event_time = case_when(BA2_event == 1 ~ as.numeric(CovidPos_1 - start_BA2),
                                      is.na(BA2_event) ~ NA_real_,
                                      TRUE ~ as.numeric(end_BA2 - start_BA2)),
           BA2_t_vax_countup = t_vax_to_BA2 + BA2_event_time,
           BA4.5_event = case_when(CovidPos_1 > start_BA4.5 &
                                       CovidPos_1 <= end_BA4.5 ~ 1,
                                   BA1_event == 1 |
                                       BA2_event == 1 ~ NA_real_,
                                   TRUE ~ 0),
           BA4.5_event_time = case_when(BA4.5_event == 1 ~ as.numeric(CovidPos_1 - start_BA4.5),
                                        is.na(BA4.5_event) ~ NA_real_,
                                        TRUE ~ as.numeric(end_BA4.5 - start_BA4.5)),
           BA4.5_t_vax_countup = t_vax_to_BA4.5 + BA4.5_event_time)


# Subset to variant periods
sdat_BA1 <- sdat %>% 
    filter(Collection.date < start_BA1 &
               Collection.date > date_vax_before_BA1_begins) %>% 
    group_by(Subject.ID) %>% 
    arrange(desc(Collection.date)) %>% 
    slice(1) %>% 
    select(c("Subject.ID", "time", "Collection.date", "conc.ug_ml", "AGE", 
             "IDorHV", "new_group", "CovidPos_1", "ab_cat", "IDorHV_cat", 
             "age_cat", "tv_cat_BA1", "date_vax_before_BA1_begins", "start_BA1", "end_BA1",
             "BA1_event", "tstart", "tstop" = "BA1_event_time", "exposure" = "BA1_t_vax_countup"))

sdat_BA2 <- sdat %>% 
    filter(Collection.date < start_BA2 &
               Collection.date > date_vax_before_BA2_begins &
               !is.na(BA2_event) &
               BA1_event == 0) %>% 
    group_by(Subject.ID) %>% 
    arrange(desc(Collection.date)) %>% 
    slice(1) %>% 
    select(c("Subject.ID", "time", "Collection.date", "conc.ug_ml", "AGE", 
             "IDorHV", "new_group", "CovidPos_1", "ab_cat", "IDorHV_cat", 
             "age_cat", "tv_cat_BA2", "date_vax_before_BA2_begins", "start_BA2", "end_BA2", 
             "BA2_event", "tstart", "tstop" = "BA2_event_time", "exposure" = "BA2_t_vax_countup"))

sdat_BA4.5 <- sdat %>% 
    filter(Collection.date < start_BA4.5 &
               Collection.date > date_vax_before_BA45_begins &
               !is.na(BA4.5_event) &
               BA1_event == 0 &
               BA2_event == 0) %>% 
    group_by(Subject.ID) %>% 
    arrange(desc(Collection.date)) %>% 
    slice(1) %>% 
    select(c("Subject.ID", "time", "Collection.date", "conc.ug_ml", "AGE", 
             "IDorHV", "new_group", "CovidPos_1", "ab_cat", "IDorHV_cat", 
             "age_cat", "tv_cat_BA4.5", "date_vax_before_BA45_begins", "start_BA4.5", "end_BA4.5",
             "BA4.5_event", "tstart", "tstop" = "BA4.5_event_time", "exposure" = "BA1_t_vax_countup"))


# BA1 KM plots
# sdat_BA1 %<>% 
#     mutate(tstart_mos = tstart / 30.5,
#            tstop_mos = tstop / 30.5)

#km_trt_fit_1 <- survfit(Surv(tstart, tstop, BA1_event) ~ ab_cat + cluster(Subject.ID), data = sdat_BA1)


# Get the pvalue
survdiff(Surv(time = (tstop - tstart), BA1_event) ~ ab_cat, data = sdat_BA1)$pvalue 

km_trt_fit_1 <- survfit(Surv(time = (tstop - tstart), BA1_event) ~ ab_cat, data = sdat_BA1)

summary(km_trt_fit_1)


# Make the plot
ggsurvplot(km_trt_fit_1,
           data = sdat_BA1,
           pval = T,
           pval.method = T,
           conf.int = T)

pal <- c("#ffb000", "#648fff", "#DC267F")
km_1 <- autoplot(km_trt_fit_1, 
               censor.shape = '*', 
               censor.size = 2, 
               fun = "event") + 
    labs(x = "Time (days)", 
         y = "Omicron BA.1\ninfection probability", 
         color = "", 
         fill = "") +
    scale_color_manual(values = pal, 
                       labels = c("Ab negative/low", "Ab medium", "Ab high")) +
    scale_fill_manual(values = pal, 
                      labels = c("Ab negative/low", "Ab medium", "Ab high")) +
    theme_fausto() + 
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 0.2)) + 
    theme(panel.border = theme_border(type = c("bottom", "left"))) + 
    annotate("text", x = 50, y = 0.1875, label = "p = 0.93", size = 6.5)

km_1$layers[[1]]$aes_params$size <- 1.25
km_1$layers[[2]]$aes_params$alpha <- 0.175

km_1


# Test the proportional hazards assumption for a Cox regression model fit (coxph).
test_BA1 <- coxph(Surv(time = (tstop - tstart), BA1_event) ~ ab_cat, data = sdat_BA1)

cox.zph(test_BA1) # no issue


# BA2 KM plots
# sdat_BA2 %<>% 
#     mutate(tstart_mos = tstart/30.5,
#            tstop_mos = tstop/30.5)

#km_trt_fit_2 <- survfit(Surv(tstart, tstop, BA2_event) ~ ab_cat + cluster(Subject.ID), data = sdat_BA2) 


# Get the pvalue
survdiff(Surv(time = (tstop - tstart), BA2_event) ~ ab_cat, data = sdat_BA2)$pvalue 

km_trt_fit_2 <- survfit(Surv(time = (tstop - tstart), BA2_event) ~ ab_cat, data = sdat_BA2)

summary(km_trt_fit_2)


# Make the plot
km_2 <- autoplot(km_trt_fit_2, 
                 censor.shape = '*', 
                 censor.size = 2, 
                 fun = "event") + 
    labs(x = "Time (days)", 
         y = "Omicron BA.2\ninfection probability", 
         color = "", 
         fill = "") +
    scale_color_manual(values = pal, 
                       labels = c("Ab negative/low", "Ab medium", "Ab high")) +
    scale_fill_manual(values = pal, 
                      labels = c("Ab negative/low", "Ab medium", "Ab high")) +
    theme_fausto() + 
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 0.4)) + 
    theme(panel.border = theme_border(type = c("bottom", "left"))) + 
    annotate("text", x = 50, y = 0.375, label = "p = 0.24", size = 6.5)

km_2$layers[[1]]$aes_params$size <- 1.25
km_2$layers[[2]]$aes_params$alpha <- 0.175

km_2

# Test the proportional hazards assumption for a Cox regression model fit (coxph).
test_BA2 <- coxph(Surv(time = (tstop - tstart), BA2_event) ~ ab_cat, data = sdat_BA2)

cox.zph(test_BA2) # 0.034


# BA4.5 KM plots
# sdat_BA4.5 %<>% 
#     mutate(tstart_mos = tstart / 30.5,
#            tstop_mos = tstop / 30.5)
# 
# km_trt_fit_3 <- survfit(Surv(tstart, tstop, BA4.5_event) ~ ab_cat + cluster(Subject.ID), data = sdat_BA4.5) 


# Get the pvalue
survdiff(Surv(time = (tstop - tstart), BA4.5_event) ~ ab_cat, data = sdat_BA4.5)$pvalue 

km_trt_fit_3 <- survfit(Surv(time = (tstop - tstart), BA4.5_event) ~ ab_cat, data = sdat_BA4.5)

summary(km_trt_fit_3)


# Make the plot
km_3 <- autoplot(km_trt_fit_3, 
                 censor.shape = '*', 
                 censor.size = 2, 
                 fun = "event") + 
    labs(x = "Time (days)", 
         y = "Omicron BA.4/5\ninfection probability", 
         color = "", 
         fill = "") +
    scale_color_manual(values = pal, 
                       labels = c("Ab negative/low", "Ab medium", "Ab high")) +
    scale_fill_manual(values = pal, 
                      labels = c("Ab negative/low", "Ab medium", "Ab high")) +
    theme_fausto() + 
    coord_cartesian(xlim = c(0, 200), ylim = c(0, 0.4)) + 
    theme(panel.border = theme_border(type = c("bottom", "left"))) + 
    annotate("text", x = 100, y = 0.375, label = "p = 0.30", size = 6.5)

km_3$layers[[1]]$aes_params$size <- 1.25
km_3$layers[[2]]$aes_params$alpha <- 0.175

km_3


# Test the proportional hazards assumption for a Cox regression model fit (coxph).
test_BA3 <- coxph(Surv(time = (tstop - tstart), BA4.5_event) ~ ab_cat, data = sdat_BA4.5)

cox.zph(test_BA3) # no issue


# Put it all together
patchwork <- km_1 / km_2 / km_3 + 
    guide_area() +
    plot_layout(guides = 'collect', 
                axis_titles = "collect",
                ncol = 1,
                heights = c(4, 4, 4, -1)) +
    theme_update(
        legend.position = "bottom",
        legend.direction = "horizontal" , 
        panel.grid.major = element_line(
            rgb(105, 105, 105, maxColorValue = 255),
            linetype = "dotted"),   
        panel.grid.minor = element_line(
            rgb(105, 105, 105, maxColorValue = 255),
            linetype = "dotted", 
            size = rel(4)),
        axis.text = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.5),
                                    hjust = 0),
        axis.title = element_text(size = rel(2)),
        plot.margin = unit(c(1,5,1,1), "mm")
    )

patchwork # 8.5 by 14 inches does the trick


# Remove useless data objects
rm(km_1, km_2, km_3, km_trt_fit_1, km_trt_fit_2, km_trt_fit_3, patchwork, pal, 
   sdat, sdat_BA1, sdat_BA2, sdat_BA4.5, test_BA1, test_BA2, test_BA3, visit)



# **** PARAGRAPH 6: BEHAVIORAL PRACTICES AND POLICY PREFERENCES **** ----



# Comparing infection-prevention attitudes and behaviors ----



# Extract needed variables and perform limited data cleaning
behavior <- behavior_raw %>%
  mutate(IDorHV = case_when(grepl("E", id, ignore.case = T) ~ "IDP",
                            grepl("H", id, ignore.case = T) ~ "HV",
                            TRUE ~ NA_character_)) %>%
  filter(!is.na(before.washed.hands.more)) %>% # Remove unanswered surveys 293 -> 153
  rename(Subject.ID = id)


# Completion rate of the behavioral survey
round(100 * length(unique(behavior$Subject.ID)) / length(unique(new_elisa$Subject.ID)), 1)


# Apply factor levels to different types of variables
# Specify the variables to apply the order
tf_variables <- c("before.washed.hands.more", "before.keep.6.feet.distance", 
                  "before.reduce.trips", "before.wore.mask", 
                  
                  "after.washed.hands.more", "after.telework.more.often", 
                  "after.keep.6.feet.distance", "after.reduce.trips", 
                  "after.wore.mask", "after.behaviors.returned")

approve_variables <- 
  c("before.limit.spread.masking", "before.limits.person.worship", "before.limits.indoor.dining",
    "before.non.essential.travel", "before.closure.schools", "before.mandatory.vaccination",
    "before.vaccine.passports", 
    
    "after.limit.spread.masking", "after.limits.person.worship", "after.limits.indoor.dining", 
    "after.non.essential.travel", "after.closure.schools", "after.increased.telework", 
    "after.mandatory.vaccination", "after.vaccine.passports")

desired_tf_order <- c("Very true for me", "Somewhat true for me", "Slightly true for me", 
                      "Slightly false for me", "Somewhat false for me", "Very false for me")

behavior <- behavior %>%
  mutate(across(all_of(tf_variables), ~ fct_relevel(., desired_tf_order)))

desired_approve_order <- c("Strongly approve", "Somewhat approve", "Slightly approve", 
                           "Slightly disapprove", "Somewhat disapprove", "Strongly disapprove")

behavior <- behavior %>%
  mutate(across(all_of(approve_variables), ~ fct_relevel(., desired_approve_order)))


# Add in the breakthrough variable
part_char_short <- new_elisa %>%
  select(Subject.ID, breakthru = fb_bt) %>%
  group_by(Subject.ID) %>%
  slice(1) %>%
  ungroup()

behavior <- behavior %>%
  filter(Subject.ID %in% part_char_short$Subject.ID) %>% # 153 -> 150
  left_join(part_char_short, by = "Subject.ID")


# Create new variables for generally pro and generally against
behavior_plus <- behavior %>%
  mutate(before.washed.hands.more_binary = 
           case_when(grepl("true", before.washed.hands.more) ~ "True", 
                     TRUE ~ "False"),
         before.keep.6.feet.distance_binary = 
           case_when(grepl("true", before.keep.6.feet.distance) ~ "True",
                     TRUE ~ "False"),
         before.reduce.trips_binary = 
           case_when(grepl("true", before.reduce.trips) ~ "True",
                     TRUE ~ "False"),
         before.wore.mask_binary = 
           case_when(grepl("true", before.wore.mask) ~ "True",
                     TRUE ~ "False"),
         before.limit.spread.masking_binary = 
           case_when(grepl("approve", before.limit.spread.masking) ~ "Approve",
                     TRUE ~ "Disapprove"),
         before.limits.person.worship_binary =
           case_when(grepl("approve", before.limits.person.worship) ~ "Approve",
                     TRUE ~ "Disapprove"),
         before.limits.indoor.dining_binary = 
           case_when(grepl("approve", before.limits.indoor.dining) ~ "Approve",
                     TRUE ~ "Disapprove"),
         before.non.essential.travel_binary = 
           case_when(grepl("approve", before.non.essential.travel) ~ "Approve",
                     TRUE ~ "Disapprove"),
         before.closure.schools_binary = 
           case_when(grepl("approve", before.closure.schools) ~ "Approve",
                     TRUE ~ "Disapprove"),
         before.increased.telework_binary = 
           case_when(grepl("approve", before.increased.telework) ~ "Approve",
                     TRUE ~ "Disapprove"),
         before.mandatory.vaccination_binary = 
           case_when(grepl("approve", before.mandatory.vaccination) ~ "Approve", 
                     TRUE ~ "Disapprove"),
         before.vaccine.passports_binary = 
           case_when(grepl("approve", before.vaccine.passports) ~ "Approve",
                     TRUE ~ "Disapprove"))

behavior_plus <- behavior_plus %>%
  mutate(after.washed.hands.more_binary = 
           case_when(grepl("true", after.washed.hands.more) ~ "True",
                     TRUE ~ "False"),
         after.keep.6.feet.distance_binary = 
           case_when(grepl("true", after.keep.6.feet.distance) ~ "True",
                     TRUE ~ "False"),
         after.reduce.trips_binary = 
           case_when(grepl("true", after.reduce.trips) ~ "True",
                     TRUE ~ "False"),
         after.telework.more.often_binary = 
           case_when(grepl("true", after.telework.more.often) ~ "True",
                     TRUE ~ "False"),
         after.wore.mask_binary = 
           case_when(grepl("true", after.wore.mask) ~ "True",
                     TRUE ~ "False"),
         after.limit.spread.masking_binary = 
           case_when(grepl("approve", after.limit.spread.masking) ~ "Approve",
                     TRUE ~ "Disapprove"),
         after.limits.person.worship_binary = 
           case_when(grepl("approve", after.limits.person.worship) ~ "Approve",
                     TRUE ~ "Disapprove"),
         after.limits.indoor.dining_binary = 
           case_when(grepl("approve", after.limits.indoor.dining) ~ "Approve",
                     TRUE ~ "Disapprove"),
         after.non.essential.travel_binary = 
           case_when(grepl("approve", after.non.essential.travel) ~ "Approve",
                     TRUE ~ "Disapprove"),
         after.closure.schools_binary = 
           case_when(grepl("approve", after.closure.schools) ~ "Approve",
                     TRUE ~ "Disapprove"),
         after.increased.telework_binary = 
           case_when(grepl("approve", after.increased.telework) ~ "Approve",
                     TRUE ~ "Disapprove"),
         after.mandatory.vaccination_binary = 
           case_when(grepl("approve", after.mandatory.vaccination) ~ "Approve",
                     TRUE ~ "Disapprove"),
         after.vaccine.passports_binary = 
           case_when(grepl("approve", after.vaccine.passports) ~ "Approve",
                     TRUE ~ "Disapprove"),
         after.behaviors.returned_binary = 
           case_when(grepl("true", after.behaviors.returned) ~ "True",
                     TRUE ~ "False"))


# Focus on the post-vaccination data and perform some data cleaning
behavior_after <- behavior_plus %>%
  select("Subject.ID", "IDorHV", "breakthru", "after.washed.hands.more_binary", 
         "after.keep.6.feet.distance_binary", "after.reduce.trips_binary", 
         "after.telework.more.often_binary", "after.wore.mask_binary", 
         "after.limit.spread.masking_binary", "after.limits.person.worship_binary", 
         "after.limits.indoor.dining_binary", "after.non.essential.travel_binary", 
         "after.closure.schools_binary", "after.increased.telework_binary", 
         "after.mandatory.vaccination_binary", "after.vaccine.passports_binary", 
         "after.behaviors.returned_binary")

behavior_after <- behavior_after %>%
  pivot_longer(cols = matches("^after"), names_to = "behavior", values_to = "perception")

behavior_after$behavior <- factor(behavior_after$behavior, 
                                  levels = c("after.washed.hands.more_binary", 
                                             "after.keep.6.feet.distance_binary", 
                                             "after.reduce.trips_binary", 
                                             "after.telework.more.often_binary", 
                                             "after.wore.mask_binary", 
                                             "after.behaviors.returned_binary",
                                             "after.limit.spread.masking_binary", 
                                             "after.limits.person.worship_binary", 
                                             "after.limits.indoor.dining_binary", 
                                             "after.non.essential.travel_binary", 
                                             "after.closure.schools_binary", 
                                             "after.increased.telework_binary", 
                                             "after.mandatory.vaccination_binary", 
                                             "after.vaccine.passports_binary"))

behavior_after <- behavior_after %>%
  group_by(IDorHV, behavior, perception) %>%
  summarise(count = n()) %>%
  mutate(total = sum(count)) %>%
  mutate(prop = count / total)

behavior_after <- behavior_after %>%
  filter(perception != "False") %>%
  filter(perception != "Disapprove")


# Make better variable names
behavior_after_fb <- behavior_after %>%
  mutate(behavior2 = 
           case_when(behavior == "after.washed.hands.more_binary" ~ "Washed hands more",
                     behavior == "after.telework.more.often_binary" ~ "WFH more often",
                     behavior == "after.keep.6.feet.distance_binary" ~ "Kept social distancing",
                     behavior == "after.reduce.trips_binary" ~ "Fewer trips outside",
                     behavior == "after.wore.mask_binary" ~ "Kept wearing masks",
                     behavior == "after.behaviors.returned_binary" ~ "Pre-pandemic behaviors returned",
                     
                     behavior == "after.limit.spread.masking_binary" ~ "Mandatory masking",
                     behavior == "after.limits.person.worship_binary" ~ "Limit in-person worship",
                     behavior == "after.limits.indoor.dining_binary" ~ "Limit indoor dining",
                     behavior == "after.non.essential.travel_binary" ~ "Limit non-essential travel",
                     behavior == "after.closure.schools_binary" ~ "Close schools",
                     behavior == "after.increased.telework_binary" ~ "Increase teleworking",
                     behavior == "after.mandatory.vaccination_binary" ~ "Mandatory vax",
                     behavior == "after.vaccine.passports_binary" ~ "Vax passports",
                     
                     TRUE ~ "oh nauuur"))


# Estimate pvalues
behavior_after_fb <- behavior_after_fb %>%
  arrange(behavior, IDorHV)

behavior_after_fb$pval <- NA
for(i in 1:nrow(behavior_after_fb)){
    
  if(i %% 2 == 1){behavior_after_fb$pval[i] <- 
    round(prop.test(x = c(behavior_after_fb$count[i], behavior_after_fb$count[i + 1]),
                    n = c(behavior_after_fb$total[i], behavior_after_fb$total[i + 1]))$p.value, 5)}
    
  if(i %% 2 != 1){behavior_after_fb$pval[i] <- behavior_after_fb$pval[i - 1]}
}

behavior_after_fb %>%
  filter(pval < 0.05)


# Factor the behavior variable
behavior_after_fb$behavior2 <- factor(behavior_after_fb$behavior2, 
                                      levels = c("Washed hands more", "WFH more often", 
                                                 "Kept social distancing", "Fewer trips outside",
                                                 "Kept wearing masks", "Pre-pandemic behaviors returned",
                                                 
                                                 "Close schools", "Limit in-person worship",
                                                 "Limit indoor dining", "Limit non-essential travel",
                                                 "Mandatory masking", "Increase teleworking",
                                                 "Mandatory vax", "Vax passports"))


# Make the heatmap
ggplot(behavior_after_fb, aes(y = IDorHV, x = behavior2, fill = prop)) + 
  theme_fausto() + 
  geom_tile() + 
  scale_y_discrete(breaks = c("HV", "IDP"),
                   labels = c("HV" = "HV",
                              "IDP" = "IDP"),
                   expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_fill_viridis(option = "A",
                     begin = 0, end = 1,
                     name = "Agreement (%) ",
                     limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "25", "50", "75", "100")) + 
  guides(fill = guide_colourbar(barheight = 1,
                                barwidth = 40)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black", 
                                   size = 25),
        axis.text.y = element_text(angle = 90, hjust = 0.5, face = "bold", 
                                   color = "black", size = 25),
        legend.text = element_text(angle = 90, hjust = 1, size = 25),
        legend.title = element_text(size = 25),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) # Export as PDF with length 19in using US letter landscape


# Remove useless data objects
rm(behavior, behavior_after, behavior_after_fb, behavior_plus, 
   approve_variables, desired_approve_order, desired_tf_order, i, tf_variables, part_char_short)



# **** PARAGRAPH 7: ADVERSE OUTCOMES OF SARS-COV-2 INFECTIONS **** ----



# N unusual outcomes ----



# Persistent infections
round(100 * sum(outcomes[outcomes$IDorHV == "PID",]$persistent_inf) / 
        nrow(outcomes[outcomes$IDorHV == "PID",]), 1) # 4.6%

round(100 * sum(outcomes[outcomes$IDorHV == "HV",]$persistent_inf) / 
        nrow(outcomes[outcomes$IDorHV == "HV",]), 1)# 3.2%


# PASC
round(100 * sum(outcomes[outcomes$IDorHV == "PID",]$Long_Covid) / 
        nrow(outcomes[outcomes$IDorHV == "PID",]), 1) # 4.6%

round(100 * sum(outcomes[outcomes$IDorHV == "HV",]$Long_Covid) / 
        nrow(outcomes[outcomes$IDorHV == "HV",]), 1) # 3.2%


# Reinfections
round(100 * sum(outcomes[outcomes$IDorHV == "PID",]$reinfected) / 
        nrow(outcomes[outcomes$IDorHV == "PID",]), 1) # 4.1%

round(100 * sum(outcomes[outcomes$IDorHV == "HV",]$reinfected) / 
        nrow(outcomes[outcomes$IDorHV == "HV",]), 1) # 1.6%


# Testing the proportions of adverse immunological outcomes be IDP and HV
prop.test(x = c(10, 2), n = c(219, 63)) # Persistent infections
prop.test(x = c(10, 2), n = c(219, 63)) # PACs
prop.test(x = c(9, 1), n = c(219, 63)) # Re-infections



# Pre-dose 1 infections ----



# How many pre-dose 1 infections where there?
pre <- new_elisa %>%
  dplyr::select(id = Subject.ID, fb_predose1_inf) %>%
  filter(fb_predose1_inf == 1) %>%
  group_by(id) %>%
  slice(1) 

nrow(pre) # 10


# How many pre-dose 1 infections had infections after the completion of their primary series?
post <- new_elisa %>%
  filter(Subject.ID %in% pre$id) %>%
  select(id = Subject.ID, vaccine_1_date, vaccine_2_date, vaccine_1, vaccine_2, 
         CovidPos_2, CovidPos_3) %>% # No Janssen so primary series is after 2nd dose
  mutate(flag = case_when(ymd(CovidPos_2) > ymd(vaccine_2_date) ~ 1,
                          ymd(CovidPos_3) > ymd(vaccine_2_date) ~ 1,
                          TRUE ~ 0)) %>%
  group_by(id) %>%
  slice(1)

sum(post$flag) # 4

rm(post, pre)



# Waffle plot of uncommon outcomes ----



# Simplify the dataset
waffle2 <- outcomes %>%
  dplyr::select(id = Subject.ID, breakthru, ever_inf = ever_infected, long_covid = Long_Covid, 
                persistent_inf, mult_bt = reinfected, predose1_inf) 


# Compare the proportions of adverse outcomes among the initial breakthrough infections, not overall
checknow <- outcomes %>% # The initial breakthroughs
    filter(breakthru == 1)

table2(checknow$Long_Covid)
table(checknow$Long_Covid, checknow$IDorHV)
prop.test(x = c(2, 10), n = c(26 + 2, 78 + 10))

table2(checknow$persistent_inf)
table(checknow$persistent_inf, checknow$IDorHV)
prop.test(x = c(2, 10), n = c(26 + 2, 78 + 10))


table2(checknow$reinfected)
table(checknow$reinfected, checknow$IDorHV)
prop.test(x = c(1, 9), n = c(1 + 27, 79 + 9))


# Add in single infected variable
waffle2$single_bt <- 0
waffle2[is.na(waffle2)] <- 0
waffle2[waffle2$breakthru == 1 & waffle2$mult_bt == 0,]$single_bt <- 1

# Reduce dataset
waffle2 <- waffle2 %>%
  dplyr::select(-ever_inf, -breakthru, -id)

# Sum all the combinations
waffle3 <- waffle2 %>%
  group_by(predose1_inf, single_bt, mult_bt, long_covid, persistent_inf) %>%
  dplyr::summarise(count = n(),
                   percent = round(100 * count / nrow(waffle2), 1))
waffle3

waffle3$name <- c("Never infected\nn=160 (56.7%)", 
                  "Multiple breakthroughs\nn=4 (1.4%)",
                  "Multiple breakthroughs & persistent infection\nn=1 (0.4%)", 
                  "Multiple breakthroughs & PASC\nn=3 (1.1%)",
                  "Multiple breakthroughs, persistent infection, & PASC\nn=1 (0.4%)", 
                  "Single breakthrough\nn=85 (30.1%)",
                  "Single breakthrough & persistent infection\nn=10 (3.5%)", 
                  "Single breakthrough & PASC\nn=8 (2.8%)",
                  "Pre-dose 1 infection\nn=6 (2.1%)",
                  "Pre-dose 1 infection & multiple breakthroughs\nn=1 (0.4%)", 
                  "Pre-dose 1 infection & single breakthrough\nn=3 (1.1%)")
waffle3$name <- as.factor(waffle3$name)

# Colors from https://mk.bcgsc.ca/colorblind/palettes.mhtml#12-color-palette-for-colorbliness
waffle3$colors <- c("#FFC33B", "#FF6E3A", "#E20134", "#A40122", "#FFB2FD",
                    "#6DE0FF", "#008DF9", "#8400CD", "#FF5AAF", "#009F81",
                    "#000000")

waffle4 <- waffle3 %>%
  ungroup() %>%
  dplyr::select(names = name, vals = count, colors)

waffle5 <- waffle(waffle4, rows = 17, colors = waffle4$colors, legend_pos = "bottom") + 
  theme(legend.position = "right",
        legend.text = element_text(size = 16),
        legend.key.size = unit(2.75, "lines"),
        legend.spacing.x = unit(0.30, "lines"),
        legend.spacing.y = unit(0.26, "lines"),
        legend.box.spacing = unit(40, "pt")) + 
  guides(fill = guide_legend(byrow = TRUE))
waffle5


# Remove useless data objects
rm(waffle2, waffle3, waffle4, waffle5, checknow)



# Analysis of persistent infections ----



# Identify the persistent positives
pp <- pp_raw %>%
    filter(persistant_inf == 1) %>%
    select(SUBJECT.ID., persistant_inf, persistant_date_began, persistant_number_days) %>%
    rename(Subject.ID = SUBJECT.ID.)

pp$persistant_date_began <- as.Date(pp$persistant_date_began, format = "%m/%d/%Y")

pp_IDs <- pp %>%
    pull(Subject.ID)


# Combine with full dataset to confirm they are not pre-dose 1 infections
infections <- inf %>%
    left_join(pp, by = c("SUBJECT.ID." = "Subject.ID"))

i_pp <- infections %>%
    filter(persistant_inf == 1) %>%
    dplyr::select(SUBJECT.ID., infection_number, infection_timing, persistant_inf,
                  persistant_date_began, persistant_date_began)


# Confirmed they are not - don't need to use this dataset ^
ELISA <- new_elisa %>%
    left_join(pp, by = "Subject.ID")


# Edit ELISA dataset
ELISA_match <- ELISA %>%
    dplyr:: select(Subject.ID, AGE, agecat, new_group, IDorHV, breakthru, CovidPos_1, conc.ug_ml,
                   Collection.date, time, COVID_pos_before_sample, persistant_date_began,
                   persistant_inf) %>%
    mutate(timepoint2 = case_when(time %in% 
                                      c("Post-dose 1", "Post-dose 2", "Post-dose 3", 
                                        "Post-dose 4", "Post-dose 5") ~ time,
                                  TRUE ~ "other timepoint 2")) %>%
    mutate(persistant_inf = case_when(Subject.ID %in% pp_IDs ~ 1,
                                      TRUE ~ 0))


# Remove people with any infection (though include people with persistent) for the match
ELISA_match <- ELISA_match %>%
    filter(persistant_inf == 1 | is.na(CovidPos_1))


# Remove bad timepoints
ELISA_match <- ELISA_match %>%
    filter(timepoint2 != "other timepoint 2")

table(ELISA_match$persistant_inf)


# Conduct the matching
Pers_M <- MatchIt::matchit(persistant_inf ~ AGE + new_group + timepoint2,
                           data = ELISA_match, 
                           replace = FALSE,
                           distance = "robust_mahalanobis",
                           exact = c("timepoint2"),
                           caliper = c(AGE = 15),
                           std.caliper = c(FALSE),
                           ratio = 2)

Pers_M <- MatchIt::match.data(Pers_M)


# Make the plot
my_comparisons <- list(c("0", "1"))
Pers_M$persistant_inf <- factor(Pers_M$persistant_inf)
Pers_M$COVID_pos_before_sample <- factor(Pers_M$COVID_pos_before_sample, levels =  c(0, 1), 
                                         labels = c("Before SARS-CoV-2 infection", 
                                                    "After SARS-CoV-2 infection"))

persistent <- Pers_M %>%
    ggplot(aes(x = persistant_inf, y = conc.ug_ml)) +
    geom_boxplot(outlier.shape = NA, 
                 show.legend = FALSE,
                 position = position_dodge(0.7, preserve = "single")) + 
    geom_jitter(aes(x = persistant_inf, 
                    y = conc.ug_ml, 
                    color = persistant_inf,
                    shape = COVID_pos_before_sample),
                width = 0.2,
                size = 3) + 
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    #  stat_compare_means(comparisons = my_comparisons, # since the comparison is non-signficnat, remove the line
    #                     method = "wilcox.test",
    #                     label = "p.signif") +
    facet_grid(. ~ COVID_pos_before_sample,
               switch = "x",
               scales = "free_x",
               space = "free_x") +
    # theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
    ggtitle("Persistent Positive Infections") +
    xlab("") +  # Remove the x-axis label
    ylab("IgG Levels") +
    scale_color_manual(values = c("#41b6c4", "plum3"), 
                       labels = c("Never had SARS-CoV-2 infection", "Persistent SARS-CoV-2 infection")) +
    #  scale_shape_manual(values = c("circle", "triangle"), labels = c("Before SARS-CoV-2 infection", "After SARS-CoV-2 infection")) +
    theme_bw() +
    theme(
        legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_text(size = 18),  # Remove the x-axis label
        axis.text.x = element_blank(),  # Remove the x-axis tick mark labels
        axis.title.y = element_text(size = 18),  
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20), 
        legend.text = element_text(size = 13),# Increase the plot title font size
        strip.text = element_text(size = 11.2)
    )  +
    guides(color = guide_legend(override.aes = list(size = 4)))

persistent

ggsave(filename = "persistent.jpg", plot = persistent, device = "jpeg", dpi = 300,
       width = 8, height = 4, units = "in")


# Set up to compare means of the before group, filter the data for the "pre-infection" group
my_comparisons <- list(c("Yes", "No"))

pre_infection_data <- Pers_M %>%
    filter(COVID_pos_before_sample == "Before SARS-CoV-2 infection")


# Check if there are at least two groups to compare in the "pre-infection" group
if (length(unique(pre_infection_data$persistant_inf)) == 2) {
    test_result <- wilcox.test(conc.ug_ml ~ persistant_inf, data = pre_infection_data)
    p_value_pre <- test_result$p.value
    print(paste("Wilcoxon test p-value (Before) for Yes vs No:", p_value_pre))
} else {
    print("Not enough groups to compare in the 'pre-infection' group")
}


# Identify the number of days of persistence for ID and HV
# Calculate the mean and standard deviation of conc.ug_ml based on persistant_inf
pp <- pp %>%
    mutate(IDorHV = case_when(grepl("E", Subject.ID) ~ "ID",
                              TRUE ~ "HV"))

stats_summary <- pp %>%
    group_by(IDorHV) %>%
    summarise(
        mean_conc = mean(persistant_number_days, na.rm = TRUE),
        sd_conc = sd(persistant_number_days, na.rm = TRUE)
    )


# Print the summary statistics
print(stats_summary)


# Remove unless data objects
rm(ELISA, ELISA_match, i_pp, my_comparisons, Pers_M, persistent, pp, pp_raw, pre_infection_data,
   stats_summary, p_value_pre, pp_IDs, test_result, infections)



# Analysis of PASC ----



# Identify PASC cases 
## Calculate illness duration
infections <- inf %>% 
    mutate(symptom_duration = as.numeric(ymd(recovery_date) - ymd(symptom_onset_date)))

# Subset to those with illness duration greater than 28 days
Long_Covid <- infections %>% 
    filter(symptom_duration > 28) %>% 
    filter(infection_timing == "Breakthrough") %>%
    select(SUBJECT.ID., covid_date, infection_number, infection_timing, symptom_onset_date,
           recovery_date, symptom_duration)

LongCovidIDs <- Long_Covid$SUBJECT.ID.

Long_Covid <- Long_Covid %>%
    mutate(Long_Covid_date_began = covid_date) %>%
    select(!c(infection_number, infection_timing, symptom_onset_date, recovery_date))

Long_Covid <- Long_Covid %>%
    mutate(Long_Covid = case_when(SUBJECT.ID. %in% LongCovidIDs ~ 1,
                                  TRUE ~ 0)) %>%
    select(!covid_date) %>%
    rename(Subject.ID = SUBJECT.ID.)

ELISA <- new_elisa %>%
    left_join(Long_Covid, by = "Subject.ID")


# Edit ELISA data and prepare it to match on age, new_group, and timepoint2
ELISA_match <- ELISA %>%
    select(Subject.ID, AGE, agecat, new_group, IDorHV, breakthru, CovidPos_1, conc.ug_ml, 
           Collection.date, time,
           COVID_pos_before_sample, Long_Covid_date_began, Long_Covid) %>%
    mutate(timepoint2 = case_when(time %in% 
                                      c("Post-dose 1", "Post-dose 2", "Post-dose 3", "Post-dose 4",
                                        "Post-dose 5") ~ time,
                                  TRUE ~ "other timepoint 2")) %>%
    mutate(Long_Covid = case_when(Subject.ID %in% LongCovidIDs ~ 1,
                                  TRUE ~ 0))

# Remove people with any infection (though include people with LC) for the match
ELISA_match <- ELISA_match %>%
    filter(is.na(CovidPos_1) | Long_Covid == 1)

table(ELISA_match$time)
table(ELISA_match$timepoint2)
table(ELISA_match$Long_Covid)


# Filter out some folks
ELISA_match <- ELISA_match %>%
    filter(timepoint2 != "other timepoint 2")

# Check the key variables
table(ELISA_match$Long_Covid, ELISA_match$timepoint2)


# Conduct the matching
LC_M <- MatchIt::matchit(Long_Covid ~ AGE + new_group + timepoint2,
                         data = ELISA_match, 
                         replace = FALSE,
                         distance = "robust_mahalanobis",
                         exact = c("timepoint2"),
                         caliper = c(AGE = 15),
                         std.caliper = c(FALSE),
                         ratio = 2)

LC_M <- MatchIt::match.data(LC_M)


# Plot the timepoints 
LC_M$COVID_pos_before_sample <- factor(LC_M$COVID_pos_before_sample)
my_comparisons <- list(c("0", "1"))
LC_M$Long_Covid <- factor(LC_M$Long_Covid, levels =  c(0, 1), 
                          labels = c("Before SARS-CoV-2 infection", "After SARS-CoV-2 infection"))

LC_M %>%
    ggplot(aes(x = Long_Covid, y = conc.ug_ml)) +
    geom_boxplot(outlier.shape = NA, 
                 show.legend = FALSE,
                 position = position_dodge(0.7, preserve = "single")) + 
    geom_jitter(aes(x = Long_Covid, 
                    y = conc.ug_ml, 
                    color = Long_Covid,
                    shape = COVID_pos_before_sample),
                width = 0.2,
                size = 3) + 
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    #  stat_compare_means(comparisons = my_comparisons, 
    #                     method = "wilcox.test",
    #                     label = "p.signif") +
    facet_grid(. ~ timepoint2,
               switch = "x",
               scales = "free_x",
               space = "free_x") +
    # theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
    ggtitle("IgG Levels in Participants with Long Covid Compared to Participants with No Infections \nMatched on Age, Immune Sub-group, and Timepoint") +
    xlab("") +  # Remove the x-axis label
    ylab("IgG Levels (conc.ug_ml)") +
    scale_color_manual(values = c("#41b6c4", "tomato1"), 
                       labels = c("Never had Covid", "Long Covid")) +
    scale_shape_manual(values = c("circle", "triangle"), 
                       labels = c("Prior to First Covid Infection", "After First Covid Infection")) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank(),  # Remove the x-axis label
        axis.text.x = element_blank(),  # Remove the x-axis tick mark labels
        axis.title.y = element_text(size = 12),  
        axis.text.y = element_text(size = 15),
        plot.title = element_text(size = 20), 
        legend.text = element_text(size = 11),# Increase the plot title font size
        strip.text = element_text(size = 11.2)
    )  


# Graph the before/after
LC_M$COVID_pos_before_sample <- factor(LC_M$COVID_pos_before_sample, levels =  c(0, 1), 
                                       c("Before SARS-CoV-2 infection", "After SARS-CoV-2 infection"))

PASC <- LC_M %>%
    ggplot(aes(x = Long_Covid, y = conc.ug_ml)) +
    geom_boxplot(outlier.shape = NA, 
                 show.legend = FALSE,
                 position = position_dodge(0.7, preserve = "single")) + 
    geom_jitter(aes(x = Long_Covid, 
                    y = conc.ug_ml, 
                    color = Long_Covid,
                    shape = COVID_pos_before_sample),
                width = 0.2,
                size = 3) + 
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    #  stat_compare_means(comparisons = my_comparisons, 
    #                     method = "wilcox.test",
    #                     label = "p.signif") +
    facet_grid(. ~ COVID_pos_before_sample,
               switch = "x",
               scales = "free_x",
               space = "free_x") +
    # theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
    ggtitle("PASC Infections") +
    xlab("") +  # Remove the x-axis label
    ylab("IgG Levels") +
    scale_color_manual(values = c("#41b6c4", "tomato1"), 
                       labels = c("Never had SARS-CoV-2 infection", "PASC infection")) +
    # scale_shape_manual(values = c("circle", "triangle"), labels = c("Prior to First Covid Infection", "After First Covid Infection")) +
    theme_bw() +
    theme(
        legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_text(size = 18),  # Remove the x-axis label
        axis.text.x = element_blank(),  # Remove the x-axis tick mark labels
        axis.title.y = element_text(size = 18),  
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20), 
        legend.text = element_text(size = 13),# Increase the plot title font size
        strip.text = element_text(size = 11.2)
    )  +
    guides(color = guide_legend(override.aes = list(size = 4)))
PASC


# Save the figure
ggsave(filename = "PASC.jpg", plot = PASC, device = "jpeg", dpi = 300, width = 8, height = 3.5, units = "in")



# Compare means of the before group
## First define comparisons
my_comparisons <- list(c("Yes", "No"))


# Filter the data for the "pre-infection" group
pre_infection_data <- LC_M %>%
    filter(COVID_pos_before_sample == "Before SARS-CoV-2 infection")


# Check if there are at least two groups to compare in the "pre-infection" group
if (length(unique(pre_infection_data$Long_Covid)) == 2) {
    test_result <- wilcox.test(conc.ug_ml ~ Long_Covid, data = pre_infection_data)
    p_value_pre <- test_result$p.value
    print(paste("Wilcoxon test p-value (Before) for Yes vs No:", p_value_pre))
} else {
    print("Not enough groups to compare in the 'pre-infection' group")
}


# Calculate the mean and standard deviation of conc.ug_ml based on PASC
Long_Covid <- Long_Covid %>%
    mutate(IDorHV = case_when(grepl("E", Subject.ID) ~ "ID",
                              TRUE ~ "HV"))

stats_summary <- Long_Covid %>%
    group_by(IDorHV) %>%
    summarise(
        mean_conc = mean(symptom_duration, na.rm = TRUE),
        sd_conc = sd(symptom_duration, na.rm = TRUE)
    )

# Print the summary statistics
print(stats_summary)


# Remove useless data objects
rm(ELISA, ELISA_match, infections, LC_M, Long_Covid, my_comparisons, PASC, pre_infection_data,
   stats_summary, test_result, LongCovidIDs, p_value_pre)



# Analysis of reinfections ----



# Define the dataset to analyze - has multiple observations for two individuals
reinfection <- reinfection_raw %>%
    filter(reinfected_post_vax == 1) %>%
    select(SUBJECT.ID., covid_date, reinfected_post_vax)

reinfection$covid_date <- as.Date(reinfection$covid_date, format = "%m/%d/%Y")

rein_IDs <- reinfection %>%
    pull(SUBJECT.ID.)

rein_IDs <- unique(rein_IDs)


# Mini analysis for Emily
reinfection_fb <- left_join(reinfection, new_elisa %>%
                                dplyr::select(SUBJECT.ID. = Subject.ID, IDorHV, vaccine_1_date,
                                              CovidPos_1, CovidPos_2, CovidPos_3) %>%
                                group_by(SUBJECT.ID.) %>%
                                slice(1) %>%
                                ungroup(), by = "SUBJECT.ID.")

reinfection_fb$CovidPos_1 <- ifelse(reinfection_fb$SUBJECT.ID. == "E0118", 
                                    NA, reinfection_fb$CovidPos_1)
reinfection_fb$CovidPos_2 <- ifelse(reinfection_fb$SUBJECT.ID. == "E0118", 
                                    NA, reinfection_fb$CovidPos_2)

reinfection_fb$time_diff <- ymd(reinfection_fb$CovidPos_2) - ymd(reinfection_fb$CovidPos_1)

mean(reinfection_fb$time_diff, na.rm = T) / 30 # 6.83 months
mean(reinfection_fb[reinfection_fb$IDorHV == "HV",]$time_diff, na.rm = T) / 30 # 8.03 months
mean(reinfection_fb[reinfection_fb$IDorHV == "PID",]$time_diff, na.rm = T) / 30  # 6.70 months


# Create the ELISA dataset with new variable on reinfections
ELISA_match <- new_elisa %>%
    dplyr::select(Subject.ID, AGE, agecat, new_group, IDorHV, breakthru, CovidPos_1, conc.ug_ml, 
                  Collection.date, time, COVID_pos_before_sample) %>%
    mutate(timepoint2 = case_when(time %in% 
                                      c("Post-dose 1", "Post-dose 2", "Post-dose 3", 
                                        "Post-dose 4", "Post-dose 5") ~ time,
                                  TRUE ~ "other timepoint 2")) %>%
    mutate(reinfection = case_when(Subject.ID %in% rein_IDs ~ 1,
                                   TRUE ~ 0))


# Remove people with any infection (though include people with reinfection) for the match
ELISA_match <- ELISA_match %>%
    filter(reinfection == 1 | is.na(CovidPos_1))

# Remove bad timepoints
ELISA_match <- ELISA_match %>%
    filter(timepoint2 != "other timepoint 2")

table(ELISA_match$reinfection)


# Conduct the matching
Rei_M <- MatchIt::matchit(reinfection ~ AGE + new_group + timepoint2,
                          data = ELISA_match, 
                          replace = FALSE,
                          distance = "robust_mahalanobis",
                          exact = c("timepoint2"),
                          caliper = c(AGE = 15),
                          std.caliper = c(FALSE),
                          ratio = 2)

Rei_M <- MatchIt::match.data(Rei_M)


# Graphing the before/after
my_comparisons <- list(c("0", "1"))
Rei_M$reinfection <- factor(Rei_M$reinfection)
Rei_M$COVID_pos_before_sample <- factor(Rei_M$COVID_pos_before_sample, levels =  c(0, 1), 
                                        labels = c("Before SARS-CoV-2 infection", "After SARS-CoV-2 infection"))

reinfections <- Rei_M %>%
    ggplot(aes(x = reinfection, y = conc.ug_ml)) +
    geom_boxplot(outlier.shape = NA, 
                 show.legend = FALSE,
                 position = position_dodge(0.7, preserve = "single")) + 
    geom_jitter(aes(x = reinfection, 
                    y = conc.ug_ml, 
                    color = reinfection,
                    shape = COVID_pos_before_sample),
                width = 0.2,
                size = 3) + 
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), limits = c(NA, 40)) +
    #  stat_compare_means(comparisons = my_comparisons, 
    #                     method = "wilcox.test",
    #                     label = "p.signif") +
    facet_grid(. ~ COVID_pos_before_sample,
               switch = "x",
               scales = "free_x",
               space = "free_x") +
    # theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
    ggtitle("Individuals with Reinfections") +
    xlab("") +  # Remove the x-axis label
    ylab("IgG Levels") +
    scale_color_manual(values = c("#41b6c4", "goldenrod"), 
                       labels = c("Never had SARS-CoV-2 infection", "Individual with SARS-CoV-2 reinfection")) +
    # scale_shape_manual(values = c("circle", "triangle"), labels = c("Prior to First Covid Infection", "After First Covid Infection")) +
    theme_bw() +
    theme(
        legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_text(size = 18),  # Remove the x-axis label
        axis.text.x = element_blank(),  # Remove the x-axis tick mark labels
        axis.title.y = element_text(size = 18),  
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 20), 
        legend.text = element_text(size = 13),# Increase the plot title font size
        strip.text = element_text(size = 11.2)
    )  +
    guides(color = guide_legend(override.aes = list(size = 4)))

reinfections

ggsave(filename = "reinfections.jpg", plot = reinfections, device = "jpeg", dpi = 300, width = 8, height = 3.5, units = "in")



# Set up to compare means of the before group, define comparisons
my_comparisons <- list(c("Yes", "No"))


# Filter the data for the "pre-infection" group
pre_infection_data <- Rei_M %>%
    filter(COVID_pos_before_sample == "Before SARS-CoV-2 infection")


# Check if there are at least two groups to compare in the "pre-infection" group
if (length(unique(pre_infection_data$reinfection)) == 2) {
    test_result <- wilcox.test(conc.ug_ml ~ reinfection, data = pre_infection_data)
    p_value_pre <- test_result$p.value
    print(paste("Wilcoxon test p-value (Before) for Yes vs No:", p_value_pre))
} else {
    print("Not enough groups to compare in the 'pre-infection' group")
}


# Remove useless data objects
rm(ELISA_match, my_comparisons, pre_infection_data, Rei_M, reinfection, reinfection_raw,
   reinfections, test_result, p_value_pre, rein_IDs)
