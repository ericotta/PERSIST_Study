

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



data_orig_batch1IgG <- read_xlsx(here("C:/Users/bustosfa/Rotation 2/Data", "Copy of 2022-04-02 PERSIST serology results RML_second run (all samples).xlsx"), 
                                 sheet = 1, skip = 1)
data_orig_batch1IgA <- read_xlsx(here("C:/Users/bustosfa/Rotation 2/Data", "Copy of 2022-04-02 PERSIST serology results RML_second run (all samples).xlsx"), 
                                 sheet = 2, skip = 1)
data_orig_batch1inhib <- read_xlsx(here("C:/Users/bustosfa/Rotation 2/Data", "Combined Batch 1, 2, 3 Meso Results.xlsx"), 
                                   sheet = 1, skip = 1)
data_orig_batch2 <- read_xlsx(here("C:/Users/bustosfa/Rotation 2/Data", "Combined Batch 1, 2, 3 Meso Results.xlsx"), 
                              sheet = 2, skip = 1)
data_orig_batch3 <- read_xlsx(here("C:/Users/bustosfa/Rotation 2/Data", "Combined Batch 1, 2, 3 Meso Results.xlsx"), 
                              sheet = 3, skip = 1)
load(here("C:/Users/bustosfa/Rotation 2/Data", "data7_6Dec2022.RData")) #From 9Dec but forgot to change name in ELISA code
load(here("C:/Users/bustosfa/Rotation 2/Data", "for_fausto.RData"))
timepoints <- read.csv(here("C:/Users/bustosfa/Rotation 2/Data", "timepoints_2022-06-10.csv"))
load(here("C:/Users/bustosfa/Rotation 2/Data", "check_bt5.RData"))



# Deal with errant timepoints in Combined Batch 1 ----




# Remove grouped tibble struture
data7_orig_vFausto <- as.data.frame(data)

# Call the dataset data6, even though it's not, to not change most of the code in there
data6 <- data
rm(data)

# Add in additional variable
data6 <- left_join(data6, check_bt5 %>%
                     dplyr::select(`Subject ID`, timepoint, sample_bad), 
                   by = c("Subject ID", "timepoint"))
data6[is.na(data6$sample_bad),]$sample_bad <- 0

# Remove people in MSD datasets that are not in ELISA datasets (v7_6Dec2022)
people <- unique(data6$`Subject ID`)
rows_to_keep <- data_orig_batch1IgG$`Subject ID` %in% people
rows_to_keep2 <- data_orig_batch1inhib$`Subject ID` %in% people
all.equal(rows_to_keep, rows_to_keep2) # gr8
data_orig_batch1IgG <- data_orig_batch1IgG[rows_to_keep,]
data_orig_batch1IgA <- data_orig_batch1IgA[rows_to_keep,]
data_orig_batch1inhib <- data_orig_batch1inhib[rows_to_keep,]

data_orig_batch2 <- data_orig_batch2[data_orig_batch2$`Subject ID` %in% people,]
data_orig_batch3 <- data_orig_batch3[data_orig_batch3$`Subject ID` %in% people,]

# Find the problematic ones
errant <- data_orig_batch1inhib %>%
  dplyr::select(`Subject ID`, `Sample #`, Timepoint, `Date of Draw`) %>%
  filter(Timepoint == "6 month")

# Get the correct date
batch1 <- data6 %>%
  dplyr::select(`Subject ID`, Timepoint = timepoint, `Date of Draw` = Collection.date) %>%
  filter(`Subject ID` %in% errant$`Subject ID`)
batch1$Timepoint <- as.character(batch1$Timepoint)

# Make an indicator for future use
errant$isOK <- 0

# Identify the problematic ones, fix the correct ones
for(i in 1:nrow(errant)){
  for(j in 1:nrow(batch1)){
    if(errant$`Subject ID`[i] == batch1$`Subject ID`[j] & 
       errant$`Date of Draw`[i] == batch1$`Date of Draw`[j]){
      errant$isOK[i] <- 1
      errant$Timepoint[i] <- batch1$Timepoint[j]
    }
  }
}

errant[errant$isOK==0,]$Timepoint <- NA

# Get the codes of interest in the right order
time <- timepoints %>%
  filter(Subject.ID %in% errant[errant$isOK==0,]$`Subject ID`) %>%
  dplyr::select(Subject.ID, timepoint, Date)

errant <- errant %>%
  arrange(isOK, `Subject ID`)

time <- time %>%
  arrange(Subject.ID) 

# Enter the correct timepoints based on the timepoints dataset
errant[errant$`Subject ID` == "E0003",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0034",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0042",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0046",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0049",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0055",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0068",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0071",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0072",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0081",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0083",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0093",]$Timepoint <- "Post-dose 3"
errant[errant$`Subject ID` == "E0106",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0125",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0132",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0159",]$Timepoint <- "Pre-dose 3"
errant[errant$`Subject ID` == "E0160",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0177",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0184",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0186",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0190",]$Timepoint <- "Pre-dose 3"
errant[errant$`Subject ID` == "E0191",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "E0193",]$Timepoint <- "Pre-dose 3"
errant[errant$`Subject ID` == "E0208",]$Timepoint <- "Pre-dose 3"
errant[errant$`Subject ID` == "E0209",]$Timepoint <- "Pre-dose 3"
errant[errant$`Subject ID` == "H0003",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "H0005",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "H0006",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "H0008",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "H0011",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "H0012",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "H0015",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "H0019",]$Timepoint <- "6 months post-dose 2"
errant[errant$`Subject ID` == "H0021",]$Timepoint <- "6 months post-dose 2"

# Replace the incorrect timepoints with the correct timepoints
check <- data_orig_batch1inhib

sum(check$Timepoint == "6 month") #41

for(i in 1:nrow(check)){
  for(j in 1:nrow(errant)){
    if(check$`Sample #`[i] == errant$`Sample #`[j]){
      check$Timepoint[i] <- errant$Timepoint[j]
    }
  }
}

sum(check$Timepoint == "6 month") #0, great!

# Check that I didn't introduce an error somewhere
check$new <- paste0(check$`Sample #`, "_", check$Timepoint)
sum(duplicated(check$new)) # 0, great!

check$new2 <- paste0(check$`Subject ID`, "_", check$Timepoint)
sum(duplicated(check$new2)) # 0, great!

rm(check)

# Now fix the underlying datasets

# First fix the inhibition dataset
for(i in 1:nrow(data_orig_batch1inhib)){
  for(j in 1:nrow(errant)){
    if(data_orig_batch1inhib$`Sample #`[i] == errant$`Sample #`[j]){
      data_orig_batch1inhib$Timepoint[i] <- errant$Timepoint[j]
    }
  }
}

data_orig_batch1inhib[data_orig_batch1inhib$Timepoint == "6 month post-dose 4",]$Timepoint <- 
  "6 months post-dose 4"

# Then fix the IgG dataset
for(i in 1:nrow(data_orig_batch1IgG)){
  for(j in 1:nrow(errant)){
    if(data_orig_batch1IgG$`Sample No`[i] == errant$`Sample #`[j]){
      data_orig_batch1IgG$Timepoint[i] <- errant$Timepoint[j]
    }
  }
}

data_orig_batch1IgG[data_orig_batch1IgG$Timepoint == "6 month post-dose 4",]$Timepoint <- 
  "6 months post-dose 4"

# Remove useless data objects before continuing to data processing
rm(batch1, errant, time, timepoints, i, j)




# Basic data cleaning of batch 1 ----




# Do basic cleaning of batch 1 IgG data
data1 <- data_orig_batch1IgG
data1 <- data1 %>%
  dplyr::select(`Subject ID`, Timepoint, "Wuhan" = `Wuhan Spike plate 23`,
                Alpha, Beta, Gamma, Delta, `Delta AY.4`, `Delta AY.4.2`, `Omicron BA.1`) %>%
  arrange(`Subject ID`)
data1 <- as.data.frame(data1)

# Deal with variable types and missing data
cols.num <- names(data1)[names(data1) %!in% c("Subject ID", "Timepoint")]
data1[cols.num] <- sapply(data1[cols.num], as.numeric)
sapply(data1, class)

data1 <- data1 %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) 

# Add in group from the ELISA data
group <- data6 %>%
  dplyr::select(`Subject ID`, group = new_group) %>%
  distinct()

data1 <- left_join(data1, group, by = "Subject ID")

data1[is.na(data1$group),]$`Subject ID` # check if E0141 is in data6
data1[is.na(data1$group),]$`Subject ID` %in% data6$`Subject ID` # no, then remove from MSD

data1 <- data1 %>%
  drop_na(group) # 157 -> 157

# Deal with the time variable
data1$time <- factor(data1$Timepoint)
data1$Timepoint <- NULL
data1 <- data1 %>%
  relocate(time, .after = `Subject ID`) %>%
  relocate(group, .after = time)
names(data1)[1:3] <- c("id", "time", "group")

data1 <- data1 %>%
  arrange(group, id, time)
data1 <- data1 %>%
  relocate(group, .before = id)

# Rename the dataset
data1_IgG <- data1
rm(data1)





# Now do basic cleaning of batch 1 IgA data
data1_IgA <- data_orig_batch1IgA
data1_IgA$`Subject ID` <- data_orig_batch1IgG$`Subject ID`
data1_IgA$Timepoint <- data_orig_batch1IgG$Timepoint
data1_IgA <- data1_IgA %>%
  dplyr::select(`Subject ID`, Timepoint, "Wuhan" = `Wuhan Spike plate 23`,
                Alpha, Beta, Gamma, Delta, `Delta AY.4`, `Delta AY.4.2`, `Omicron BA.1`) %>%
  arrange(`Subject ID`)
data1_IgA <- as.data.frame(data1_IgA)

# Deal with variable types and missing data
cols.num <- names(data1_IgA)[names(data1_IgA) %!in% c("Subject ID", "Timepoint")]
data1_IgA[cols.num] <- sapply(data1_IgA[cols.num], as.numeric)
sapply(data1_IgA, class)

data1_IgA <- data1_IgA %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) 

# Add in group from the ELISA data
data1_IgA <- left_join(data1_IgA, group, by = "Subject ID")

data1_IgA[is.na(data1_IgA$group),]$`Subject ID` # check if E0141 is in data6
data1_IgA[is.na(data1_IgA$group),]$`Subject ID` %in% data6$`Subject ID` # no, then remove from MSD

data1_IgA <- data1_IgA %>%
  drop_na(group) # 157 -> 157

# Deal with the time variable
data1_IgA$time <- factor(data1_IgA$Timepoint)
data1_IgA$Timepoint <- NULL
data1_IgA <- data1_IgA %>%
  relocate(time, .after = `Subject ID`) %>%
  relocate(group, .after = time)
names(data1_IgA)[1:3] <- c("id", "time", "group")

data1_IgA <- data1_IgA %>%
  arrange(group, id, time)
data1_IgA <- data1_IgA %>%
  relocate(group, .before = id)








# Now do basic cleaning of batch 1 inhibition data
data1_inhib <- as.data.frame(data_orig_batch1inhib)
data1_inhib <- data1_inhib %>%
  dplyr::select(`Subject ID`, Timepoint, `Wuhan` = `CoV-2 Spike`, 
                Alpha = B.1.1.7, Beta = B.1.351, Gamma = P.1, Delta = B.1.617.2,
                `Delta AY.4` = AY.4, `Delta AY.4.2` = AY.4.2, `Omicron BA.1` = B.1.1.529) %>%
  arrange(`Subject ID`) %>%
  drop_na(Beta) # Remove 3 with no inhib data # 157 -> 154

# Deal with those below the detection limit
data1_inhib[data1_inhib <= 5] <- 5 / sqrt(2)

# Deal with variable types and missing data
cols.num <- names(data1_inhib)[names(data1_inhib) %!in% c("Subject ID", "Timepoint")]
data1_inhib[cols.num] <- sapply(data1_inhib[cols.num], as.numeric)
sapply(data1_inhib, class)

data1_inhib <- data1_inhib %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) 

# Add in group from the ELISA data
data1_inhib <- left_join(data1_inhib, group, by = "Subject ID")

data1_inhib[is.na(data1_inhib$group),]$`Subject ID` # check if E0141 is in data6
data1_inhib[is.na(data1_inhib$group),]$`Subject ID` %in% data6$`Subject ID` # no, then remove from MSD

data1_inhib <- data1_inhib %>%
  drop_na(group) # 154 -> 154

# Deal with the time variable
data1_inhib$time <- factor(data1_inhib$Timepoint)
data1_inhib$Timepoint <- NULL
table2(data1_inhib$time)

# Move things around
data1_inhib <- data1_inhib %>%
  relocate(time, .after = `Subject ID`) %>%
  relocate(group, .after = time)
names(data1_inhib)[1:3] <- c("id", "time", "group")

data1_inhib <- data1_inhib %>%
  arrange(group, id, time)
data1_inhib <- data1_inhib %>%
  relocate(group, .before = id)




# Basic data cleaning of batch 2 ----




# Cleaning of batch 2 IgG data
data2_IgG <- as.data.frame(data_orig_batch2)
data2_IgG <- data2_IgG %>%
  dplyr::select(`Subject ID`, Timepoint, `Wuhan` = `CoV-2 Spike...6`, 
                Alpha = B.1.1.7...11, Beta = B.1.351...12, Gamma = P.1...10, 
                Delta = B.1.617.2...13, `Delta AY.4` = AY.4...9, 
                `Delta AY.4.2` = AY.4.2...8, `Omicron BA.1` = B.1.1.529...7) %>%
  arrange(`Subject ID`)

# Deal with variable types and missing data
cols.num <- names(data2_IgG)[names(data2_IgG) %!in% c("Subject ID", "Timepoint")]
data2_IgG[cols.num] <- sapply(data2_IgG[cols.num], as.numeric)
sapply(data2_IgG, class)

data2_IgG <- data2_IgG %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) 

# Add in group from the ELISA data
data2_IgG <- left_join(data2_IgG, group, by = "Subject ID")

data2_IgG[is.na(data2_IgG$group),]$`Subject ID` # check if H0026 is in data6
data2_IgG[is.na(data2_IgG$group),]$`Subject ID` %in% data6$`Subject ID` # no, then remove from MSD

data2_IgG <- data2_IgG %>%
  drop_na(group) # 145 -> 145

# Deal with the time variable
data2_IgG$time <- factor(data2_IgG$Timepoint)
data2_IgG$Timepoint <- NULL

# Move things around
data2_IgG <- data2_IgG %>%
  relocate(time, .after = `Subject ID`) %>%
  relocate(group, .after = time)
names(data2_IgG)[1:3] <- c("id", "time", "group")

data2_IgG <- data2_IgG %>%
  arrange(group, id, time)
data2_IgG <- data2_IgG %>%
  relocate(group, .before = id)








# Cleaning of batch 2 IgA data
data2_IgA <- as.data.frame(data_orig_batch2)
data2_IgA <- data2_IgA %>%
  dplyr::select(`Subject ID`, Timepoint, `Wuhan` = `CoV-2 Spike...14`, 
                Alpha = B.1.1.7...19, Beta = B.1.351...20, Gamma = P.1...18, 
                Delta = B.1.617.2...21, `Delta AY.4` = AY.4...17, 
                `Delta AY.4.2` = AY.4.2...16, `Omicron BA.1` = B.1.1.529...15) %>%
  arrange(`Subject ID`)

# Deal with variable types and missing data
cols.num <- names(data2_IgA)[names(data2_IgA) %!in% c("Subject ID", "Timepoint")]
data2_IgA[cols.num] <- sapply(data2_IgA[cols.num], as.numeric)
sapply(data2_IgA, class)

data2_IgA <- data2_IgA %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) 

# Add in group from the ELISA data
data2_IgA <- left_join(data2_IgA, group, by = "Subject ID")

data2_IgA[is.na(data2_IgA$group),]$`Subject ID` # check if H0026 is in data6
data2_IgA[is.na(data2_IgA$group),]$`Subject ID` %in% data6$`Subject ID` # no, then remove from MSD

data2_IgA <- data2_IgA %>%
  drop_na(group) # 145 -> 145

# Deal with the time variable
data2_IgA$time <- factor(data2_IgA$Timepoint)
data2_IgA$Timepoint <- NULL

# Move things around
data2_IgA <- data2_IgA %>%
  relocate(time, .after = `Subject ID`) %>%
  relocate(group, .after = time)
names(data2_IgA)[1:3] <- c("id", "time", "group")

data2_IgA <- data2_IgA %>%
  arrange(group, id, time)
data2_IgA <- data2_IgA %>%
  relocate(group, .before = id)






# Now do basic cleaning of batch 2 inhibition data
data2_inhib <- as.data.frame(data_orig_batch2)
data2_inhib <- data2_inhib %>%
  dplyr::select(`Subject ID`, Timepoint, `Wuhan` = `CoV-2 Spike...22`, 
                Alpha = B.1.1.7...27, Beta = B.1.351...28, Gamma = P.1...26, 
                Delta = B.1.617.2...29, `Delta AY.4` = AY.4...25, 
                `Delta AY.4.2` = AY.4.2...24, `Omicron BA.1` = B.1.1.529...23) %>%
  arrange(`Subject ID`) %>%
  drop_na(Beta) # 145 -> 102 Remove those without inhib data

# Deal with those below the detection limit
data2_inhib[data2_inhib <= 5] <- 5 / sqrt(2)

# Deal with variable types and missing data
cols.num <- names(data2_inhib)[names(data2_inhib) %!in% c("Subject ID", "Timepoint")]
data2_inhib[cols.num] <- sapply(data2_inhib[cols.num], as.numeric)
sapply(data2_inhib, class)

data2_inhib <- data2_inhib %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) 

# Add in group from the ELISA data
data2_inhib <- left_join(data2_inhib, group, by = "Subject ID")

data2_inhib[is.na(data2_inhib$group),]$`Subject ID` 
data2_inhib[is.na(data2_inhib$group),]$`Subject ID` %in% data6$`Subject ID` 

data2_inhib <- data2_inhib %>%
  drop_na(group) # 102 -> 102

# Deal with the time variable
data2_inhib$time <- factor(data2_inhib$Timepoint)
data2_inhib$Timepoint <- NULL
table2(data2_inhib$time)

# Move things around
data2_inhib <- data2_inhib %>%
  relocate(time, .after = `Subject ID`) %>%
  relocate(group, .after = time)
names(data2_inhib)[1:3] <- c("id", "time", "group")

data2_inhib <- data2_inhib %>%
  arrange(group, id, time)
data2_inhib <- data2_inhib %>%
  relocate(group, .before = id)







# Basic data cleaning of batch 3 ----


# Cleaning of batch 3 IgG data
data3_IgG <- as.data.frame(data_orig_batch3)
data3_IgG <- data3_IgG %>%
  dplyr::select(`Subject ID`, Timepoint, `Wuhan` = `CoV-2 Spike...6`, 
                Alpha = B.1.1.7...11, Beta = B.1.351...12, Gamma = P.1...10, 
                Delta = B.1.617.2...13, `Delta AY.4` = AY.4...9, 
                `Delta AY.4.2` = AY.4.2...8, `Omicron BA.1` = B.1.1.529...7) %>%
  arrange(`Subject ID`)

# Deal with variable types and missing data
cols.num <- names(data3_IgG)[names(data3_IgG) %!in% c("Subject ID", "Timepoint")]
data3_IgG[cols.num] <- sapply(data3_IgG[cols.num], as.numeric)
sapply(data3_IgG, class)

data3_IgG <- data3_IgG %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) 

# Add in group from the ELISA data
data3_IgG <- left_join(data3_IgG, group, by = "Subject ID")

data3_IgG[is.na(data3_IgG$group),]$`Subject ID` # check if the following are in data6
# "E0122" "E0141" "E0223" "E0224" "E0225" "H0048"
data3_IgG[is.na(data3_IgG$group),]$`Subject ID` %in% data6$`Subject ID` # no, then remove from MSD

data3_IgG <- data3_IgG %>%
  drop_na(group) # 136 -> 136

# Deal with the time variable
data3_IgG$time <- factor(data3_IgG$Timepoint)
data3_IgG$Timepoint <- NULL

# Move things around
data3_IgG <- data3_IgG %>%
  relocate(time, .after = `Subject ID`) %>%
  relocate(group, .after = time)
names(data3_IgG)[1:3] <- c("id", "time", "group")

data3_IgG <- data3_IgG %>%
  arrange(group, id, time)
data3_IgG <- data3_IgG %>%
  relocate(group, .before = id)








# Cleaning of batch 3 IgA data
data3_IgA <- as.data.frame(data_orig_batch3)
data3_IgA <- data3_IgA %>%
  dplyr::select(`Subject ID`, Timepoint, `Wuhan` = `CoV-2 Spike...14`, 
                Alpha = B.1.1.7...19, Beta = B.1.351...20, Gamma = P.1...18, 
                Delta = B.1.617.2...21, `Delta AY.4` = AY.4...17, 
                `Delta AY.4.2` = AY.4.2...16, `Omicron BA.1` = B.1.1.529...15) %>%
  arrange(`Subject ID`)

# Deal with variable types and missing data
cols.num <- names(data3_IgA)[names(data3_IgA) %!in% c("Subject ID", "Timepoint")]
data3_IgA[cols.num] <- sapply(data3_IgA[cols.num], as.numeric)
sapply(data3_IgA, class)

data3_IgA <- data3_IgA %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) 

# Add in group from the ELISA data
data3_IgA <- left_join(data3_IgA, group, by = "Subject ID")

data3_IgA[is.na(data3_IgA$group),]$`Subject ID` # check if the following are in data6
#  "E0122" "E0141" "E0223" "E0224" "E0225" "H0048"
data3_IgA[is.na(data3_IgA$group),]$`Subject ID` %in% data6$`Subject ID` # no, then remove from MSD

data3_IgA <- data3_IgA %>%
  drop_na(group) # 136 -> 136

# Deal with the time variable
data3_IgA$time <- factor(data3_IgA$Timepoint)
data3_IgA$Timepoint <- NULL

# Move things around
data3_IgA <- data3_IgA %>%
  relocate(time, .after = `Subject ID`) %>%
  relocate(group, .after = time)
names(data3_IgA)[1:3] <- c("id", "time", "group")

data3_IgA <- data3_IgA %>%
  arrange(group, id, time)
data3_IgA <- data3_IgA %>%
  relocate(group, .before = id)









# Now do basic cleaning of batch 3 inhibition data
data3_inhib <- as.data.frame(data_orig_batch3)
data3_inhib <- data3_inhib %>%
  dplyr::select(`Subject ID`, Timepoint, `Wuhan` = `CoV-2 Spike...22`, 
                Alpha = B.1.1.7...27, Beta = B.1.351...28, Gamma = P.1...26, 
                Delta = B.1.617.2...29, `Delta AY.4` = AY.4...25, 
                `Delta AY.4.2` = AY.4.2...24, `Omicron BA.1` = B.1.1.529...23) %>%
  arrange(`Subject ID`) %>%
  drop_na(Beta) # 136 -> 119 drop people without inhib values

# Deal with those below the detection limit
data3_inhib[data3_inhib <= 5] <- 5 / sqrt(2)

# Deal with variable types and missing data
cols.num <- names(data3_inhib)[names(data3_inhib) %!in% c("Subject ID", "Timepoint")]
data3_inhib[cols.num] <- sapply(data3_inhib[cols.num], as.numeric)
sapply(data3_inhib, class)

data3_inhib <- data3_inhib %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) 

# Add in group from the ELISA data
data3_inhib <- left_join(data3_inhib, group, by = "Subject ID")

data3_inhib[is.na(data3_inhib$group),]$`Subject ID` # check if the following are in data6
# E0122, E0141, E0223, E0224, E0225, H0048
data3_inhib[is.na(data3_inhib$group),]$`Subject ID` %in% data6$`Subject ID` # no, then remove from MSD

data3_inhib <- data3_inhib %>%
  drop_na(group) # 119 -> 119

# Deal with the time variable
data3_inhib$time <- factor(data3_inhib$Timepoint)
data3_inhib$Timepoint <- NULL
table2(data3_inhib$time)

# Move things around
data3_inhib <- data3_inhib %>%
  relocate(time, .after = `Subject ID`) %>%
  relocate(group, .after = time)
names(data3_inhib)[1:3] <- c("id", "time", "group")

data3_inhib <- data3_inhib %>%
  arrange(group, id, time)
data3_inhib <- data3_inhib %>%
  relocate(group, .before = id)








# Fix the immunological group ----
data6$group <- NULL
unique(data6$new_group)

small <- data6 %>%
  dplyr::select(id = `Subject ID`, new_group) %>%
  distinct()
unique(small$group)

data1_IgG <- left_join(data1_IgG, small, by = "id") %>%
  dplyr::select(new_group, id, time, Wuhan, Alpha, Beta, Gamma, Delta,
                `Delta AY.4`, `Delta AY.4.2`, `Omicron BA.1`)

data2_IgG <- left_join(data2_IgG, small, by = "id") %>%
  dplyr::select(new_group, id, time, Wuhan, Alpha, Beta, Gamma, Delta,
                `Delta AY.4`, `Delta AY.4.2`, `Omicron BA.1`)

data3_IgG <- left_join(data3_IgG, small, by = "id") %>%
  dplyr::select(new_group, id, time, Wuhan, Alpha, Beta, Gamma, Delta,
                `Delta AY.4`, `Delta AY.4.2`, `Omicron BA.1`)



data1_inhib <- left_join(data1_inhib, small, by = "id") %>%
  dplyr::select(new_group, id, time, Wuhan, Alpha, Beta, Gamma, Delta,
                `Delta AY.4`, `Delta AY.4.2`, `Omicron BA.1`)

data2_inhib <- left_join(data2_inhib, small, by = "id") %>%
  dplyr::select(new_group, id, time, Wuhan, Alpha, Beta, Gamma, Delta,
                `Delta AY.4`, `Delta AY.4.2`, `Omicron BA.1`)

data3_inhib <- left_join(data3_inhib, small, by = "id") %>%
  dplyr::select(new_group, id, time, Wuhan, Alpha, Beta, Gamma, Delta,
                `Delta AY.4`, `Delta AY.4.2`, `Omicron BA.1`)




# Put the IgG data together ----
all_IgG <- rbind(data1_IgG, data2_IgG, data3_IgG)

# Combine 6 months post-dose 2 and pre-dose 3
all_IgG$timepoint <- as.character(all_IgG$time)
check <- all_IgG %>%
  dplyr::select(id, timepoint) %>%
  filter(timepoint %in% c("6 months post-dose 2", "Pre-dose 3")) %>%
  distinct() %>%
  arrange(id) %>%
  group_by(id) %>% 
  filter(n() > 1) 
check <- as.vector(unique(check$id)) #ID with both timepoints: E0111

all_IgG$timepoint2 <- all_IgG$timepoint
for(i in 1:nrow(all_IgG)){
  if(all_IgG$id[i] %!in% check & all_IgG$timepoint[i] == "Pre-dose 3"){
    all_IgG$timepoint2[i] <- "6 months post-dose 2"
  }
  if(all_IgG$id[i] %in% check & all_IgG$timepoint[i] == "Pre-dose 3"){
    all_IgG$timepoint2[i] <- "Remove"
  }
}
check <- all_IgG %>%
  dplyr::select(id, time, timepoint, timepoint2)

all_IgG <- all_IgG[all_IgG$timepoint2 != "Remove",] #438 -> 437

all_IgG$timepoint <- NULL
all_IgG$time <- all_IgG$timepoint2
all_IgG$timepoint2 <- NULL

all_IgG$time <- factor(all_IgG$time) # To get rid of factor levels with 0 entries

# Restrict to the time points of interest
all_IgG <- all_IgG %>%
  filter(time %in% c("Post-dose 2","6 months post-dose 2",
                     "Post-dose 3", "6 months post-dose 3"))

# Add in the pos_before variable from data6
data6 %>%
  filter(analyte == "Spike IgG" &
           timepoint %in% c("Post-dose 2","6 months post-dose 2",
                            "Post-dose 3", "6 months post-dose 3")) %>%
  dplyr::select(id = `Subject ID`, time = timepoint, evusheld_before_sample,
                COVID_pos_before_sample, N_pos_before_sample,
                sample_bad) -> restriction
restriction$pos_before <- ifelse(restriction$evusheld_before_sample == 1 | 
                                   restriction$COVID_pos_before_sample == 1 | 
                                   restriction$N_pos_before_sample == 1, 1, 0)
round(prop.table(table(restriction$pos_before)), 2)

restriction <- restriction %>%
  dplyr::select(id, time, pos_before, sample_bad)

all_IgG <- left_join(all_IgG, restriction, by = c("id" = "id", "time" = "time"))
all_IgG[is.na(all_IgG$sample_bad),]$sample_bad <- 0

table2(all_IgG$sample_bad)

# Identify any samples in MESO data that are not in data6
table(all_IgG$pos_before, useNA = "ifany") # 18 NA

for_emily <- all_IgG %>%
  filter(is.na(pos_before)) %>%
  dplyr::select(id, time)
#save(for_emily, file = "for_emily.RData")

# Add in additional info on previously positive individuals that aren't in data6 
# Decision: Only remove those samples that show evidence of previous infection
# Nothing else suggests the others are previously positive
for_fausto2 <- for_fausto %>%
  filter(evusheld_before_sample == 1 | Covid_pos_before_sample == 1)

all_IgG <- arrange(all_IgG, id, time) # Arrange dataset by time to make for loop work
all_IgG <- distinct(all_IgG, id, time, .keep_all = TRUE) # Remove dupes, 328 -> 271

for(i in 1:(nrow(all_IgG) - 1)){
  if(all_IgG$id[i] == all_IgG$id[i + 1] &
     all_IgG$sample_bad[i] == 1){
    all_IgG$sample_bad[i + 1] <- 1
  }
}


table(all_IgG$pos_before, useNA = "ifany") # 17 NA

check <- all_IgG %>%
  filter(id == "E0159")

for(i in 1:nrow(all_IgG)){ # Add in Emily's data points
  for(j in 1:nrow(for_fausto2)){
    if(all_IgG$id[i] == for_fausto2$id[j] & 
       all_IgG$time[i] == for_fausto2$time[j]){
      all_IgG$pos_before[i] <- 1
    }
  }
}

table(all_IgG$pos_before, useNA = "ifany") # recheck missing, 15 missing

for(i in 1:(nrow(all_IgG) - 1)){ # Deploy the for loop above
  if(all_IgG$id[i] == all_IgG$id[i + 1] & 
     !is.na(all_IgG$pos_before[i]) & all_IgG$pos_before[i] == 1){
    all_IgG$pos_before[i + 1] <- 1
  }
}

table(all_IgG$pos_before, useNA = "ifany") # recheck missing, 11 missing

all_IgG[is.na(all_IgG$pos_before),]$pos_before <- 0

table(all_IgG$pos_before, useNA = "ifany") # recheck missing, 0 missing
# 107 pos_before, 164 not pos_before


# Additional data cleaning
all_IgG2 <- all_IgG %>%
  pivot_longer(!c(new_group, id, time, pos_before, sample_bad), 
               names_to = "antigen", values_to = "binding")

all_IgG2$antigen <- factor(all_IgG2$antigen,
                           levels = c("Wuhan", "Alpha", 
                                      "Beta", "Gamma", "Delta",
                                      "Delta AY.4", "Delta AY.4.2", "Omicron BA.1"))

all_IgG2$new_group <- factor(all_IgG2$new_group,
                             levels = c("Healthy volunteer", "Antibody deficiency", 
                                        "PIRD", "Combined immunodeficiency",
                                        "Other IEI", "Other immune disorder"))

all_IgG2$time <- factor(all_IgG2$time,
                        levels = c("Post-dose 2","6 months post-dose 2",
                                   "Post-dose 3", "6 months post-dose 3"))

all_IgG2 %>%
  filter(antigen %!in% c("Delta AY.4", "Delta AY.4.2")) -> all_IgG2

# Deal with conversion factor
all_IgG2$binding <- all_IgG2$binding * 0.00901

# Deal with limit of detection issues for IgG 
all_IgG2[all_IgG2$binding <= 21.05 * 0.00901 & 
           !is.na(all_IgG2$binding),]$binding <- 21.05 * 0.00901 / sqrt(2)






# Put the inhibition data together ----
all_inhib <- rbind(data1_inhib, data2_inhib, data3_inhib)

# Combine 6 months post-dose 2 and pre-dose 3
all_inhib$timepoint <- as.character(all_inhib$time)
check <- all_inhib %>%
  dplyr::select(id, timepoint) %>%
  filter(timepoint %in% c("6 months post-dose 2", "Pre-dose 3")) %>%
  distinct() %>%
  arrange(id) %>%
  group_by(id) %>% 
  filter(n() > 1) 
check <- as.vector(unique(check$id)) #ID with both timepoints: E0111

all_inhib$timepoint2 <- all_inhib$timepoint
for(i in 1:nrow(all_inhib)){
  if(all_inhib$id[i] %!in% check & all_inhib$timepoint[i] == "Pre-dose 3"){
    all_inhib$timepoint2[i] <- "6 months post-dose 2"
  }
  if(all_inhib$id[i] %in% check & all_inhib$timepoint[i] == "Pre-dose 3"){
    all_inhib$timepoint2[i] <- "Remove"
  }
}
check <- all_inhib %>%
  dplyr::select(id, time, timepoint, timepoint2)

all_inhib <- all_inhib[all_inhib$timepoint2 != "Remove",] #375 -> 374

all_inhib$timepoint <- NULL
all_inhib$time <- all_inhib$timepoint2
all_inhib$timepoint2 <- NULL

all_inhib$time <- factor(all_inhib$time) # To get rid of factor levels with 0 entries

# Restrict to the timepoints of interest
all_inhib <- all_inhib %>%
  filter(time %in% c("Post-dose 2","6 months post-dose 2",
                     "Post-dose 3", "6 months post-dose 3"))

# Add in the pos_before variable
all_inhib <- left_join(all_inhib, restriction, by = c("id" = "id", "time" = "time"))


all_inhib[is.na(all_inhib$sample_bad),]$sample_bad <- 0

all_inhib <- arrange(all_inhib, id, time) # Arrange dataset by time to make for loop work
all_inhib <- distinct(all_inhib, id, time, .keep_all = TRUE) # Remove dupes, 328 -> 271

for(i in 1:(nrow(all_inhib) - 1)){
  if(all_inhib$id[i] == all_inhib$id[i + 1] &
     all_inhib$sample_bad[i] == 1){
    all_inhib$sample_bad[i + 1] <- 1
  }
}

table2(all_inhib$sample_bad)

table(all_inhib$pos_before, useNA = "ifany") #17 missing

for(i in 1:nrow(all_inhib)){ # Add in Emily's data points
  for(j in 1:nrow(for_fausto2)){
    if(all_inhib$id[i] == for_fausto2$id[j] & 
       all_inhib$time[i] == for_fausto2$time[j]){
      all_inhib$pos_before[i] <- 1
    }
  }
}

table(all_inhib$pos_before, useNA = "ifany") # recheck missing, 15 missing

for(i in 1:(nrow(all_inhib) - 1)){ # Deploy the for loop above
  if(all_inhib$id[i] == all_inhib$id[i + 1] & 
     !is.na(all_inhib$pos_before[i]) & all_inhib$pos_before[i] == 1){
    all_inhib$pos_before[i + 1] <- 1
  }
}

table(all_inhib$pos_before, useNA = "ifany") # recheck missing, 11 missing

all_inhib[is.na(all_inhib$pos_before),]$pos_before <- 0

table(all_inhib$pos_before, useNA = "ifany") # recheck missing, 0 missing
# 107 pos_before, 164 not pos_before

# Additional data cleaning
all_inhib2 <- all_inhib %>%
  pivot_longer(!c(new_group, id, time, pos_before, sample_bad),
               names_to = "antigen", values_to = "binding")

all_inhib2$antigen <- factor(all_inhib2$antigen,
                             levels = c("Wuhan", "Alpha", 
                                        "Beta", "Gamma", "Delta",
                                        "Delta AY.4", "Delta AY.4.2", "Omicron BA.1"))

all_inhib2$new_group <- factor(all_inhib2$new_group,
                               levels = c("Healthy volunteer", "Antibody deficiency", 
                                          "PIRD", "Combined immunodeficiency",
                                          "Other IEI", "Other immune disorder"))

all_inhib2$time <- factor(all_inhib2$time,
                          levels = c("Post-dose 2","6 months post-dose 2",
                                     "Post-dose 3", "6 months post-dose 3"))

all_inhib2 %>%
  filter(antigen %!in% c("Delta AY.4", "Delta AY.4.2")) -> all_inhib2

# Remove unnecessary data objects
rm(check_bt5, data_orig_batch1IgA, data1_IgA, data1_IgG, data1_inhib,
   data2_IgA, data2_IgG, data2_inhib, data3_IgA, data3_IgG, data3_inhib, 
   data6, for_emily, for_fausto, for_fausto2, group, restriction, small,
   i, j, people, rows_to_keep, rows_to_keep2, cols.num, check)

# Set LOD variables
LOD_IgG <- 21.05 * 0.00901
LOD_inhib <- 5

# For both IgG and inhib, change Wuhan to Ancestral
names(all_IgG)[names(all_IgG) == "Wuhan"] <- "Ancestral"

names(all_inhib)[names(all_inhib) == "Wuhan"] <- "Ancestral"

all_IgG2$antigen <- as.character(all_IgG2$antigen)
all_IgG2$antigen[all_IgG2$antigen == "Wuhan"] <- "Ancestral"
all_IgG2$antigen <- factor(all_IgG2$antigen, 
                           levels = c("Ancestral", "Alpha", "Beta",
                                      "Gamma", "Delta", "Omicron BA.1"))

all_inhib2$antigen <- as.character(all_inhib2$antigen)
all_inhib2$antigen[all_inhib2$antigen == "Wuhan"] <- "Ancestral"
all_inhib2$antigen <- factor(all_inhib2$antigen, 
                             levels = c("Ancestral", "Alpha", "Beta",
                                        "Gamma", "Delta", "Omicron BA.1"))



#### Figure S5 - IgG concentration over time by subgroup ----



# Make the labels
fig_s5_labels <- c(expression("10"^-1), expression("10"^1), 
                   expression("10"^3), expression("10"^5))

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
               fun = "median", size = 0.5, geom = "crossbar", show.legend = FALSE) + 
  
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



#### Compare IgG concentration and ACE2 neutralization at pd2 and 6mpd2 ---- 



# Decrease from pd2 to 6mpds: IgG concentration
all_IgG2 %>%
  group_by(time) %>%
  summarise(titers = median(binding, na.rm = T)) %>% 
  filter(time %in% c("Post-dose 2", "6 months post-dose 2")) %>%
  pivot_wider(names_from = time, 
              values_from = titers) %>%
  mutate(decrease = paste0(round(100 * (1 - `6 months post-dose 2`/`Post-dose 2`), 1), "%"),
         fold = `Post-dose 2` / `6 months post-dose 2`)

# Decrease from pd2 to 6mpds: ACE2 neutralization
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
corr_df <- corr_df[corr_df$binding < 4e+06 * 0.00901,] # 981, to drop biologically implausible outliers
corr_df <- na.omit(corr_df) # 975, can't do correlation with missing values



#### Overall correlation ----



# Note: This section is commented out because the code takes a long time to run. Some results and given in the comments. Uncomment the code to run the analysis.



# # Assuming independence at the row level
# ci_cor(corr_df$binding, corr_df$inhib, seed = 1989,
#        method = "kendall", type = "bootstrap", boot_type = "bca") # 0.53, (0.49, 0.56)
# 
# 	Two-sided 95% bootstrap confidence interval for the true Kendall correlation
# 	coefficient based on 9999 bootstrap replications and the bca method
# 
# Sample estimate: 0.5360391
# Confidence interval:
#      2.5%      97.5%
# 0.4920645    0.5585753
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


# # Now repeat for the Wuhan variant
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



# Binding efficiency and ACE2 neutralization by time
dose.labs <- c("Post-dose 2", "6 months post-dose 2",
               "Post-dose 3", "6 months post-dose 3")
names(dose.labs) <- c("Post-dose 2", "6 months post-dose 2", 
                      "Post-dose 3", "6 months post-dose 3")

# Binding efficiency and ACE2 neutralization by time, with random effects
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
fisher.test(df_time) # p-value = 7.431e-10

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
# need a reasonably regular grid in order to get sensibal contour plots
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



#### Fig S16 ----



# Comparing ELISA IgG to MSD
powers_of_ten_labels3 <- c(expression("10"^-3), expression("10"^-2), 
                           expression("10"^-1), expression("10"^0),
                           expression("10"^1), expression("10"^2), 
                           expression("10"^3), expression("10"^4), 
                           expression("10"^5), expression("10"^6))

elisa <- data7_orig_vFausto
elisa <- elisa %>%
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
