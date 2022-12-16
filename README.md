# PERSIST_Study

# Title: PERSIST study - Creation and application of 4 paramater logistic model

# Persistent Identifier: [A unique persistent identified (PID) such as a digital object identifier (DOI) or accession number supports data discovery, reporting and assessment.]

# Author(s): Emily Ricotta

# Grant Number: [AI001335]

This is code to create the functions and then calculate ELISA antibody concentrations from optical density using a 4 parameter logistic model. 

ELISAs were conducted using the methods described in Klumpp-Thomas, C., Kalish, H., Drew, M. et al. Standardization of ELISA protocols for serosurveys of the SARS-CoV-2 pandemic using clinical and at-home blood sampling. Nat Commun 12, 113 (2021). https://doi.org/10.1038/s41467-020-20383-x. Concentrations were calculated using the equations in supplemental figure 10 in the supplementary material. 

This code was created using R version 4.2.1 (2022-06-23). Packages required include the following and all dependencies:

tidyverse_1.3.2
 Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, 
  Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller
  K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C,
  Woo K, Yutani H (2019). “Welcome to the tidyverse.” _Journal of Open Source
  Software_, *4*(43), 1686. doi:10.21105/joss.01686
  <https://doi.org/10.21105/joss.01686>
  
magrittr_2.0.3
  Bache S, Wickham H (2022). _magrittr: A Forward-Pipe Operator for R_. R
    package version 2.0.3, <https://CRAN.R-project.org/package=magrittr>
    
drc_3.0-1
  Ritz, C., Baty, F., Streibig, J. C., Gerhard, D. (2015) Dose-Response
    Analysis Using R PLOS ONE, 10(12), e0146021
 
dr4pl_2.0.0
  Landis J, An H, Bailey A, Dittmer D, Marron J (2021). _dr4pl: Dose Response
    Data Analysis using the 4 Parameter Logistic (4pl) Model_. R package version
    2.0.0, <https://CRAN.R-project.org/package=dr4pl>
    
    
Simulated data has been provided as an example.

Please direct any questions or comments to emily.ricotta@nih.gov

