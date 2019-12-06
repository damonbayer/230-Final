library(tidyverse)
library(here)

cancer.dat <- read_csv(here("Data", "breast-cancer-wisconsin.data"),
                       col_names = c('Sample code number',
                                     'Clump Thickness',
                                     'Uniformity of Cell Size',
                                     'Uniformity of Cell Shape',
                                     'Marginal Adhesion',
                                     'Single Epithelial Cell Size',
                                     'Bare Nuclei',
                                     'Bland Chromatin',
                                     'Normal Nucleoli',
                                     'Mitoses',
                                     'Class'), na = '?') %>% 
  na.omit()