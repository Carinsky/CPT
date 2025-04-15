# Preparation -------------------------------------------------------------

# load packages
pacman::p_load(tidyverse, R2jags, abind, patchwork)
library(dplyr)
library(abind)
library(ggplot2)
library(dplyr)
library(patchwork)
library(purrr)

# load data
choices <- read_rds("data/trial_summaries.rds.bz2")

unique(choices$paper)
unique(choices$dom)
# Data Selection --------------------------------------------------------------------


dat_cpt <- choices %>% 
  
  # get problems where sampling is autonomous, outcomes are gains only, and with up to 2 outcomes per option
  
  filter(type == "free", 
         dom == "Gain" | dom == "Loss", 
         probA3 == 0 & probA4 == 0 & probA5 == 0, 
         probB3 == 0 & probB4 == 0 & probB5 == 0) %>% 
  
  # identify higher and lower outcome for each option (as CPT uses rank-dependent transformations)
  
  mutate(
    
    HA = if_else(outA1 > outA2, outA1, outA2) , 
    LA = if_else(outA1 < outA2, outA1, outA2) ,
    HB = if_else(outB1 > outB2, outB1, outB2) , 
    LB = if_else(outB1 < outB2, outB1, outB2) ,
    sprobHA = if_else(outA1 > outA2, sprobA1, sprobA2) ,
    sprobLA = if_else(outA1 < outA2, sprobA1, sprobA2) , 
    sprobHB = if_else(outB1 > outB2, sprobB1, sprobB2) ,
    sprobLB = if_else(outB1 < outB2, sprobB1, sprobB2)
    
  ) %>%  
  
  # select only relevant variables
  
  select(index:subject , 
         HA, sprobHA, LA, sprobLA , 
         HB, sprobHB, LB, sprobLB , 
         choice,r_switch, dom) 

unique(dat_cpt$paper)

# compute number of problems and subjects in each data set 

nproblem <- dat_cpt %>% 
  group_by(paper) %>% 
  distinct(problem) %>% 
  summarise(nproblem = n()) 

nsubjects <- dat_cpt %>% 
  group_by(paper) %>% 
  distinct(subject) %>% 
  summarise(nsubject = n()) 

# overview of max problems and subjects per paper 
overview <- left_join(nproblem, nsubjects, by=join_by(paper)) %>% 
  mutate(nchoice = nproblem * nsubject)

# check for paper picks
table(dat_cpt$paper)

papernames <- unique(dat_cpt$paper)


#loop across paper to save the data separetely 

for (i in papernames) {
  # pick paper
  papername <- i
  paper <- dat_cpt %>% filter(paper == papername)
  
  # check for NAs in sprob and remove the respective rows
  NAcolumns <- paper  %>% filter(is.na(sprobHA) == TRUE | is.na(sprobLA) == TRUE |is.na(sprobHB) == TRUE|is.na(sprobLB) == TRUE)
  
  paper <- setdiff(paper,NAcolumns)
  
  # replace NAs with 0
  paper[is.na(paper)] <-0
  
  # get vector with solved problems per participant
  nsolvedprob <- paper %>% 
    group_by(subject) %>% 
    summarise(nsolvedprob = n()) 
  
  # merge nsolvedprob in original working table
  paper <- merge(paper,nsolvedprob, by = "subject")
  
  # vector of solved problems per participant
  nsolvedprob <- nsolvedprob$nsolvedprob
  
  # check duplicates (add only for Frey, Hau, Erev, Wullff)
  paper_unique <- paper %>%
    group_by(subject) %>%
    mutate(unique_subject_problem = row_number()) %>%
    ungroup() %>%
    select(-problem) %>% # replace old problem column with unique ID
    rename(problem = unique_subject_problem)
  
  # ensure right order of the columns ((add only for Frey, Hau, Erev, Wullff))
  paper <- paper_unique
  paper <- paper %>%
    select(1:4, problem, everything()) 
  
  # calculate average switchrate per subject
  paper <- paper %>%
    group_by(subject) %>% 
    mutate(avg_switchrate = mean(r_switch, na.rm = TRUE), min_switchrate = min(r_switch, na.rm = TRUE), max_switchrate = max(r_switch, na.rm = TRUE)) %>% 
    ungroup() 
  
  #Get median
  median_switch <- median(paper$avg_switchrate)
  paper$cat_switch = ifelse(paper$avg_switchrate>= median_switch,1,0)
  
  ## data that is the same for all subjects (but see ISSUES above)
  
  problems <- paper %>% distinct(problem) 
  
  ## data that differs between subjects 
  ## create data frames that store the trial/problem-specific choices and sampled probabilities for each subject
  ## in model code referred to as [j=subject=rows, i=problem=columns]
  
  choices <- paper %>% select(subject, problem, choice) %>% pivot_wider(names_from = problem, names_prefix = "p", values_from = choice) %>% select(-subject)
  sprobHA <- paper %>% select(subject, problem, sprobHA) %>% pivot_wider(names_from = problem, names_prefix = "p", values_from = sprobHA) %>% select(-subject)
  sprobLA <- paper %>% select(subject, problem, sprobLA) %>% pivot_wider(names_from = problem, names_prefix = "p", values_from = sprobLA) %>% select(-subject)
  sprobHB <- paper %>% select(subject, problem, sprobHB) %>% pivot_wider(names_from = problem, names_prefix = "p", values_from = sprobHB) %>% select(-subject)
  sprobLB <- paper %>% select(subject, problem, sprobLB) %>% pivot_wider(names_from = problem, names_prefix = "p", values_from = sprobLB) %>% select(-subject)
  
  nprob <- nrow(problems)
  nsub <- nrow(choices)
  
  # sanity checks (should be true)
  sprobHA+sprobLA == 1
  sprobHA+sprobLA
  sprobHB+sprobLB == 1
  nrow(choices) == nrow(sprobHA) 
  nrow(sprobHA) ==  nrow(sprobLA)
  nrow(sprobLA) ==  nrow(sprobHB)
  nrow(sprobHB) ==  nrow(sprobLB)
  
  write_rds(paper, paste("data/PreprocessedPaperData/cpt_", tolower(papername), ".rds.bz2", sep =""), compress = "bz2")
}


# Dataset: XY -------------------------------------------------------------

