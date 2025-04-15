# README --------------------------------------------------------------------

'The dataset of Wulff et al. (2018) can be retrieved from: https://www.dirkwulff.org/#data)'

# Preparation -------------------------------------------------------------

# load packages 
pacman::p_load(tidyverse)

# read data (can be retrieved from: ) 
data <- read.table("data/exp.txt") %>% as_tibble()

unique(data$paper)
num_part <- data  %>% group_by(paper) %>% summarise(sum_sub = n_distinct(subject))
num_dec <- data  %>% group_by(paper, subject) %>% summarise(sum_problem = n_distinct(problem)) 
num_dec <- num_dec  %>% group_by(paper) %>% summarise(num_problem = sum(sum_problem))

table <-merge(num_dec,num_part, by = "paper")
sum_row <- data.frame(
  paper = "Total",
  num_problem = sum(table$num_problem, na.rm = TRUE),
  sum_sub = sum(table$sum_sub, na.rm = TRUE)
)

# Append sum row to the dataframe
table <- rbind(table, sum_row)

# View the updated dataframe
print(table)
####

# Preprocessing -----------------------------------------------------------

dat <- data %>% 
  
  rename(
    
    sample = trial, # sample number (trial is confusing as it is typically used to refer to a problem that a person solved)
    attended = option # sampled option
    
  ) %>% 
  
  mutate( 
    
    ## ISSUE: some papers treat the same outcome as different outcome or did not specify the second probability; needs to be fixed
    
    # determine which outcome from which option was sampled to compute sampled probabilities (see below)      
    seenA1 = if_else(attended == 0 & outcome == outA1, 1, 0) ,
    seenA2 = if_else(attended == 0 & probA2 > 0 & outcome == outA2, 1, 0) ,
    seenA3 = if_else(attended == 0 & probA3 > 0 & outcome == outA3, 1, 0) ,
    seenA4 = if_else(attended == 0 & probA4 > 0 & outcome == outA4, 1, 0) ,
    seenA5 = if_else(attended == 0 & probA5 > 0 & outcome == outA5, 1, 0) ,
    
    seenB1 = if_else(attended == 1 & outcome == outB1, 1, 0) ,
    seenB2 = if_else(attended == 1 & probB2 > 0 & outcome == outB2, 1, 0) , 
    seenB3 = if_else(attended == 0 & probB3 > 0 & outcome == outB3, 1, 0) ,
    seenB4 = if_else(attended == 0 & probB4 > 0 & outcome == outB4, 1, 0) ,
    seenB5 = if_else(attended == 0 & probB5 > 0 & outcome == outB5, 1, 0) ,
    
    #Making sure that if two cases are the same within either option A or B it is treated as a safe option
    
    probB1 = ifelse(probB1 == probB2 & seenB1 == seenB2 & outB1 == outB2, 1, probB1),
    probB2 = ifelse(probB1 == 1 & seenB1 == seenB2 & outB1 == outB2, 0, probB2),
    seenB1 = ifelse(probB1 == probB2 & seenB1 == seenB2 & outB1 == outB2, 1, seenB1),
    seenB2 = ifelse(probB1 == 1 & seenB1 == seenB2 & outB1 == outB2, 0, seenB2)
    
    ) %>% 
  
  group_by(paper, id, subject, problem) %>% # to compute variables on trial level 
  
  mutate(
    
    # compute sample sizes 
    
    n_sample = max(sample) , # sample size
    n_sample_1 = sum(attended) , # sample size option 1/B
    n_sample_0 = n_sample - n_sample_1 , # sample size option 0/A
    
    # computes sampled probabilities 
    
    sprobA1 = sum(seenA1)/n_sample_0 , # sampled probability of outcome A1
    sprobA2 = sum(seenA2)/n_sample_0 , 
    sprobA3 = sum(seenA3)/n_sample_0 , 
    sprobA4 = sum(seenA4)/n_sample_0 ,
    sprobA5 = sum(seenA5)/n_sample_0 , 

    sprobB1 = sum(seenB1)/n_sample_1 , 
    sprobB2 = sum(seenB2)/n_sample_1 , 
    sprobB3 = sum(seenB3)/n_sample_1 , 
    sprobB4 = sum(seenB4)/n_sample_1 ,
    sprobB5 = sum(seenB5)/n_sample_1 ,
    
    # compute switch rate 
    
    switch = ifelse(attended != lag(attended), 1, 0) , # did switch occur
    n_switch = sum(switch, na.rm = TRUE) , # number of switches
    r_switch = round(n_switch/((n_sample - 1)), 2) , # observed switch rate
         
  ) %>% 
  
  ungroup()

unique(dat$paper)
#Prefilter events that don't add up to a prob of 1 or sprob of 1
dat <- dat %>% filter((probA1 + probA2 + probA3 + probA4 + probA5 == 1) & 
                        (probB1 + probB2 + probB3 + probB4 + probB5 ==1) & 
                        (sprobA1 + sprobA2 + sprobA3 + sprobA4 + sprobA5 ==1) & 
                        (sprobB1 + sprobB2 + sprobB3 + sprobB4 + sprobB5 ==1))



## ISSUE: Not all sampled probabilities add up to 1 (see ISSUE comment above); paper are excluded for now, but issue should be fixed

excl <-  dat %>% 
  
  mutate(
  
  checkA = round(sprobA1 + sprobA2 + sprobA3 + sprobA4 + sprobA5, 2) , 
  checkB = round(sprobB1 + sprobB2 + sprobB3 + sprobB4 + sprobB5, 2)
  
  ) %>%
  
  filter(checkA != 1 | checkB != 1) %>% 
  
  distinct(paper)

dat <- dat %>% filter(! paper %in% excl$paper)

# remove redundant rows (only one row per trial)

choices <- dat %>% 
  group_by(paper, id, subject, problem) %>% 
  mutate(stop = ifelse(sample == n_sample, 1, 0) ) %>% 
  ungroup() %>% 
  filter(stop == 1)

unique(choices$paper)
# Storing -----------------------------------------------------------------

write_rds(choices, "data/trial_summaries.rds.bz2", compress = "bz2")
