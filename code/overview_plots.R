# Preparation -------------------------------------------------------------

# load packages
pacman::p_load(tidyverse, R2jags)
library(dplyr)
library(abind)
library(ggplot2)
library(dplyr)
library(patchwork)
library(purrr)


#load data
choices <- read_rds("data/trial_summaries.rds.bz2")

choices <- choices %>% 
  
  filter(paper == "Camilleri09b" |  
           paper == "Erev10" |  
           paper == "Frey15" |  
           paper == "Gloeckner12" |  
           paper == "Gloeckner16" |  
           paper == "Hau08" |  
           paper == "Hertwig04" |  
           paper == "Hertwig10" |  
           paper == "Kellen16" |  
           paper == "Noguchi15" |  
           paper == "Rakow08" |  
           paper == "Ungemach09" |  
           paper == "Wulff12" |  
           paper == "Wulff15a"   )
           

# prepare data for CPT --------------------------------------------------------------------

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
         choice, r_switch, dom) 

#Adjust for Losses
dat_cpt2 <- dat_cpt %>% mutate(
  sprobHA = if_else(HA == 0 & sprobHA == 0, sprobLA, sprobHA) , 
  sprobLA = 1-sprobHA,
  sprobHB = if_else(HB == 0 & sprobHB == 0, sprobLB, sprobHB) , 
  sprobLB = 1-sprobHB)
  
  
  

# compute number of problems and subjects in each data set 

nproblem <- dat_cpt %>% 
  group_by(paper) %>% 
  distinct(problem) %>% 
  summarise(nproblem = n()) 

nsubjects <- dat_cpt %>% 
  group_by(paper) %>% 
  distinct(subject) %>% 
  summarise(nsubject = n()) 

overview <- left_join(nproblem, nsubjects, by=join_by(paper)) %>% 
  mutate(nchoice = nproblem * nsubject)

dat_cpt <- dat_cpt %>%  
  group_by(paper, subject) %>% 
  mutate(avg_switchrate = mean(r_switch, na.rm = TRUE)) %>% 
  ungroup() 

#Calculate Median per Subject
median_sub <- dat_cpt %>% group_by(paper, subject) %>% summarise(MedianSR = median(avg_switchrate))

#Calculate Median per Paper
median_paper <- median_sub %>% summarise(Median_SR = median(MedianSR))

#Merge Median per Paper in working datasets
dat_cpt <- merge(dat_cpt, median_paper, by = "paper")

# function to create individual plot pair of sprobHA and sprobHB for each paper
create_plot_sprob <- function(df, paper_name, bg_color) {
  p1 <- ggplot(df, aes(x = sprobHA)) +
    geom_histogram(bins = 20, fill = "skyblue", color = "black", alpha = 0.7) +
    labs(x = "sprobHA", y = "Frequency") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_text(size = 7),
      axis.title.y = element_text(size = 7),
      panel.background = element_rect(fill = bg_color, color = NA),
      plot.background = element_rect(fill = bg_color, color = NA)
    ) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.25),
      labels = function(x) ifelse(x %in% c(0, 1), as.character(x), sub("\\.", ",", as.character(x)))
    )
  
  p2 <- ggplot(df, aes(x =sprobHB)) +
    geom_histogram(bins = 20, fill = "orange", color = "black", alpha = 0.7) +
    labs(x = "sprobHB", y = "Frequency") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_text(size = 7),
      axis.title.y = element_text(size = 7),
      panel.background = element_rect(fill = bg_color, color = NA),
      plot.background = element_rect(fill = bg_color, color = NA)
    ) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.25),
      labels = function(x) ifelse(x %in% c(0, 1), as.character(x), sub("\\.", ",", as.character(x)))
    )
  
  # side by side combination of the two plots for the same paper
  combined_plot <- p1 | p2
  
  title <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = paper_name, size = 3, hjust = 0.5) +
    theme_void() 
  
  title / combined_plot + plot_layout(heights = c(0.15, 0.85)) 
}


# number of columns in grid
ncol_grid <- 4

# create the plots of sprob
plots <- dat_cpt2 %>%
  group_split(paper) %>%
  imap(~ {
    row <- ceiling(.y / ncol_grid)   
    col <- (.y - 1) %% ncol_grid + 1 
    bg_color <- ifelse((row + col) %% 2 == 0, "#F0F0F0", "white")
    create_plot_sprob(.x, as.character(.x$paper[1]), bg_color)
  })

# Showing plot as 4x4 matrix
overview_plot_sprob <- wrap_plots(plots, ncol = ncol_grid) & 
  theme(plot.margin = margin(5, 5, 5, 5)) 
overview_plot_sprob

ggsave("plots/sprob_14.jpg", plot = overview_plot_sprob)


# Function to create r_switch plot for each paper (with mean)
create_plot_r_switch <- function(df, paper_name, bg_color) {
  median_avg_switchrate <- unique(df$Median_SR)

  p <- ggplot(df, aes(x = r_switch)) +
    geom_histogram(
      bins = 20,
      fill = "skyblue",
      color = "black",
      alpha = 0.7
    ) +
    geom_vline(
      xintercept = median_avg_switchrate,
      color = "red",
      linetype = "dashed",
      size = 0.7
    ) +
    annotate(
      "text",
      x = 0.5,
      y = Inf, 
      label = paste0("Median: ", round(median_avg_switchrate, 2)),
      color = "red",
      size = 3,
      vjust = 1 
    ) +
    labs(x = "r_switch", y = "Frequency") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_text(size = 7),
      axis.title.y = element_text(size = 7),
      panel.background = element_rect(fill = bg_color, color = NA),
      plot.background = element_rect(fill = bg_color, color = NA)
    ) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.25),
      labels = function(x) ifelse(x %in% c(0, 1), as.character(x), sub("\\.", ",", as.character(x)))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) 
  
  title <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = paper_name, size = 3, hjust = 0.5) +
    theme_void()
  
  title / p + plot_layout(heights = c(0.15, 0.85))
}

# number column in grid
ncol_grid <- 4



# Generating plots
plots <- dat_cpt %>%
  group_split(paper) %>%
  map2(1:length(.), ~ {
    row <- ceiling(.y / ncol_grid)
    col <- (.y - 1) %% ncol_grid + 1
    bg_color <- ifelse((row + col) %% 2 == 0, "#F0F0F0", "white") # Schachbrettmuster
    create_plot_r_switch(.x, as.character(.x$paper[1]), bg_color)
  })

# Showing plot as 4x4 matrix
overview_plot_r_switch <- wrap_plots(plots, ncol = ncol_grid) & 
  theme(plot.margin = margin(5, 5, 5, 5))

overview_plot_r_switch

ggsave("plots/r_switch_14.jpg", plot = overview_plot_r_switch)



