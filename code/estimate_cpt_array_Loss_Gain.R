# Preparation -------------------------------------------------------------

# load pkgs
pacman::p_load(tidyverse , 
               R2jags , 
               viridis , 
               glue ,
               tools ,
               brms ,
               abind , 
               coda)

# Initialization of regression and parameter plot data frames ------------------------------------------------------------

reg_results_list <-  data.frame()  
  
parameters_plot_data <- data.frame(
    subject = integer(0),      
    alpha = numeric(0),        
    delta = numeric(0),        
    deviance = numeric(0),     
    gamma = numeric(0),        
    mu.alpha = numeric(0),     
    mu.delta = numeric(0),     
    mu.gamma = numeric(0),    
    mu.rho = numeric(0),       
    rho = numeric(0),          
    paper = character(0),      
    cat_switch = numeric(0),   
    avg_switchrate = numeric(0)
)

parameters_df <- data.frame(
  paper = character(),
  subject = numeric(),
  alpha = numeric(),
  gamma = numeric(),
  rho = numeric(),
  delta = numeric()
)


# Constitute valid papers ------------------------------------------------------------
valid_papers <- c("Gloeckner12", "Ungemach09", "Rakow08", "noguchi15","hertwig04")

#Single Test 
#valid_papers <- c("Gloeckner12")
#valid_papers <- c("Ungemach09")
#paper <- read_rds(glue("data/PreprocessedPaperData/cpt_rakow08.rds.bz2"))

#Loop across all papers using JAGS + Generate Estimates + Regression
for (p in valid_papers) {
  #paper <- read_rds(glue("data/PreprocessedPaperData/cpt_{p}.rds.bz2"))
  paper <- read_rds("data/PreprocessedPaperData/cpt_gloeckner16.rds.bz2")
  papername <- unique(paper$paper)

  #sort Subject then Gain then loss problems
  paper_sorted <- paper[order(paper$subject, paper$dom), ] %>%
  group_by(subject) %>%
    mutate(
      nGain = sum(dom == "Gain"),
      nLoss = sum(dom == "Loss"))
  
  #nGain Vector to iterate through model code
  nGain_vector <- paper_sorted %>%
    distinct(subject, nGain) %>%  
    pull(nGain) 
  
  #Get vector of the number of solved problems per participant
  nsolvedprob <- paper_sorted %>% 
    group_by(subject) %>% 
    summarise(nsolvedprob = n()) 

  nsolvedprob <- nsolvedprob$nsolvedprob

  #Create a list of subject and the respective switch rate
  switch_rate_data <- paper_sorted %>%
    select(paper, problem,subject, cat_switch, avg_switchrate,dom) %>% 
    distinct() 
  switch_rate_data <- switch_rate_data %>% distinct(subject, .keep_all = T) %>% select(subject, avg_switchrate)
    
  #switch_rate_data$subject <- seq_len(nrow(switch_rate_data))
  
  
  ## create data object for JAGS
  # create Array 
  preArrayData <- subset(paper_sorted, select = - c(index, problem,id))

  preArrayData <- preArrayData %>%
    select(-paper, paper)

  split_data <- split(preArrayData, preArrayData$subject)

  #get max number of problems 
  problems <- paper %>% distinct(problem) 
  nprob <- nrow(problems) 

  #get number of subjects
  nsub <-length(unique(paper_sorted$subject)) 

  #number of columns without subject
  num_cols <- ncol(preArrayData) - 1 

  #initiate Array
  data_array <- array(NA, dim = c(nprob, num_cols, nsub))

  # fill 3D-Array
  for (i in 1:nsub) {
    subject_data <- as.matrix(split_data[[i]][, -1])  # Data without subject ID
    data_array[1:nrow(subject_data), 1:ncol(subject_data), i] <- subject_data  # Fill the data
  }

  # check dimension
  dim(data_array)

  
  # Converge array in list again to have NAs in the end of each element vector
  # Due to control variable in model code the NAs won't be taken into account

  ## load and prepare data ---------------------------------------------------

  data_list <- list(
    nprob = nsolvedprob,          
    nsub =nsub,   
    nGain = nGain_vector,
    HA = apply(t(data_array[, 1, ]),2, as.numeric), 
    LA = apply(t(data_array[, 3, ]),2, as.numeric), 
    HB = apply(t(data_array[, 5, ]),2, as.numeric),
    LB = apply(t(data_array[, 7, ]),2, as.numeric),
    sprobHA = apply(t(data_array[, 2, ]),2, as.numeric),
    sprobLA = apply(t(data_array[, 4, ]),2, as.numeric),
    sprobHB = apply(t(data_array[, 6, ]),2, as.numeric),
    sprobLB = apply(t(data_array[, 8, ]),2, as.numeric),  
    choice = apply(t(data_array[, 9, ]) ,2, as.numeric),
    dom = apply(t(data_array[, 11, ]), 2, as.character)
  )
 

  ## fitting ------------------------------------------------------------

  params <- c("mu.alpha", "alpha", "mu.gamma", "gamma", "mu.delta", "delta", "mu.rho", "rho", "mu.lambda", "lambda") # fitted parameters that should be shown in the results

  # function for creating parameter values to initialize the MCMC chains with
  # when rnorm(1, .4, .1), then .4 is the initial value for the parameter on the desired scale

  params_init <- function(){
    list("mu.probit.alpha" = 0 , # hyper parameters (mu. prefix refers to  group level)
         "mu.probit.gamma" = 0 ,
         "mu.probit.delta" = 0 , 
         "mu.probit.rho" = -5 , 
        "mu.probit.lambda" = 0,
        "probit.alpha" = rep(0, nsub) , # individual level parameters
        "probit.gamma" = rep(0, nsub) ,
        "probit.delta" = rep(0, nsub) , 
        "probit.rho" = rep(-5, nsub) ,
        "probit.lambda" = rep(0, nsub)) 
  } 


  # run JAGS model

  mfit <- jags.parallel(
    data = data_list , # data list
    inits = params_init , # creates list of initial values for each chain
    parameters.to.save = params , 
    model.file = "code/models/CPT_hierarchical_array_Loss_Gain.txt" , # model code, see file
    n.chains = 6 , # number of MCMC chains
    n.iter = 10000 , # number of iterations (should be set much higher once it's clear that the model works)
    n.burnin = 5000 , # first 1000 samples of each chain are discarded
    n.thin = 5 , # with 1, every sample is stored, with 2, every 2nd sample is stored, ... to reduce autocorrelatons, use higher values. however, higher values require more iterations
    n.cluster = 6 , # compute MCMC chains in parallel (on different cores of the computer)
    DIC = TRUE # store all posterior samples
  )


  # sanity check 2
  data_list$sprobHA+data_list$sprobLA
  data_list$sprobHB+data_list$sprobLB

  ## store results -----------------------------------------------------------------


  # summary of posterior distributions

  estimates <- mfit$BUGSoutput$summary %>% round(4) %>%  as_tibble(rownames = "parameter")
  posteriors <- mfit$BUGSoutput$sims.matrix %>% as_tibble()

  papername <- unique(paper$paper)
  write_rds(estimates, paste("data/PostEst/cpt_", tolower(papername), "_estimates.rds.bz2", sep =""), compress = "bz2")
  write_rds(posteriors, paste("data/PostEst/cpt_", tolower(papername), "_posteriors.rds.bz2", sep =""), compress = "bz2")


  # Create summary data frame with estimates per subject
  parameter_data <- estimates %>%
   mutate(subject = as.numeric(gsub(".*\\[(\\d+)\\]", "\\1", parameter)),  
          param_name = gsub("\\[.*", "", parameter)) %>%                  
   select(subject, param_name, mean) %>%                                
    pivot_wider(names_from = param_name, values_from = mean)              

  # summary of posterior distributions
  fits <- mfit$BUGSoutput$summary %>% round(4) %>%  as_tibble(rownames = "parameter")


  # plot weighting function

  ## weighting function for group level
  mu.fits <- fits %>% 
    filter(parameter %in% c("mu.alpha", "mu.gamma", "mu.delta", "mu.rho", "mu.lambda")) %>% 
    mutate(subject = 0, 
         parameter = c("alpha", "delta", "gamma", "rho", "lambda")) 

  mu.weights <- mu.fits %>% 
   select(subject, parameter, mean) %>% 
    pivot_wider(names_from = parameter, values_from = mean) %>% 
 
    expand_grid(p = seq(0, 1, .01)) %>% 
    mutate(w = round(  (delta * p^gamma)/ ((delta * p^gamma)+(1-p)^gamma), 4)) 

  ## weighting function for individual level
    ind.fits <- fits %>% 
    filter(! parameter %in% c("mu.alpha", "mu.gamma", "mu.delta", "mu.rho", "mu.lambda", "deviance")) %>% 
    mutate(subject = rep(1:nsub, 5),
           parameter = c(rep("alpha", nsub), rep("delta", nsub),  rep("gamma", nsub), rep("rho", nsub), rep("lambda", nsub))
    )

  ind.weights <- ind.fits %>%
    select(subject, parameter, mean) %>% 
    pivot_wider(names_from = parameter, values_from = mean) %>% 
    select(-c(alpha, rho)) %>%
    expand_grid(p = seq(0, 1, .01)) %>% 
    mutate(w = round(  (delta * p^gamma)/ ((delta * p^gamma)+(1-p)^gamma), 4)) 
  
  ind.weights <- left_join(ind.weights, switch_rate_data, by = "subject")

  # Calculation median switching rate
  median_switch_rate <- median(ind.weights$avg_switchrate, na.rm = TRUE)

  # Plot with median switching rate
  prob_weight_plot <- mu.weights %>% 
    ggplot(aes(p, w)) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.25), 
      labels = scales::number_format(accuracy = 0.01)
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, by = 0.25),
      labels = scales::number_format(accuracy = 0.01) 
    ) +
    labs(
      title = papername,
      x = "p",
      y = "w(p)",
      caption = paste("Median Switching Rate:", round(median_switch_rate,2))
    ) +
    geom_line(data = ind.weights, aes(group = subject, color = avg_switchrate), alpha = 0.8) +
    scale_color_viridis_c(
      option = "plasma",  
      name = "Average Switching Rate",  
      breaks = c(0, 0.25, 0.5, 0.75, 1),  
      labels = scales::number_format(accuracy = 0.01),
      limits = c(0, 1) 
    ) +
    geom_abline(intercept = 0, slope = 1, linewidth = 1, color = "black", linetype = "dashed") +
    geom_line(linewidth = 1.2, color = "black"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      legend.position = "bottom",  
      legend.direction = "horizontal",  
      legend.box = "horizontal",  
      legend.title.align = 0.5, 
      legend.text = element_text(size = 10, color = "black"),  
      legend.title = element_text(size = 12, face = "plain", hjust = 0.5, vjust = 0.8),  
      legend.key.width = unit(2, "cm"),  
      legend.spacing.y = unit(2.5, "cm"),  
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA),
      plot.caption = element_text(size = 12, hjust = 0.5, color = "black")  
    )


  prob_weight_plot

  ggsave(filename = paste0("plots/ProbWeighting/ProbWeighting_", papername, ".png"), plot = prob_weight_plot)

  #regression data preparations
  data_regression <- fits %>%
    mutate(
      subject = as.integer(str_extract(parameter, "\\d+")), 
      parameter = str_remove(parameter, "\\[\\d+\\]")        
    ) %>%
    select(subject, parameter, mean) %>%                     
    pivot_wider(
      names_from = parameter,                                
      values_from = mean                                     
    )

  data_regression <- left_join(data_regression, switch_rate_data, by = "subject")
  parameters_plot_data <- rbind(parameters_plot_data, data_regression)
  
  parameters_df_paper<- data_regression %>%
    select(subject, alpha, gamma, rho, delta) %>%
    mutate(paper = p)  
  
  # Saving parameters for each paper in data frame
  parameters_df <- rbind( parameters_df, parameters_df_paper)
  
  #Linear regression Gamma 
  reg_model_gamma <- brm(gamma ~ avg_switchrate, 
                           data = data_regression, 
                           family = gaussian(), 
                           prior = c(set_prior("normal(0,1)", class = "b"),
                                     set_prior("normal(0,10)", class = "Intercept")),
                           chains = 4, 
                           iter = 2000, 
                           cores = parallel::detectCores())
  
  summary(reg_model_gamma)
  
  #Linear regression Delta
  reg_model_delta <- brm(delta ~ avg_switchrate, 
                           data = data_regression, 
                           family = gaussian(), 
                           prior = c(set_prior("normal(0,1)", class = "b"),
                                     set_prior("normal(0,10)", class = "Intercept")),
                           chains = 4, 
                           iter = 2000, 
                           cores = parallel::detectCores())
  
  summary(reg_model_delta)
  
  #Linear regression Alpha
  
  reg_model_alpha <- brm(alpha ~ avg_switchrate, 
                         data = data_regression, 
                         family = gaussian(), 
                         prior = c(set_prior("normal(0,1)", class = "b"),
                                   set_prior("normal(0,10)", class = "Intercept")),
                         chains = 4, 
                         iter = 2000, 
                         cores = parallel::detectCores())
  
  summary(reg_model_alpha)
  
  #Linear regression Rho 
  reg_model_rho <- brm(rho ~ avg_switchrate, 
                         data = data_regression, 
                         family = gaussian(), 
                         prior = c(set_prior("normal(0,1)", class = "b"),
                                   set_prior("normal(0,10)", class = "Intercept")),
                         chains = 4, 
                         iter = 2000, 
                         cores = parallel::detectCores())
  
  summary(reg_model_rho)
  
  
  #Store Model Summaries and create dataframe with all regression data
  model_summary_gamma <-summary(reg_model_gamma)
  gamma_coefficients <- as.data.frame(model_summary_gamma$fixed)
  gamma_distributional_params <- as.data.frame(model_summary_gamma$spec_pars)
  model_summary_delta <-summary(reg_model_delta)
  delta_coefficients <- as.data.frame(model_summary_delta$fixed)
  delta_distributional_params <- as.data.frame(model_summary_delta$spec_pars)
  model_summary_alpha<-summary(reg_model_alpha)
  alpha_coefficients <- as.data.frame(model_summary_alpha$fixed)
  alpha_distributional_params <- as.data.frame(model_summary_alpha$spec_pars)
  model_summary_rho <-summary(reg_model_rho)
  rho_coefficients <- as.data.frame(model_summary_rho$fixed)
  rho_distributional_params <- as.data.frame(model_summary_rho$spec_pars)
  
  
  # Regression Data Gamma
  gamma_coefficients <- data.frame(
    Parameter = rownames(gamma_coefficients),
    gamma_coefficients
  )
  gamma_distributional <- data.frame(
    Parameter = "sigma",
    gamma_distributional_params
  )
  gamma_data <- rbind(gamma_coefficients, gamma_distributional)
  gamma_data$Model <- "Gamma"
  
  # Regression Data Delta
  delta_coefficients <- data.frame(
    Parameter = rownames(delta_coefficients),
    delta_coefficients
  )
  delta_distributional <- data.frame(
    Parameter = "sigma",
    delta_distributional_params
  )
  delta_data <- rbind(delta_coefficients, delta_distributional)
  delta_data$Model <- "Delta"
  
  # # Regression Data Alpha
  alpha_coefficients <- data.frame(
    Parameter = rownames(alpha_coefficients),
    alpha_coefficients
  )
  alpha_distributional <- data.frame(
    Parameter = "sigma",
   alpha_distributional_params
  )
  alpha_data <- rbind(alpha_coefficients, alpha_distributional)
  alpha_data$Model <- "Alpha"
  
  # # Regression Data Rho
  rho_coefficients <- data.frame(
    Parameter = rownames(rho_coefficients),
    rho_coefficients
  )
  rho_distributional <- data.frame(
    Parameter = "sigma",
    rho_distributional_params
  )
  rho_data <- rbind(rho_coefficients, rho_distributional)
  rho_data$Model <- "Rho"
  
  # Combining all parameter data 
  combined_data <- rbind(gamma_data, delta_data, alpha_data, rho_data)
  
  
  # add papername
  combined_data$Paper <- papername
  
  # Fill each paper data in one data frame
  reg_results_list <- rbind(reg_results_list, combined_data)
}
  
#Plot Regression Avg Switchrate and Gamma
gamma_reg_results_list <- reg_results_list %>% filter(Parameter=="avg_switchrate" & Model == "Gamma")


reg_plot_gamma <- ggplot(gamma_reg_results_list, aes(x = Estimate, y = Paper)) +
  geom_errorbarh(aes(xmin = l.95..CI, xmax = u.95..CI), 
                 height = 0.2, color = "darkgreen") +
  geom_point(size = 3, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Linear Regression Average Switching Rate and Gamma",
    x = "Estimate",
    y = "Paper"
  ) +
  theme_minimal() +
  theme(
    axis.title.y  = element_text(size = 12),
    axis.text.y   = element_text(size = 10),
    axis.title.x  = element_text(size = 12),
    axis.text.x   = element_text(size = 10),
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey70", linetype = "dotted"),
    panel.grid.major.y = element_blank()
  )
ggsave(filename = paste0("plots/LinearRegression/Gamma_Regression_", papername, ".png"), plot = reg_plot_gamma ,width = 8, height = 6, dpi = 300)

#Plot Regression Avg Switchrate and Alpha
alpha_reg_results_list <- reg_results_list %>% filter(Parameter=="avg_switchrate" & Model == "Alpha")

reg_plot_alpha <- ggplot(alpha_reg_results_list, aes(x = Estimate, y = Paper)) +
  geom_errorbarh(aes(xmin = l.95..CI, xmax = u.95..CI), 
                 height = 0.2, color = "orange") +
  geom_point(size = 3, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Linear Regression Average Switching Rate and Alpha",
    x = "Estimate",
    y = "Paper"
  ) +
  theme_minimal() +
  theme(
    axis.title.y  = element_text(size = 12),
    axis.text.y   = element_text(size = 10),
    axis.title.x  = element_text(size = 12),
    axis.text.x   = element_text(size = 10),
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey70", linetype = "dotted"),
    panel.grid.major.y = element_blank()
  )
ggsave(filename = paste0("plots/LinearRegression/Alpha_Regression_", papername, ".png"), plot = reg_plot_alpha ,width = 8, height = 6, dpi = 300)

#Plot Regression Avg Switchrate and Delta
delta_reg_results_list <- reg_results_list %>% filter(Parameter=="avg_switchrate" & Model == "Delta")

reg_plot_delta <- ggplot(delta_reg_results_list, aes(x = Estimate, y = Paper)) +
  geom_errorbarh(aes(xmin = l.95..CI, xmax = u.95..CI), 
                 height = 0.2) +
  geom_point(size = 3, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Linear Regression Average Switching Rate and Delta",
    x = "Estimate",
    y = "Paper"
  ) +
  theme_minimal() +
  theme(
    axis.title.y  = element_text(size = 12),
    axis.text.y   = element_text(size = 10),
    axis.title.x  = element_text(size = 12),
    axis.text.x   = element_text(size = 10),
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey70", linetype = "dotted"),
    panel.grid.major.y = element_blank()
  )
ggsave(filename = paste0("plots/LinearRegression/Delta_Regression_", papername, ".png"), plot = reg_plot_delta,width = 8, height = 6, dpi = 300)

#Plot Regression Avg Switchrate and Rho
rho_reg_results_list <- reg_results_list %>% filter(Parameter=="avg_switchrate" & Model == "Rho")

reg_plot_rho <- ggplot(rho_reg_results_list, aes(x = Estimate, y = Paper)) +
  geom_errorbarh(aes(xmin = l.95..CI, xmax = u.95..CI), 
                 height = 0.2, color = "purple") +
  geom_point(size = 3, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Linear Regression Average Switching Rate and Rho",
    x = "Estimate",
    y = "Paper"
  ) +
  theme_minimal() +
  theme(
    axis.title.y  = element_text(size = 12),
    axis.text.y   = element_text(size = 10),
    axis.title.x  = element_text(size = 12),
    axis.text.x   = element_text(size = 10),
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey70", linetype = "dotted"),
    panel.grid.major.y = element_blank()
  )
ggsave(filename = paste0("plots/LinearRegression/Rho_Regression_", papername, ".png"), plot = reg_plot_rho,width = 8, height = 6, dpi = 300)


#Loop parameter/avg switching rate for each paper
for (p in valid_papers) {
  papername <- toTitleCase(p)
  
  plot_data <- subset(parameters_plot_data, paper == papername)
  
  
  plot <- ggplot(plot_data, aes(x = avg_switchrate, y = gamma)) +
    geom_point(size = 4, alpha = 0.8, color = "darkgreen",na.rm = TRUE) + 
    labs(
      title = paste("Gamma-Distribution:", papername),
      subtitle = "Relationship between Gamma-Values and Switching Rate",  
      x = "Average Switchrate",
      y = "Gamma-Values"
    ) +
    theme_minimal(base_size = 14) + 
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray50"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    ) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "darkred")  
  
 
  
 ggsave(filename = paste0("plots/ParameterDistribution/Gamma_Distribution_", papername, ".png"), plot = plot,width = 8, height = 6, dpi = 300)
  
  plot <- ggplot(plot_data, aes(x = avg_switchrate, y = alpha)) +
    geom_point(size = 4, alpha = 0.8, color = "orange", na.rm = TRUE) +  
    labs(
      title = paste("Alpha-Distribution:", papername),
      subtitle = "Relationship between Alpha-Values and Switching Rate", 
      x = "Average Switchrate",
      y = "Alpha-Values"
    ) +
    theme_minimal(base_size = 14) +  
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray50"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    ) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "darkred")  
  
  ggsave(filename = paste0("plots/ParameterDistribution/Alpha_Distribution_", papername, ".png"), plot = plot,width = 8, height = 6, dpi = 300)
  
  plot <- ggplot(plot_data, aes(x = avg_switchrate, y = delta)) +
    geom_point(size = 4, alpha = 0.8, na.rm = TRUE) +  
    labs(
      title = paste("Delta-Distribution:", papername),
      subtitle = "Relationship between Delta-Values and Switching Rate",  
      x = "Average Switchrate",
      y = "Delta-Values"
    ) +
    theme_minimal(base_size = 14) +  
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray50"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    ) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "darkred")  
  
  ggsave(filename = paste0("plots/ParameterDistribution/Delta_Distribution_", papername, ".png"), plot = plot,width = 8, height = 6, dpi = 300)
  
  plot <- ggplot(plot_data, aes(x = avg_switchrate, y = rho)) +
    geom_point(size = 4, alpha = 0.8, color = "purple",na.rm = TRUE) + 
    labs(
      title = paste("Rho-Distribution:", papername),
      subtitle = "Relationship between Rho-Values and Switching Rate",  
      x = "Average Switchrate",
      y = "Rho-Values"
    ) +
    theme_minimal(base_size = 14) + 
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray50"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    ) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dotted", color = "darkred")  
  
  
  
    ggsave(filename = paste0("plots/ParameterDistribution/Rho_Distribution_", papername, ".png"), plot = plot,width = 8, height = 6, dpi = 300)
  
}

