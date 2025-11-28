# MERMAID fish data download and clean

rm(list=ls())

# Set up ====

library(mermaidr)
library(tidyverse)
library(janitor)
library(here)

time_it <- function(expr, with_sound = F, print_header = '') {
  start_time <- Sys.time()
  result <- eval(expr, envir = parent.frame())
  end_time <- Sys.time() - start_time
  
  print(print_header)
  cat("Elapsed Time: ", end_time, "\n")
  if(with_sound) { 
    beepr::beep(sound = 2) 
    Sys.sleep(2)
  }
  
  invisible(result)  # Return the result without printing it
}

out_dir <- "Z:/Dropbox/Lenfest_kenya/data/ecological/marine/coral/cover/field_observations/global/2001_to_ONGOING/mermaid/"

out_name <- "mermaid_all_public_summaries_site_level_fish"


target.methods <- c("beltfish" #,
                    # "benthiclit", 
                    # "benthicpit",
                    # "benthicpqt" #, 
                    #"habitatcomplexity", 
                    #"bleachingqc"
                    )

target.policies <- c(#"Private",
                     #"Public",
                     "Public Summary"
                     )

target.aggregation_levels <- c(#"observations", 
                               #"sampleunits", 
                               "sampleevents"
                               )


id_vars.projects <- c('id', 
                      'name', 
                      'countries', 
                      'num_sites')


# Main ====

if(!dir.exists(out_dir)) dir.create(out_dir, recursive = T)


# Get projects
projects <- mermaid_get_projects()


# Check policy distribution
projects.long <- projects %>% 
  pivot_longer(cols = c("data_policy_beltfish",
                        "data_policy_benthiclit", 
                        "data_policy_benthicpit",
                        "data_policy_benthicpqt", 
                        "data_policy_habitatcomplexity", 
                        "data_policy_bleachingqc"),
               names_to = 'method',
               values_to = 'policy'
               )

ggplot(data = projects.long, aes(x = method, fill= policy)) +
  geom_bar() +
  
  scale_fill_manual(values = c("#AA2222", "#22AA22", "#AA6622")) +
  
  theme_bw(base_size = 16)



# Get and filter data
set_ups <- expand_grid(method = target.methods,
                       policy = target.policies)

project_data <- lapply(1:nrow(set_ups), function(x) {
  
  method_x <- set_ups$method[x]
  policy_x <- set_ups$policy[x]
  
  
  projects_x <- projects %>% 
    filter(get(paste0('data_policy_',method_x)) %in% c(policy_x)) %>% 
    select(any_of(c(id_vars.projects, paste0('data_policy_',method_x)))) %>% 
    filter(num_sites > 0) %>% 
    arrange(num_sites)
  
 
  # Test getting one observation
  
  method_x2 <- ifelse(method_x == "beltfish", 'fishbelt', method_x)
  
  data_x <- projects_x %>% 
    mermaid_get_project_data(method_x2, target.aggregation_levels) %>%
    time_it(T)
  
  data_x <- data_x %>%
    mutate(method = method_x,
           policy = policy_x) %>%
    relocate(method, policy)
  
}) %>% time_it(T)

data.all <- bind_rows(project_data)

write_csv(x = data.all,
          file = paste0(out_dir,out_name,'_up_to_',most_recent_data,'.csv'))
