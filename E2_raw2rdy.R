
# this script transforms the raw datafiles from Experiment 2 (top vs. bottom) into a csv file
# written by Vincent Ott


###########################################################################################
# import packages
###########################################################################################

library(tidyverse)
library(psych)


###########################################################################################
# import individual txt-files as a whole dataframe
###########################################################################################

l_all_datafiles <- list.files(path = "./E2_raw", recursive = TRUE, # specify path to raw data
                          pattern = "\\.txt$", 
                          full.names = TRUE)

list_of_trials = list()
list_of_dem = list()

for (i in l_all_datafiles) {
  
  if (grepl("trials", i, fixed = TRUE)) {
    len = length(list_of_trials)
    list_of_trials[[len+1]] = i
    
  } else {
    len = length(list_of_dem)
    list_of_dem[[len+1]] = i
  }

}


df_trials <- read_delim(list_of_trials[[1]], delim = "\t", quote = "")
for (i in 1:length(list_of_trials)) {
  if (i > 1) {
    frame <- read_delim(list_of_trials[[i]], delim = "\t", quote = "")
    df_trials <- bind_rows(df_trials, frame)
  }
}

df_dems <- read_delim(list_of_dem[[1]], delim = "\t", quote = "")
for (i in 1:length(list_of_dem)) {
  if (i > 1) {
    frame <- read_delim(list_of_dem[[i]], delim = "\t", quote = "")
    df_dems <- bind_rows(df_dems, frame)
  }
}

View(df_trials)
View(df_dems)

# there are now two dataframes. one contains the trial data for all participants
# the other one contains the demographic information for all participants
# now they need to be amalmagated


###########################################################################################
# check that the dfs are not messed up and prepare for amalmagation
###########################################################################################

df_trials %>% describe()
df_dems %>% describe()

# amalmagate
df = cbind(df_trials, df_dems)
View(df)

# renaming columns and variables
df = rename(df, id = subj_id)
df = df %>% mutate(
  item_type = case_when(
    item_type == "*EEE*" ~ "R item",
    item_type == "*VVV*" ~ "F item",
    item_type == "*nnn*" ~ "distractor"
  )
)

# remove unnecessary columns
df <- df %>% select(-"subj_id_no")


View(df)


###########################################################################################
# export DF
###########################################################################################

write.csv(df, "./E2_data.csv", row.names = FALSE)

