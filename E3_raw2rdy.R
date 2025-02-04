#
#
# transform the raw datafiles from Experiment 3 (left vs. right) into a csv file
# written by Vincent Ott
#



# Import packages --------------------------------------------------------------
library(tidyverse)
library(psych)


# Create two dataframes - one for trial data and one for demographics ----------

l_all_datafiles <- list.files(path = "./E3_raw", # specify path to raw data
                              recursive = TRUE,
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



# Check format of dfs and prepare for fusion -----------------------------------

# both n columns should be multiples of 160
df_trials %>% describe()
df_dems %>% describe()

# fuse
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


# Data Audit -------------------------------------------------------------------

# while there are 81 different participants, participant 114 mistakenly
# has a subject id 113. but the trial and dem file for 114 are labeled
# correctly.
# we will thus assign the correct id to participant 114

df %>% filter(id == 113) %>% nrow() / 160

fine_df = df %>% filter(id != 113)
nrow(fine_df) / 160

fix_df = df %>% filter(id == 113)
nrow(fix_df) / 160

new_ids = c(rep(113, 160), rep(114, 160))  # for participants 113 and 114
                                           # 113 stays 113 and 114 gets 114

fix_df$id = new_ids

unique(fix_df$id)  # successful


# now we merge fix_df and fine_df again
df = rbind(fix_df, fine_df)

# check if everything worked
nrow(df) / 160
length(unique(df$id))


df = df %>% arrange(-desc(id))


# Rename and add columns to reuse code from E2 ---------------------------------

df <- rename(df, ori_source = source)

# left -> top
# right -> bottom
df = df %>% mutate(
  source = case_when(
    ori_source == "nosou" ~ "nosou",
    ori_source == "left" ~ "top",
    ori_source == "right" ~ "bottom"
  )
)
df = df %>% relocate(source, .after = ori_source)


# assign
df = df %>% mutate(
  assign = case_when(
    moni_response == "N" ~ "N",
    moni_response == "L" ~ "T",
    moni_response == "R" ~ "B"
  )
)
df = df %>% relocate(assign, .after = moni_response)


# Export -----------------------------------------------------------------------

write.csv(df, "./E3_data.csv", row.names = FALSE)

