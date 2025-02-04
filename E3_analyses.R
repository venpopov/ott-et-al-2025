###
###
###
###       Better Source Memory for Recognized To-Be-Forgotten Items
###       than for Recognized To-Be-Remembered Items
###
###       Experiment 3: locations left vs. right
###
###       Data Analysis
###       written by Vincent Ott
###
###
###



###
###
### IMPORTANT NOTE:
###
###
### In order to be able to reuse code from Experiment 2,
### Source A (left) is referred to as "top" (Source A in Experiment 2) and
### Source B (right) is referred to as "bottom" (Source B in Experiment 2)
###
###



# Libraries --------------------------------------------------------------------

library(tidyverse)

library(psych)
library(moments)

library(afex)
library(emmeans)
library(jmv)



### Load dataframe ---------------------------------------------------------------

df <- read_csv("E3_data.csv")


#
# description of the dataframe
#

# id: participant id

# trial_no: trial number in the actual test phase (160 per participant)

# item_type: R vs. F vs. new distractor

# ori_source: the original source in the experiment, left vs. right vs. nosou
# ("no source) in case of a new distractor

# source: top vs. bottom vs. nosou
# CREATED TO BE ABLE TO RE-USE THE CODE FROM EXPERIMENT 2

# moni_response: participant´s response during source monitoring
# L (left) vs. R (right) vs. N (new)

# assign: T (top) vs. B (bottom) vs. N (new)
# CREATED TO BE ABLE TO RE-USE THE CODE FROM EXPERIMENT 2

# rt_moni: response time in seconds until moni_response was given

# word: stimulus word for that trial

# study_pos: position at which the word appeared during the study phase.
# -95 in case of new distractor.

# age: participant age

# gender: m = male, w = female, d = diverse

# occupation: non_stud (not a student) vs.
# stud_non_psy (not a psychology student)
# vs. psy_stud (psychology student)

# strategy1: did they use a strategy? 0 = "no", 1 = "yes"

# strategy2: description of strategy, in case they used one



##### Absolute number of stimuli in E3 ---------------------------------------------

no_RnF_items = 40 # i.e. there were 40 R and 40 F items in E2
no_distractors = 80



##### pp_df ------------------------------------------------------------------------
# Dataframe to store values per participant; thus nrow() = N

pp_df = df %>% filter(trial_no == 1)



##### Hit rates, False Alarm rates, corrected Hit Rates ----------------------------
# token1 dataframe

# recog column
df <- df %>% mutate(
  recog = case_when(
    
    item_type == "R item" & moni_response == "L" ~ "old",
    item_type == "R item" & moni_response == "R" ~ "old",
    
    item_type == "F item" & moni_response == "L" ~ "old",
    item_type == "F item" & moni_response == "R" ~ "old",
    
    item_type == "distractor" & moni_response == "L" ~ "old",
    item_type == "distractor" & moni_response == "R" ~ "old",
    
    TRUE ~ "new"
    
  )
)
df = df %>% relocate(recog, .after = moni_response)


# sigdec_event
df = df %>% mutate(
  sigdec_event = case_when(
    item_type == "R item" & recog == "old" ~ "hit",
    item_type == "R item" & recog == "new" ~ "miss",
    item_type == "F item" & recog == "old" ~ "hit",
    item_type == "F item" & recog == "new" ~ "miss",
    item_type == "distractor" & recog == "old" ~ "false alarm",
    TRUE ~ "corr_reject"
  ))
df = df %>% relocate(sigdec_event, .after = recog)



# token1 dataframe: absolute hits and false alarms
token1 = df %>% group_by(id, item_type, sigdec_event) %>% summarize (amount = n())

# Participants 110 and 140 have no hits for F items!
# Perfect directed forgetting!
token1 %>% filter(id %in% c(110, 140))

# This should mean that nrow(df) is not divisible by 6 * 81
nrow(df) / (6 * 81)

token1 = token1 %>% filter(sigdec_event == "hit" | sigdec_event == "false alarm")

# hit rates and false alarm rates
token1 = token1 %>% mutate(
  rate = case_when(
    item_type == "distractor" ~ amount/no_distractors,
    TRUE ~ amount/no_RnF_items
  )
)

token1 %>% filter(id %in% c(110, 140))  # No hits & rates for F items!


# discriminability rates
# they used to be named "corrected Hit rates",
# therefore the columns are called "cHR_R" and "cHR_F"
token1 = token1 %>% select(-amount)
token1 = token1 %>% unite("eventdepitem", item_type, sigdec_event) %>%
  spread(eventdepitem, rate)
token1 = token1 %>% rename(FAR = "distractor_false alarm")
token1[is.na(token1)] = 0 # in case that there are no false alarms, then FAR = 0



token1 %>% filter(id %in% c(110, 140))  # Now correctly set to 0
sum(is.na(token1))



# rename and append to pp_df
pp_df = pp_df %>% mutate(FAR = token1$FAR)

pp_df = pp_df %>% mutate(HR_R = token1$"R item_hit")
pp_df = pp_df %>% mutate(cHR_R = HR_R - FAR)

pp_df = pp_df %>% mutate(HR_F = token1$"F item_hit")
pp_df = pp_df %>% mutate(cHR_F = HR_F - FAR)

hist(pp_df$HR_R, breaks = seq(0, 1, by = 0.05))

hist(pp_df$HR_F, breaks = c(0, 0.01, 0.05, seq(0.1, 1, by = 0.05)))
pp_df %>% filter(HR_F <= 0.05) %>% select(id, c(FAR:cHR_F))

hist(pp_df$FAR, breaks = c(0, 0.001, 0.05, seq(0.1, 1, by = 0.01)))

##### Exclusion of participants ----------------------------------------------------

# Exclude those with HR_R <= FAR
exclusions = pp_df %>% filter(!HR_R > FAR) %>% select(id) %>% pull()
no_excl = length(exclusions)

pp_df = pp_df %>% filter(!id %in% exclusions)
df = df %>% filter(!id %in% exclusions)
(final_sample_size = pp_df %>% nrow())


##### Demographics -----------------------------------------------------------------

boxplot(pp_df$age, horizontal = TRUE)
pp_df %>% filter(age < 18)  # remove for age descriptives
pp_df %>% filter(age >= 18) %>% select(age) %>%
  pull() %>% as.numeric() %>% describe()

occupation_distribution = pp_df %>% group_by(occupation) %>%
  summarize(absolute = n(), relative = n()/final_sample_size)

gender_distribution = pp_df %>% group_by(gender) %>% 
  summarize(absolute = n(), relative = n()/final_sample_size)


##### Source Monitoring Events -----------------------------------------------

df = df %>% mutate(
  sm_event = case_when(
    item_type == "R item" & source == "top" & assign == "T" ~ "RC", # "C" stands for "correct"
    item_type == "R item" & source == "top" & assign == "B" ~ "RI", # "I" for "incorrect"
    item_type == "R item" & source == "top" & assign == "N" ~ "RN", # "N" for "new"
    item_type == "R item" & source == "bottom" & assign == "B" ~ "RC",
    item_type == "R item" & source == "bottom" & assign == "T" ~ "RI",
    item_type == "R item" & source == "bottom" & assign == "N" ~ "RN",
    
    item_type == "F item" & source == "top" & assign == "T" ~ "FC",
    item_type == "F item" & source == "top" & assign == "B" ~ "FI",
    item_type == "F item" & source == "top" & assign == "N" ~ "FN",
    item_type == "F item" & source == "bottom" & assign == "B" ~ "FC",
    item_type == "F item" & source == "bottom" & assign == "T" ~ "FI",
    item_type == "F item" & source == "bottom" & assign == "N" ~ "FN",
    
    item_type == "distractor" & assign == "N" ~ "NN",
    item_type == "distractor" & assign == "T" ~ "NT",
    item_type == "distractor" & assign == "B" ~ "NB"
  )
)
df = relocate(df, sm_event, .after = "assign")



##### Source Identification Measures -----------------------------------------------

# token3 dataframe
# (token2 is used below for corrected hit rates ANOVA)

token3 = df %>% filter(item_type != "distractor")%>%
  filter(sm_event != "FN", sm_event != "RN") %>%
  group_by(id, sm_event, source) %>% summarize(frequency = n())

token3 = token3 %>% unite("placeholder", sm_event, source) %>%
  spread(placeholder, frequency)


# Remove participants 110 and 114 as they have no F item hits.
# So, we cannot calculate CSIM and ACSIM values for them.
token3 = token3 %>% filter(!id %in% c(110, 140))
nrow(token3)


# in case there is not a single event of correct identification
# for a combination of item type and source
# for a participant the event frequency is 0

# however, this is swallowed by the group_by function and
# we need to add the following line of code:
token3[is.na(token3)] = 0


# CSIM for all 4 factor combinations of item type and source
token3 = token3 %>% mutate(CSIM_top_R =
                             RC_top/(RC_top + RI_top))

token3 = token3 %>% mutate(CSIM_bottom_R =
                             RC_bottom/(RC_bottom + RI_bottom))

token3 = token3 %>% mutate(CSIM_top_F =
                             FC_top/(FC_top + FI_top))

token3 = token3 %>% mutate(CSIM_bottom_F =
                             FC_bottom/(FC_bottom + FI_bottom))


##### Effect size function ---------------------------------------------------------
# adjusted partial eta squared
# see Mordkoff (2019), doi: doi.org/10.1177/2515245919855053

adparn2 = function(dfeffect,dferror,Fvalue) {
  numerator = (Fvalue - 1) * dfeffect
  denominator = Fvalue * dfeffect + dferror
  print(numerator/denominator)
}

# dfeffect for ttest is = 1
# Fvalue for ttest is = t^2



##### FAR descriptives (for Appendix) ----------------------------------------------

describe(pp_df$FAR)


##### ANOVA: UN corrected hit rates ------------------------------------------------

# Step 1
# get UNcorrected hit rates for each participant depending on item type and location
# N participants -> 4 combinations -> dataframe should have N*4 rows
tokenRates = df %>% filter(item_type != "distractor", sigdec_event == "hit") %>%
  group_by(id, item_type, source) %>% summarise(absolute_hits = n())

tokenRates$item_type = substr(tokenRates$item_type, 1, 1)


# Some participants might have no counts for some of the
# item_type x source combinations
tokenRates = tokenRates %>% unite("placeholder", item_type, source) %>%
  spread(placeholder, absolute_hits)
sum(is.na(tokenRates))  # checks out! due to participants 110 and 140
tokenRates[is.na(tokenRates)] = 0

# Revert to long format
tokenRates = tokenRates %>% gather(item_type, absolute_hits, 2:5)
tokenRates = tokenRates %>% separate(item_type, into = c("item_type", "source"))


# get hit rates from the absolute hits we just calculated
tokenRates = tokenRates %>% mutate(hitrate = absolute_hits/(no_RnF_items/2))
# we divide no_RnF_items by 2 as we grouped_by source


# now we should just be able to feed this into aov_ez
# check whether ANOVA requirements are met: descriptive statistics
tokenRates %>% group_by(item_type, source) %>%
  summarize(m = mean(hitrate), sd = sd(hitrate), n = n(),
            skw = skewness(hitrate), kurt = kurtosis(hitrate))

# ANOVA
(Rates_ANOVA = aov_ez(id = "id",
                      dv = "hitrate",
                      within = c("item_type", "source"),
                      data = tokenRates))
# calculate effect sizes with adparn2()
adparn2(1, 80, 205.14)

# plot (just to get a feeling for the data, this is not reported in the manuscript)
emmip(Rates_ANOVA, item_type~source, CIs = T, ylab = "Hit Rate") + theme_classic()

# contrast item type: get the descriptive statistics
tokenRates %>% group_by(id, item_type) %>% summarize(tmean = mean(hitrate)) %>% 
  group_by(item_type) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())




##### ANOVA: corrected hit rates ---------------------------------------------------

# calculate corrected Hitrates depending on item type and location

token2 = tokenRates


# the false alarm rate for each participant does not and can not depend on item location
# we can thus use FAR from pp_df for the correction of hitrates
FARvector = pp_df %>% select(FAR) %>% slice(rep(1:n(), each = 4)) # FAR is matched to long format of token2
token2 = cbind(token2, FARvector)
token2 = token2 %>% mutate(cHR = hitrate - FAR)

hist(token2$cHR, breaks = 20)


token2 = token2 %>% relocate(cHR, .after = source)
token2 = token2 %>% select(-c(absolute_hits:FAR))


# check whether ANOVA requirements are met: descriptive statistics
token2 %>% group_by(item_type, source) %>%
  summarize(m = mean(cHR), sd = sd(cHR), n = n(),
            skw = skewness(cHR), kurt = kurtosis(cHR))

# calculating the cHR_ANOVA
(cHR_ANOVA = aov_ez(id = "id",
                    dv = "cHR",
                    within = c("item_type", "source"),
                    data = token2))
# calculate effect sizes with adparn2()
adparn2(1, 80, 176.05)

# plot
emmip(cHR_ANOVA, item_type~source, CIs = T, ylab = "corrected Hit Rate") + theme_classic()

# contrast item type: simply get descriptive stats
token2 %>% group_by(id,item_type) %>% summarize(tmean = mean(cHR)) %>% 
  group_by(item_type) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())




##### ANOVA: CSIM ------------------------------------------------------------------

df_for_CSIM_ANOVA = token3 %>% select(id, CSIM_top_R:CSIM_bottom_F) %>%
  gather(factor_combination, CSIM, 2:5) %>%
  separate(factor_combination,into=c("placeholder","source","item_type")) %>%
  select(-placeholder) %>%
  arrange(-desc(id))

sum(is.na(df_for_CSIM_ANOVA))


# ANOVA requirements?
df_for_CSIM_ANOVA %>% group_by(item_type,source) %>%
  summarize(m = mean(CSIM), sd = sd(CSIM), N = n(),
            skw = skewness(CSIM), kurt = kurtosis(CSIM))

# calculating the CSIM_ANOVA
(CSIM_ANOVA = aov_ez(id = "id",
                     dv = "CSIM",
                     within = c("item_type", "source"),
                     data = df_for_CSIM_ANOVA))
# no effects, so no effect sizes to calculate

# plot
emmip(CSIM_ANOVA, item_type~source,CIs = T, ylab = "CSIM") + theme_classic()

# contrast item type: descriptive stats
df_for_CSIM_ANOVA %>% group_by(id,item_type) %>% summarize(tmean = mean(CSIM)) %>% 
  group_by(item_type) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())



##### Frequencies Source Monitoring Events -----------------------------------------

five_df = df %>% group_by(item_type, source, assign) %>% summarize(event_frequency = n())
View(five_df)


##### Strategies -------------------------------------------------------------------

strat_usage = pp_df %>% group_by(strategy1) %>%
  summarize(absolute = n(), relative = n()/final_sample_size)

View(strat_usage)

# One participant did not answser whether they used a strategy.
# But they are someone different from the person who did not answer age.
pp_df %>% filter(age == -1)
pp_df %>% filter(strategy1 == -1)
# Did 111 still describe their strategy?
pp_df %>% filter(strategy1 == -1) %>%
  select(id, age, strategy1, strategy2) # yes, they did!


# Are there people who denied to have used a strat but described something
# nonetheless?
pp_df %>% filter(strategy1 == 0) %>% select(strategy2)


# export those participants described a strategy
strat_description = pp_df %>% filter(strategy1 != 0) %>% select(id, strategy2)
nrow(strat_description) / 81  # used a strat


write.csv(strat_description, "./E3_strats.csv", row.names = FALSE)
