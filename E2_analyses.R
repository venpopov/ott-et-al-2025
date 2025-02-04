###
###
###
###       Better Source Memory for Recognized To-Be-Forgotten Items
###       than for Recognized To-Be-Remembered Items
###
###       Experiment 2: locations top vs. bottom
###
###       Data Analysis
###       written by Vincent Ott
###
###
###


###########################################################################################
# importing packages
###########################################################################################

library(tidyverse)

library(psych)
library(moments)

library(afex)
library(emmeans)
library(jmv)


###########################################################################################
# importing the data from Experiment 2: Locations top vs. bottom
###########################################################################################

df <- read_csv("E2_data.csv")
df = df %>% arrange(-desc(id))

#
# description of the dataframe
#

# id: participant id

# trial_no: number of the trial in the actual test phase (160 per participant)

# item_type: R vs. F vs. new distractor

# source: top vs. bottom vs. nosou ("no source) in case of a new distractor

# recog: participant´s response during recognition; either "old" or "new".

# rt_recog: response time in seconds until recognition response was given

# assign: source identification event. item was either labeled as stemming from top "T", bottom "B", or being new "N"

# rt_assign: response time in seconds for source identification prompt.
# either -66 when item was regarded as new or a positive value in case of items regarded as old.

# word: stimulus word for that trial

# study_pos: position at which the word appear during the study phase. -95 in case of new distractor.

# age: participant age

# gender: m = male, w = female, d = diverse

# occupation: 0 = non-student, 1 = non-psychology student, 2 = psychology student

# strategy1: did they use a strategy? 0 = "no", 1 = "yes"

# strategy2: description of strategy, in case they used one


###########################################################################################
# defining absolute number of stimuli for E2
###########################################################################################

no_RnF_items = 40 # i.e. there were 40 R and 40 F items in E2
no_distractors = 80


###########################################################################################
# creating a new dataframe: values per participant will be stored here ===> nrow() = N
###########################################################################################

pp_df = df %>% filter(trial_no == 1)


###########################################################################################
# calculating Hit Rates, False Alarm Rates, and corrected Hit Rates: token1 data frame
###########################################################################################

# creating new column: the signal detection event of each trial (sigdec_event)
# "hit" vs. "false alarm" (there are also the two events of "miss" and "correct rejection")
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
token1 = df %>% group_by(id,item_type,sigdec_event) %>% summarize (amount = n())
token1 = token1 %>% filter(sigdec_event == "hit" | sigdec_event == "false alarm")

# hit rates and false alarm rates
token1 = token1 %>% mutate(
  rate = case_when(
    item_type == "distractor" ~ amount/no_distractors,
    TRUE ~ amount/no_RnF_items
  )
)

# discriminability rates
# they used to be named "corrected Hit rates", therefore the columns are called "cHR_R" and "cHR_F"
token1 = token1 %>% select(-amount)
token1 = token1 %>% unite("eventdepitem", item_type, sigdec_event) %>% spread(eventdepitem, rate)
token1 = token1 %>% rename(FAR = "distractor_false alarm")
token1[is.na(token1)] = 0 # in case that there are no false alarms, then FAR = 0

# rename and append to pp_df
pp_df = pp_df %>% mutate(FAR = token1$FAR)

pp_df = pp_df %>% mutate(HR_R = token1$"R item_hit")
pp_df = pp_df %>% mutate(cHR_R = HR_R - FAR)

pp_df = pp_df %>% mutate(HR_F = token1$"F item_hit")
pp_df = pp_df %>% mutate(cHR_F = HR_F - FAR)

# descriptives for FAR (for Appendix)
describe(pp_df$FAR)


###########################################################################################
# exclusion of participants:
# are there any participants, where the hit Rates for R or F items are equal to or smaller than FAR?
# this would indicate, that they did not achieve those hits due to recognition, but due their lax response-criterion
# as indicated by their FAR
###########################################################################################

# the following people remain
exclusions = pp_df %>% filter(!HR_R > FAR | !HR_F > FAR) %>% select(id) %>% pull()

pp_df = pp_df %>% filter(!id %in% exclusions)
df = df %>% filter(!id %in% exclusions)

no_excl = length(exclusions)
final_sample_size = pp_df %>% nrow() # no participant has to be excluded


###########################################################################################
# demographics
###########################################################################################
describe(pp_df$age)

gender_distribution = pp_df %>% group_by(gender) %>% summarize(absolute = n(), relative = n()/final_sample_size)

occupation_distribution = pp_df %>% group_by(occupation) %>% summarize(absolute = n(), relative = n()/final_sample_size)


###########################################################################################
# column for source monitoring event: sm_event
###########################################################################################

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


###########################################################################################
# calculating CSIM and ACSIM in token3 dataframe - token 2 is used later on
###########################################################################################

token3 = df %>% filter(item_type != "distractor")%>%
  filter(sm_event != "FN", sm_event != "RN") %>%
  group_by(id, sm_event, source) %>% summarize(frequency = n())
token3 = token3 %>% unite("placeholder", sm_event, source) %>% spread(placeholder, frequency)

# in case there is not a single event of correct identification for a combination of item type and source
# for a participant the event frequency is 0
# however, this is swallowed by the group_by function and we need to add the following line of code:
token3[is.na(token3)] = 0


# CSIM for all 4 factor combinations
token3 = token3 %>% mutate(CSIM_top_R =
                             RC_top/(RC_top + RI_top))

token3 = token3 %>% mutate(CSIM_bottom_R =
                             RC_bottom/(RC_bottom + RI_bottom))

token3 = token3 %>% mutate(CSIM_top_F =
                             FC_top/(FC_top + FI_top))

token3 = token3 %>% mutate(CSIM_bottom_F =
                             FC_bottom/(FC_bottom + FI_bottom))


###########################################################################################
# effect size function:
# adjusted partial eta squared - see Mordkoff (2019), doi: doi.org/10.1177/2515245919855053
###########################################################################################

adparn2 = function(dfeffect,dferror,Fvalue) {
  numerator = (Fvalue - 1) * dfeffect
  denominator = Fvalue * dfeffect + dferror
  print(numerator/denominator)
}

# dfeffect for ttest is = 1
# Fvalue for ttest is = t^2




                                            #######

                                    # 2 x 2 within ANOVAs #

                                            #######

###########################################################################################
# ANOVA for UNcorrected hit rates: tokenRates dataframe
###########################################################################################

# step 1: get UNcorrected hit rates for each participant depending on item type and location
# 64 participants -> 4 combinations -> dataframe should have 256 rows
tokenRates = df %>% filter(item_type != "distractor", sigdec_event == "hit") %>%
  group_by(id, item_type, source) %>% summarise(absolute_hits = n())

# get hit rates from the absolute hits we just calculated
tokenRates = tokenRates %>% mutate(hitrate = absolute_hits/(no_RnF_items/2)) # we divide no_RnF_items by 2 as we grouped_by source

# now we should just be able to feed this into aov_ez
# check whether ANOVA requirements are met: descriptive statistics
tokenRates %>% group_by(item_type, source) %>%
  summarize(m = mean(hitrate), sd = sd(hitrate), n = n(), skw = skewness(hitrate), kurt = kurtosis(hitrate))

# ANOVA
(Rates_ANOVA = aov_ez(id = "id",
                    dv = "hitrate",
                    within = c("item_type", "source"),
                    data = tokenRates))
adparn2(1,63,223.47)

# plot (just to get a feeling for the data, this is not reported in the manuscript)
emmip(Rates_ANOVA, item_type~source, CIs = T, ylab = "Hit Rate") + theme_classic()

# contrast item type: get the descriptive statistics
tokenRates %>% group_by(id,item_type) %>% summarize(tmean = mean(hitrate)) %>% 
  group_by(item_type) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())


###########################################################################################
# corrected Hit Rate ANOVA; token2 dataframe
###########################################################################################

# calculate corrected Hitrates depending on item type and location
token2 = df %>% filter(item_type != "distractor", sigdec_event == "hit") %>%
  group_by(id, item_type, source) %>% summarise(absolute_hits = n())
token2 = token2 %>% mutate(hitrate = absolute_hits/(no_RnF_items/2))
# the false alarm rate for each participant does not and can not depend on item location
# we can thus use FAR from pp_df for the correction of hitrates
FARvector = pp_df %>% select(FAR) %>% slice(rep(1:n(), each = 4)) # FAR is matched to long format of token2
token2 = cbind(token2, FARvector)
token2 = token2 %>% mutate(cHR = hitrate - FAR)

token2 = token2 %>% relocate(cHR, .after = source)
token2 = token2 %>% select(-c(absolute_hits:FAR))


# check whether ANOVA requirements are met: descriptive statistics
token2 %>% group_by(item_type, source) %>%
  summarize(m = mean(cHR), sd = sd(cHR), n = n(), skw = skewness(cHR), kurt = kurtosis(cHR))

# calculating the cHR_ANOVA
(cHR_ANOVA = aov_ez(id = "id",
                    dv = "cHR",
                    within = c("item_type", "source"),
                    data = token2))
adparn2(1,63,223.47)

# plot
emmip(cHR_ANOVA, item_type~source, CIs = T, ylab = "corrected Hit Rate") + theme_classic()

# contrast item type: simply get descriptive stats
token2 %>% group_by(id,item_type) %>% summarize(tmean = mean(cHR)) %>% 
  group_by(item_type) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())


###########################################################################################
# Recog Hit RT ANOVA, token4 dataframe
###########################################################################################

# outlier analysis
token4 = df %>% filter(sigdec_event == "hit") %>% select(-recog,-sigdec_event,-c(assign:strategy2))

boxplot(data = token4, rt_recog~id*item_type*source, ylab = "recog hit RT")
# we want to remove all circles (outliers)

l_id = pp_df %>% select(id) %>% pull()
l_item_types = c("R item", "F item")
l_sources = c("top","bottom")
out_recog_rts = c()

for (subject in l_id) {
  for (type in l_item_types) {
    for (location in l_sources) {
      
      dummytoken = token4 %>% filter(id == subject, item_type == type, source == location)
      
      outliers = boxplot.stats(dummytoken$rt_recog)$out
      
      if (length(outliers) > 0) {
        out_recog_rts = append(out_recog_rts, outliers)  
      }
      
    }
  }
}
token4 %>% nrow() # 3264 hit trials across participants

token4 = token4 %>% filter(!rt_recog %in% out_recog_rts)
token4 %>% nrow() # 226 excluded outlier trials

boxplot(data = token4, rt_recog~id*item_type*source, ylab = "recog hit RT")
# now all previous outliers should have been removed. of course there are now new outliers, since its new data

token4 = token4 %>% group_by(id,item_type,source) %>% summarize(mean_hit_rt = mean(rt_recog))

# check whether ANOVA requirements are met: descriptive statistics
token4 %>% group_by(item_type, source) %>% summarize(m = mean(mean_hit_rt), sd = sd(mean_hit_rt), n = n(),
                                                     skw = skewness(mean_hit_rt), kurt = kurtosis(mean_hit_rt))

# ANOVA
(RT_RECOG_ANOVA = aov_ez(id = "id",
                         dv = "mean_hit_rt",
                         within = c("item_type", "source"),
                         data = token4))
adparn2(1,63,48.14)
adparn2(1,63,4.29)

# plot
emmip(RT_RECOG_ANOVA, item_type~source, CIs = T, ylab = "Recognition Hit RT") + theme_classic()

# contrast for item type: descriptive stats
token4 %>% group_by(id,item_type) %>% summarize(tmean =mean(mean_hit_rt)) %>%
  group_by(item_type) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())

# contrast for source: descriptives stats
token4 %>% group_by(id,source) %>% summarize(tmean = mean(mean_hit_rt)) %>%
  group_by(source) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())


###########################################################################################
# CSIM ANOVA; token3 dataframe
###########################################################################################

df_for_CSIM_ANOVA = token3 %>% select(id, CSIM_top_R:CSIM_bottom_F) %>%
  gather(factor_combination, CSIM, 2:5) %>%
  separate(factor_combination,into=c("placeholder","source","item_type")) %>%
  select(-placeholder) %>%
  arrange(-desc(id))


# ANOVA requirements?
df_for_CSIM_ANOVA %>% group_by(item_type,source) %>% summarize(m = mean(CSIM), sd = sd(CSIM), N = n(),
                                                               skw = skewness(CSIM), kurt = kurtosis(CSIM))

# calculating the CSIM_ANOVA
(CSIM_ANOVA = aov_ez(id = "id",
                    dv = "CSIM",
                    within = c("item_type", "source"),
                    data = df_for_CSIM_ANOVA))
adparn2(1,63,4.11)

# plot
emmip(CSIM_ANOVA, item_type~source,CIs = T, ylab = "CSIM") + theme_classic()

# contrast item type: descriptive stats
df_for_CSIM_ANOVA %>% group_by(id,item_type) %>% summarize(tmean = mean(CSIM)) %>% 
  group_by(item_type) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())

# post hoc comparisons to explore interaction
(m3a = emmeans(CSIM_ANOVA, ~ item_type + source))
#Pairewise Comparisions
contrast(m3a, method="tukey")


###########################################################################################
# RT for correct source identification ANOVA
###########################################################################################

# any participant for which there are no source monitoring events of correct idenfication
# for any of the four possible combinations of item type and item location must be excluded from the SM Hit RT ANOVA.
# This is because there is no data for this participant in the corresponding ANOVA cell.

excl_from_SM_RT_ANOVA = df_for_CSIM_ANOVA %>% filter(CSIM == 0) %>% select(id) %>% pull
excl_from_SM_RT_ANOVA # four participants

token3 %>% filter(id %in% excl_from_SM_RT_ANOVA)

token5 = df %>% filter(!id %in% excl_from_SM_RT_ANOVA) %>%
  filter(sm_event == "RC" | sm_event == "FC") %>% select(-c(recog:assign)) %>% select(-sm_event,-c(word:strategy2))

# outlier analysis
boxplot(data = token5, rt_assign~id*item_type*source, ylab = "correct assign RT")
# all circles should be removed

l_id = pp_df %>% select(id) %>% filter(!id %in% excl_from_SM_RT_ANOVA) %>% pull()
l_item_types = c("R item", "F item")
l_sources = c("top","bottom")
out_assign_rts = c()

for (subject in l_id) {
  for (type in l_item_types) {
    for (location in l_sources) {
      
      dummytoken = token5 %>% filter(id == subject, item_type == type, source == location)
      
      outliers = boxplot.stats(dummytoken$rt_assign)$out
      
      if (length(outliers) > 0) {
        out_assign_rts = append(out_assign_rts, outliers)  
      }
      
    }
  }
}

token5 %>% nrow() # 1878 trials of correct source identification across participants

token5 = token5 %>% filter(!rt_assign %in% out_assign_rts)
boxplot(data = token5, rt_assign~id*item_type*source, ylab = "recog hit RT")
# now all previous outliers should have been removed. of course there are now new outliers, since its new data

token5 = token5 %>% group_by(id,item_type,source) %>% summarize(mean_assign_rt = mean(rt_assign))


# ANOVA requirements
token5 %>% group_by(item_type, source) %>% summarize(mean = mean(mean_assign_rt), sd = sd(mean_assign_rt),
                                                     skw = skewness(mean_assign_rt), kurt = kurtosis(mean_assign_rt))

# calculating ANOVA
(RT_ASSIGN_ANOVA = aov_ez(id = "id",
                         dv = "mean_assign_rt",
                         within = c("item_type", "source"),
                         data = token5))
adparn2(1,59,4.42)

# plot
emmip(RT_ASSIGN_ANOVA, item_type~source, CIs = T, ylab = "RT for correct source identification") + theme_classic()

# contrast source: descriptive stats
token5 %>% group_by(id, source) %>% summarize(tmean = mean(mean_assign_rt)) %>% 
  group_by(source) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())

# contrast item type: descriptive stats
token5 %>% group_by(id, item_type) %>% summarize(tmean = mean(mean_assign_rt)) %>% 
  group_by(item_type) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())


                                    #######

                            # Trees and Strategies #

                                    #######

###########################################################################################
# 2HT model
# count occurences of source monitoring events for subsequent analysis in multiTree
###########################################################################################

# there are fifteen different source monitoring events in the 5 tree 2HT model
five_df = df %>% group_by(item_type, source, assign) %>% summarize(event_frequency = n())
View(five_df)


###########################################################################################
# exploratory analyses: use of strategies and description of strategies
###########################################################################################
strat_usage = pp_df %>% group_by(strategy1) %>% summarize(absolute = n(), relative = n()/final_sample_size)
View(strat_usage)

# only export those participants who stated they used a strategy
strat_description = pp_df %>% filter(strategy1 == 1) %>% select(id,strategy2)

write.csv(strat_description, "./E2_strats.csv", row.names = FALSE)
