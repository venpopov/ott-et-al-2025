###
###
###
###       Better Source Memory for Recognized To-Be-Forgotten Items
###       than for Recognized To-Be-Remembered Items
###
###       Experiment 1: colors cyan vs. yellow
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
### Source A (cyan) is referred to as "top" (Source A in Experiment 2) and
### Source B (yellow) is referred to as "bottom" (Source B in Experiment 2)
###
###





###########################################################################################
# importing packages
###########################################################################################

library(tidyverse)

library(psych)
library(moments)

library(car) # for leveneTest()

library(afex)
library(emmeans)
library(jmv)


###########################################################################################
# importing data from Experiment 1
###########################################################################################

df <- read_csv("data/E1_data.csv")
df = df %>% rename(age = Age, trial_no = Trial, id = Subject, gender = Sex) %>%
  mutate(PresCon = Gruppe) %>%
  mutate(PresCon = if_else(PresCon == "IN", "intrinsic", "extrinsic"))

#
# description of the dataframe. only for the relevant columns used in the conducted analyses.
#

# id: participant id

# Gruppe: Presentation Condition. IN = intrinsic, EX = extrinsic
# PresCon: Presentation Contdition. intrinsic vs. extrinsic.        YES! these two are redundant! only PresCon is used here
# PresCon is at the very end of the dataframe, last column

# trial_no: number of the trial in the actual test phase (104 per participant)

# item_type: R vs. F vs. new distractor

# color: cyan vs. yellow vs. new distractor
# source: top vs. bottom vs. nosou ("no source) in case of a new distractor   only the source column is used here

# recog: participant´s response during recognition; either "old" or "new".

# ONTest.RT: response time in milliseconds until recognition response was given

# RKResponse: "remember" vs. "know" response to the corresponding prompt after an "old" reponse was given

# tag: tag that participants assigned to the testword -> "R" vs. "F" vs. "N" for new distractor

# ori_assign: the original source assignment participants made: either cyan "C" or yellow "Y" or new distractor "N"
# assign: source identification event. item was either labeled as stemming from top "T", bottom "B", or being new "N"
# only the assign column is used here.

# age: participant age

# gender: männlich = male, weiblich = female

# studies: whether they study psychology "Psychologie" or something different "Anderes"

# Strategieabfrage.RESP: whether they used a strategy "j" or didnt "n"


###########################################################################################
# defining absolute number of stimuli in Experiment 1
###########################################################################################

no_RnF_items = 32 # there were 32 R items and 32 F items
no_distractors = 40


###########################################################################################
# creating a new dataframe: values per participant will be stored here ===> nrow() = N
###########################################################################################

pp_df = df %>% filter(trial_no == 1)


###########################################################################################
# calculating Hit Rates, False Alarm Rates, and corrected Hit Rates (discriminability rates): token1 data frame
###########################################################################################

# creating a column: the signal detection event of each trial (sigdec_event)
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


# calculating the UNcorrected rates for each participant in token1 dataframe
token1 = df %>% group_by(id,item_type,sigdec_event) %>% summarize (amount = n())
token1 = token1 %>% filter(sigdec_event == "hit" | sigdec_event == "false alarm")
token1 = token1 %>% mutate(
  rate = case_when(
    item_type == "distractor" ~ amount/no_distractors,
    TRUE ~ amount/no_RnF_items
  )
)

# calculating corrected hit rates (discriminability rates)
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


# descriptives FAR for Appendix
pp_df %>% group_by(PresCon) %>% summarize(describe(FAR))
# significant difference? (not reported)
ttestIS(data = pp_df, vars = "FAR", group = "PresCon")


# Comparison of correct rejections between the two presentation conditions
# out of interest
correj_df = df %>% filter(sigdec_event == "corr_reject") %>% 
  group_by(id, PresCon) %>% summarize(amount = n())
correj_df = correj_df %>% mutate(
  rate = amount/no_distractors,
)
# significant difference? (not reported)
ttestIS(data = correj_df, vars = "rate", group = "PresCon", desc = TRUE)
# turns out it is the exact same as for False Alarm rate!


###########################################################################################
# exclusion of participants:
# are there any participants, where the uncorrected Hit Rate for R items is equal to or smaller than FAR?
# this would indicate, that they did not achieve those hits due to recognition, but due their lax response-criterion
# as indicated by their FAR
###########################################################################################

exclusions = pp_df %>% filter(HR_R <= FAR) %>% select(id) %>% pull()
# this criterion applies to no participant

# but what about F items? (perfect directed forgetting?)
perfect = pp_df %>% filter(HR_F <= FAR) %>% select(id) %>% pull()
# participant 48 shows perfect directed forgetting as indicated by their HR_F which is equal to or even lower than FAR
pp_df %>% filter(id == 48) %>% select(FAR, HR_R, HR_F)


pp_df = pp_df %>% filter(!id %in% exclusions)
df = df %>% filter(!id %in% exclusions)

no_excl = length(exclusions)
final_sample_size = pp_df %>% nrow() # no participant has to be excluded


###########################################################################################
# demographics
###########################################################################################
describe(pp_df$age)

gender_distribution = pp_df %>% group_by(gender) %>% summarize(absolute = n(), relative = n()/final_sample_size)
View(gender_distribution)

occupation_distribution = pp_df %>% group_by(studies) %>%
    summarize(absolute = n(), relative = n()/final_sample_size)

# PresCon distribution
pp_df %>% group_by(PresCon) %>% summarize(n = n())


###########################################################################################
# column for tagging event: tag_event
# token0 dataframe

# analyses refer to the conditional correct tagging measure (CCTM)
# in analogy to the CSIM measure, we only analyse correct tagging among recognized R and F items

###########################################################################################

token0 = df %>% filter(sigdec_event == "hit") %>%
  mutate(tag_event = case_when(
           item_type == "R item" & tag == "R" ~ "RCT", # CT stands for correct tag
           item_type == "R item" & tag == "F" ~ "RIT", # IT stands for incorrect tag
           item_type == "F item" & tag == "F" ~ "FCT",
           item_type == "F item" & tag == "R" ~ "FIT"
           ))
token0 = relocate(token0, tag_event, .after = "tag")

token0 = token0 %>% group_by(id, PresCon,source,tag_event) %>% summarize(frequency = n())
token0 = token0 %>% unite("placeholder", source, tag_event) %>% spread(placeholder, frequency)
token0[is.na(token0)] = 0

# calculate Conditional Correct Tagging Measure CCTM for all 4 within-factor combinations
token0 = token0 %>% mutate(CCTM_top_R =
                             top_RCT/(top_RCT + top_RIT))

token0 = token0 %>% mutate(CCTM_bottom_R =
                             bottom_RCT/(bottom_RCT + bottom_RIT))

token0 = token0 %>% mutate(CCTM_top_F =
                             top_FCT/(top_FCT + top_FIT))

token0 = token0 %>% mutate(CCTM_bottom_F =
                             bottom_FCT/(bottom_FCT + bottom_FIT))


###########################################################################################
# creating column for source monitoring event: sm_event
###########################################################################################

df = df %>% mutate(
  sm_event = case_when(
    item_type == "R item" & source == "top" & assign == "T" ~ "RC",
    item_type == "R item" & source == "top" & assign == "B" ~ "RI",
    item_type == "R item" & source == "top" & assign == "N" ~ "RN",
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
# calculating CSIM and ACSIM, token3 dataframe - token 2 is used later on
###########################################################################################

token3 = df %>% filter(item_type != "distractor")%>%
  filter(sm_event != "FN", sm_event != "RN") %>%
  group_by(id, PresCon, sm_event, source) %>% summarize(frequency = n())
token3 = token3 %>% unite("placeholder", sm_event, source) %>% spread(placeholder, frequency)

# in case there is not a single event of correct identification for a combination of item type and source
# for a participant, the event frequency is thus 0
# however, this is swallowed by the group_by function and we need the following line of code:
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
  print(numerator/denominator) }

# dfeffect for ttest is = 1
# Fvalue for ttest is = t^2





                                            #######

                                      # Behavioral Analyses #

                   # Mixed ANOVAs (2 within factors and one between factor) #

                                            #######

###########################################################################################
# ANOVA for UNcorrected hit rates: tokenRates dataframe
###########################################################################################

tokenRates = df %>% filter(item_type != "distractor", sigdec_event == "hit") %>%
  group_by(id, PresCon, item_type, source) %>% summarise(absolute_hits = n())

# get hit rates from the absolute hits we just calculated
tokenRates = tokenRates %>% mutate(hitrate = absolute_hits/(no_RnF_items/2))

# now we should just be able to feed this into aov_ez

# check whether ANOVA requirements are met: descriptive statistics
tokenRates %>% group_by(PresCon, item_type, source) %>%
  summarize(m = mean(hitrate), sd = sd(hitrate), n = n(), skw = skewness(hitrate), kurt = kurtosis(hitrate))

# Levene-Test
leveneTest(hitrate ~ PresCon, data = tokenRates)

# ANOVA
(Rates_ANOVA = aov_ez(id = "id",
                      dv = "hitrate",
                      within = c("item_type", "source"),
                      between = c("PresCon"),
                      data = tokenRates))
adparn2(1,99,427.66)

# item type contrast: descriptive stats
tokenRates %>% group_by(id,item_type) %>% summarize(tmean = mean(hitrate)) %>% 
  group_by(item_type) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())


###########################################################################################
# corrected Hit Rate ANOVA; token2 dataframe
###########################################################################################

# calculate corrected Hitrates depending on item type and source
token2 = df %>% filter(item_type != "distractor", sigdec_event == "hit") %>%
  group_by(id, PresCon, item_type, source) %>% summarise(absolute_hits = n())
token2 = token2 %>% mutate(hitrate = absolute_hits/(no_RnF_items/2))

# the false alarm rate for each participant does not and can not depend on item location
# we can thus use FAR from pp_df for the correction of hitrates
FARvector = pp_df %>% select(FAR) %>% slice(rep(1:n(), each = 4)) # FAR is matched to long format of token2
token2 = cbind(token2, FARvector)
token2 = token2 %>% mutate(cHR = hitrate - FAR)

# check ANOVA requirements
token2 %>% group_by(PresCon, item_type, source) %>%
  summarize(m = mean(cHR), sd = sd(cHR), n = n(), skw = skewness(cHR), kurt = kurtosis(cHR))

# Levene-Test
leveneTest(cHR ~ PresCon, data = token2)

# calculating the cHR_ANOVA
(cHR_ANOVA = aov_ez(id = "id",
                    dv = "cHR",
                    within = c("item_type", "source"),
                    between = c("PresCon"),
                    data = token2))
adparn2(1,99,427.66)
adparn2(1,99,6.40)

# item type contrast: descriptive stats
token2 %>% group_by(id,item_type) %>% summarize(tmean = mean(cHR)) %>% 
  group_by(item_type) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())

# presentation condition contrast: descriptive stats
token2 %>% group_by(id,PresCon) %>% summarize(tmean = mean(cHR)) %>% 
  group_by(PresCon) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())


###########################################################################################
# Recog Hit RT ANOVA
###########################################################################################

# outlier analysis
token4 = df %>% filter(sigdec_event == "hit")
boxplot(data = token4, ONTest.RT~id*item_type*source, ylab = "recog hit RT")
# all circles (outliers) should be removed

l_id = pp_df %>% select(id) %>% pull()
l_item_types = c("R item", "F item")
l_sources = c("top","bottom")
out_recog_rts = c()

for (subject in l_id) {
  for (type in l_item_types) {
    for (location in l_sources) {
      
      dummytoken = token4 %>% filter(id == subject, item_type == type, source == location)
      
      outliers = boxplot.stats(dummytoken$ONTest.RT)$out
      
      if (length(outliers) > 0) {
        out_recog_rts = append(out_recog_rts, outliers)  
      }
      
    }
  }
}
token4 %>% nrow() #4324 #308 (7%)

token4 = token4 %>% filter(!ONTest.RT %in% out_recog_rts)
boxplot(data = token4, ONTest.RT~id*item_type*source, ylab = "recog hit RT")
# now all previous outliers should have been removed. of course there are now new outliers, since its new data

token4 = token4 %>% group_by(id, PresCon, item_type,source) %>% summarize(mean_hit_rt = mean(ONTest.RT))
token4 = token4 %>% filter(!id %in% exclusions)


# check ANOVA requirements
token4 %>% group_by(PresCon, item_type, source) %>% summarize(m = mean(mean_hit_rt), sd = sd(mean_hit_rt),
                                                              n = n(), skw = skewness(mean_hit_rt), kurt = kurtosis(mean_hit_rt))

# Levene-Test
leveneTest(mean_hit_rt ~ PresCon, data = token4)

# calculating ANOVA
(RT_RECOG_ANOVA = aov_ez(id = "id",
                         dv = "mean_hit_rt",
                         within = c("item_type", "source"),
                         between = c("PresCon"),
                         data = token4))
adparn2(1,99,79.47)
adparn2(1,99,6.11)

# plot
emmip(RT_RECOG_ANOVA, PresCon~item_type~source, CIs = T, ylab = "Recognition Hit RT") + theme_classic()
# it looks like that the interaction of presentation condition and item type pertains to extrinsic vs intrinsic F items

# contrast item type: descriptive stats
token4 %>% group_by(id,item_type) %>% summarize(tmean = mean(mean_hit_rt)) %>% 
  group_by(item_type) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())



# investigating interaction of item type and presentation condition
# post hoc comparisons

# if we collapse over the factor of color, then there are four combinations/cells
# here, the Ms and SDs for the four comparison cells are calculated:
token4 %>% group_by(id, item_type, PresCon) %>% summarize(tmean = mean(mean_hit_rt)) %>%
  group_by(item_type, PresCon) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())

# but how do we calculate the actual comparisons?
# i dont think that we can rely on emmeans here:
em = emmeans(RT_RECOG_ANOVA, specs = c("item_type", "PresCon"))
contrast(em, method = "tukey")
# the degrees of freedom always look like a between-participant comparison was conducted

# two of the post hocs happen within each presentation condition
# R item intrinsic vs. R item extrinsic
# and F item intrinsic vs. F item extrinsic

# one cannot arrive at df = 99 here even if we were to calculate a between comparison,
# then the dfs should add up to:
# t(n1 + n2 - 2) = t(100) respectively t(98)
# since intrinsic N = 51 and extrinsic N = 50

# but if the data allows for a within comparison, then why not do just that?
# i.e. t(n - 1) = t(50) respectively t(49)?

# accordingly, we will calculate the four between-participant comparisons
# and the two within-participant comparisons manually with the jamovi package

# to hedge against alpha-error inflation, we then simply employ a bonferroni-correction
# each result will only be interpreted as significant, if the p-value
# is below .05/6 = .008

# we start with token4 which stores the mean hit RT for each participant in dependence of source and item type
View(token4)
# but for the post-hocs we collapse over the source factor, thus:
ph_token4 = token4 %>% group_by(id, PresCon, item_type) %>% summarize(meanRT = mean(mean_hit_rt))

# now we split this dataframe into four cells (combinations of PresCon and item type):
ph_Rin = ph_token4 %>% filter(PresCon == "intrinsic", item_type == "R item")
ph_Rex = ph_token4 %>% filter(PresCon == "extrinsic", item_type == "R item")
ph_Fin = ph_token4 %>% filter(PresCon == "intrinsic", item_type == "F item")
ph_Fex = ph_token4 %>% filter(PresCon == "extrinsic", item_type == "F item")


# BETWEEN-participant comparisons:

# R item Intrinsic - R Item Extrinsic       ### congruent with emmeans
ph1 = rbind(ph_Rin, ph_Rex)
ttestIS(data = ph1,
        vars = "meanRT",
        group = "PresCon",
        hypothesis = "different")

# R item Intrinsic - F item Extrinsic       ### incongruent
ph2 = rbind(ph_Rin, ph_Fex)
ttestIS(data = ph2,
        vars = "meanRT",
        group = "PresCon",
        hypothesis = "different")
adparn2(1,99,4.41^2)

# R item Extrinsic - F item Intrinsic       ### incongruent
ph3 = rbind(ph_Rex, ph_Fin)
ttestIS(data = ph3,
        vars = "meanRT",
        group = "PresCon",
        hypothesis = "different")
adparn2(1,99,5.20^2)

# F item Intrinsic - F item Extrinsic         ### congruent with emmeans
ph4 = rbind(ph_Fin, ph_Fex)
ttestIS(data = ph4,
        vars = "meanRT",
        group = "PresCon",
        hypothesis = "different")


# WITHIN-comparisons:

# R item intrinsic - F item Intrinsic
t.test(ph_Rin$meanRT, y = ph_Fin$meanRT, paired = TRUE)
adparn2(1,50,6.73^2)

# R item Extrinsic - F item Extrinsic
t.test(ph_Rex$meanRT, y = ph_Fex$meanRT, paired = TRUE)
adparn2(1,49,6.14^2)


###########################################################################################
# ANOVA for Remember/Know Response
###########################################################################################

# compare the proportions of recognized R items (vs. recognized F items) with "remember" responses
# take into account, that F items receive fewer hits and thus there are fewer absolute remember responses possible

tokenRK = df %>% filter(item_type != "distractor", sigdec_event == "hit") %>% group_by(id, PresCon, item_type, source) %>% 
  summarize(absolute_hits = n()) #404 columns

help_vector = df %>% filter(item_type != "distractor", RKResponse == "remember") %>% 
  group_by(id, PresCon, item_type, source) %>% summarize(R_responses = n())

help_vector = help_vector %>% mutate(item_type = if_else(item_type == "F item", "F", "R"))

help_vector = help_vector %>% unite("placeholder",item_type,source) %>% spread(placeholder,R_responses)
help_vector[is.na(help_vector)] = 0
help_vector = help_vector %>% gather(item_type, r_responses, 3:6) %>%
  separate(item_type, into=c("item_type", "source")) %>% arrange(-desc(id))

tokenRK$r_responses = help_vector$r_responses
tokenRK = tokenRK %>% mutate(account = r_responses/absolute_hits)

tokenRK = tokenRK %>% filter(!id %in% exclusions)

# check ANOVA requirements
tokenRK %>% group_by(PresCon, item_type, source) %>% summarize(m = mean(account), sd = sd(account), n = n(),
                                                               skw = skewness(account), kurt = kurtosis(account))

# levene test
leveneTest(account ~ PresCon, data = tokenRK)

(RK_ANOVA = aov_ez(id = "id",
                   dv = "account",
                   within = c("item_type", "source"),
                   between = c("PresCon"),
                   data = tokenRK))
adparn2(1,99,282.33)

# plot
emmip(RK_ANOVA, item_type~source~PresCon, CIs = T, ylab = "account R reponses") + theme_classic()

# contrast item type: descriptive stats
tokenRK %>% group_by(id,item_type) %>% summarize(tmean = mean(account)) %>% 
  group_by(item_type) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())


###########################################################################################
# ANOVA for tagging: token0 dataframe (CCTM: conditional correct tagging measure)
###########################################################################################

df_for_CCTM_ANOVA = token0 %>% select(id, PresCon, CCTM_top_R:CCTM_bottom_F) %>%
  gather(factor_combination, CCTM, 3:6) %>% 
  separate(factor_combination, into=c("placeholder", "source", "item_type")) %>% 
  select(-placeholder) %>% 
  arrange(-desc(id))

df_for_CCTM_ANOVA = df_for_CCTM_ANOVA %>% filter(!id %in% exclusions)

# check ANOVA requirements
df_for_CCTM_ANOVA %>% group_by(PresCon, item_type, source) %>% summarize(m = mean(CCTM), sd = sd(CCTM), n = n(),
                                                                         skw = skewness(CCTM), kurt = kurtosis(CCTM))

# levene
leveneTest(CCTM ~ PresCon, data = df_for_CCTM_ANOVA)
adparn2(1,402,6.76)

# calculating the CCTM_ANOVA
(CCTM_ANOVA = aov_ez(id = "id",
                     dv = "CCTM",
                     within = c("item_type", "source"),
                     between = c("PresCon"),
                     data = df_for_CCTM_ANOVA))
adparn2(1,99,4.31)
adparn2(1,99,0)

# plot
emmip(CCTM_ANOVA, item_type~source~PresCon, CIs = T, ylab = "CCTM") + theme_classic()

# contrast prescon: descriptive stats
df_for_CCTM_ANOVA %>% group_by(id,PresCon) %>% summarize(tmean = mean(CCTM)) %>% 
  group_by(PresCon) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())

# contrast item type: descriptive stats
df_for_CCTM_ANOVA %>% group_by(id,item_type) %>% summarize(tmean = mean(CCTM)) %>% 
  group_by(item_type) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())


###########################################################################################
# CSIM color ANOVA; token3 dataframe
###########################################################################################

df_for_CSIM_ANOVA = token3 %>% select(id,PresCon,CSIM_top_R:CSIM_bottom_F) %>%
  gather(factor_combination, CSIM, 3:6) %>%
  separate(factor_combination,into=c("placeholder","source","item_type")) %>%
  select(-placeholder) %>%
  arrange(-desc(id))

# check ANOVA requirements
df_for_CSIM_ANOVA %>% group_by(PresCon, item_type, source) %>% summarize(m = mean (CSIM), sd = sd(CSIM),
                                                                         skw = skewness(CSIM), kurt = kurtosis(CSIM))

# levene
leveneTest(CSIM ~ PresCon, data = df_for_CSIM_ANOVA)

# calculating the CSIM_ANOVA
(CSIM_ANOVA = aov_ez(id = "id",
                    dv = "CSIM",
                    within = c("item_type", "source"),
                    between = c("PresCon"),
                    data = df_for_CSIM_ANOVA))
adparn2(1,99,32.04)
adparn2(1,99,14.57)

# plot to investigate whether the interaction effect is not completely crossed
emmip(CSIM_ANOVA, item_type~source, CIs = T, ylab = "CSIM") + theme_classic()
emmip(CSIM_ANOVA, source~item_type, CIs = T, ylab = "CSIM") + theme_classic()

# contrast color: descriptive stats
df_for_CSIM_ANOVA %>% group_by(id,source) %>% summarize(tmean = mean(CSIM)) %>% 
  group_by(source) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())

# contrast item type: descriptive stats
df_for_CSIM_ANOVA %>% group_by(id,item_type) %>% summarize(tmean = mean(CSIM)) %>% 
  group_by(item_type) %>% summarize(mean = mean(tmean), sd = sd(tmean), n = n())


# post hoc comparisons to investigate interaction
# cells averaged across presentation condition
df_for_CSIM_ANOVA %>% group_by(item_type, source) %>% summarize(m = mean(CSIM), sd = sd(CSIM), n = n())

# post hoc contrasts
(m3a = emmeans(CSIM_ANOVA, ~ item_type + source))
#Pairewise Comparisions
contrast(m3a, method="tukey")

# effect sizes
adparn2(1,99,7.26^2)
adparn2(1,99,4.36^2)
adparn2(1,99,5.54^2)
adparn2(1,99,3.71^2)
adparn2(1,99,1.76^2)
adparn2(1,99,1.72^2)





                                           #######

          # Extraction of Frequencies for Model-Based Analyses in multiTree (Moshagen) #

                                           #######


###########################################################################################
# 5 tree 2HT model: intrinsic PresCon only
###########################################################################################

in_five_df = df %>% filter(PresCon == "intrinsic") %>% group_by(item_type, source, assign) %>% 
  summarize(event_frequency = n())

###########################################################################################
# 5 tree 2HT model: extrinsic PresCon only
###########################################################################################

ex_five_df = df %>% filter(PresCon == "extrinsic") %>% group_by(item_type, source, assign) %>% 
  summarize(event_frequency = n())

###########################################################################################
# 5 tree 2HT model: aggregated across both PresCons // collapsing between-factor
###########################################################################################

five_df = df %>% group_by(item_type, source, assign) %>% summarize(event_frequency = n())




                                              #######

                                        # Use of Strategies #

                                              #######

pp_df %>% group_by(Strategieabfrage.RESP) %>% summarize(n = n(), r = n/final_sample_size)
