# this script transforms the raw excel file from Experiment 1 (color) into a cleaner csv file
# written by Vincent Ott

library(tidyverse)
library(readxl)

df <- read_excel("data-raw/E1_raw.xlsx")

# cleaning up a lot
# note: i created a superfluous dataframe "rdf" for "clean dataframe"
rdf = df
rdf = rdf %>% select(-ExperimentName)
rdf = rdf %>% select(-Session)

rdf = rdf %>% relocate(Gruppe,.after = Subject)

rdf = rdf %>% select(-c(`Clock.Information`:ExperimentVersion))
rdf = rdf %>% select(-c(RandomSeed:RuntimeVersionExpected))
rdf = rdf %>% select(-c(StudioVersion:welche))

rdf = rdf %>% relocate(Trial, .after = Subject)
rdf = rdf %>% relocate(c(col,font_source), .before = Age)
rdf = rdf %>% relocate(FontResponse, .before = Age)
rdf = rdf %>% relocate(itemtype, .before = Age)
rdf = rdf %>% relocate(`Procedure[Trial]`, .before = Age)
rdf = rdf %>% relocate(studyword, .before = Age)
rdf = rdf %>% relocate(testword, .before = Age)

rdf = rdf %>% filter(!is.na(Trial))

rdf = rdf %>% filter(`Procedure[Trial]` == "testprocnoda")
rdf = rdf %>% select(-c(`Procedure[Trial]`:studyword))

rdf = rdf %>% relocate(ONResponse, .before = Age)

rdf = rdf %>% arrange(-desc(Subject))
rdf = rdf %>% select(-col)

rdf = rdf %>% relocate(ONTest.RT, .before = Age)
rdf = rdf %>% relocate(position_source, .before = Age)
rdf = rdf %>% relocate(PosResponse, .before = Age)
rdf = rdf %>% relocate(RKResponse, .before = Age)
rdf = rdf %>% relocate(`RKTest.RT`, .before = Age)
rdf = rdf %>% relocate(`SMTest.RT`, .before = Age)

rdf = rdf %>% select(-c(UebungsFarbe:Uebungswort))


# now we can arrange the important columns in the correct order at the beginning (first few columns)
rdf = rdf %>% relocate(Gruppe, .before = Trial)
rdf = rdf %>% relocate(position_source, .after = Trial)

rdf = rdf %>% select(-itemtype)

rdf = rdf %>% relocate(c(ONResponse,ONTest.RT), .after = font_source)
rdf = rdf %>% relocate(PosResponse, .before = FontResponse)
rdf = rdf %>% relocate(RKResponse, .before = PosResponse)
rdf = rdf %>% relocate(`RKTest.RT`, .after = RKResponse)
rdf = rdf %>% relocate(`SMTest.RT`, .after = FontResponse)

# renaming and recoding varibles
rdf = rdf %>% rename(item_type = position_source)
rdf = rdf %>% rename(color = font_source)
rdf = rdf %>% rename(recog = ONResponse)
rdf = rdf %>% rename(tag = PosResponse)
rdf = rdf %>% rename(ori_assign = FontResponse)

rdf = rdf %>% mutate(
  item_type = recode(item_type, "new" = "distractor", "EEE" = "R item", "VVV" = "F item"))

rdf$tag = rdf$tag %>% replace_na("nnn")
rdf$RKResponse = rdf$RKResponse %>% replace_na("new")
rdf$ori_assign = rdf$ori_assign %>% replace_na("N")

rdf = rdf %>% mutate(
  ori_assign = recode(ori_assign, "cyan" = "C", "yellow" = "Y"),
  tag = recode(tag, "EEE" = "R", "VVV" = "F", "nnn" = "N"))

# add columns with no additional information, which will allow us to reuse a lot the code from Experiment 2 (location)
# cyan -> top
# yellow -> bottom
rdf = rdf %>% mutate(
  source = case_when(
    color == "new" ~ "nosou",
    color == "cyan" ~ "top",
    color == "yellow" ~ "bottom"
  )
)
rdf = rdf %>% relocate(source, .after = color)

rdf = rdf %>% mutate(
  assign = case_when(
    ori_assign == "N" ~ "N",
    ori_assign == "C" ~ "T",
    ori_assign == "Y" ~ "B"
  )
)
rdf = rdf %>% relocate(assign, .after = ori_assign)


# export
write.csv(rdf, "data/E1_data.csv", row.names = FALSE)
