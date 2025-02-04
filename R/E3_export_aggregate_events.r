### VEN: export aggregated data for MPT

library(dplyr)
library(tidyr)
library(readr)

df <- read_csv("data/E3_data.csv")

events_df <- df %>% 
  mutate(
    source = case_when(
      ori_source == "left" ~ "A",
      ori_source == "right" ~ "B",
      ori_source == "nosou" ~ "N"
    ),
    resp = case_when(
      assign == "T" ~ "A",
      assign == "B" ~ "B",
      assign == "N" ~ "N"
    )
  ) |> 
  count(id, item_type, source, resp) |> 
  complete(id, nesting(item_type, source, resp), fill = list(n = 0))

write_csv(events_df, "data/E3_events.csv")
