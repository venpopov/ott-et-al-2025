### VEN: export aggregated data for MPT
aggregate_event_counts <- function(data = data.frame()) {
  withr::local_package("dplyr")
  data %>% 
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
    tidyr::complete(id, tidyr::nesting(item_type, source, resp), fill = list(n = 0))
}

