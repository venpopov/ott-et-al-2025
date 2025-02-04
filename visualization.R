###
###
###
###       Better Source Memory for Recognized To-Be-Forgotten Items than for
###       Recognized To-Be-Remembered Items
###
###       Visualisation of Parameter Estimates
###
###       written by Vincent Ott
###       
###       a tutorial that was very helpful can be found here:
###       https://www.youtube.com/watch?v=AmeVAVxvSqg
###



### Libraries --------------------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(ggpubr)


### Item memory parameters -------------------------------------------------------

# Estimates, SEs, and CIs are from the multiTree analyses

df_item <- data.frame(
  
  Experiment = c(1, 1, 2, 2, 3, 3),
  Parameter = c("DR", "DF", "DR", "DF", "DR", "DF"),
  
  #             DR   DF
  Estimate = c(.76, .26,  # E1
               .73, .28,  # E2
               .69, .20   # E3
              ),
  
  
  #       DR   DF
  SE = c(.01, .02,  # E1
         .01, .02,  # E2
         .01, .02   # E3
        )

)

df_item <- df_item %>% mutate(
  CI_lower = Estimate - SE * 1.96,
  CI_upper = Estimate + SE * 1.96,
  )

# Reorder factors
df_item$Experiment <- factor(df_item$Experiment, levels = c(1, 2, 3))
df_item$Parameter <- factor(df_item$Parameter, levels = c("DR", "DF"))



### Source memory parameters -----------------------------------------------------

# Estimates, SEs, and CIs are from the multiTree analyses

df_source <- data.frame(
  
  Experiment = c(1, 1, 2, 2, 3, 3),
  Parameter = c("dR", "dF", "dR", "dF", "dR", "dF"),
  
  #             dR   dF
  Estimate = c(.02, .12,  # E1
               .20, .44,  # E2
               .40, .88   # E3
  ),
  
  
  #       dR   dF
  SE = c(.02, .05,  # E1
         .02, .05,  # E2
         .02, .10   # E3
  )
  
)

df_source <- df_source %>% mutate(
  CI_lower = Estimate - SE * 1.96,
  CI_upper = Estimate + SE * 1.96,
  )

# Cap CI at 1.00
df_source[6, "CI_upper"] = 0.999  # instead of 1 so that ggplot
                                  # does not cut off half of
                                  # the horizontal end of
                                  # the error bar

# Reorder factors
df_source$Experiment <- factor(df_source$Experiment, levels = c(1, 2, 3))
df_source$Parameter <- factor(df_source$Parameter, levels = c("dR", "dF"))





### create_plot() ----------------------------------------------------------------

create_plot <- function(
    df, Experiment, Estimate, Parameter, CI_lower, CI_upper
  ) {
  
  ggplot(
    df, aes(x = Experiment, y = Estimate, fill = Parameter)
  ) +
    
  geom_col(
    width = .65,
    position = position_dodge(.8),
    colour = "black",
  ) +
  
  geom_errorbar(
    position = position_dodge(.8),
    width = .175,
    aes(
      ymin = CI_lower,
      ymax = CI_upper
    )
  ) +
    
  coord_cartesian(ylim = c(0, 1)) +
    
  scale_y_continuous(
    breaks = seq(from = 0, to = 1, by = .20),
    labels = c("0", ".20", ".40", ".60", ".80", "1"),
    expand = c(0, 0)  # Remove gap between x bars and x-axis
  ) +
    
  theme(
    
    plot.title = element_text(
      face = "bold", size = 14, hjust = 0.5,
      margin = margin(b = 15)
    ),
    
    # Remove grid
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    
    # Add axes
    axis.line = element_line(color = 'black'),
    
    # Y ticks inwards
    axis.ticks.length.y = unit(-0.20, "cm"),
    
    # Remove x ticks
    axis.ticks.x = element_blank(),
      
    # Labels
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x = element_text(
      size = 11, color = "black",
      margin = margin(t = 6, b = 5)  # distance from axis and axis title
    ),
      
    axis.title.y = element_text(size = 13, margin = margin(r = 9)),  # distance from axis
    axis.text.y = element_text(size = 11, color = "black"),
      
    legend.position = c(1.065, 0.93),
    legend.text = element_text(size = 12),
    legend.spacing.y = unit(0.1, 'cm'),  # note also guides() below
      
  ) + # End theme
    
  guides(
    fill = guide_legend(
      override.aes = list(size = 0.9),
      byrow = TRUE
    )
  )

}  # End create_plot()

### Item and source plots --------------------------------------------------------

item_plot <- create_plot(
  df_item,
  df_item$Experiment,
  df_item$Estimate,
  df_item$Parameter,
  df_item$CI_lower,
  df_item$CI_upper
  ) +
  scale_fill_manual(
    labels = c(
      expression(italic(D)[R]),
      expression(italic(D)[F])
    ),
    values = c("white","grey60")
  ) +
  labs(
    fill = element_blank(),  # Remove legend header
    title = "Item Memory"
  ) +
  theme(  # c(top, right, bottom, left)
    plot.margin = unit(c(0.75, 1.25, 0.75, 0.80), "cm")  # unit(c(0.75, 0.85, 0.75, 0.75), "cm")
  )
# End item_plot

source_plot <- create_plot(
  df_source,
  df_source$Experiment,
  df_source$Estimate,
  df_source$Parameter,
  df_source$CI_lower,
  df_source$CI_upper
  ) +
  scale_fill_manual(
    labels = c(
      expression(italic(d)[R]),
      expression(italic(d)[F])
    ),
    values = c("white","grey60")
  ) +
  labs(
    fill = element_blank(),  # Remove legend header
    title = "Source Memory"
  ) +
  theme(  # # c(top, right, bottom, left)
    plot.margin = unit(c(0.75, 1.80, 0.75, 0.50), "cm"),
    legend.text = element_text(size = 14)
  )
# End source plot



### Combine to one figure --------------------------------------------------------

# item_plot
# source_plot

ggarrange(item_plot, source_plot, nrow = 1, ncol = 2) # +
  # theme(plot.margin = unit(rep(1, 4), "cm"))

# Now simply re-size in plot window of R studio
# and then save.
# Chose width 916, height 473

### Appendix C: -----------------------------------------------------------------

# C: Item memory parameters ----
C_df_item <- data.frame(
  
  Condition = c(rep("intrinsic", 2), rep("extrinsic", 2)),
  Parameter = c("DR", "DF", "DR", "DF"),
  
  #             DR   DF
  Estimate = c(.80, .29,  # in
               .73, .22   # ex
  ),
  
  
  #       DR   DF
  SE = c(.01, .02,  # in
         .01, .02   # ex
  )
  
)

C_df_item <- C_df_item %>% mutate(
  CI_lower = Estimate - SE * 1.96,
  CI_upper = Estimate + SE * 1.96
)

# Reorder factors
C_df_item$Condition <- factor(C_df_item$Condition, levels = c("intrinsic", "extrinsic"))
C_df_item$Parameter <- factor(C_df_item$Parameter, levels = c("DR", "DF"))



# C: Source memory parameters did not differ from collapsed estimates, no plot! ----



# C: Item plot ----

C_item_plot <- create_plot(
  C_df_item,
  C_df_item$Condition,
  C_df_item$Estimate,
  C_df_item$Parameter,
  C_df_item$CI_lower,
  C_df_item$CI_upper
) +
  scale_fill_manual(
    labels = c(
      expression(italic(D)[R]),
      expression(italic(D)[F])
    ),
    values = c("white","grey60")
  ) +
  labs(
    fill = element_blank(),  # Remove legend header
    title = "Item Memory"
  ) +
  theme(  # c(top, right, bottom, left)
    legend.position = c(1.01, 0.93),
    plot.margin = unit(c(0.75, 1.25, 0.75, 0.80), "cm")  # unit(c(0.75, 0.85, 0.75, 0.75), "cm")
  ) +
  labs(x = "Presentation Condition")
# End item_plot


# C: Export item plot ----

C_item_plot


# Now simply re-size in plot window of R studio
# and then save.
# Chose width 916/2 = 458, height 473


