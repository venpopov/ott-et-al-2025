Better Source Memory for Recognized To-Be-Forgtten Items than for Recognized To-Be-Remembered Items

This readme.txt explains the files inside the R-Project folder.


Each file (except for visualization.R, see below) either pertains to Experiment 1 (color), Experiment 2 (top vs. bottom),
or Experiment 3 (left vs. right).
Thus, each file either starts with "E1", "E2", or "E3".


For all experiments there is raw data.
E1: raw data is stored in an Excel file.
E2: a folder with two txt-files for each participants 1) trial-data and 2) demographics
E3: as for E2.


The first step is get the raw data ready for use in R.
For instance read/import all the individual txt-files or re-structure the dataframes a tiny bit.
This is done with the R-Scripts named "raw2rdy".


The outcome of the "raw2rdy" files are the csv files:
"E1_data.csv"
"E2_data.csv"
"E3_data.csv"
These are now ready for statistical analyses.


The corresponding R-Scripts for the analyses contain "analyses" in their file names.


For E2 and E3 there are two additional files.
"E2_strats.csv" and "E3_strats.csv" are created at the end of the analyses files.
They store the strategy descriptions that were given by participants.

To code these descriptions, files "E2_stratsCoded.xlsx" and "E2_stratsCoded.xlsx" were created
(same for E3).
The coding/classification happened in these very files.
Therefore, anyone interested can review this coding process, too.

The file visualization.R creates a ggplot of the important parameter estimates from all three experiments (see Figure_2 file). It also creates Figure_C1.