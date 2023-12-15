# TemporalOriginNestedness

Please cite: 
Staniczenko, P.P.A. & Panja, D. (2023). Temporal origin of nestedness in interaction networks. PNAS Nexus, 2, pgad412

Author: Phillip P.A. Staniczenko
pstaniczenko@brooklyn.cuny.edu
https://staniczenkoresearch.net

version 1 -- 2023-09-21

Code file is TemporalOriginNestedness.R

Example data in folder ExampleData.

DATA DESCRIPTION

ExampleData folder contains daily interaction matrices and metadata files listing how they are organized into weeks, months, and years, adapted from:

CaraDonna PJ. Temporal variation in plant-pollinator interactions, Rocky Mountain Biological Laboratory, CO, USA, 2013–2015, ver 1. Environmental data initiative, 2020 [accessed October 2022]. https://doi.org/10.6073/pasta/27dc02fe1655e3896f20326fed5cb95f

There are 106 daily interaction matrices. Rows represent plants (46 species) and columns represent pollinators (93 species); entries represent the number of interactions (i.e., visits) recorded between each pair of plants and pollinators.

days_in_weeks.csv lists the week number (first column; e.g., week 1) and the number of daily interaction matrices associated with that week (second column; e.g., 4 daily interaction matrices are associated with week 1). The 106 daily interaction matrices are ordered by time, so the matrices 1, 2, 3, and 4 are associated with week 1, matrices 5, 6, 7, 8, and 9 with week 2, and so on.

The same format is used for days_in_months.csv and days_in_years, but the first column lists month and year number, respectively, rather than week number.

CODE DESCRIPTION

INPUT: daily interaction matrices (e.g., daily_adjacency_1.csv) and metadata file of matrices partitioned by week/month/year (e.g., days_in_weeks.csv)

Number of synthetic matrices for nestedness analysis can be changed on L285.
Parameters set by data are hardcoded on L288–L209.
Level of temporal aggregation (weeks, months, years) for analysis can be changed on L291; if a different data set is used, the total number of weeks/months/years must be changed on L309, L312, L315, respectively.

OUTPUT: Results.csv lists mean (row 3) and standard deviation (row 4) of six metrics (nestedness, sensitivity, specificity, F-score, Informedness, phi coefficient; row 2) for three models (random graph, degree distribution, phenology; row 1).