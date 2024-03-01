library(sleuth)
library(data.table)
library(dplyr)

# Read in the sample table
stab <- read.table("cov.txt", header=TRUE, stringsAsFactors = FALSE, sep='\t')

# Create a Sleuth object
so <- sleuth_prep(stab)

# Fit a model comparing conditions
so <- sleuth_fit(so, ~condition, 'full')

# Fit a reduced model
so <- sleuth_fit(so, ~1, 'reduced')

# Likelihood ratio test for differential expression
so <- sleuth_lrt(so, 'reduced','full')

# Extract results
sleuth_table <- sleuth_results(so, 'reduced:full','lrt',show_all=FALSE)

# Filter significant results
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

# Select relevant columns
sig_sleuth <- sleuth_significant %>% select(target_id)

# Write target IDs to a file
write.table(sig_sleuth, file="significant_transcripts.txt", quote=FALSE, row.names=FALSE)