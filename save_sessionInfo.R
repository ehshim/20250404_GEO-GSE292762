# save_sessionInfo.R
# Save current R session information to a text file

sink(file = "./250404_sessionInfo.txt")
sessionInfo()
sink()