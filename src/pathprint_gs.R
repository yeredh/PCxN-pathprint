# 06/18/2016
#
# Retrieve and save the gene sets as an RDS object

# retrieve pathprint gene sets
library(pathprint)
data("pathprint.Hs.gs")

# save gene sets from pathprint v.1.2.3
saveRDS(pathprint.Hs.gs,"data/pathprint_v1_2_3.RDS")
