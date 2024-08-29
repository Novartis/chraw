

chromHMMDict <- rbind(
    read.delim("inst/files/chromHMM_18states_hg38.txt", check.names=FALSE)[,c("chromHmmLabs", "simpleLabs")],
    read.delim("inst/files/chromHMM_18states_mm10.txt", check.names=FALSE)[,c("chromHmmLabs", "simpleLabs")])

chromHMMDict <- unique(chromHMMDict)

save(chromHMMDict, file="data/chromHMMDict.RData")
