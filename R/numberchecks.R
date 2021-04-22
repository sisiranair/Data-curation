library(PharmacoGx)
library(SummarizedExperiment)

gcsi <- readRDS("gCSI.rds")
gcsi.cell <- as.data.frame(gcsi@cell)
gcsi.cell.nsclc <- subset(gcsi.cell, gcsi.cell$Cellosaurus.Disease.Type == "Non-small cell lung carcinoma", drop = F)
gcsi.dr.no <- subset(gcsi@sensitivity$info, gcsi@sensitivity$info$cellid %in% gcsi.cell.nsclc$cellid, drop = F)
length(unique(gcsi.dr.no$drugid))

ccle <- readRDS("CCLE.rds")
ccle.cell <- as.data.frame(ccle@cell)
ccle.cell.nsclc <- subset(ccle.cell, ccle.cell$Cellosaurus.Disease.Type == "Non-small cell lung carcinoma", drop = F)
ccle.dr.no <- subset(ccle@sensitivity$info, ccle@sensitivity$info$cellid %in% ccle.cell.nsclc$cellid, drop = F)
length(unique(ccle.dr.no$drugid))

ctrpv2 <- readRDS("CTRPv2.rds")
ctrpv2.cell <- as.data.frame(ctrpv2@cell)
ctrpv2.cell.nsclc <- subset(ctrpv2.cell, ctrpv2.cell$Cellosaurus.Disease.Type == "Non-small cell lung carcinoma", drop = F)
ctrpv2.dr.no <- subset(ctrpv2@sensitivity$info, ctrpv2@sensitivity$info$cellid %in% ctrpv2.cell.nsclc$cellid, drop = F)
length(unique(ctrpv2.dr.no$drugid))

fimm <- readRDS("FIMM.rds")
fimm.cell <- as.data.frame(fimm@cell)
fimm.cell.nsclc <- subset(fimm.cell, fimm.cell$tissueid == "lung", drop = F)#no lung tissue

gdscv2 <- readRDS("GDSC2.rds")
gdscv2.cell <- as.data.frame(gdscv2@cell)
gdscv2.cell.nsclc <- subset(gdscv2.cell, gdscv2.cell$Cellosaurus.Disease.Type == "Non-small cell lung carcinoma", drop = F)
gdscv2.dr.no <- subset(gdscv2@sensitivity$info, gdscv2@sensitivity$info$cellid %in% gdscv2.cell.nsclc$cellid, drop = F)
length(unique(gdscv2.dr.no$drugid))

gdscv1 <- readRDS("GDSC1.rds")
gdscv1.cell <- as.data.frame(gdscv1@cell)
gdscv1.cell.nsclc <- subset(gdscv1.cell, gdscv1.cell$Cellosaurus.Disease.Type == "Non-small cell lung carcinoma", drop = F)
gdscv1.dr.no <- subset(gdscv1@sensitivity$info, gdscv1@sensitivity$info$cellid %in% gdscv1.cell.nsclc$cellid, drop = F)
length(unique(gdscv1.dr.no$drugid))
