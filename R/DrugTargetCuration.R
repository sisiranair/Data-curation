root_path <- "/Users/sisira/"
badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]|[(]|[)]"
#drugs to annotate with FDA approval status

fdaList <- read.csv(paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/FDAdrugNames_final.csv", sep = ""), stringsAsFactors = FALSE)
fdaList1 <- read.csv(paste(root_path,"Google Drive/BHKLab/Drug_revisit/FDAStatus/FDAdrugNames.csv", sep = ""), stringsAsFactors = FALSE)

#nci60filtered <- read.csv(paste(root_path,"Google Drive/BHKLab/Drug_revisit/NCI60_filtered.csv", sep = ""),stringsAsFactors = FALSE)
#length(intersect(fdaList$drugNames, nci60filtered$unique.drugid))
#setdiff(nci60filtered$unique.drugid,fdaList$drugNames)
#which(duplicated(fdaList$drugNames))
#extracting drugbank data (version 5.1.2, 20th dec 18)

drugbankxml <- read.delim(paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/Heewon/Drugbank/drugbank_5.1.2.txt", sep = ""),sep="\t",stringsAsFactors = FALSE)

drugbankxml$cleannames <- toupper(gsub(badchars,"",drugbankxml$DrugName))

fdaList$cleannames <- toupper(gsub(badchars, "", fdaList$drugNames))

combinedlist.drugbank <- drugbankxml[drugbankxml$cleannames %in% fdaList$cleannames,,drop = FALSE]

#filter out unique drug and target combo  
drug.bank.targets <- combinedlist.drugbank[!duplicated(combinedlist.drugbank[c("DrugName","GeneSymbol")]),,drop = FALSE] 

#subsetDrugbank[which(subsetDrugbank$TARGET_NAME == ""),] - to check for blank spaces
subsetDrugbank <- subset(x = drug.bank.targets, select = c("DrugBankID","DrugName","GeneSymbol"))

subsetDrugbank$DATABASE <- "DrugBank"

colnames(subsetDrugbank) <- c("ID","MOLECULE_NAME", "TARGET_NAME", "DATABASE")
write.csv(subsetDrugbank, paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/Final_Drugbank.csv", sep = ""))
#subsetDrugbank <- read.csv(paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/Final_Drugbank.csv", sep = ""))
#######################################################################################################################################################################################

#extracting ChEMBL data (chembl_drugtargets-19_15_27_38.txt)

chembllist <- read.delim(paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/chembl_drugtargets-19_15_27_38.txt", sep = ""), sep = "\t",stringsAsFactors = FALSE)

#subsetchembl <- subset(x = chembllist, select = c("MOLECULE_NAME", "TARGET_NAME", "MOLECULE_TYPE", "MECHANISM_OF_ACTION","FIRST_APPROVAL","ACTION_TYPE", "TARGET_TYPE"))

subsetchembl <- subset(x = chembllist, select = c("MOLECULE_CHEMBL_ID","MOLECULE_NAME", "TARGET_NAME"))
#filter out drugs without target info
subsetchembl <- subsetchembl[!(subsetchembl$TARGET_NAME == ""),,drop = FALSE]

subsetchembl$cleannames <- toupper(gsub(badchars,"",subsetchembl$MOLECULE_NAME))

combinedlist.chembl <- subsetchembl[subsetchembl$cleannames %in% fdaList$cleannames,,drop = FALSE]

combinedlist.chembl$DATABASE <- "ChEMBL"

combinedlist.chembl <- subset(x = combinedlist.chembl, select = -cleannames)

colnames(combinedlist.chembl)[1] <- "ID"
write.csv(combinedlist.chembl, paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/Final_Chembl.csv", sep = ""))

#combinedlist.chembl <- read.csv(paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/Final_Chembl.csv", sep = ""))
#######################################################################################################################################################################################

#extracting CTRPv2 data (from PSet, downloaded on 19th March 2019)

ctrpv2target <- CTRPv2@drug

subsetCtrpv2 <- subset(x = ctrpv2target, select = c("broad_cpd_id","drugid","gene_symbol_of_protein_target"))

subsetCtrpv2 <- subsetCtrpv2[!(subsetCtrpv2$gene_symbol_of_protein_target == ""),,drop = FALSE]

subsetCtrpv2$cleannames <- toupper(gsub(badchars, "", subsetCtrpv2$drugid))

combinedlist.ctrpv2 <- subsetCtrpv2[subsetCtrpv2$cleannames %in% fdaList$cleannames,,drop = FALSE]

combinedlist.ctrpv2$DATABASE <- "CTRPv2"

combinedlist.ctrpv2 <- subset(x = combinedlist.ctrpv2, select = -cleannames)

rownames(combinedlist.ctrpv2) <- NULL

colnames(combinedlist.ctrpv2) <- c("ID","MOLECULE_NAME", "TARGET_NAME", "DATABASE")

write.csv(combinedlist.ctrpv2, paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/Final_CTRPv2.csv", sep = ""))

#combinedlist.ctrpv2 <- read.csv(paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/Final_CTRPv2.csv", sep = ""))
#######################################################################################################################################################################################

#extracting Drug target commons data
#DOENLOADED DtcDrugTargetInteractions.csv ON March 25, 2019
#targetDTC <- read.csv("/Users/sisira/Downloads/DtcDrugTargetInteractions.csv")
targetDTC <- readRDS(paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/DTC.rds", sep =""))

subsetDTC <- subset(x = targetDTC, select = c("compound_id","compound_name","gene_names"))

subsetDTC <- subsetDTC[!(subsetDTC$gene_names == ""),,drop = FALSE]

subsetDTC$cleannames <- toupper(gsub(badchars,"", subsetDTC$compound_name))

combinedlist.DTC <- subsetDTC[subsetDTC$cleannames %in% fdaList$cleannames,,drop = FALSE]

combinedlist.DTC <- combinedlist.DTC[!duplicated(combinedlist.DTC[c("compound_name","gene_names")]),,drop = FALSE] 

combinedlist.DTC$DATABASE <- "DTC"

combinedlist.DTC <- subset(x = combinedlist.DTC, select = -cleannames)

colnames(combinedlist.DTC) <- c("ID","MOLECULE_NAME", "TARGET_NAME", "DATABASE")
write.csv(combinedlist.DTC, paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/Final_DTC.csv", sep = ""))

#combinedlist.DTC <- read.csv(paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/Final_DTC.csv", sep = ""))
#No specific ID . Uses chembl, drugbank IDS along with other databases
#######################################################################################################################################################################################
#extracting clue.io data
#IDs are BROAD IDs present in ctrpv2
targetClueio <- read.delim(paste(root_path, "Google Drive/BHKLab/Drug_revisit/Drug_targets/repurposing_drugs_20180907_clue.txt", sep =""), sep = "\t", stringsAsFactors = FALSE)
#delete first 9 rows which contain descriptions, not data
targetClueioSubset <- targetClueio[-c(1:9),]

targetClueioSubset <- subset(x = targetClueioSubset, select = c("X.Source", "X.1"))

targetClueioSubset <- targetClueioSubset[!(targetClueioSubset$X.1 == ""),,drop = FALSE]

colnames(targetClueioSubset) <- c("MOLECULE_NAME", "TARGET_NAME")
#split "|"
edittarget <- strsplit(targetClueioSubset$TARGET_NAME, split="|", fixed=TRUE)
#assign the split target names to drugs
edittarget <- data.frame(MOLECULE_NAME = rep(targetClueioSubset$MOLECULE_NAME,
                                             sapply(edittarget, length)), TARGET_NAME = unlist(edittarget),stringsAsFactors = FALSE)

edittarget$DATABASE <- "CLUE.IO"

edittarget$cleannames <- toupper(gsub(badchars,"",edittarget$MOLECULE_NAME))

combinedlist.clue.io <- edittarget[edittarget$cleannames %in% fdaList$cleannames,,drop = FALSE]

combinedlist.clue.io <- subset(x = combinedlist.clue.io, select = -cleannames)

combinedlist.clue.io$ID <- "NA"

write.csv(combinedlist.clue.io, paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/Final_CLUE.csv", sep = ""))

#combinedlist.clue.io <- read.csv(paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/Final_CLUE.csv", sep = ""))
#######################################################################################################################################################################################
#Trends in the exploitation of novel drug targets [TEND] Version: 01-August-2011
readTEND <- read.csv(paste(root_path,"/Google Drive/BHKLab/Drug_revisit/Drug_targets/TEND download.csv", sep = ""), stringsAsFactors = FALSE)

subsetTEND <- subset(x = readTEND, select = c("Drug","Target.gene."))

subsetTEND$DATABASE <- "TEND"

colnames(subsetTEND) <- c("MOLECULE_NAME", "TARGET_NAME", "DATABASE")

subsetTEND$cleannames <- toupper(gsub(badchars,"",subsetTEND$MOLECULE_NAME))

combinedlist.TEND <- subsetTEND[subsetTEND$cleannames %in% fdaList$cleannames,,drop = FALSE]

combinedlist.TEND <- subset(x = combinedlist.TEND, select = -cleannames)

combinedlist.TEND$ID <- "NA"

write.csv(combinedlist.TEND, paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/Final_TEND.csv", sep = ""))

#combinedlist.TEND <- read.csv(paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/Final_TEND.csv", sep = ""))
#######################################################################################################################################################################################

#combining all drug data 
rowbindtargets <- rbind(subsetDrugbank, combinedlist.chembl, combinedlist.ctrpv2,combinedlist.DTC,combinedlist.clue.io,combinedlist.TEND)
resultfile <- rowbindtargets[!duplicated(rowbindtargets[c(2:4)]),,drop = FALSE]
#merge fdalist with resultfile
resultfile$cleannames <- toupper(gsub(badchars, "", resultfile$MOLECULE_NAME))
mergefdaresultfile <- merge(fdaList, resultfile, by.x = "cleannames", by.y = "cleannames", all.x = T)

drugtargetCompilation <- mergefdaresultfile[,-c(1,2)]
colnames(drugtargetCompilation) <- c("BHK_drugname", "FDA_APPROVAL_STATUS","ID" ,"MOLECULE_NAME", "TARGET_NAME", "DATABASE")
#drugnames changed after Anthony updated teh duplicates.
drugtargetCompilation$BHK_drugname <- gsub("Tretinoin","ATRA",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Lestaurtinib(CEP-701)","CEP-701",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Vismodegib","GDC-0449",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Velcade","Bortezomib",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("BEZ-235","NVP-BEZ235",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("serdemetan","JNJ-26854165",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("rucaparib","AG-014699",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Obatoclax","Obatoclax Mesylate",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("linsitinib","OSI-906",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("NVP-AUY922","AUY922",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("sirolimus","Rapamycin",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("NVP-TAE684","TAE684",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Sunitinib Malate","Sunitinib",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("CAY10576","KIN001-135",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Dovitinib","TKI258",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Sodium dichloroacetate","Acetic acid, 2,2-dichloro-",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("JNK Inhibitor II","pyrazolanthrone",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Homoharringtonine","omacetaxine mepesuccinate",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("AZD1775","MK 1775",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("pyrvinium","pyrvinium pamoate",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("17-DMAG","alvespimycin",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Quinacrine (Acrichine)","mepacrine",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("(-)-MK-801","dizocilpine",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("ASN 05257430","ML083",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("NVP-LDE225","erismodegib",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("AC220","quizartinib",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("I-BET-762","GSK525762A",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("CYT-387","momelotinib",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("XL-184","cabozantinib",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Flavopiridol","alvocidib",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Enzastaurin","LY317615",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("10-methoxyharmalan","6-methoxyharmalan",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Galanthamine","galantamine",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Prestwick-642","isotretinoin",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Prestwick-665","salbutamol",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Prestwick-685","clofazimine",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Prestwick-689","androsterone",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Prestwick-692","conessine",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Prestwick-984","butacaine",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("PD-0325901","PD325901",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("Ispinesib Mesylate","Ispinesib",drugtargetCompilation$BHK_drugname)
drugtargetCompilation$BHK_drugname <- gsub("GSK1904529","GSK-1904529A",drugtargetCompilation$BHK_drugname)

#drugtargetCompilation$cleannames <- toupper(gsub(badchars,"", drugtargetCompilation$MOLECULE_NAME))
#length(unique(drugtargetCompilation$cleannames))
drugtargetCompilation <- drugtargetCompilation[,c("BHK_drugname","MOLECULE_NAME","TARGET_NAME","ID","DATABASE","FDA_APPROVAL_STATUS")]
saveRDS(drugtargetCompilation, paste(root_path, "/Google Drive/BHKLab/Drug_revisit/Drug_targets/DrugtargetCompilation.rds",sep = ""))
write.csv(DrugTargetCompilation, paste(root_path, "/Google Drive/BHKLab/Drug_revisit/Drug_targets/DrugtargetCompilation.csv",sep = ""))
#removing duplicate drug-target pairs 
# resultfile2 <- resultfile
# resultfile2$cleannames <- toupper(gsub(badchars,"",resultfile$MOLECULE_NAME))
#length(unique(resultfile$cleannames))
# which(subsetTEND$MOLECULE_NAME %in% resultfile$MOLECULE_NAME)
x <- unique(pharmacodbtargets$chembl.targets$MOLECULE_NAME)
y <- unique(pharmacodbtargets$clue.io.targets$MOLECULE_NAME)
z <- unique(pharmacodbtargets$ctrp.targets$MOLECULE_NAME)
a <- unique(pharmacodbtargets$drug.bank.targets$MOLECULE_NAME)

ss <- append(x,y)
ss1 <- append(ss,z)
ss2 <- append(ss1,a)
length(ss2)


#combining
ochembl <- combinedlist.chembl
odtc <- combinedlist.DTC
oclue <- combinedlist.clue.io
octrpv2 <- combinedlist.ctrpv2
odrugbank <- subsetDrugbank
otend <- combinedlist.TEND


