#to map genes to drug from all databases, ultimately

#filter drugs from unique list (only human readable names). One-time run.
library(readxl)

#input list of drugs
druglist <- read.csv("data/drugs_with_ids.csv", stringsAsFactors = FALSE)


source("R/filter_drugs.R")

outputDf <- getFilteredDrugList(druglist, 2, 2)

#write.csv(outputDf, "outputDf.csv")
###########################################################################################################################################################################################################
#Step 1 - run the code to pull all target genes from all four databases - drugbank, Chembl, clue, ctrpv2

#drugs (#2654)
drugNames <- outputDf

badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]|[(]|[)]"

#extracting drugbank data (version 5.1.2, 20th dec 18)
drugbankxml <- read.delim("data/drugbank_Heewon/drugbank_5.1.4_target_info.txt", sep = "\t", stringsAsFactors = FALSE)

drugbankxml$cleannames <- toupper(gsub(badchars,"",drugbankxml$DrugName))

drugNames$cleannames <- toupper(gsub(badchars, "", drugNames$unique.drugid))

commonNames <- drugbankxml[drugbankxml$cleannames %in% drugNames$cleannames,,drop = FALSE]

# #filter out unique drug and target combo
drug.bank.targets <- commonNames[!duplicated(commonNames[c("DrugName","GeneSymbol")]),,drop = FALSE]

#combinedlist.drugbank[which(combinedlist.drugbank$TARGET_NAME == ""),] #- to check for blank spaces

combinedlist.drugbank <- subset(x = drug.bank.targets, select = c("DrugBankID","DrugName","GeneSymbol"))

combinedlist.drugbank$DATABASE <- "DrugBank"

colnames(combinedlist.drugbank) <- c("ID","MOLECULE_NAME", "TARGET_NAME", "DATABASE")
#######################################################################################################################################################################################

#extracting ChEMBL data (chembl_drugtargets-19_21_23_20.txt)

chembllist <- read.delim("data/chembl_drugtargets-19_21_23_20.txt",sep = "\t", stringsAsFactors = FALSE)

subsetchembl <- subset(x = chembllist, select = c("MOLECULE_CHEMBL_ID","MOLECULE_NAME", "TARGET_NAME","TARGET_CHEMBL_ID"))

#filter out drugs without target info
subsetchembl <- subsetchembl[!(subsetchembl$TARGET_NAME == ""),,drop = FALSE]

#map target chembl target name to gene symbol
##downloaded txt file with mapping from target chembl id to uniprot gene id.

chembl_uniprot <- read.delim("data/chembl_uniprot_mapping.txt", sep = "\t", stringsAsFactors = F, skip = 1)

#map uniprot ids to gene symbols
##download xls containing mapping after inptting all uniprot ids
library(readxl)
uniprot_sym <- read_excel("data/uniprot-yourlist_M20200106E5A08BB0B2D1C45B0C7BC3B55FD265566EF5A3M-filtered-org--.xlsx")

#extract only first word in the column gene names
library(stringr)
uniprot_sym$Symbol <- word(uniprot_sym$`Gene names`, 1)

#subset necessary cols

sub_uniprot_sym <- subset(uniprot_sym, select = c("Entry", "Symbol"), drop = F)

#merge gene symbols to chembl target ids using uniprot ids

mer_gs_ch <- merge(chembl_uniprot, sub_uniprot_sym, by.x = "Uniprot", by.y = "Entry", all.x = T)

#merge the above mapping to chembl target file

mer_ch_tr <- merge(subsetchembl, mer_gs_ch, by.x = "TARGET_CHEMBL_ID", by.y = "chemblid", all.x = T)

mer_ch_tr$cleannames <- toupper(gsub(badchars,"",mer_ch_tr$MOLECULE_NAME))

for(dr in 1:nrow(mer_ch_tr)){

  if(is.na(mer_ch_tr$Symbol[dr])){
    mer_ch_tr$Symbol[dr] <- mer_ch_tr$TARGET_NAME[dr]
  }
}


combinedlist.chembl <- mer_ch_tr[mer_ch_tr$cleannames %in% drugNames$cleannames,,drop = FALSE]

combinedlist.chembl <- subset(combinedlist.chembl, select = c("MOLECULE_CHEMBL_ID", "MOLECULE_NAME", "Symbol"), drop = F)
combinedlist.chembl$DATABASE <- "ChEMBL"

colnames(combinedlist.chembl) <- c("ID", "MOLECULE_NAME", "TARGET_NAME",  "DATABASE")
#write.csv(combinedlist.chembl, paste(root_path,"Google Drive/BHKLab/Drug_revisit/Drug_targets/Final_Chembl.csv", sep = ""))

#######################################################################################################################################################################################

#extracting CTRPv2 data ( ftp://anonymous:guest@caftpd.nci.nih.gov/pub/OCG-DCC/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/)

ctrpv2target <- read.delim("data/v20.meta.per_compound.txt",sep = "\t", stringsAsFactors = FALSE)

subsetCtrpv2 <- subset(x = ctrpv2target, select = c("broad_cpd_id","cpd_name","gene_symbol_of_protein_target"))

subsetCtrpv2 <- subsetCtrpv2[!(subsetCtrpv2$gene_symbol_of_protein_target == ""),,drop = FALSE]

subsetCtrpv2$cleannames <- toupper(gsub(badchars, "", subsetCtrpv2$cpd_name))

combinedlist.ctrpv2 <- subsetCtrpv2[subsetCtrpv2$cleannames %in% drugNames$cleannames,,drop = FALSE]

combinedlist.ctrpv2 <- subset(x = combinedlist.ctrpv2, select = -cleannames)

rownames(combinedlist.ctrpv2) <- NULL

#split ";"
edittarget_ctrpv2 <- strsplit(combinedlist.ctrpv2$gene_symbol_of_protein_target, split=";", fixed=TRUE)
#assign the split target names to drugs
edittarget_ctrpv2 <- data.frame(cpd_name = rep(combinedlist.ctrpv2$cpd_name,
                                                    sapply(edittarget_ctrpv2, length)), gene_symbol_of_protein_target = unlist(edittarget_ctrpv2),stringsAsFactors = FALSE)

pre_combinedlist.ctrpv2 <- merge(edittarget_ctrpv2, combinedlist.ctrpv2, by = "cpd_name", all.x = T)

pre_combinedlist.ctrpv2 <- pre_combinedlist.ctrpv2[,-4]

pre_combinedlist.ctrpv2$DATABASE <- "CTRPv2"


colnames(pre_combinedlist.ctrpv2) <- c("MOLECULE_NAME", "TARGET_NAME", "ID","DATABASE")

combinedlist.ctrpv2 <- pre_combinedlist.ctrpv2[c(3,1,2,4)]

#######################################################################################################################################################################################
#extracting clue.io data
#IDs are BROAD IDs present in ctrpv2
targetClueio <- read.delim("data/repurposing_drugs_20180907.txt", sep = "\t", stringsAsFactors = FALSE)
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

combinedlist.clue.io <- edittarget[edittarget$cleannames %in% drugNames$cleannames,,drop = FALSE]

combinedlist.clue.io <- subset(x = combinedlist.clue.io, select = -cleannames)

combinedlist.clue.io$ID <- "NA"

#######################################################################################################################################################################################

#combining all drug data
rowbindtargets <- rbind(combinedlist.drugbank, combinedlist.chembl, combinedlist.ctrpv2,combinedlist.clue.io)

#merge with BHK names

rowbindtargets$cleannames <- toupper(gsub(badchars, "", rowbindtargets$MOLECULE_NAME))

drugNames$cleannames <- toupper(gsub(badchars, "", drugNames$unique.drugid))

mergerowbindtargets <- merge(drugNames, rowbindtargets, by.x = "cleannames", by.y = "cleannames", all.x = T)

#######################################################################################################################################################################################
