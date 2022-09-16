#' To compile and re-arrange Tier 1 and 2 genes from Cancer Gene Census 
#' Remote file download using basic authentication
#' Data accessed from @COSMIC 
#' Get credentials by registering on COSMIC.
#' Follow instructions - Cancer Gene Census section in link - https://cancer.sanger.ac.uk/cosmic/download 

#required libraries
library(RCurl)
library(httr)

#' @example 
uname = "name@uhnresearch.ca"
pwd = ""

#' @param username or @param uname - email id after registering on COSMIC
#' @param password or @param pwd - password after registering on COSMIC

get_remoteFile <- function(username, password){
  unam.pwd <- paste(username, password, sep = ":")
  benc <- RCurl::base64(unam.pwd)
  #gurl <- httr::GET("https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v96/cancer_gene_census.csv", config(add_headers(Authorization = paste("Basic ", benc, sep = ""))))
  gurl <- httr::GET("https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v96/cancer_gene_census.csv", config(add_headers(Authorization = paste("Basic ", benc, sep = ""))))
  URL <- unlist(content(gurl)$url)
  fetch.data <- RCurl::getURL(URL)
  fetch.data.out <- read.csv(textConnection(fetch.data), stringsAsFactors = F)
#write.csv(fetch.data.out, "../cgc.genes.op.csv")
  write.csv(fetch.data.out, "Desktop/icb_data_exploration_and_analysis/Data/cgc.genes.op.csv")

return(fetch.data.out)
}
#username and password input
fetch.op <- get_remoteFile(username = uname, password = pwd) #SSL certificate error, downloading manually
fetch.op <- read.csv("~/Downloads/cancer_gene_census.csv")
#clean data
ensmbl.split <- strsplit(fetch.op$Synonyms, split=",", fixed=TRUE)
ensmbl.split <- data.frame(Gene.Symbol = rep(fetch.op$Gene.Symbol,
                                             sapply(ensmbl.split, length)), Ensembl.id = unlist(ensmbl.split),stringsAsFactors = FALSE)
ensmbl.split.sub <- ensmbl.split[grep("ENSG000", ensmbl.split$Ensembl.id),]
ensmbl.split.sub$Ensembl.id.clean <- gsub("\\..*", "", ensmbl.split.sub$Ensembl.id)
#merge downloaded data with ensembl ids
mer.data <- merge(fetch.op, ensmbl.split.sub, by = "Gene.Symbol", all.x = T)
col_order <- c("Gene.Symbol","Ensembl.id","Ensembl.id.clean", "Entrez.GeneId","Name","Tier","Genome.Location","Hallmark","Chr.Band","Somatic","Germline","Tumour.Types.Somatic.","Tumour.Types.Germline.","Cancer.Syndrome","Tissue.Type","Molecular.Genetics","Role.in.Cancer","Mutation.Types","Translocation.Partner","Other.Germline.Mut","Other.Syndrome", "COSMIC.ID", "cosmic.gene.name", "Synonyms")
mer.data <- mer.data[,col_order]
#write.csv(mer.data, "../cgc.genes.csv")
write.csv(mer.data, "Desktop/icb_data_exploration_and_analysis/Data/FS_cgc_genes_v96.csv")

  
  
  
  
  
