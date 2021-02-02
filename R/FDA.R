fdaFileLink <- "https://www.fda.gov/media/76860/download"
#download and unzip location
downloadLocation <- paste("~/","data/",sep="")
downloadedZipFile <- paste(downloadLocation,"/fdazip.zip", sep = "")
download.file(fdaFileLink,destfile = downloadedZipFile)

unzip(downloadedZipFile,exdir = downloadLocation)
#downloading and reading the products file 
productsFileLocation <- paste(downloadLocation,"/products.txt", sep = "")

productsFile <- read.delim(file = productsFileLocation,sep = "~", stringsAsFactors = FALSE)
#reading the drug names input
#drugNames <- as.vector(unlist(as.list(read.csv(paste(root_path,"Google Drive/BHKLab/Drug_revisit/FDAStatus/drugNames.csv", sep ="")))))

drugNames <- read.delim("data/compounds.txt", sep = "\t", stringsAsFactors = F)
colnames(drugNames) <- "drugNames"
#adding two new columns
drugNames$FDA <- FALSE
drugNames$ActiveIngredient <- ""

badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]|[(]|[)]"

productsFile$modDrugname <- toupper(gsub(badchars, "", productsFile$Trade_Name))
productsFile$modActiveIngre <- toupper(gsub(badchars, "", productsFile$Ingredient))

for(i in 1:nrow(drugNames)){
  newdrugname <- gsub(pattern = badchars,replacement = "",x = drugNames$drugNames[i])
  newdrugname <- toupper(newdrugname)
  #print(newdrugname)
  #x <- productsFile[productsFile$DrugName == newdrugname,]
  #productdf <- productsFile[newdrugname == productsFile$modDrugname | newdrugname == productsFile$modActiveIngre,]
  
  productdf <- productsFile[grepl(pattern = newdrugname,x = productsFile$modDrugname)
                            | grepl(pattern = newdrugname, x=productsFile$modActiveIngre),]
  #if there are more tahn 1 row, means that there is a match.
  if(nrow(productdf) > 0){
    drugNames$FDA[i] = TRUE
    drugNames$ActiveIngredient[i] <- as.character(productdf$Ingredient[1])
    
  }
}


write.csv(drugNames, "data/FDA_output.csv")
#THE LIST WAS MANUALLY CHECKED AND FILTERED FOR VALID FDA STATUS

