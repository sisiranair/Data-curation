#to retrieve Pubchem names
library("httr")
library("stringi")
library("XML")
getPubChemIUPAC <- function(drugList){
  require("httr")
  require("stringi")
  require("XML")
  baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/DRUG_NAME/property/IUPACName/XML"
  baseMatrix <- matrix("", nrow = length(drugList), ncol = 2)
  
  for(i in 1:length(drugList)){
    if(i %% 5 == 0){
      #Wait a second before next execution as PubChem restricts the requests at 5 requests per second
      Sys.sleep(1)
    }
    currentDrugName <- drugList[i]
    requestURL <- gsub("DRUG_NAME",currentDrugName, baseURL)
    #fetching drug data from PubChem using HTTP GET request
    requestURL <- URLencode(requestURL)
    
    response <- GET(requestURL)
    
    responseText <- content(response,"text")
    
    responseXML <- XML::xmlParse(responseText)
    
    pubChemList <- XML::xmlToList(responseXML)
    
    
    baseMatrix[i, 1] <- drugList[i]
    
    if(!is.null(pubChemList$Properties$IUPACName[1])){
      baseMatrix[i, 2] <- pubChemList$Properties$IUPACName[1]
      
    }else{
      print(paste(drugList[i]," ","fetched null values"))
    }
    
  }
  
  drugDataFrame <- as.data.frame(baseMatrix)
  colnames(drugDataFrame) <- c("Drug", "IUPACNAme")
  return(drugDataFrame)
  
  
}


