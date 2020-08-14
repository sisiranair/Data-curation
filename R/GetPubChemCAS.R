#to retrieve Pubchem names
library("httr")
library("stringi")
library("XML")
getPubChemCAS <- function(CIDlist){
  require("httr")
  require("stringi")
  require("XML")
  baseURL <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/CID/XML"
  
  baseMatrix <- matrix("", nrow = length(CIDlist), ncol = 4)
  
  for(i in 1:length(CIDlist)){
    if(i %% 5 == 0){
      #Wait a second before next execution as PubChem restricts the requests at 5 requests per second
      Sys.sleep(1)
    }
    currentCID <- CIDlist[i]
    #print(paste("cid found ", currentCID))
    requestURL <- gsub("CID",currentCID, baseURL)
    #fetching drug data from PubChem using HTTP GET request
    requestURL <- URLencode(requestURL)
    
    response <- GET(requestURL)
    
    responseText <- content(response,"text")
    
    if(is.na(responseText)){
      baseMatrix[i, 1] <- CIDlist[i]
      baseMatrix[i, 2] <- NA
      next
    }
    
    responseXML <- XML::xmlParse(responseText)
    #----xPath code
    nsDefs <- xmlNamespaceDefinitions(responseXML)
    ns <- structure(sapply(nsDefs, function(x) x$uri), names = names(nsDefs))
    names(ns)[1] <- "x"
    casnm <- xpathSApply(responseXML, "//x:Section/x:TOCHeading[text()='CAS']/parent::x:Section/x:Information/x:Value/x:StringWithMarkup/x:String", namespaces = ns, fun = xmlValue)
    #--xPath code
    
    
    baseMatrix[i, 1] <- CIDlist[i]
    
    if(!is.null(casnm)){
      baseMatrix[i, 2] <- casnm[1]
      
    }else{
      print(paste(CIDlist[i]," ","fetched null values"))
   
    }
  }
  
  drugDataFrame <- as.data.frame(baseMatrix)
  colnames(drugDataFrame) <- c("cid", "cas")
  return(drugDataFrame)
  
}