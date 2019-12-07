#to filter out drugs with non human readable names 
library(stringr)

#counts the occurence of a character 
getCount <- function(name, charList){
  require(stringr)
  count <- 0
  for(i in 1 : length(charList)){
    count <- count + str_count(name, charList[i])
  }
  return(count)
}

#function to filter only drug names with specified character occurence
#dataFrame - input data
#colIndex - which column in the data
#threshold - no. of occurence of the characters

getFilteredDrugList <- function(dataFrame, colIndex, threshold){
  #charList <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]|[(]|[)]"
  targetDataFrame <- data.frame(matrix(ncol = ncol(dataFrame), nrow = 0))
  colnames(targetDataFrame) <- colnames(dataFrame)
  charList <- c("[-]","[\\[]", "[\\]]", "[{]", "[}]", "[(]", "[)]","[~]", "[:]", "[']")
  targetRow <- 1
  for(row in 1 : nrow(dataFrame)){
     dName <- dataFrame[row, colIndex]
     count <- getCount(dName, charList)
     if(count <= threshold){
       targetDataFrame[targetRow, ] <- dataFrame[row,]
       targetRow <- targetRow + 1
     }
  }

return(targetDataFrame)
}

