#Plot word cloud of BHK's papers

#Step 1 - Retrieve data from PubMed
install.packages("easyPubMed")
library(easyPubMed)

# Cast PubMed record info of BHK into a data.frame

bhk_query <- "Benjamin Haibe-Kains[AU]"
bhk_on_pubmed <- get_pubmed_ids(bhk_query)
bhk_abstracts_xml <- fetch_pubmed_data(bhk_on_pubmed, encoding = "ASCII")
output <- table_articles_byAuth(pubmed_data = bhk_abstracts_xml,
                            included_authors = "last",
                            max_chars = 100,
                            autofill = TRUE)
#list can be further filtered using last name

#Step 2 - create a word cloud of journal names
install.packages("wordcloud")
library(wordcloud)
install.packages("RColorBrewer")
library(RColorBrewer)
install.packages("wordcloud2")
library("wordcloud2")

#Create a vector containing only the journal abbreviation
jabbr <- output$jabbrv
df.jabbr <- as.data.frame(table(jabbr))
df.jabbr <- df.jabbr[order(df.jabbr$Freq, decreasing = T),]

pdf("bhk.journal.pdf")
set.seed(234) #change for a different layout of words

# more color palettes for brewer.pal here - https://www.datanovia.com/en/blog/the-a-z-of-rcolorbrewer-palette/
wordcloud(words = df.jabbr$jabbr, freq = df.jabbr$Freq, min.freq = 1, random.order=FALSE, 
                              rot.per=0.35,colors=brewer.pal(8, "Dark2")) 
dev.off()
