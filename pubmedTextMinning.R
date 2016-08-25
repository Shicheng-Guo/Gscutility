# build my research word cloud

library(RCurl)
library(XML)
library(snippets)
library(tm)
# esearch
url <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
q   <- "db=pubmed&term=Shicheng+Guo[au]&usehistory=y"
esearch <- xmlTreeParse(getURL(paste(url, q, sep="")), useInternal = T)
webenv  <- xmlValue(getNodeSet(esearch, "//WebEnv")[[1]])
key     <- xmlValue(getNodeSet(esearch, "//QueryKey")[[1]])
# efetch
url <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
q   <- "db=pubmed&retmode=xml&rettype=abstract"
efetch <- xmlTreeParse(getURL(paste(url, q, "&WebEnv=", webenv, "&query_key=", key, sep="")), useInternal = T)
abstracts <- getNodeSet(efetch, "//AbstractText")
# words
abstracts <- sapply(abstracts, function(x) { xmlValue(x) } )
words <- tolower(unlist(lapply(abstracts, function(x) strsplit(x, " "))))
# remove parentheses, comma, [semi-]colon, period, quotation marks
words <- words[-grep("[\\)\\(,;:\\.\'\"]", words)]
words <- words[-grep("^\\d+$", words)]
words <- words[!words %in% stopwords()]
wt <- table(words)
wt <- (wt[wt > 5])
cloud(wt, col = col.br(wt, fit=TRUE))
