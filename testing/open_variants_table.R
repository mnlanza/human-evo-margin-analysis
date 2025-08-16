v <- read.delim("input/variants.txt", sep="\t", check.names=FALSE)
write.csv(v, "testing/variants.csv", row.names=FALSE)
try(system2("open", shQuote(normalizePath("testing/variants.csv")), wait=FALSE), silent=TRUE)
