args <- commandArgs(TRUE)
samout <- args[1]
spl <- args[2]

input <- read.table(samout,sep = "\t", dec = ".", header = F, stringsAsFactors = F, 
                    col.names = c("QNAME","FLAG","RNAME","POS","MAPQ","NM","MD"))
#remove nucleotides from MD tag
md <- strsplit(input$MD,"[ACGT\\^]")
md <- lapply(md, function(x){x[!x %in% ""]}) #remove empty cells
#Sum contiguous match length to ref
md.len <- unlist(lapply(md, function(x){sum(as.numeric(x))}))

#Sum NM + contigous match length to get aligned to ref length
TLEN <- input$NM + md.len

#build output dataframe

output <- data.frame(input, MDLEN=md.len, TLEN=TLEN)

write.table(output, file = spl,quote = F, col.names = T, row.names = F, sep = "\t")
