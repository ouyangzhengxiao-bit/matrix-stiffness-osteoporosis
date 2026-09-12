.libPaths(c('work/high_impact/r_libs',.libPaths()))
suppressPackageStartupMessages({library(AnnotationDbi);library(hugene20sttranscriptcluster.db)})
read_series <- function(path){
 x <- readLines(gzfile(path),warn=FALSE)
 i <- which(x=='!series_matrix_table_begin')+1
 j <- which(x=='!series_matrix_table_end')-1
 z <- read.delim(textConnection(x[i:j]),check.names=FALSE,quote='"',row.names=1)
 colnames(z) <- gsub('"','',colnames(z)); z
}
x <- read_series('work/high_impact/data/GSE156508_series_matrix.txt.gz')
mp <- AnnotationDbi::select(hugene20sttranscriptcluster.db,keys=rownames(x),keytype='PROBEID',columns='SYMBOL')
mp <- unique(mp[!is.na(mp$SYMBOL),c('PROBEID','SYMBOL')])
nmap <- table(mp$PROBEID)
mp <- mp[mp$PROBEID %in% names(nmap)[nmap==1],]
x$gene <- mp$SYMBOL[match(rownames(x),mp$PROBEID)]
x <- x[!is.na(x$gene),]
g <- aggregate(x[,setdiff(colnames(x),'gene')],list(gene=x$gene),median)
dir.create('outputs/研究升级/high_impact_v4/osteoporosis_osteoblast_validation',recursive=TRUE,showWarnings=FALSE)
write.csv(g,'outputs/研究升级/high_impact_v4/osteoporosis_osteoblast_validation/GSE156508_gene_expression.csv',row.names=FALSE)
write.csv(mp,'outputs/研究升级/high_impact_v4/osteoporosis_osteoblast_validation/GPL16686_probe_symbol_mapping.csv',row.names=FALSE)
cat(nrow(x),'mapped probes;',nrow(g),'genes\n')
