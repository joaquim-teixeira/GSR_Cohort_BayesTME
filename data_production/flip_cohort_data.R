library(readr)
ring.genes <- c("FTH1", "EEF2", "BEST1", "LRRC59", "PRDX1", "CD63", "DYNC1H1", "ENO1",
                "PSMB3", "RNF187", "RNASE1", "CFL1", "GRN", "UBC", "TAX1BP3", "COX4I1",
                "CUTA", "NME1", "H3F3B", "AKR7A2", "IMPDH2")

flip_tsv=function (fls,flout)
{
  tst=read.table (fls)
  gn=colnames(tst)
  ln=rownames(tst)
  vals=as.matrix (tst)
  df=data.frame (cbind (gn, t(vals)))
  tdf=df[setdiff(rownames(df),ring.genes), ]
  ndf=c('gene',ln)
  names (tdf)=ndf
  write_tsv(tdf,flout)

}
fls=list.files ('../her2st-master/data/ST-cnts/',full.names = T)
fls=fls[seq(1,length(fls),by=2)]
fouts=paste('flipped_cohort_data',substr(fls,32,37),sep='/')
for (i in 1:length(fls))
{
  flip_tsv(fls[i],fouts[i])
  print (fls[i])
}