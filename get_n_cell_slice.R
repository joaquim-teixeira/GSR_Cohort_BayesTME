library(readr)

fls=list.files ('GSR_22_Tansey_Banerjee/her2st-master/res/ST-cluster/lbl/',full.names = T)[seq(1,72,by=2)]
val=vector()
for (f in fls)
{
  df=read_tsv (f)
  val=c(val, length(unique(df$label)))
}
write.csv (val,'GSR_22_Tansey_Banerjee/teixeira_code/n_cell_type_cluster_expression.csv')
