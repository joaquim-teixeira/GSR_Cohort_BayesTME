libra


get_nc_list=function (ncf)
{
  nc=nc_open (ncf)
  x=list()
  for (n in names (nc$var))
  {
    x[[n]]=ncvar_get(nc,n)
  }
  return (x)
  nc_close (nc)
}



analysis_slice=function (f_coords, f_btme, f_lda)
{
  dat=read.table(f_coords)
  cds=get_coords_cdata(dat[1,-1])
  lda=get_nc_list(f_lda)
  btme=get_nc_list(f_btme)
  lda$gene_probs=lda$components_/rowSums(lda$components_)
  #match=match_clusters(cds, apply  (t(btme$cell_num_trace)[,-1],1,which.max), apply  (t(lda$transform),1,which.max))
  probs_lda=data.frame (cds, t(lda$transform))
  names (probs_lda)=c('x','y',paste ('cluster_',1:6,sep=''))
  probs_btme=data.frame (cds, t(btme$cell_prob_trace[-1,]))
  names (probs_btme)=c('x','y',paste ('cluster_',1:6,sep=''))
  probs_lda=pivot_longer(probs_lda,cols=starts_with('cluster'))
  probs_btme=pivot_longer(probs_btme,cols=starts_with('cluster'))
  probs_lda$type='LDA'
  probs_btme$type='BayesTME'
  probs=rbind (probs_lda, probs_btme)
  g1= ggplot (probs, aes (x,y, fill=value))+geom_tile(color='black')+theme_bw()+facet_grid(name~type)+scale_fill_gradientn(colors=jet.colors (11),breaks =seq(0,1,length.out=11))+labs (x='',y='',fill='Spot Probability per Cell Type')  
  dfl=list()
  for (i in 1:6)
  {
    dfl[[i]]=data.frame (lda$gene_probs[btme$marker_genes[,i]+1,],cell_type=i,gene=lda$gene_names[btme$marker_genes[,i]+1])
  }
  df=rbindlist(dfl)
  names (df)=c(1:6, 'cell_type','gene')
  df=pivot_longer(df, cols=1:6)
  df$name=as.numeric(as.character(df$name))
  df$cell_type=paste('BayesTME_cell_type_',df$cell_type,sep='')
  g2=ggplot (df, aes (name,value, color=gene))+geom_line()+facet_wrap(~cell_type)+theme_bw()+theme(legend.position = 'bottom')+scale_x_continuous(breaks=seq(1,6,1))+labs (x='LDA Cell Type',y='LDA Gene Expression Probability\n for BayesTME Marker Genes',color='BayesTME Marker Genes')
  return (list(g1,g2,lda,btme))
}

analysis_lda_patient=function (f_coords, f_btme, f_lda_slice,f_lda_patient, f_clust,indvals)
{
  dat=read.table(f_coords)
  cds=get_coords_cdata(dat[1,-1])
  lda=get_nc_list(f_lda_slice)
  btme=get_nc_list(f_btme)
  lda$gene_probs=lda$components_/rowSums(lda$components_)
  probs_lda=data.frame (cds, t(lda$transform))
  names (probs_lda)=c('x','y',paste ('cluster_',1:6,sep=''))
  probs_btme=data.frame (cds, t(btme$cell_prob_trace[-1,]))
  names (probs_btme)=c('x','y',paste ('cluster_',1:6,sep=''))
  probs_lda=pivot_longer(probs_lda,cols=starts_with('cluster'))
  probs_btme=pivot_longer(probs_btme,cols=starts_with('cluster'))
  probs_lda$type='2-LDA_slice'
  probs_btme$type='1-BayesTME'
  hcldf=read.table (f_clust)
  probs_hclust=data.frame (x=hcldf$x,y=hcldf$y, name=paste ('cluster_',hcldf$label+1,sep=''))  
  probs_hclust$value=1
  probs_hclust$type='4-Hierarchical'
  lda_pat=get_nc_list(f_lda_patient)
  probs_lda_patient=data.frame (cds, t(lda_pat$transform)[(indvals[1]:indvals[2]),])
  names (probs_lda_patient)=c('x','y',paste ('cluster_',1:6,sep=''))
  probs_lda_patient=pivot_longer(probs_lda_patient,cols=starts_with('cluster'))
  probs_lda_patient$type='3-LDA_patient'
  probs=rbind (probs_lda, probs_btme,probs_hclust,probs_lda_patient)
  
  
  g1= ggplot (probs, aes (x,y, fill=value))+geom_tile(color='black')+theme_bw()+facet_grid(name~type)+scale_fill_gradientn(colors=jet.colors (11),breaks =seq(0,1,length.out=11))+labs (x='',y='',fill='Spot Probability by\n Cell Type')  
  dfl=list()
  for (i in 1:6)
  {
    dfl[[i]]=data.frame (lda$gene_probs[btme$marker_genes[,i]+1,],cell_type=i,gene=lda$gene_names[btme$marker_genes[,i]+1])
  }
  df=rbindlist(dfl)
  names (df)=c(1:6, 'cell_type','gene')
  df=pivot_longer(df, cols=1:6)
  df$name=as.numeric(as.character(df$name))
  df$cell_type=paste('BayesTME_cell_type_',df$cell_type,sep='')
  g2=ggplot (df, aes (name,value, color=gene))+geom_line()+facet_wrap(~cell_type)+theme_bw()+theme(legend.position = 'bottom')+scale_x_continuous(breaks=seq(1,6,1))+labs (x='LDA Cell Type',y='LDA Gene Expression Probability\n for BayesTME Marker Genes',color='BayesTME Marker Genes')
  return (list(g1,g2,lda,btme))
}



flc=list.files ('GSR_22_Tansey_Banerjee/teixeira_code/flipped_cohort_data',full.names = T)
flb=list.files('GSR_22_Tansey_Banerjee/teixeira_code/deconvolution_results_1000_sims',full.names = T)
fll=list.files('GSR_22_Tansey_Banerjee/teixeira_code/ldas_by_slice',full.names = T)
fllp=list.files('GSR_22_Tansey_Banerjee/teixeira_code/ldas_by_patient',full.names = T)
fllp=c(rep(fllp[1],6),rep(fllp[2],6),rep(fllp[3],6),rep(fllp[4],6),rep(fllp[5],3),rep(fllp[6],3),rep(fllp[7],3))

flh=list.files ('GSR_22_Tansey_Banerjee/her2st-master/res/ST-cluster/lbl/',full.names = T)[seq(1,72,by=2)]
inds=read.csv ('GSR_22_Tansey_Banerjee/teixeira_code/dim_per_slice.csv')
inds=as.matrix(inds[,-1])
vals=list()
for (i in 1:5)
{
  vals[[i]]=analysis_lda_patient(flc[i],flb[i],fll[i],fllp[i],flh[i],inds[i,])
}

