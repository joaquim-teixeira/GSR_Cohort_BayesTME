library(readr)

fls=list.files ('GSR_22_Tansey_Banerjee/teixeira_code/flipped_cohort_data/',full.names = T)
patients=substr(fls,59,59)
pt='A'
vale=vector()
vp=0
vals=vector()
i=1
for (f in fls)
{
  if (patients[i]==pt)
  {
    print (f)
    df=read_tsv (f)[,-1]
    vale=c(vale, dim(df)[2]+vp)
    vp=dim(df)[2]+vp
    pt=patients[i]
  }
  else
  {
    print ('test')
    print (f)
    df=read_tsv (f)[,-1]
    vale=c(vale, dim(df)[2])
    vp=dim(df)[2]
    pt=patients[i]
  }
  i=i+1
}
vals=c(1,vale+1)[-37]
df=data.frame (vals, vale)
write.csv (df,'GSR_22_Tansey_Banerjee/teixeira_code/dim_per_slice.csv')
