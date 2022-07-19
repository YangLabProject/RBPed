##for calculating fisher test pvalue and delta editing levels


#data is the editing levels for specific editing site from eCLIP and RNA-seq dataset


Pvalue=list()
i=1
while(i<=dim(data)[1]){
data2=data[i,]
#V7 is number of base I in eclip, V8 is total base covered this site in ECLIP, V11 is number of base I in RNA-seq ,V10 is total base covered this site in RNA-seq.
Pvalue[i]=fisher.test(matrix(c(data2$V7,data2$V11-data2$V7,data2$V8,data2$V10-data2$V8),nrow=2,ncol=2))$p.value
i=i+1 
}

delta=data$V7/data$V8-data$V11/data$V10
m2=cbind(data,unlist(Pvalue),delta,eclipEL=data$V7/data$V8,RNAEL=data$V11/data$V10)

#m2 is the results for potential binding preference.