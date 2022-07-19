# RBPed
RBPed is a bioinformatics pipeline to identify differential RNA editing events between eCLIP-seq and RNA-seq

## Authors
Xiaolin Hu, Qin Zou, Li Yao and Xuerui Yang

## 1. python dependency
RBPed was compromised by samtools, REDItools V1, and R. Among those softwares, REDItools V1 depends on python2.7, besides, it depends on pysam 0.7.6. You can install with code: 
* python -m pip install pysam==0.7.6

Also, it depends on some bioinformatics packages
* samtools

## 2. reference file
the reference genome (hg19) downloaded from GENCODE website. Known editing events downloaded from REDIportal database.
Those data can be downloaded within this pipeline.

## 3. Run
* ### A. Downloading mapping files (bam) from ENCODE websites.
```
wget http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg19.txt.gz
gunzip TABLE1_hg19.txt.gz
grep -v 'Region' TABLE1_hg19.txt |awk '{print $1,$2,$5}' | sed 's/ /\t/g' | sed 's/chr//g' | sort -k1,1 -k2,2n  > REDIportals.forREDItools.txt

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
gunzip GRCh37.primary_assembly.genome.fa.gz
sed 's/chr//g' GRCh37.primary_assembly.genome.fa > hg19.ref.fa
```
* ### B. Sort, index bam files
```
samtools  sort --threads 16  -o XXX.sorted.bam XXX.eCLIP.bam
samtools index XXX.sorted.ba
```
* ### C. Detection RNA editing events with REDItools.
Notice: if this code not work out properly, chromosome name (with chr or not) should be checked in bam, reference and REDIportals data.
```
python2.7 REDItoolKnown.py  -i XXX.sorted.bam -f hg19.ref.fa -l REDIportals.forREDItools.txt -t 16 -c 10 -T 6-0  -p -e -d -u -m20  -v 0 -n 0.0 -o outputDir
```
* ### D. Merge RNA editing events in different replicate and between different eCLIP and RNA-seq
* ### E. Using R to calculate delta editing levels and fisher test P-values. core code listed:
```
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

```


## License
RBPed is licensed under the MIT license
