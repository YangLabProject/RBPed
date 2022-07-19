#
#this script is for calling editing events from eclip/RNA-seq dataset

#for reference downloading
wget http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg19.txt.gz
gunzip TABLE1_hg19.txt.gz
grep -v 'Region' TABLE1_hg19.txt |awk '{print $1,$2,$5}' | sed 's/ /\t/g' | sed 's/chr//g' | sort -k1,1 -k2,2n  > REDIportals.forREDItools.txt

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
gunzip GRCh37.primary_assembly.genome.fa.gz
sed 's/chr//g' GRCh37.primary_assembly.genome.fa > hg19.ref.fa

#sorting bam downloaded from ENCODE for eCLIP or RNA-seq
samtools  sort --threads 16  -o XXX.sorted.bam XXX.eCLIP.bam

#for REDIKnown.py detection, if this code not work out properly, chromosome name (with chr or not) should be checked in bam, reference and REDIportals data.
python2.7 REDItoolKnown.py  -i XXX.sorted.bam -f hg19.ref.fa -l REDIportals.forREDItools.txt -t 16 -c 10 -T 6-0  -p -e -d -u -m20  -v 0 -n 0.0 -o outputDir

