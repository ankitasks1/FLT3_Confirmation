# FLT3 Confirmation
## Repository for "Tools and materials to confirm FLT3"

### USAGE
### check the usage and parameters
python3 ~/Documents/scripts/confirm_flt3_v6.py -h

### with ITDs detection and annotation:
## python3 ~/Documents/scripts/confirm_flt3_v6.py --lspan 0 --rspan 2 --alt_type dup --bamfile sampleid.Aligned.sortedByCoord.out.bam --config config.txt --sid sampleid --id anyid --vcftype f --annotation y --itdetect y

### without ITDs detection and annotation
## python3 ~/Documents/scripts/confirm_flt3_v6.py --lspan 0 --rspan 2 --alt_type dup --bamfile sampleid.Aligned.sortedByCoord.out.bam --config config.txt --sid sampleid --id anyid --vcftype f --annotation n --itdetect n

### with ITDs detection but No annotation
```bash
python3 ~/Documents/scripts/confirm_flt3_v6.py --lspan 0 --rspan 2 --alt_type dup --bamfile sampleid.Aligned.sortedByCoord.out.bam --config config.txt --sid sampleid --id anyid --vcftype f --annotation n --itdetect y </code>
```

### CONFIG (Path to softwares and databases)
samtools	/home/tools_av/samtools-1.20/samtools
bedtools	/usr/local/bin/bedtools
vep	/home/Documents/checks/flt3/ensembl-vep/vep
itdetectpy	ITDetect.py
genome_fasta	/home/Documents/checks/picard/hg38.fa
genomic_portion_file	/home/Documents/checks/flt3/FLT3.exon.14.15.hg38.txt
path_to_itdetect	/home/Documents/checks/flt3/itdetect/
sub_region	FLT3.exon.14.15.hg38.bed
cachedir	/home/Documents/checks/flt3/.vep/

### Note:
#### For FLT3 duplication this script should be run with --alt_type dup
#### For large deletion this should be run --alt_type del
Download: Get the ITdetect from https://github.com/cpmsnuh/ITDetect

### References:
1. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-023-05173-8
2. https://asia.ensembl.org/info/docs/tools/vep/index.html
3. https://bedtools.readthedocs.io/en/latest/
4. https://www.htslib.org/
