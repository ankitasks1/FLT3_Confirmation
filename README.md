# FLT3 Confirmation
## Repository for "Tools and materials to confirm FLT3"

### check the usage and parameters
python3 ~/Documents/scripts/confirm_flt3_v5.py -h

### with ITDs detection:
python3 ~/Documents/scripts/confirm_flt3_v5.py --span 6 --alt_type dup --bamfile sample1.Aligned.sortedByCoord.out.bam --config config.txt --sid sample1 --id rxt59 --vcftype f --itdetect y

### without ITDs detection
python3 ~/Documents/scripts/confirm_flt3_v5.py --span 6 --alt_type dup --bamfile sample1.Aligned.sortedByCoord.out.bam --config config.txt --sid sample1 --id rxt59 --vcftype f 

### Note:
#### For FLT3 duplication this script should be run with --alt_type dup
#### For large deletion this should be run --alt_type del
Download: Get the ITdetect from https://github.com/cpmsnuh/ITDetect
