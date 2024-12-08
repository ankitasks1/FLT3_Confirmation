import os
import sys
import re
import subprocess
import argparse
import time

## python3 ~/Documents/scripts/confirm_flt3_v4.py --span 6 --alt_type dup --bamfile 8807602.rmdup.realign.recal.bam --config config.txt --sid 8807602 --id x46 --vcftype f --itdetect y

#""""Ankit Verma (ankitverma9079@gmail.com)"""

def parse_arguments():
    parser = argparse.ArgumentParser(prog='finditds',description='Find ITD events and calculate allele burden for large mutated genomic regions eg. del and dup')
    parser.add_argument('--span', help='number of bases to cover inside and outside from both sides of ref and alt allele')
    parser.add_argument('--id', help='any identifier', default='x1234')
    parser.add_argument('--alt_type', help='eg. dup or del',  default='none')
    parser.add_argument('--bamfile', help='sorted bam file',  default='none')
    parser.add_argument('--itdetect', help='whether to run itdetect or not, if yes --itdetect y else do not give this option',  default='none')
    parser.add_argument('--sid', help='123456', default='123456')
    parser.add_argument('--vcftype', help='f=filtered, o=original', default='f')
    parser.add_argument('--config', help='path to config file which contains path to all the softwares', default='config.txt')

    # Print help message if no arguments
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    return args

def config_file(file):
    config = {}
    with open(file, "r") as configfile:
        configfile = configfile.readlines()
        for n in configfile:
            n = n.strip().split('\t')
            config[n[0]] = n[1]
    return config

def subset_bam(path_to_bedtools, samtools, bamfile, sub_region, sid):
    try:
        print("bam filtering started...")
        subbam = [path_to_bedtools,"intersect", "-a", bamfile, "-b", sub_region]
        with open(sid+"_sub.bam", 'wb') as out_file:
            subprocess.run(subbam, stdout=out_file)
        subbami = [samtools, "index", f"{sid}_sub.bam"]
        print("generating bam index file...")
        subprocess.run(subbami, capture_output=False, text=True)
        print("bam subsetted...")

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing subprocess: {e}")
    except Exception as e:
        print(f"Unexpected error in subset_bam: {e}")

def itd_detect(bamfile, itdetectpy,path_to_itdetect,genome_fasta, genomic_portion_file, sid):
    try:
        print("itd detection started...")
        itdetect = ["python3",itdetectpy, "-p", path_to_itdetect, "-r", genome_fasta, "-t", genomic_portion_file, "-b", bamfile, "-o", f"{sid}_itdetect.vcf"]
        subprocess.run(itdetect, capture_output=False, text=True)
        print("itd detected...")
    except FileNotFoundError as e:
        print(f"Error: {e}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing subprocess: {e}")
    except Exception as e:
        print(f"Unexpected error in itd_detect: {e}")


def parse_vcf(vcf):
    #don't run the command if the vcf file is empty
    if os.stat(vcf).st_size == 0:
        print(f"File {vcf} is empty. Skipping parsing.")
        raise ValueError(f"The vcf file {file} is empty. Cannot proceed.")
    variants_coord = []
    variants_info = []
    variants_freq = []
    print("parsing vcf file...")
    with open (vcf, 'r') as myfile:
        for i in myfile:
            # skip any line that start with "#"
            if not i.startswith('#'):
                line = i.strip().split('\t')
                chrm = line[0]
                start = int(line[1])
                end = start + (len(line[4]) - 1)
                ref= line[3]
                alt = line[4]
                freq = line[7].strip().split(';')
                variants_coord.append(f"{chrm}:{start}-{end}")
                variants_info.append(f"{chrm}:{start}{ref}>{alt}")
                variants_freq.append([f"{chrm}:{start}{ref}>{alt}",freq])
                # print(chr ,start ,end )
    return(variants_coord, variants_info, variants_freq)

def prep_regions(var_regions_list1, span, id, out_file, in_file):
    count = 1
    out_files_list = []
    in_files_list = []
    indiv_regions = [] # useful for getting variant coordinates
    for c in var_regions_list1:
        # print(c)
        region_out = []
        region_ins = []
        chrm = c.strip().split(':')[0]
        start = int(c.strip().split(':')[1].split('-')[0])
        end = int(c.strip().split(':')[1].split('-')[1])
        
        for n_span in range(0,(int(span)+1)):
            # print(n_span)
            region_out.append("".join(chrm +'\t' + str(start-n_span)))
            region_out.append("".join(chrm +'\t' + str(end+n_span)))
            region_ins.append("".join(chrm +'\t' + str(start+n_span)))
            region_ins.append("".join(chrm +'\t' + str(end-n_span)))
        # print(region_out, region_ins)
        """"Capture the files name"""
        out_files_list.append(out_file+ '_' + str(count) + '.txt')
        in_files_list.append(in_file+ '_' + str(count) + '.txt')
        indiv_regions.append(c)
        """Remove if any existing files with similar name"""
        if os.path.exists(out_file+ '_' + str(count) + '.txt'):
            os.remove(out_file+ '_' + str(count) + '.txt')
        if os.path.exists(in_file + '_' + str(count) + '.txt'):
            os.remove(in_file + '_' + str(count) + '.txt')
        """"Capture the output of splitted regions in individual file"""
        with open(out_file+ '_' + str(count) + '.txt', 'a') as onefile:
            for lines1 in region_out:
                onefile.write(lines1+'\n')
            onefile.close

        with open(in_file + '_' + str(count) + '.txt', 'a') as secfile:
            for lines2 in region_ins:
                secfile.write(lines2+'\n')
            secfile.close()
        count += 1
    print(f"generated coordinate files = {out_file} 1 to {int(count)-1}, {in_file} 1 to {int(count)-1}\n")
    return (in_files_list, out_files_list, indiv_regions)

def get_depth(samtools, in_files_list, out_files_list,indiv_regions, bamfile, sid):
    print('depth calculation started...')
    out_depth = []
    in_depth = []
    var = []
    count  = 1
    for pair in zip(in_files_list, out_files_list, indiv_regions):
        print(f"variant files {count}  {pair[0]} {pair[1]}")
        in_sam = [samtools, "depth", "-b", pair[0], f"{sid}_sub.bam"]
        out_sam = [samtools, "depth", "-b", pair[1], f"{sid}_sub.bam"]
        result_in_sam = subprocess.run(in_sam, capture_output=True, text=True)
        # print(result_in_sam)
        result_out_sam = subprocess.run(out_sam, capture_output=True, text=True)
        # print(result_out_sam)
        in_depth.append(result_in_sam)
        out_depth.append(result_out_sam)
        var.append(pair[2])
        # print(in_depth, out_depth)
        count  += 1
    print('depth calculaton done...')
    return in_depth, out_depth, var

def annotation(vep, vcf, sid, cachedir):
    print('annotation started...')
    try:
    #     command=[vep, "-i",vcf,"-o",f"{sid}_anno.txt","--cache", "--dir_cache", cachedir, "--force_overwrite","--hgvsg","--hgvs"]
    #     subprocess.run(command, capture_output=False, text=True)
        print('annotation completed...')
        with open(f"{sid}_anno.txt", "r") as myfile:
            for i in myfile:
                if not i.startswith('#'):
                    line = i.strip().split('\t')
                    variant = line[1] + ">" + line[2]
                    for annotype in line[13].split(';'):
                        if re.match("HGVSg", annotype):
                            genomic = annotype.split('=')[1]
                        if re.match("HGVSc", annotype):
                            cdna = annotype.split('=')[1]
                        if re.match("HGVSp", annotype):
                            protein = annotype.split('=')[1] 
                    consequences = line[6]
                    print('\nAnnotaion: ')
                    print(variant)
                    print(f"HGVSg: {genomic}")
                    print(f"HGVSc: {cdna} ")
                    print(f"HGVSp: {protein}")
                    print(f"VarType: {consequences} \n ")
    except FileNotFoundError as e:
        print(f"Error: {e}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing subprocess: {e}")

def manually_calculate_allele_burden(in_depth, out_depth, var, alt_attribute):
    count = 1
    for res in zip(in_depth, out_depth, var):
        out_num = []
        in_num = []
        for res1 in res[0].stdout.strip().split('\n'):
            res1 = res1.split('\t')
            in_num.append(int(res1[2]))
        for res2 in res[1].stdout.strip().split('\n'):
            res2 = res2.split('\t')
            out_num.append(int(res2[2]))
        res3 = res[2]
        # print(sum(out_num), sum(in_num))
        in_avg = round(sum(in_num)/len(in_num),2)
        out_avg = round(sum(out_num)/len(out_num),2)
        print(f"variant file {count} is for {res3}")
        if alt_attribute == "dup":
            alt_burden = round(((in_avg - out_avg)/in_avg)*100,2)
            print(f"Total depth (Alt+Ref) = {in_avg}\n Alt depth = {out_avg}\n Allele Burden = {alt_burden}%\n")
        if alt_attribute == "del":
            alt_burden = round(((out_avg - in_avg)/out_avg)*100,2)
            print(f"Total depth (Alt+Ref) = {out_avg}\n Alt depth = {in_avg}\n Allele Burden = {alt_burden}%\n")
        count += 1 

def software_calculate_allele_burden(freqs, alt_attribute):
    for lines in freqs:
        variant_desc = lines[0]
        for desc  in lines[1]:
            if re.match("RawReadCount",desc):
                total_depth = desc.split('=')[1]
            elif re.match("AltReadCount",desc):
                alt_depth = desc.split('=')[1]
            elif re.match("VAF",desc):
                alt_burden_per = float(desc.split('=')[1]) * 100
        print(variant_desc)
        if alt_attribute == "dup":
            print(f"Total depth (Alt+Ref) = {total_depth}\n Alt depth = {alt_depth}\n Allele Burden = {alt_burden_per}%\n")
        if alt_attribute == "del":
            print(f"Total depth (Alt+Ref) = {total_depth}\n Alt depth = {alt_depth}\n Allele Burden = {alt_burden_per}%\n")
        
        
    
def main():
    start_time = time.time()

    args = parse_arguments()
    
    # --Execute commands-- #
    print('\nAnalysis Started\n-------------------\nReading config file...\n')
    config = config_file(args.config)
    # print(config)

    print('\nSubsetting BAM...')
    subset_bam(config['bedtools'], config['samtools'], args.bamfile, config['sub_region'], args.sid)

    print('\nDetecting ITDs...')
    if args.itdetect.lower() == 'y':
        print("ITD detect is enabled.")
        itd_detect(args.bamfile, config['itdetectpy'], config['path_to_itdetect'], config['genome_fasta'], config['genomic_portion_file'], args.sid)
    else:
        print("ITD detect is not enabled or 'itdetect' argument is missing.")
        pass

    print('\nParsing VCF...')
    if args.vcftype == 'f':
        vcf_info = parse_vcf(args.sid + "_itdetect.vcf.filtered.vcf")
    elif args.vcftype == 'o':
        vcf_info = parse_vcf(args.sid + "_itdetect.vcf")

    print('\nPreparing Variant Region files...')
    var_regions_files = prep_regions(vcf_info[0], args.span, args.id, f"region_out_{args.id}", f"region_in_{args.id}")

    print('\nGetting Depth for...')
    res1,res2,res3 = get_depth(config['samtools'],var_regions_files[0], var_regions_files[1], var_regions_files[2], args.bamfile, args.sid)

    print("\nGetting annotation...")
    if args.vcftype == 'f':
        annotation(config['vep'], args.sid + "_itdetect.vcf.filtered.vcf", args.sid, config['cachedir'])
    elif args.vcftype == 'o':
        annotation(config['vep'], args.sid + "_itdetect.vcf", args.sid, config['cachedir'])

    print("\nManually calculated: Allele burden")
    print("-----------------------------------")
    manually_calculate_allele_burden(res1, res2, res3, args.alt_type)

    print("\nSoftware calculated: Allele burden")
    print("-----------------------------------")
    software_calculate_allele_burden(vcf_info[2], args.alt_type)

    print('Thanks')
    end_time = time.time()
    seconds = end_time - start_time
    minutes = seconds / 60
    print('Time taken in minutes:', round(minutes,2)) 


if __name__ == '__main__':
    main()
