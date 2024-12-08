import sys, getopt, os
import collections
####ITDetect v.1.4 (2021.04.30)####
####Heejun Jang (starz77@snu.ac.kr)####
####Sungyoung Lee (me@lsy.io)####

program = 'ITDetect'
version = 'v.1.4'
lastupdate = '2021.04.19'

print(program)
print('Version: '+version)
print('Last Update: '+lastupdate)
###Command: python run_ITDetect.py \
#                -p 'path/to/ITDetect/' \
#                -b 'path/to/input.bam' \
#                -o 'outputDir' \
#                -r 'reference.file' \
#                -c 'minimum number of support reads for the ITD called (default: 3)' \
#                -p 'minimum number of repeats in the region where ITD is called (default: 10)' \
#                -t 'target gene file (default: FLT3.hg19.exon14.exon15.txt)' \
#                -m 'temporary directory (default: .)' \
#                -v 'Run ITDetect in verbose mode (default: disabled)' \
#                -h 'Print help message'
#                -e 'run with example input files'

path = '.'
bam = ''
output = ''
reference = ''
rcnt = 3
ecnt = 6
repeat = 10
target = path+'/data/FLT3.exon.hg19.txt'
tmp = '.'
verbose = ''
help1 = ''
exam = ''

try:
    opts, args = getopt.getopt(sys.argv[1:], "p:b:o:r:c:n:q:t:m:vhe",["path=", "bam=", "output=", "reference=", "rcnt=", "ecnt=", "repeat=", "target=", "tmp=", "verbose", "help", "example"])
except getopt.GetoptError as err:
    print(str(err))
    sys.exit(1)
for opt, arg in opts:
    if opt in ('-p', '--path'):
        path = arg
    elif opt in ('-b', '--bam'):
        bam = arg
    elif opt in ('-o', '--output'):
        output = arg
    elif opt in ('-r', '--reference'):
        reference = arg
    elif opt in ('-c', '--rcnt'):
        rcnt = arg
    elif opt in ('-n', '--ecnt'):
        ecnt = arg
    elif opt in ('-q', '--repeat'):
        repeat = arg
    elif opt in ('-t', '--target'):
        target = arg
    elif opt in ('-m', '--tmp'):
        tmp = arg
    elif opt in ('-v', '--verbose'):
        verbose = '-verbose'
    elif opt in ('-h', '--help'):
        help1 = '-help'
    elif opt in ('-e', '--example'):
        bam = ''
        output = ''
        reference = ''
        target = ''

if help1 == '':
    err = False
    if len(bam) == 0:
        print('\nERROR : -b/--bam is required')
        err = True
    if len(reference) == 0:
        print('\nERROR : -r/--reference is required')
        err = True
    if len(output) == 0:
        print('\nERROR : -o/--output is required')
        err = True

    if err:
        exit()
    if not os.path.exists(path+'/bin/ITDetect.jar'):
        print('\nERROR : jar path ['+path+'/bin/ITDetect.jar] does not exist!')
        exit()
    cmd = 'java -jar '+path+'/bin/ITDetect.jar'+' -b '+bam+' -o '+output+' -r '+reference+' -t '+target+' -rcnt '+str(rcnt)+' -repeat '+str(repeat)+' -tmp '+tmp+' '+verbose
    print(cmd)
    os.system(cmd)
    if os.path.exists(output):
        cmd2 = 'python '+path+'/bin/ITDetect_Filter.py '+output+' '+target+ ' '+str(ecnt)
        os.system(cmd2)
else:
    cmd = 'java -jar '+path+'/bin/ITDetect.jar'+' '+help1
    print(cmd)
    os.system(cmd)
    
if tmp == '.':
    cmd = 'mv '+tmp+'/*.candipos ' +output+'.candipos.txt  >/dev/null 2>&1'
    os.system(cmd)
else:
    pass
