import os
import argparse
import gzip

# create variables that can be entered in the command line
parser = argparse.ArgumentParser()
parser.add_argument('-v', type=str, metavar='vcf_path', required=True, help='path to vcfs')
parser.add_argument('-w', type=str, metavar='Window_Size', required=True, help='size of scaffold window')
parser.add_argument('-dp', type=int, metavar='min_depth', required=True, help='minimum genotype depth to be considered non-missing')
parser.add_argument('-o', type=str, metavar='output_path', required=True, help='full path of output file')
parser.add_argument('-gz', type=str, metavar='gzipped?', required=True, help='are vcfs gzipped (true) or not (false)')

args = parser.parse_args()

if args.v.endswith("/") is False:
    args.v += "/"
vcf_list = []

for file in os.listdir(args.v):  # get names of vcf files in args.v directory
    if args.gz == 'true':
        if file[-3:] == '.gz':
            vcf_list.append(file)

    elif args.gz == 'false':
        if file[-3:] == 'vcf':
            vcf_list.append(file)

    else:
        print('error')



out1 = open(args.o, 'w')
out1.write("scaff\tstart\tend\tnumSites\tmissing\n")
for iii, vcf in enumerate(vcf_list):
    if args.gz == 'true':
        src = gzip.open(args.v + vcf)
    elif args.gz == 'false':
        src = open(args.v + vcf)
    M = 0.0
    tot_count = 0
    site_count = 0
    start = 0
    end = int(args.w)
    window_size = int(args.w)
    # evaluate contents of each line of input file
    for line_idx, line in enumerate(src):  # Cycle over lines in the VCF file
        cols = line.replace('\n', '').split('\t')  # Split each line of vcf
        if line_idx % 10000 == 0:
            print(line_idx)
        if len(cols) < 2:          # This should be info just before header
            pass
        elif cols[0] == "#CHROM":  # This should be header
            pass
        else:
            genos = []
            scaff = cols[0]               # parse important info from each line
            pos = int(cols[1])
            info = cols[7].split(";")
            AN = float(info[2].split("=")[1])
            AC = float(info[0].split("=")[1])
            num_missing = 0
            num_ind = float(len(cols[9:]))
            m = 0.0
            if pos > start and pos <= end:
                site_count += 1
                for ind in cols[9:]:
                    ind = ind.split(":")
                    geno = ind[0].split("/")
                    if geno[0] == ".":
                        m += 1.0 / num_ind
                    else:
                        try:
                            if int(ind[2]) < args.dp:
                                m += 1.0 / num_ind
                        except (IndexError, ValueError):
                            pass
                M += m

            elif pos > end:
                M = M / float(site_count)
                out1.write(scaff + '\t' +
                           str(start) + '\t' +
                           str(end) + '\t' +
                           str(args.w) + '\t' +
                           str(site_count) + '\t' +
                           str(M) + '\n')
                site_count = 0
                M = 0.0

                while pos > end:
                    end += window_size
                    start = end - window_size
                    if pos > end:
                        out1.write(scaff + '\t' +
                                   str(start) + '\t' +
                                   str(end) + '\t' +
                                   str(args.w) + '\t' +
                                   str(site_count) + '\t' +
                                   str(-99) + '\n')

                if int(pos) > start and int(pos) <= end:
                    site_count += 1
                    for ind in cols[9:]:
                        ind = ind.split(":")
                        geno = ind[0].split("/")
                        if geno[0] == ".":
                            m += 1.0 / num_ind
                        else:
                            try:
                                if int(ind[2]) < args.dp:
                                    m += 1.0 / num_ind
                            except (IndexError, ValueError):
                                pass
                    M += m


    out1.write(scaff + '\t' +
               str(start) + '\t' +
               str(end) + '\t' +
               str(args.w) + '\t' +
               str(site_count) + '\t' +
               str(M) + '\n')


out1.close()
