"""Calculate within-population metrics for diversity and selection.  Input is a tab delimited file with no header, containing Scaff, pos, an, dp, and genotypes coded as 0-4.
Input file custom: Filename should end with _XXX_raw.table.recode.txt, where XXX is three-letter population abbreviation"""

import argparse
import numpy
import math
from random import randint

# export PYTHONPATH="$PYTHONPATH:/Users/monnahap/Documents/Research/code/GenomeScan/"

# Substitute missingness for number of individuals to downsample to.  That way you can enter one number for all populations to be downsampled to.


def calcAFS(input_file, output, prefix, window_size, num_bootstraps, sampind=5):

    scaff_lengths = [33684748, 19642879, 24872290, 23717143, 21575646, 25532148, 25060017, 23333815]  # Scaffold lengths determined from alygenomes.fasta
    num_winds = sum([float(x) / float(window_size) for x in scaff_lengths])
    per_scaff = [int(math.ceil(float(x) / float(window_size))) for x in scaff_lengths]
    wind_counts = []
    for j in range(0, num_bootstraps):
        wind_counts.append([0 for x in range(0, int(num_winds))])
        for x in range(0, int(num_winds)):
            wind_counts[j][randint(0, int(num_winds) - 1)] += 1


    with open(input_file, 'r') as infile:
        for i, line in enumerate(infile):
            # print(line)

            line = line.strip("\n")
            line = line.strip("\t")
            line = line.split("\t")

            pop, ploidy, scaff, pos, an, dp = line[:6]

            gt = line[6:]

            if i % 100000 == 0:
                print(i)
            if i == 0:
                outfile = output + prefix + "_AFS.txt"
                out1 = open(outfile, 'w')
                AN = int(sampind * int(float(ploidy)))
                AFS = [0 for cat in range(0, AN + 1)]
                BSFS = {}
                for rep in range(0, num_bootstraps):
                    BSFS[str(rep)] = [0 for cat in range(0, AN + 1)]

            gt = [x for x in gt if x != "-9"]
            if len(gt) >= sampind:
                sgt = numpy.random.choice(gt, size=sampind, replace=False)
                sac = sum([int(x) for x in sgt])
                AFS[sac] += 1
                scaff_num = int(scaff.split("_")[1])
                cur_pos = float(pos)
                wind_num = sum(per_scaff[:scaff_num - 1]) + int(math.ceil(cur_pos / float(window_size))) - 1
                for j, rep in enumerate(wind_counts):
                    for samp_num in range(0, rep[wind_num]):
                        BSFS[str(j)][sac] += 1
        print(AFS)
        print(BSFS['1'])
        for k, j in enumerate(AFS):
            out1.write(prefix + '\t-9\t' + str(k) + "\t" + str(j) + '\n')
        for g in range(0, num_bootstraps):
            for k, j in enumerate(BSFS[str(g)]):
                out1.write(prefix + '\t' + str(g + 1) + '\t' + str(k) + "\t" + str(j) + '\n')


    out1.close()
    return AFS

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar='input_file', required=True, help='input file created with recode012.py')
    parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Output Directory')
    parser.add_argument('-p', type=str, metavar='prefix', required=True, help='prefix for file name of output')
    parser.add_argument('-sampind', type=int, metavar='DownSampled_individuals', required=False, default='5', help='Number of individuals to downsample the data to')
    parser.add_argument('-ws', type=int, metavar='window_size', required=False, default='50000', help='Size of windows for bootstrap replicates')
    parser.add_argument('-nrep', type=int, metavar='number_of_bootstrap_reps', required=False, default='200', help='Number of bootstrap replicates of the AFS')

    args = parser.parse_args()

    j1 = calcAFS(args.i, args.o, args.p, args.ws, args.nrep, args.sampind)
