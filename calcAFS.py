"""Calculate within-population metrics for diversity and selection.  Input is a tab delimited file with no header, containing Scaff, pos, an, dp, and genotypes coded as 0-4.
Input file custom: Filename should end with _XXX_raw.table.recode.txt, where XXX is three-letter population abbreviation"""

import argparse
import numpy

# export PYTHONPATH="$PYTHONPATH:/Users/monnahap/Documents/Research/code/GenomeScan/"

# Substitute missingness for number of individuals to downsample to.  That way you can enter one number for all populations to be downsampled to.


def calcAFS(input_file, output, prefix, sampind=5):

    with open(input_file, 'r') as infile:
        for i, line in enumerate(infile):

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
            else:
                gt = [x for x in gt if x != "-9"]
                if len(gt) >= sampind:
                    sgt = numpy.random.choice(gt, size=sampind, replace=False)
                    sac = sum([int(x) for x in sgt])
                    AFS[sac] += 1

        for k, j in enumerate(AFS):
            out1.write(prefix + '\t' + str(k) + "\t" + str(j) + '\n')

    out1.close()
    return AFS

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar='input_file', required=True, help='input file created with recode012.py')
    parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Output Directory')
    parser.add_argument('-p', type=str, metavar='prefix', required=True, help='prefix for file name of output')
    parser.add_argument('-sampind', type=int, metavar='DownSampled_individuals', required=False, default='5', help='Number of individuals to downsample the data to')

    args = parser.parse_args()

    j1 = calcAFS(args.i, args.o, args.p, args.sampind)
