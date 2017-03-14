
import argparse

parser = argparse.ArgumentParser(description='This program takes a merged vcf with highest-coverage individuals from specific pops as input, and outputs a lookup key from the merged vcf')
parser.add_argument('-i', type=str, metavar='input_file', required=True, help='full path to input file')
parser.add_argument('-r', type=str, metavar='repolarization_key', required=True, help='key containing scaffold and position of sites that need to be flipped')
parser.add_argument('-o', type=str, metavar='Output_Prefix', required=True, help='should be absolute path to output directory plus population name.')

args = parser.parse_args()

new = open(args.o + ".table.repol.txt", 'w')

key = open(args.r, 'r')
key_list = []

for i, line in enumerate(key):
    line = line.strip("\n").split("\t")
    scaff = line[0]
    pos = line[1]
    if i == 0:
        old_scaff = scaff
        Scaff = [pos]
    elif scaff == old_scaff:
        Scaff.append(pos)
    else:
        key_list.append(Scaff)
        Scaff = [pos]
        old_scaff = scaff
    key_list.append(Scaff)

with open(args.i) as inf:
    for j, line in enumerate(inf):
        cols = line.strip("\n").strip("\t").split("\t")
        pop, ploidy, scaff, pos, an, dp = cols[:6]
        scaff = int(scaff.split("_")[1])
        gt = ""
        if pos in key_list[scaff - 1]:
            for geno in cols[6:]:
                if geno != "-9":
                    gt += str(int(float(ploidy)) - int(geno)) + "\t"
                else:
                    gt += "-9\t"
            gt.strip("\t")
            newline = pop + "\t" + str(ploidy) + "\tscaffold_" + str(scaff) + "\t" + pos + "\t" + an + "\t" + dp + "\t" + gt + "\n"
            new.write(newline)
        else:
            new.write(line)
