import argparse

# create variables that can be entered in the command line
parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, metavar='input_table_file', required=True, help='path to input table file (from -VariantsToTable) of scaffolds containing just scaff, pos, ac, an, dp and GT fields in that order')
parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='')
parser.add_argument('-pop', type=str, metavar='population_name', required=True, help='')

args = parser.parse_args()

filename = args.i.split('/')[-1]

if args.o.endswith("/") is False:
    args.o += "/"

with open(args.i, "rU") as table:
    prefix = args.i[:-6]
    GTfile = open(args.o + filename + ".recode.txt", "w")
    for i, line in enumerate(table):
        line = line.strip("\n")
        line = line.split("\t")
        if line[0].split("_")[0] == 'scaffold':
            ref = '0'
            numind = len(line[5:])
            scaff = line[0]
            pos = line[1]
            ac = line[2]
            an = line[3]
            dp = int(line[4])
            ploidy = float(len(line[5].split("/")))
            GT = args.pop + '\t' + str(ploidy) + '\t' + scaff + '\t' + pos + '\t' + ac + '\t' + an + '\t' + str(dp) + '\t'
            alt = False
            numobs = 0
            for j, gt in enumerate(line[5:]):
                gt = gt.split("/")
                if gt[0] == '.':
                    GT += '-9\t'
                else:
                    gtc = int(len(gt))
                    if ref == '0':
                        ref = gt[0]
                    for g in gt:
                        if g == ref:
                            gtc -= 1
                        else:
                            alt = True
                    GT += str(gtc) + "\t"
                    numobs += 1
            if i % 100000 == 0:
                print(i)
            GT.strip("\t")
            GT += "\n"
            GTfile.write(GT)

    print("number of individuals = ", numind)
