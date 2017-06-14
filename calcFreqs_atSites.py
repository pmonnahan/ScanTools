import argparse


def calcFreqs(input_directory, output_directory, output_filename, sites_file, pops, suffix):
    if input_directory.endswith("/") is False:
        input_directory += "/"
    if output_directory.endswith("/") is False:
        output_directory += "/"
    if output_filename.endswith(".txt") is False:
        output_filename += ".txt"

    out1 = open(output_directory + output_filename, 'w')

    site_list = [[] for i in range(0, 9)]
    with open(sites_file) as sf:
        for site in sf:
            site = site.split("\t")
            site_list[int(site[0])].append(site[1])

    pop_list = pops.split(",")
    print(pop_list)

    for pop in pop_list:
        with open(input_directory + pop + suffix) as pop_data:
            for line in pop_data:
                line = line.strip("\n").strip("\t").split("\t")
                scaff = int(line[2].split("_")[1])
                pos = line[3]
                an = float(line[4])
                if pos in site_list[scaff]:
                    genos = line[6:]
                    ac = float(sum([int(x) for x in genos if x != "-9"]))
                    freq = ac / an
                    if freq != 0.0:
                        print(line)
                        print(ac, an, freq, genos)
                    out1.write(pop + "\t" + "scaffold_" + str(scaff) + "\t" + pos + "\t" + str(freq) + "\n")

    return site_list


if __name__ == '__main__':  # Used to run code from command line

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar='input_directory', required=True, help='input directory containing files created with recode012.py')
    parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Output Directory')
    parser.add_argument('-of', type=str, metavar='output_filename', required=True, help='Output Filename')
    parser.add_argument('-s', type=str, metavar='sites_file', required=True, help='full absolute path to file containing scaffold and position for each site')
    parser.add_argument('-pops', type=str, metavar='populations', required=True, help='Size of windows in bp')
    parser.add_argument('-suffix', type=str, metavar='pop_data_suffx', required=True, help='suffix of file containing population data')
    args = parser.parse_args()

    j1 = calcFreqs(args.i, args.o, args.of, args.s, args.pops, args.suffix)
