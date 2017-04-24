"""Calculate within-population metrics for diversity and selection.  Input is a tab delimited file with no header, containing Scaff, pos, an, dp, and genotypes coded as 0-4.
Input file custom: Filename should end with _XXX_raw.table.recode.txt, where XXX is three-letter population abbreviation"""

import argparse
import numpy

# export PYTHONPATH="$PYTHONPATH:/Users/monnahap/Documents/Research/code/GenomeScan/"

# Substitute missingness for number of individuals to downsample to.  That way you can enter one number for all populations to be downsampled to.
def calcwpm(input_file, output, prefix, sampind=5, window_size=50000, minimum_snps=2):

    snp_count = 0
    Snp_count = 0
    tot_count = 0
    Tot_count = 0
    start = 0.0
    end = window_size
    winexclcount = 0
    num_wind = 0

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
                outfile = output + prefix + "_WPM.txt"
                out1 = open(outfile, 'w')
                out1.write("pop\tploidy\tsampind\tscaff\tstart\tend\twin_size\tnum_snps\tnum_sites\tnum_singletons\tavg_freq\tavg_Ehet\tDiversity\tThetaW\tPi\tThetaH\tThetaL\tD\tH\tE\n")
                oldscaff = scaff
                AN = int(sampind * int(float(ploidy)))
                n = float(AN)
                p = []
                Ehet = []
                num_sing = 0
                Num_sing = 0
                afs = [0 for cat in range(0, AN + 1)]
                AFS = [0 for cat in range(0, AN + 1)]
                aw = 0.0
                bw = 0.0  # This is b sub n+1 in Zeng. a just goes to n but b goes to n+1
                a2 = 0.0
                for j in range(1, AN):
                    aw += 1.0 / float(j)  # a1 in Tajima 1989 and an in Zeng 2006
                    a2 += 1.0 / float(j**2)  # This is bn according to Zeng 2006
                    bw += 1.0 / float(j**2)
                bw += 1.0 / float(n**2)  # This is the n+1 part from Zeng
                b1 = (n + 1) / (3 * (n - 1))  # b1 through e2 are taken directly rom Tajima 1989
                b2 = (2 * (n**2 + n + 3)) / (9 * (n * (n - 1)))
                c1 = b1 - (1 / aw)
                c2 = b2 - (n + 2) / (aw * n) + a2 / aw**2
                e1 = c1 / aw
                e2 = c2 / (aw**2 + a2)

            if int(pos) > start and int(pos) <= end and int(an) >= AN and scaff == oldscaff:
                gt = [x for x in gt if x != "-9"]
                if len(gt) >= sampind:
                    sgt = numpy.random.choice(gt, size=sampind, replace=False)
                    sac = sum([int(x) for x in sgt])
                    tot_count += 1
                    Tot_count += 1
                    if sac != 0 and sac != AN:
                        if sac == 1:
                            num_sing += 1
                            Num_sing += 1
                        p1 = float(sac) / float(AN)
                        p.append(p1)
                        Ehet.append(p1 * (1 - p1))
                        afs[sac] += 1
                        AFS[sac] += 1
                        snp_count += 1
                        Snp_count += 1

            elif int(pos) > end or scaff != oldscaff:
                Pi = 0.0
                h = 0.0
                L = 0.0
                n = float(AN)

                if snp_count >= minimum_snps:

                    num_wind += 1
                    S = float(sum(afs[1:-1]))
                    W = S / aw
                    W2 = S * (S - 1) / ((aw**2) + a2)
                    for j in range(1, AN):

                        Pi += afs[j] * j * (AN - j)
                        h += afs[j] * (j**2)
                        L += afs[j] * j

                    Pi = 2 * Pi / (n * (n - 1))
                    h = 2 * h / (n * (n - 1))
                    L = L / (n - 1)
                    div = Pi / float(tot_count)


                    # When calculating variance of D, H, and E should we use W value for window or for entire genome.  Probably for window??

                    varPi_W = e1 * S + e2 * S * (S - 1)
                    varPi_L1 = (((n - 2) / (6 * n - 6)) * W)
                    varPi_L2 = (((18 * n**2 * (3 * n + 2) * bw) - (88 * n**3 + 9 * n**2 - 13 * n + 6)) / (9 * n * (n - 1)**2) * W2)
                    varL_W = (((n / (2 * n - 2)) - 1 / aw) * W) + ((a2 / aw**2 + (2 * (n / (n - 1))**2) * a2 - 2 * (n * a2 - n + 1) / ((n - 1) * aw) - (3 * n + 1) / (n - 1)) * W2)
                    varPi_L = varPi_L1 + varPi_L2
                    try:
                        D = (Pi - W) / (varPi_W**0.5)
                        H = (Pi - L) / (varPi_L**0.5)  # Normalized H according to Zeng 2006
                        E = (L - W) / (varL_W**0.5)
                        out1.write(pop + '\t' +
                                   str(ploidy) + '\t' +
                                   str(sampind) + '\t' +
                                   scaff + '\t' +
                                   str(start) + '\t' +
                                   str(end) + '\t' +
                                   str(window_size) + '\t' +
                                   str(snp_count) + '\t' +
                                   str(tot_count) + '\t' +
                                   str(num_sing) + '\t' +
                                   str(numpy.mean(p)) + '\t' +
                                   str(numpy.mean(Ehet)) + '\t' +
                                   str(div) + '\t' +
                                   str(W) + '\t' +
                                   str(Pi) + '\t' +
                                   str(h) + '\t' +
                                   str(L) + '\t' +
                                   str(D) + '\t' +
                                   str(H) + '\t' +
                                   str(E) + '\n')
                    except ZeroDivisionError:
                        print('ZeroDivisionError:', varPi_W, aw, bw, a2, b1, b2, c1, c2, e1, e2, S)


                else:
                    winexclcount += 1

                snp_count = 0
                tot_count = 0
                p = []
                Ehet = []
                num_sing = 0
                afs = [0 for cat in range(0, AN + 1)]

                if float(pos) > end:
                    while float(pos) > end:
                        end += window_size
                    start = end - window_size
                elif scaff != oldscaff:
                    oldscaff = scaff

                    start = 0.0
                    end = window_size

                if int(pos) > start and int(pos) <= end and int(an) >= AN and scaff == oldscaff:
                    gt = [x for x in gt if x != "-9"]
                    if len(gt) >= sampind:
                        sgt = numpy.random.choice(gt, size=sampind, replace=False)
                        sac = sum([int(x) for x in sgt])
                        tot_count += 1
                        Tot_count += 1
                        if sac != 0 and sac != AN:
                            if sac == 1:
                                num_sing += 1
                                Num_sing += 1
                            p1 = float(sac) / float(AN)
                            p.append(p1)
                            Ehet.append(p1 * (1 - p1))
                            afs[sac] += 1
                            AFS[sac] += 1
                            snp_count += 1
                            Snp_count += 1

        if snp_count >= minimum_snps:
            Pi = 0.0
            h = 0.0
            L = 0.0
            n = float(AN)
            num_wind += 1
            S = float(sum(afs[1:-1]))
            W = S / aw
            W2 = S * (S - 1) / ((aw**2) + a2)
            for j in range(1, AN):

                Pi += afs[j] * j * (AN - j)
                h += afs[j] * (j**2)
                L += afs[j] * j

            Pi = 2 * Pi / (n * (n - 1))
            h = 2 * h / (n * (n - 1))
            L = L / (n - 1)

            div = Pi / float(tot_count)

            # When calculating variance of D, H, and E should we use W value for window or for entire genome.  Probably for window??

            varPi_W = e1 * S + e2 * S * (S - 1)
            varPi_L1 = (((n - 2) / (6 * n - 6)) * W)
            varPi_L2 = (((18 * n**2 * (3 * n + 2) * bw) - (88 * n**3 + 9 * n**2 - 13 * n + 6)) / (9 * n * (n - 1)**2) * W2)
            varL_W = (((n / (2 * n - 2)) - 1 / aw) * W) + ((a2 / aw**2 + (2 * (n / (n - 1))**2) * a2 - 2 * (n * a2 - n + 1) / ((n - 1) * aw) - (3 * n + 1) / (n - 1)) * W2)
            varPi_L = varPi_L1 + varPi_L2
            try:
                D = (Pi - W) / (varPi_W**0.5)
                H = (Pi - L) / (varPi_L**0.5)  # Normalized H according to Zeng 2006
                E = (L - W) / (varL_W**0.5)
                out1.write(pop + '\t' +
                           str(ploidy) + '\t' +
                           str(sampind) + '\t' +
                           scaff + '\t' +
                           str(start) + '\t' +
                           str(end) + '\t' +
                           str(window_size) + '\t' +
                           str(snp_count) + '\t' +
                           str(tot_count) + '\t' +
                           str(num_sing) + '\t' +
                           str(numpy.mean(p)) + '\t' +
                           str(numpy.mean(Ehet)) + '\t' +
                           str(div) + '\t' +
                           str(W) + '\t' +
                           str(Pi) + '\t' +
                           str(h) + '\t' +
                           str(L) + '\t' +
                           str(D) + '\t' +
                           str(H) + '\t' +
                           str(E) + '\n')
            except ZeroDivisionError:
                print('ZeroDivisionError:', varPi_W, aw, bw, a2, b1, b2, c1, c2, e1, e2)

    Pi = 0.0
    h = 0.0
    L = 0.0
    n = float(AN)
    S = float(sum(AFS[1:-1]))
    W = S / aw
    W2 = S * (S - 1) / ((aw**2) + a2)
    for j in range(1, AN):
        Pi += AFS[j] * j * (AN - j)
        h += AFS[j] * (j**2)
        L += AFS[j] * j

    Pi = 2 * Pi / (n * (n - 1))
    h = 2 * h / (n * (n - 1))
    L = L / (n - 1)

    div = Pi / float(Tot_count)

    # When calculating variance of D, H, and E should we use W value for window or for entire genome.  Probably for window??

    varPi_W = e1 * S + e2 * S * (S - 1)
    varPi_L1 = (((n - 2) / (6 * n - 6)) * W)
    varPi_L2 = (((18 * n**2 * (3 * n + 2) * bw) - (88 * n**3 + 9 * n**2 - 13 * n + 6)) / (9 * n * (n - 1)**2) * W2)
    varL_W = (((n / (2 * n - 2)) - 1 / aw) * W) + ((a2 / aw**2 + (2 * (n / (n - 1))**2) * a2 - 2 * (n * a2 - n + 1) / ((n - 1) * aw) - (3 * n + 1) / (n - 1)) * W2)
    varPi_L = varPi_L1 + varPi_L2
    try:
        D = (Pi - W) / (varPi_W**0.5)
        H = (Pi - L) / (varPi_L**0.5)  # Normalized H according to Zeng 2006
        E = (L - W) / (varL_W**0.5)
        out1.write(pop + '\t' +
                   str(ploidy) + '\t' +
                   str(sampind) + '\t' +
                   "Genome" + '\t' +
                   "-99" + '\t' +
                   "-99" + '\t' +
                   "-99" + '\t' +
                   str(Snp_count) + '\t' +
                   str(Tot_count) + '\t' +
                   str(Num_sing) + '\t' +
                   "-99" + '\t' +
                   "-99" + '\t' +
                   str(div) + '\t' +
                   str(W) + '\t' +
                   str(Pi) + '\t' +
                   str(h) + '\t' +
                   str(L) + '\t' +
                   str(D) + '\t' +
                   str(H) + '\t' +
                   str(E) + '\n')
    except ZeroDivisionError:
        print('ZeroDivisionError:', varPi_W, aw, bw, a2, b1, b2, c1, c2, e1, e2)


    out1.close()

    return num_wind, winexclcount, W, Pi, h, L, D, H, E


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar='input_file', required=True, help='input file created with recode012.py')
    parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Output Directory')
    parser.add_argument('-p', type=str, metavar='prefix', required=True, help='prefix for file name of output')
    parser.add_argument('-sampind', type=int, metavar='DownSampled_individuals', required=False, default='5', help='Number of individuals to downsample the data to')
    parser.add_argument('-ws', type=float, metavar='window_size', required=False, default='10000.0', help='Size of windows in bp')
    parser.add_argument('-ms', type=int, metavar='minimum_snps', required=False, default='2', help='minimum number of snps in a window')

    args = parser.parse_args()

    j1, j2, j3, j4, j5, j6, j7, j8, j9 = calcwpm(args.i, args.o, args.p, args.sampind, args.ws, args.ms)
