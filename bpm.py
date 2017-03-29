"""Calculate between-population divergence metrics.  Input is a tab delimited file with no header, containing Scaff, pos, ac, an, dp, and genotypes coded as 0-4.
Input file custom: Filename should end with _XXX_raw.table.recode.txt, where XXX is three-letter population abbreviation
Calculations are roughly equivalent to F. Kolar's calculations for toy dataset.  
His rho's:  0.176 and 0.226.  Mine: 0.169 and 0.197
His Fst's:  0.05 and 0.099.  Mine: 0.093 and 0.131"""

import argparse

# export PYTHONPATH="$PYTHONPATH:/Users/monnahap/Documents/Research/code/GenomeScan/"



def calcBPM(input_file, output, outname, window_size, minimum_snps, num_pops):

    def NestedAnova(locus_list):
        locus = []
        for l in locus_list:  # Remove missing data from the site
            locus.append([x for x in l if x != "-9"])
        r = float(len(locus))  # Number of populations:  Should always be 2 for pairwise analysis
        p_i = []  # Allele frequency of each population
        n_i = []  # Number of individuals in each population
        ploidy_list = []  # Ploidy of each population
        Fssg = 0.0  # Among population sum of squares for Fst
        Fssi = 0.0  # Among individuals within populations sum of squares for Fst
        Fssw = 0.0  # Within individuals sum of squares for Fst
        Fsst = 0.0  # Total sum of squares:  confirmed that this equals sum of components within rounding error
        Rssg = 0.0  # Among population sum of squares for Rho
        Rssi = 0.0  # Among individual within populations sum of squares for Rho
        Rsst = 0.0  # Total sum of squres:  confirmed that this equals sum of components within rounding error
        p_ij = []  # Within-individual allele frequency for each individual in each population
        tac = 0  # Total alternative allele count used for determining allele frequency
        tan = 0  # Total allele number
        ac_ij = []  # Alternative allele count per individual
        df_w = 0.0  # degrees of freedom for Fssw
        FS_Nij2 = 0.0  # Attempting to mirror spagedi code; Used for calculation of sample size coefficients in variance components
        FS_SNij2 = 0.0
        RS_SNij2 = 0.0
        FS_SNij2_over_SNij = 0.0
        for pop_site in locus:
            FSNij2temp = 0.0
            p = []  # List containing individual allele frequencies for a population
            ac_i = []  # List containing individual ac counts for a population
            ploidy = float(pop_site[1])
            nnn = float(len(pop_site[6:]))  # sample size for thispopulation
            an = ploidy * nnn  # Calculation an for the population
            ac = sum([float(geno) for geno in pop_site[6:]])  # Calculate ac for a population
            n_i.append(nnn)
            p_i.append(ac / an)
            ploidy_list.append(ploidy)
            pop_an = 0
            for ind in pop_site[6:]:
                p_ind = float(ind) / ploidy  # Calculate individual's allele frequency
                p.append(p_ind)
                ac_i.append(float(ind))
                num_alleles = 0
                for c in range(0, int(ploidy) - int(ind)):  # Calculation of tan and tac for p_bar calculation
                    num_alleles += 1
                    tan += 1
                for c in range(0, int(ind)):
                    num_alleles += 1
                    tac += 1
                    tan += 1
                FS_Nij2 += ploidy**2  # Squared sample size per individual
                FSNij2temp += ploidy**2
                pop_an += ploidy
                assert num_alleles == int(ploidy)
                df_w += float(num_alleles) - 1.0
            FS_SNij2_over_SNij += FSNij2temp / pop_an
            FS_SNij2 += pop_an**2  # Squared total number of observations per population
            RS_SNij2 += nnn**2  # Squared total number of observations per population
            p_ij.append(p)
            ac_ij.append(ac_i)

        p_bar = float(tac) / float(tan)

        # Calculate degrees of freedom and sample size coefficients
        df_g = r - 1.0
        df_i = sum([x - 1.0 for x in n_i])
        df_t = df_g + df_i + df_w
        fn0bis = (FS_SNij2_over_SNij - (FS_Nij2 / tan)) / df_g
        fn0 = (tan - FS_SNij2_over_SNij) / df_i
        fnb0 = (tan - (FS_SNij2 / tan)) / df_g
        rnb0 = (float(sum(n_i)) - (RS_SNij2 / float(sum(n_i)))) / df_g

        assert df_t == float(tan) - 1.0
        for i, pop in enumerate(ac_ij):
            for j, ind in enumerate(pop):
                for ref in range(0, int(int(ploidy_list[i]) - int(ind))):  # Loop over number of reference alleles
                    Fssg += (p_i[i] - p_bar)**2
                    Fssi += (p_ij[i][j] - p_i[i])**2
                    Fssw += (0 - p_ij[i][j])**2
                    Fsst += (0 - p_bar)**2
                for alt in range(0, int(ind)):  # Loop over number of alternative alleles
                    Fssg += (p_i[i] - p_bar)**2
                    Fssi += (p_ij[i][j] - p_i[i])**2
                    Fssw += (1 - p_ij[i][j])**2
                    Fsst += (1 - p_bar)**2
                Rssi += (p_ij[i][j] - p_i[i])**2
                Rssg += (p_i[i] - p_bar)**2
                Rsst += (p_ij[i][j] - p_bar)**2
        # Calculate Mean Squares
        FMS_g = Fssg / df_g
        FMS_i = Fssi / df_i
        FMS_w = Fssw / df_w
        RMS_g = Rssg / df_g
        RMS_i = Rssi / df_i
        # Calculate variance components
        fs2_w = FMS_w
        fs2_i = (FMS_i - fs2_w) / fn0
        fs2_g = (FMS_g - fs2_w - fn0bis * fs2_i) / fnb0
        rs2_i = RMS_i
        rs2_g = (RMS_g - rs2_i) / rnb0
        # Calculate numerator and denominators of rho and fst.
        rnum = rs2_g
        rden = rs2_i + rs2_g
        fnum = fs2_g
        fden = fs2_w + fs2_g + fs2_i

        if all(p == 0.0 for p in p_i) or all(p == 1.0 for p in p_i):
            poly = False
        else:
            poly = True

        return rnum, rden, fnum, fden, poly


    def calcDxy(locus_list):
        locus = []
        for l in locus_list:  # Remove missing data from the site
            locus.append([x for x in l if x != "-9"])
        for i, pop_site in enumerate(locus):
            ploidy = float(pop_site[1])
            nnn = float(len(pop_site[6:]))
            an = ploidy * nnn
            ac = sum([float(geno) for geno in pop_site[6:]])
            if i == 0:
                p1 = (ac / an)
            if i == 1:
                p2 = (ac / an)

        dxy = (p1 * (1.0 - p2)) + (p2 * (1.0 - p1))
        afd = abs(p1 - p2)
        return dxy, afd

    # Prepare output file
    outfile = output + outname + '_BPM.txt'
    out1 = open(outfile, 'w')
    out1.write("outname\tscaff\tstart\tend\twin_size\tnum_snps\tRho\tFst\tdxy\tFixedDiff\tAFD\n")

    # Begin loop over data file
    snp_count = 0
    Snp_count = 0
    start = 0.0
    end = window_size
    winexclcount = 0
    num_wind = 0
    with open(input_file, 'r') as data:
        for i, line in enumerate(data):
            line = line.strip("\n").strip("\t").split("\t")
            pop, ploidy, scaff, pos, an, dp = line[:6]
            pos = float(pos)

            if i % 100000 == 0:
                print(i)
            if i == 0:
                old_pos = pos
                Locus = []  # Hold information of single locus across populations
                oldscaff = scaff
                Fst = [0.0, 0.0]  # Fst[0] is numerator for genome, [1] is denominator for genome
                Rho = [0.0, 0.0]  # Rho[0] is numerator for genome, [1] is denominator for genome
                fst = [0.0, 0.0]  # Fst[0] is numerator for window, [1] is denominator for window
                rho = [0.0, 0.0]  # Rho[0] is numerator for window, [1] is denominator for window
                dxy = 0.0  # Dxy for window
                Dxy = 0.0
                dn = 0
                Dn = 0
                afd = 0.0
                AFD = 0.0

            if pos > start and pos <= end and scaff == oldscaff:
                if pos == old_pos:  # Accruing information from multiple populations but same locus
                    Locus.append(line)
                elif len(Locus) == num_pops:  # Within current window but have moved on from previous locus
                    rnum, rden, fnum, fden, poly = NestedAnova(Locus)
                    if poly is True:
                        if num_pops == 2:
                            d, a = calcDxy(Locus)
                            dxy += d
                            Dxy += d
                            afd += a
                            AFD += a
                            if a == 1.0:
                                dn += 1
                                Dn += 1
                        snp_count += 1
                        fst = [sum(x) for x in zip(fst, [fnum, fden])]  # Adds numerator and denominators from current site to running sums for window and genome respectively
                        rho = [sum(x) for x in zip(rho, [rnum, rden])]
                        Fst = [sum(x) for x in zip(Fst, [fnum, fden])]
                        Rho = [sum(x) for x in zip(Rho, [rnum, rden])]
                    Locus = []  # Clear site information
                    Locus.append(line)  # Append current site to site info
                    old_pos = pos

                else:  # Previous site contained data from only one population, so skip calculations
                    Locus = []
                    Locus.append(line)
                    old_pos = pos


            elif int(pos) > end or scaff != oldscaff:   # Current snp is onto next window, but must do calculation for previous locus before moving on
                if len(Locus) == num_pops:  # Skip locus calc if data not available for all populations
                    rnum, rden, fnum, fden, poly = NestedAnova(Locus)
                    if poly is True:
                        snp_count += 1
                        fst = [sum(x) for x in zip(fst, [fnum, fden])]
                        rho = [sum(x) for x in zip(rho, [rnum, rden])]
                        Fst = [sum(x) for x in zip(Fst, [fnum, fden])]
                        Rho = [sum(x) for x in zip(Rho, [rnum, rden])]
                        if num_pops == 2:
                            d, a = calcDxy(Locus)
                            dxy += d
                            Dxy += d
                            afd += a
                            AFD += a
                            if a == 1.0:
                                dn += 1
                                Dn += 1

                if snp_count >= minimum_snps:  # Report or exclude window
                    Snp_count += snp_count
                    num_wind += 1
                    try:
                        fst = fst[0] / fst[1]  # fst and rho and ratios of running sums of numerator and denominator calculations
                        fac = rho[0] / rho[1]

                        rho_i = fac / (1 + fac)
                        if num_pops == 2:
                            dxy = dxy / float(snp_count)
                            afd = afd / float(snp_count)
                        else:
                            dxy = "-9"
                            dn = "-9"
                            afd = "-9"
                        out1.write(outname + '\t' + scaff + '\t' +  # Write calculations to file
                                   str(start) + '\t' +
                                   str(end) + '\t' +
                                   str(window_size) + '\t' +
                                   str(snp_count) + '\t' +
                                   str(rho_i) + '\t' +
                                   str(fst) + '\t' +
                                   str(dxy) + '\t' +
                                   str(afd) + '\t' +
                                   str(dn) + '\n')
                    except ZeroDivisionError:
                        print('ZeroDivisionError:', snp_count, fst, rho, Fst, Rho)

                else:
                    winexclcount += 1
                # Reset running sums for current window
                snp_count = 0
                fst = [0.0, 0.0]
                rho = [0.0, 0.0]
                dxy = 0.0
                dn = 0
                afd = 0.0

                # Moving on to deal with current SNP.  Must reset window boundaries based on current position
                if float(pos) > end:
                    while float(pos) > end:
                        end += window_size  # Increment the end boundary by units of window-size until end is greater than current position
                    start = end - window_size
                elif scaff != oldscaff:
                    oldscaff = scaff

                    start = 0.0
                    end = window_size

                if int(pos) > start and int(pos) <= end and scaff == oldscaff:  # Ensures that current snp falls within window bounds we just determined
                    Locus = []
                    Locus.append(line)
                    old_pos = pos

    if len(Locus) == num_pops:  # Final window calculations
        rnum, rden, fnum, fden, poly = NestedAnova(Locus)
        if poly is True:
            snp_count += 1
            Snp_count += 1
            fst = [sum(x) for x in zip(fst, [fnum, fden])]
            rho = [sum(x) for x in zip(rho, [rnum, rden])]
            Fst = [sum(x) for x in zip(Fst, [fnum, fden])]
            Rho = [sum(x) for x in zip(Rho, [rnum, rden])]
            if num_pops == 2:
                d, a = calcDxy(Locus)
                dxy += d
                Dxy += d
                afd += a
                AFD += a
                if a == 1.0:
                    dn += 1
                    Dn += 1

    if snp_count >= minimum_snps:  # Use or exclude window
        num_wind += 1
        try:
            fst = fst[0] / fst[1]  # fst and rho and ratios of running sums of numerator and denominator calculations
            fac = rho[0] / rho[1]
            rho_i = fac / (1 + fac)
            if num_pops == 2:
                dxy = dxy / float(snp_count)
                afd = afd / float(snp_count)
            else:
                dxy = "-9"
                dn = "-9"
                afd = "-9"
            out1.write(outname + '\t' + scaff + '\t' +  # Write calculations to file
                       str(start) + '\t' +
                       str(end) + '\t' +
                       str(window_size) + '\t' +
                       str(snp_count) + '\t' +
                       str(rho_i) + '\t' +
                       str(fst) + '\t' +
                       str(dxy) + '\t' +
                       str(afd) + '\t' +
                       str(dn) + '\n')
        except ZeroDivisionError:
            print('ZeroDivisionError:', snp_count, fst, rho, Fst, Rho)

    else:
        winexclcount += 1
    # Calculate genome-wide metrics
    FAC = Rho[0] / Rho[1]
    rho_G = FAC / (1 + FAC)
    Fst_G = Fst[0] / Fst[1]
    if num_pops == 2:
        Dxy = Dxy / float(Snp_count)
        AFD = AFD / float(Snp_count)
    else:
        Dxy = "-9"
        Dn = "-9"
        AFD = "-9"
    out1.write(outname + '\t' + "Genome" + '\t' +
               "-9" + '\t' +
               "-9" + '\t' +
               str(window_size) + '\t' +
               str(Snp_count) + '\t' +
               str(rho_G) + '\t' +
               str(Fst_G) + '\t' +
               str(AFD) + '\t' +
               str(Dn) + '\n')

    out1.close()

    return num_wind, winexclcount, rho


if __name__ == '__main__':  # Used to run code from command line

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar='input_file', required=True, help='input file created with recode012.py')
    parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Output Directory')
    parser.add_argument('-prefix', type=str, metavar='output_file_prefix', required=True, help='Name indicating populations in input file')
    parser.add_argument('-ws', type=float, metavar='window_size', required=False, default='10000.0', help='Size of windows in bp')
    parser.add_argument('-ms', type=int, metavar='minimum_snps', required=False, default='2', help='minimum number of snps in a window')
    parser.add_argument('-np', type=int, metavar='number_populations', required=False, default='2', help='Number of populations')

    args = parser.parse_args()

    j1, j2, j3 = calcBPM(args.i, args.o, args.prefix, args.ws, args.ms, args.np)
