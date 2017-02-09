#!/usr/bin/env python35

import os
import subprocess
import argparse
import pandas
import math
import statistics
import numpy as np

# Directions:
# import WPM_Wrapper2 into python console
# and pass path of Individual_Key file as object to class PopGen

# USE BELOW COMMANDS TO CONCAT POP FILES after all are finished running
# head -1 BEL_WPM.txt > WPM_All.txt
# tail -n+2 -q *_WPM.txt >> WPM_All.txt

# create variables that can be entered in the command line


class scantools:

    def __init__(self, WorkingDir):
        if WorkingDir.endswith("/") is False:
            WorkingDir += "/"
        if os.path.exists(WorkingDir) is False:
            os.mkdir(WorkingDir)
            os.mkdir(WorkingDir + "OandE/")
        if os.path.exists(WorkingDir + "OandE/") is False:
            os.mkdir(WorkingDir + "OandE/")
        POP_file = pandas.read_csv(WorkingDir + "PopKey.csv", header=0)
        POP_names = list(POP_file.columns.values)[1:]
        sample_names = list(POP_file['Samples'])
        samps = {}
        samp_nums = {}
        for pop in POP_names:
            pop_list = []
            include_index = list(POP_file[pop])
            for i, sample in enumerate(sample_names):
                if include_index[i] == 1:
                    pop_list.append(sample)
            samps[pop] = pop_list
            samp_nums[pop] = len(pop_list)


        # Determine number of individuals to downsample all populations to
        min_ind = min([sum(list(POP_file[pop])) for pop in POP_names])
        self.pops = POP_names
        self.samps = samps
        self.samp_nums = samp_nums
        self.min_ind = min_ind
        self.dir = WorkingDir
        self.oande = WorkingDir + "OandE/"
        self.code_dir = os.getcwd()

        self.split_dirs = []
        for path in os.listdir(self.dir):
            if path.split("_")[0] == "VCF":
                self.split_dirs.append(self.dir + path)


    def removePop(self, popname):
        '''Purpose: remove population from all object and recalculate min_ind'''
        popname = str(popname)
        if popname in self.pops:
            self.pops.remove(popname)
            self.samps.pop(popname, None)
            self.samp_nums.pop(popname, None)
            # Recalculate min_ind
            min_ind = min([self.samp_nums[pop] for pop in self.pops])
            self.min_ind = min_ind
        else:
            print("Population does not exist")


    def combinePops(self, pops, popname):
        new_samps = []
        for pop in pops:
            for samp in self.samps[pop]:
                new_samps.append(samp)
        self.pops.append(popname)
        self.samps[popname] = new_samps
        self.samp_nums[popname] = len(new_samps)


    def splitVCFs(self, vcf_dir, ref_path, min_dp, mffg, repolarization_key="-99", pops='all', mem=16000, time='0-04:00', numcores=1, print1=False, overwrite=False, partition="long", keep_intermediates=False):
        '''Purpose:  Find all vcfs in vcf_dir and split them by population according to samples associated with said population.
                    Then, take only biallelic snps and convert vcf to table containing scaff, pos, ac, an, dp, and genotype fields.
                    Finally, concatenate all per-scaffold tables to one giant table. Resulting files will be put into ~/Working_Dir/VCFs/
            Notes: mffg is maximum fraction of filtered genotypes.  Number of actual genotypes allowed will be rounded up.
                    If you want to print the batch scripts, you must set print1 and overwrite to True'''

        if vcf_dir.endswith("/") is False:
            vcf_dir += "/"
        vcf_dir_name = vcf_dir.split("/")[-2]
        outdir = self.dir + "VCF_" + str(vcf_dir_name) + "_DP" + str(min_dp) + ".M" + str(mffg) + "/"
        self.vcf_dir = vcf_dir
        if outdir not in self.split_dirs:
            self.split_dirs.append(outdir)

        mem1 = int(mem / 1000)

        if os.path.exists(outdir) is True and overwrite is False:
            print("VCF directory already exists.  Set 'overwrite = True' if you want to overwrite existing files")
        else:
            if os.path.exists(outdir) is False:
                os.mkdir(outdir)
            elif overwrite is True:
                print("Overwriting files in existing VCF directory")

            if pops == 'all':
                pops = self.pops

            for pop in pops:
                # Add samples to list for each population according to PF file
                sample_string1 = ""
                for samp in self.samps[pop]:
                    sample_string1 += " -sn " + samp
                joblist = []

                mfg = int(math.ceil(float(self.samp_nums[pop]) * float(mffg)))

                vcf_list = []
                vcf_basenames = []
                for file in os.listdir(vcf_dir):
                    if file[-6:] == 'vcf.gz':
                        vcf_list.append(file)
                        vcf_basenames.append(file[:-7])
                    elif file[-3:] == 'vcf':
                        vcf_list.append(file)
                        vcf_basenames.append(file[:-4])
                for v, vcf in enumerate(vcf_list):
                    # Select single population and biallelic SNPs for each scaffold and convert to variants table
                    shfile1 = open(pop + '.sh', 'w')
                    shfile1.write('#!/bin/bash\n' +
                                  '#SBATCH -J ' + pop + '.sh' + '\n' +
                                  '#SBATCH -e ' + self.oande + pop + vcf + '.gatk.err' + '\n' +
                                  '#SBATCH -o ' + self.oande + pop + vcf + '.gatk.out' + '\n' +
                                  '#SBATCH -p nbi-' + str(partition) + '\n' +
                                  '#SBATCH -n ' + str(numcores) + '\n' +
                                  '#SBATCH -t ' + str(time) + '\n' +
                                  '#SBATCH --mem=' + str(mem) + '\n' +
                                  'source GATK-nightly.2016.09.26\n' +
                                  'java -Xmx' + str(mem1) + 'g -jar /nbi/software/testing/GATK/nightly.2016.09.26/x86_64/jars/GenomeAnalysisTK.jar -T SelectVariants -R ' + ref_path + ' -V ' + vcf_dir + vcf + sample_string1 + ' -o ' + outdir + vcf_basenames[v] + '.' + pop + '.vcf\n' +
                                  'gunzip ' + outdir + vcf_basenames[v] + '.' + pop + '.vcf.gz\n' +
                                  'java -Xmx' + str(mem1) + 'g -jar /nbi/software/testing/GATK/nightly.2016.09.26/x86_64/jars/GenomeAnalysisTK.jar -T VariantFiltration -R ' + ref_path + ' -V ' + outdir + vcf_basenames[v] + '.' + pop + '.vcf --genotypeFilterExpression "DP < ' + str(min_dp) + '" --genotypeFilterName "DP" -o ' + outdir + vcf_basenames[v] + '.' + pop + '.dp' + str(min_dp) + '.1.vcf\n' +
                                  'java -Xmx' + str(mem1) + 'g -jar /nbi/software/testing/GATK/nightly.2016.09.26/x86_64/jars/GenomeAnalysisTK.jar -T VariantFiltration -R ' + ref_path + ' -V ' + outdir + vcf_basenames[v] + '.' + pop + '.dp' + str(min_dp) + '.1.vcf --setFilteredGtToNocall -o ' + outdir + vcf_basenames[v] + '.' + pop + '.dp' + str(min_dp) + '.vcf\n' +
                                  'java -Xmx' + str(mem1) + 'g -jar /nbi/software/testing/GATK/nightly.2016.09.26/x86_64/jars/GenomeAnalysisTK.jar -T SelectVariants -R ' + ref_path + ' -V ' + outdir + vcf_basenames[v] + '.' + pop + '.dp' + str(min_dp) + '.vcf --maxNOCALLnumber ' + str(mfg) + ' -o ' + outdir + vcf_basenames[v] + '.' + pop + '.m' + str(mffg) + '.dp' + str(min_dp) + '.bi.vcf\n' +
                                  'java -Xmx' + str(mem1) + 'g -jar /nbi/software/testing/GATK/nightly.2016.09.26/x86_64/jars/GenomeAnalysisTK.jar -T VariantsToTable -R ' + ref_path + ' -V ' + outdir + vcf_basenames[v] + '.' + pop + '.m' + str(mffg) + '.dp' + str(min_dp) + '.bi.vcf -F CHROM -F POS -F AC -F AN -F DP -GF GT -o ' + outdir + vcf_basenames[v] + '.' + pop + '_raw.table\n')

                    shfile1.close()

                    if print1 is False:  # send slurm job to NBI SLURM cluster
                        cmd1 = ('sbatch ' + pop + '.sh')
                        p1 = subprocess.Popen(cmd1, shell=True)
                        sts1 = os.waitpid(p1.pid, 0)[1]
                        joblist.append(p1.pid)

                    else:
                        file1 = open(pop + '.sh', 'r')
                        data1 = file1.read()
                        print(data1)

                    os.remove(pop + '.sh')

                # combine all variants table for each scaffold within a population

                shfile3 = open(pop + '.sh', 'w')

                shfile3.write('#!/bin/bash\n' +
                              '#SBATCH -J ' + pop + '.sh' + '\n' +
                              '#SBATCH -e ' + self.oande + pop + '.cat.err' + '\n' +
                              '#SBATCH -o ' + self.oande + pop + '.cat.out' + '\n' +
                              '#SBATCH -p nbi-medium\n' +
                              '#SBATCH -n ' + str(numcores) + '\n' +
                              '#SBATCH -t 0-12:00\n' +
                              '#SBATCH --mem=' + str(mem) + '\n' +
                              'source python-3.5.1\n' +
                              'cat ' + outdir + '*' + pop + '_raw.table | tail -n+2 > ' + outdir + pop + '.table\n' +
                              'python3 ' + self.code_dir + '/recode012.py -i ' + outdir + pop + '.table -pop ' + pop + ' -o ' + outdir + '\n')
                if repolarization_key == "-99":
                    print("No repolarization key provided.  Repolarized input files will not be produced.  Must set 'use_repol' to True in subsequent steps")
                else:
                    shfile3.write('python3 ' + self.code_dir + '/repol.py -i ' + outdir + pop + '.table.recode.txt -o ' + outdir + pop + ' -r ' + repolarization_key + '\n')
                    if keep_intermediates is False:
                        shfile3.write('rm ' + outdir + pop + '.table.recode.txt\n')
                if keep_intermediates is False:
                    shfile3.write('rm ' + outdir + '*.' + pop + '.dp' + str(min_dp) + '.vcf\n')
                    shfile3.write('rm ' + outdir + '*.' + pop + '.dp' + str(min_dp) + '.vcf.idx\n')
                    shfile3.write('rm ' + outdir + '*.' + pop + '.m' + str(mffg) + '.dp' + str(min_dp) + '.bi.vcf\n')
                    shfile3.write('rm ' + outdir + '*.' + pop + '.m' + str(mffg) + '.dp' + str(min_dp) + '.bi.vcf.idx\n')
                    shfile3.write('rm ' + outdir + '*.' + pop + '.vcf\n')
                    shfile3.write('rm ' + outdir + '*.' + pop + '.vcf.idx\n')
                    shfile3.write('rm ' + outdir + '*.' + pop + '.dp' + str(min_dp) + '.1.vcf\n')
                    shfile3.write('rm ' + outdir + '*.' + pop + '.dp' + str(min_dp) + '.1.vcf.idx\n')
                    shfile3.write('rm ' + outdir + '*' + pop + '_raw.table\n' +
                                  'rm ' + outdir + pop + '.table\n')
                shfile3.close()

                if print1 is False:
                    cmd3 = ('sbatch -d singleton ' + pop + '.sh')
                    p3 = subprocess.Popen(cmd3, shell=True)
                    sts3 = os.waitpid(p3.pid, 0)[1]
                else:
                    file3 = open(pop + '.sh', 'r')
                    data3 = file3.read()
                    print(data3)

                os.remove(pop + '.sh')

    def recode(self, split_dir, pops="all", print1=False, mem=4000, numcores=1, partition="medium"):
        '''Call: recode(self, min_avg_dp, missingness, print1=False, mem=16000, numcores=1)
           Purpose: Take the concatenated table files in ~/Working_Dir/VCFs/ and recode them so that genotypes are represented as number of alternative alleles
           Notes: split_dir is the directory of split vcfs generated by splitVCFs().  Also, adds population name and ploidy as columns in file.  Files are output to a folder named ~/Working_Dir/Recoded.DPXX.MX.X/ where the Xs are user specified
                    These directories will be primary depositories for all subsequent commands.'''

        # Determine if the recoded vcf files already exist and if so, set VCF_Parse to False
        if split_dir.endswith("/"):
            recode_dir = split_dir[:-1]
            recode_dir += ".Recoded/"
        else:
            recode_dir = split_dir + ".Recoded/"
            split_dir += "/"

        self.recode_dirs.append(recode_dir)

        if os.path.exists(split_dir) is True:
            existing_files = []
            if os.path.exists(recode_dir) is False:
                os.mkdir(recode_dir)
            else:
                for file in os.listdir(recode_dir):
                    if file.endswith('.table.recode.txt'):
                        existing_files.append(file.split('.')[0])
            if set(self.pops).issuperset(set(existing_files)) is True and len(existing_files) != 0:
                print("Recoded vcf files already exist.  Delete folder or change parameters")
            # Look for '.recode' table files corresponding to POP_names...if they exist...skip all but wpm step and print message to output.

            else:
                if pops == "all":
                    pops = self.pops
                for pop in pops:
                    # FORMAT TABLE FOR WPM FILE.

                    shfile3 = open(pop + '.sh', 'w')
                    shfile3.write('#!/bin/bash\n' +
                                  '#SBATCH -J ' + pop + '.sh' + '\n' +
                                  '#SBATCH -e ' + self.oande + pop + '.recode012.err' + '\n' +
                                  '#SBATCH -o ' + self.oande + pop + '.recode012.out' + '\n' +
                                  '#SBATCH -p nbi-' + str(partition) + '\n' +
                                  '#SBATCH -n ' + str(numcores) + '\n' +
                                  '#SBATCH -t 1-00:00\n' +
                                  '#SBATCH --mem=' + str(mem) + '\n' +
                                  'source python-3.5.1\n' +
                                  'source env/bin/activate\n' +
                                  'python3 ' + self.code_dir + '/recode012.py -i ' + split_dir + pop + '.table -pop ' + pop + ' -o ' + recode_dir + '\n')
                    shfile3.close()

                    if print1 is False:
                        cmd3 = ('sbatch -d singleton ' + pop + '.sh')
                        p3 = subprocess.Popen(cmd3, shell=True)
                        sts3 = os.waitpid(p3.pid, 0)[1]
                    else:
                        file3 = open(pop + '.sh', 'r')
                        data3 = file3.read()
                        print(data3)

                    os.remove(pop + '.sh')
        else:
            print("Must run splitVCFs prior to running recode")


    def repolarize(self, recode_dir, repolarization_key, pops="all", time='0-02:00', mem='8000', numcores=1):
        '''Call: repolarize(self, recode_dir, in_file, repolarization_key)
           Purpose: repolarize the .recode.txt files in recode_dir according to a key generated by repolarization_lookupKey.py
           Notes: This is necessary for generateFSC2input as this method looks for files named .repol.txt'''
        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if pops == 'all':
            pops = self.pops

        for pop in pops:
            if os.path.exists(recode_dir + pop + '.table.recode.txt') is True:
                shfile4 = open(pop + '.sh', 'w')
                shfile4.write('#!/bin/bash\n' +
                              '#SBATCH -J ' + pop + '.repol.sh' + '\n' +
                              '#SBATCH -e ' + self.oande + pop + '.repol.err' + '\n' +
                              '#SBATCH -o ' + self.oande + pop + '.repol.out' + '\n' +
                              '#SBATCH -p nbi-short\n' +
                              '#SBATCH -n ' + str(numcores) + '\n' +
                              '#SBATCH -t ' + str(time) + '\n' +
                              '#SBATCH --mem=' + str(mem) + '\n' +
                              'source python-3.5.1\n' +
                              'source env/bin/activate\n' +
                              'python3 ' + self.code_dir + '/repol.py -i ' + recode_dir + pop + '.table.recode.txt -o ' + recode_dir + pop + ' -r ' + repolarization_key + '\n')
                shfile4.close()

                cmd1 = ('sbatch -d singleton ' + pop + '.sh')
                p1 = subprocess.Popen(cmd1, shell=True)
                sts1 = os.waitpid(p1.pid, 0)[1]

                os.remove(pop + '.sh')

            else:
                print("Did not find .table.recode.txt file for population: ", pop)



    def getPloidies(self, recode_dir):
        '''Purpose: Create new methods of scantools object containing ploidy of each population (.ploidies) as well as a list of dips (.dips) and tetraploid populations (.tets)
           Notes: Can only be executed after recode has been executed on vcfs'''

        print("Be sure that 'recode' scripts have all finished")
        if recode_dir.endswith("/") is False:
            recode_dir += "/"
        ploidies = {}
        dips = []
        tets = []
        if os.path.exists(recode_dir) is True:
            for pop in self.pops:
                try:
                    tmp = open(recode_dir + pop + '.table.recode.txt', 'r')
                    line = tmp.readline()
                    ploidy = line.split("\t")[1]
                    ploidies[pop] = ploidy
                    if ploidy == "4.0":
                        tets.append(pop)
                    elif ploidy == "2.0":
                        dips.append(pop)
                    else:
                        print("Ploidy level not recognized")
                except (FileNotFoundError, IndexError):
                    print("Error determining ploidy for population: ", pop)
            self.ploidies = ploidies
            self.dips = dips
            self.tets = tets
        else:
            print("recode_dir does not exist")


    # CALCULATE WITHIN POPULATION METRICS
    def calcwpm(self, recode_dir, window_size, min_snps, pops="all", print1=False, mem=16000, numcores=1, sampind="-99", partition="medium", use_repol=True):
        '''Purpose: Calculate within population metrics including: allele frequency, expected heterozygosity, Wattersons theta, Pi, ThetaH, ThetaL and neutrality tests: D, normalized H, E
           Notes:  Currently, all populations are downsampled to same number of individuals.  By default, this minimum individuals across populations minus 1 to allow for some missing data
                    It is worth considering whether downsampling should be based on number of individuals or number of alleles.
                    Results are held ~/Working_Dir/Recoded/ in series of files ending in _WPM.txt.  These can be concatenated using concatWPM'''

        if use_repol is True:
            suffix = '.table.repol.txt'
        else:
            suffix = '.table.recode.txt'

        if sampind == "-99":
            sind = self.min_ind - 1
        else:
            sind = sampind

        if pops == "all":
            pops = self.pops
        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True and sind > 3:

            for pop in pops:
                shfile3 = open(pop + '.sh', 'w')

                shfile3.write('#!/bin/bash\n' +
                              '#SBATCH -J ' + pop + '.sh' + '\n' +
                              '#SBATCH -e ' + self.oande + pop + '.wpm.err' + '\n' +
                              '#SBATCH -o ' + self.oande + pop + '.wpm.out' + '\n' +
                              '#SBATCH -p nbi-' + str(partition) + '\n' +
                              '#SBATCH -n ' + str(numcores) + '\n' +
                              '#SBATCH -t 1-00:00\n' +
                              '#SBATCH --mem=' + str(mem) + '\n' +
                              'source python-3.5.1\n' +
                              'python3 ' + self.code_dir + '/wpm.py -i ' + recode_dir + pop + suffix + ' -o ' + recode_dir + ' -sampind ' + str(sind) + ' -ws ' + str(window_size) + ' -ms ' + str(min_snps) + '\n')
                shfile3.close()

                if print1 is False:
                    cmd3 = ('sbatch -d singleton ' + pop + '.sh')
                    p3 = subprocess.Popen(cmd3, shell=True)
                    sts3 = os.waitpid(p3.pid, 0)[1]
                else:
                    file3 = open(pop + '.sh', 'r')
                    data3 = file3.read()
                    print(data3)

                os.remove(pop + '.sh')
                pops.remove(pop)
            for pop in pops:
                print("Did not find input files for: ", pop)

        elif sind <= 3:
            print("Number of individuals to be used/downsampled to is <= 3.  Unable to calculate within-population-metrics on so few individuals.")
        else:
            print("Did not find recode_dir.  Must run splitVCFs followed by recode before able to calculate within population metrics")

    def concatWPM(self, recode_dir, pops='all'):
        '''Purpose:  Concatenate _WPM.txt files corresponding to populations indicated in pops parameter.'''

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True:
            new = open(recode_dir + "All_WPM.txt", 'w')
            if pops == 'all':
                pops = self.pops
            for i, pop in enumerate(pops):
                try:
                    with open(recode_dir + pop + "_WPM.txt", 'r') as inf:
                        for j, line in enumerate(inf):
                            if j == 0 and i == 0:
                                new.write(line)
                            else:
                                new.write(line)
                except FileNotFoundError:
                    print("Did not find _WPM.txt file for population: ", pop)


    def calcbpm(self, recode_dir, pops, output_name, window_size, minimum_snps, print1=False, mem=16000, numcores=1, partition="medium", use_repol=True, keep_intermediates=False):
        '''Purpose:  Calculate between population metrics including: Dxy, Fst (using Weir and Cockerham 1984), and Rho (Ronfort et al. 1998)
           Notes: User provides a list of populations to be included.  For pairwise estimates, simply provide two populations
                    Calculations are done for windows of a given bp size.  User also must specify the minimum number of snps in a window
                    for calculations to be made'''

        if use_repol is True:
            suffix = '.table.repol.txt'
        else:
            suffix = '.table.recode.txt'

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True and len(pops) > 1:

            # Concatenate input files and sort them
            print("Concatenating input files")
            concat_file = open(recode_dir + output_name + '.concat.txt', 'w')
            pop_num = 0
            for pop in pops:  # Add data from all populations to single, huge listg
                try:
                    with open(recode_dir + pop + suffix, 'r') as in1:
                        for line in in1:
                            concat_file.write(line)
                    pop_num += 1
                except IOError:
                    print("Did not find input file for pop ", pop)
            print("Finished preparing input data")
            if len(pops) != pop_num:
                print("Did not find all input files!!  Aborting.")
                os.remove(recode_dir + output_name + '.concat.txt')
            else:
                shfile3 = open(output_name + '.bpm.sh', 'w')

                shfile3.write('#!/bin/bash\n' +
                              '#SBATCH -J ' + output_name + '.bpm.sh' + '\n' +
                              '#SBATCH -e ' + self.oande + output_name + '.bpm.err' + '\n' +
                              '#SBATCH -o ' + self.oande + output_name + '.bpm.out' + '\n' +
                              '#SBATCH -p nbi-' + str(partition) + '\n' +
                              '#SBATCH -n ' + str(numcores) + '\n' +
                              '#SBATCH -t 1-00:00\n' +
                              '#SBATCH --mem=' + str(mem) + '\n' +
                              'source python-3.5.1\n' +
                              'python3 ' + self.code_dir + '/bpm.py -i ' + recode_dir + output_name + '.concat.txt' + ' -o ' + recode_dir + ' -prefix ' + output_name + ' -ws ' + str(window_size) + ' -ms ' + str(minimum_snps) + ' -np ' + str(pop_num) + '\n')
                if keep_intermediates is False:
                    shfile3.write('rm ' + recode_dir + output_name + '.concat.txt')
                shfile3.close()

                if print1 is False:
                    cmd3 = ('sbatch -d singleton ' + output_name + '.bpm.sh')
                    p3 = subprocess.Popen(cmd3, shell=True)
                    sts3 = os.waitpid(p3.pid, 0)[1]
                else:
                    file3 = open(output_name + '.bpm.sh', 'r')
                    data3 = file3.read()
                    print(data3)

                os.remove(output_name + '.bpm.sh')
        elif len(pops) < 2:
            print("'pops' argument must be a list of strings specifiying two or more population names as they appear in input file prefixes.  len(pops) was < 2")
        else:
            print("Did not find recode_dir.  Must run splitVCFs followed by recode before able to calculate between population metrics")

    def calcPairwisebpm(self, recode_dir, pops, window_size, minimum_snps, print1=False, mem=16000, numcores=1, partition="medium", use_repol=True, keep_intermediates=False):
        '''Purpose:  Calculate between population metrics including: Dxy, Fst (using Weir and Cockerham 1984), and Rho (Ronfort et al. 1998)
           Notes: User provides a list of populations to be included.  For pairwise estimates, simply provide two populations
                    Calculations are done for windows of a given bp size.  User also must specify the minimum number of snps in a window
                    for calculations to be made'''

        if use_repol is True:
            suffix = '.table.repol.txt'
        else:
            suffix = '.table.recode.txt'

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True and len(pops) > 1:
            # Concatenate input files and sort them
            for i, pop1 in enumerate(pops):  # Add data from all populations to single, huge listg
                for pop2 in pops[i + 1:]:
                    print(pop1,pop2)
                    output_name = pop1 + pop2
                    concat_file = open(recode_dir + output_name + '.concat.txt', 'w')
                    skip = False
                    try:
                        with open(recode_dir + pop1 + suffix, 'r') as in1:
                            for line in in1:
                                concat_file.write(line)
                    except IOError:
                        print("Did not find input file for pop ", pop1)
                        skip = True
                    try:
                        with open(recode_dir + pop2 + suffix, 'r') as in1:
                            for line in in1:
                                concat_file.write(line)
                    except IOError:
                        print("Did not find input file for pop ", pop2)
                        skip = True
                    if skip is True:
                        print("Did not find all input files!!  Aborting pairwise bpm for contrast: ", output_name)
                        os.remove(recode_dir + output_name + '.concat.txt')
                    else:
                        shfile3 = open(output_name + '.bpm.sh', 'w')

                        shfile3.write('#!/bin/bash\n' +
                                      '#SBATCH -J ' + output_name + '.bpm.sh' + '\n' +
                                      '#SBATCH -e ' + self.oande + output_name + '.bpm.err' + '\n' +
                                      '#SBATCH -o ' + self.oande + output_name + '.bpm.out' + '\n' +
                                      '#SBATCH -p nbi-' + str(partition) + '\n' +
                                      '#SBATCH -n ' + str(numcores) + '\n' +
                                      '#SBATCH -t 1-00:00\n' +
                                      '#SBATCH --mem=' + str(mem) + '\n' +
                                      'source python-3.5.1\n' +
                                      'python3 ' + self.code_dir + '/bpm.py -i ' + recode_dir + output_name + '.concat.txt' + ' -o ' + recode_dir + ' -prefix ' + output_name + ' -ws ' + str(window_size) + ' -ms ' + str(minimum_snps) + ' -np 2\n')
                        if keep_intermediates is False:
                            shfile3.write('rm ' + recode_dir + output_name + '.concat.txt')
                        shfile3.close()

                        if print1 is False:
                            cmd3 = ('sbatch -d singleton ' + output_name + '.bpm.sh')
                            p3 = subprocess.Popen(cmd3, shell=True)
                            sts3 = os.waitpid(p3.pid, 0)[1]
                        else:
                            file3 = open(output_name + '.bpm.sh', 'r')
                            data3 = file3.read()
                            print(data3)

                        os.remove(output_name + '.bpm.sh')
        elif len(pops) < 2:
            print("'pops' argument must be a list of strings specifiying two or more population names as they appear in input file prefixes.  len(pops) was < 2")
        else:
            print("Did not find recode_dir.  Must run splitVCFs followed by recode before able to calculate between population metrics")


    def findOutliers(self, recode_dir, in_file, column_index_list, percentile, tails='upper'):
        '''Purpose:  Take output from either calcwpm or calcbpm and determine outlier metrics for a given percentile.
           Notes: Output will be two csv files (one containing all sites with outliers indicated by 0 or 1 and another containing just outliers)
                  as well as a bed file to be used in annotateOutliers'''

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if percentile > 1.0:
            print("!!percentile parameter should be coded as a proportion, not as a percentage!!")
            return False

        if os.path.exists(recode_dir) is True:

            data = pandas.read_table(recode_dir + in_file, header=0)
            metrics = []
            for i in column_index_list:
                metrics.append(list(data.columns.values)[i])

            for metric in metrics:
                data[metric + '.out'] = 0
                if tails == 'both':
                    data[metric + '.out'].loc[(data[metric] > data[metric].quantile(q=percentile))] = 1
                    data[metric + '.out'].loc[(data[metric] < data[metric].quantile(q=1.0 - percentile))] = 1
                elif tails == 'lower':
                    data[metric + '.out'].loc[(data[metric] < data[metric].quantile(q=1.0 - percentile))] = 1
                elif tails == 'upper':
                    data[metric + '.out'].loc[(data[metric] > data[metric].quantile(q=percentile))] = 1
                else:
                    print("Did not specify tails option correctly.  Options are: both, upper, and lower")
            data['num_outliers'] = data.iloc[:, -len(metrics):].sum(1)
            data.to_csv(recode_dir + in_file.replace(".txt", "") + '_' + str(percentile) + 'tile_OutLabelled.csv', index=False)
            # select all windows that are outliers for at least one metric
            df_outlier = data[(data.num_outliers != 0)]
            df_outlier.to_csv(recode_dir + in_file.replace(".txt", "") + '_' + str(percentile) + 'tile_OutOnly.csv', index=False)
            df_outlier.to_csv(recode_dir + in_file.replace(".txt", "") + '_' + str(percentile) + 'tile_OutOnly.bed', index=False, sep='\t', columns=["scaffold", "start", "end"], header=False)
        return data

    def annotateOutliers(self, recode_dir, in_file, basename, annotation_file, overlap_proportion=0.000001):
        '''Purpose: annotate bed file from findOutliers using information in annotation_file
           Notes: The output (suffix ol_genes.gff) only contains the window locations along with annotation info and does not contain
                    the original metric information used to determine outliers.  Use mergeAnnotation to merge original outlier file with annotation info'''

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True:
            shfile1 = open(recode_dir + 'bedtools_gff.sh', 'w')
            shfile1.write('#!/bin/bash\n' +
                          '#SBATCH -J GS.bedtools.sh' + '\n' +
                          '#SBATCH -e GS.bedtools.err\n' +
                          '#SBATCH -o GS.bedtools.out\n' +
                          '#SBATCH -p nbi-short\n' +
                          '#SBATCH -n 1\n' +
                          '#SBATCH -t 0-02:00\n' +
                          '#SBATCH --mem=16000\n' +
                          'source bedtools-2.17.0\n' +
                          'bedtools intersect -a ' + recode_dir + in_file + ' -b ' + annotation_file + ' -f ' + str(overlap_proportion) + ' -wb | ' +
                          """awk '{$1=$2=$3=""; print $4,$5,$6,$7,$8,$9,$10,$11,$12}'""" +
                          '| grep transcript | grep -v transcription | sort -u | ' +
                          """tr ' ' '\t' """
                          '> ' + recode_dir + basename + '_' + str(overlap_proportion * 100) + 'ol_genes.gff')
            shfile1.close()

            cmd1 = ('sbatch ' + recode_dir + 'bedtools_gff.sh')
            p1 = subprocess.Popen(cmd1, shell=True)
            sts1 = os.waitpid(p1.pid, 0)[1]

            os.remove(recode_dir + 'bedtools_gff.sh')

        else:
            print("recode_dir not found")


    def mergeAnnotation(self, recode_dir, outlier_file, annotated_outlier_file):
        '''Purpose: Merge the annotation information with the original outlier file results from findOutliers.'''

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True:
            try:
                outliers = pandas.read_table(recode_dir + outlier_file, header=0)
                annotation = pandas.read_table(recode_dir + annotated_outlier_file, names=["scaffold", "start", "end", "info1", "info2", "info3", "info4", "info5", "info6", "info7", "info8", "info9"])
            except IOError:
                print("Did not find either original outlier file or the annotated outlier file")
            merged = pandas.merge(outliers, annotation, ["scaffold", "start", "end"],)
            merged.to_csv(recode_dir + outlier_file.replace(".txt", "") + '_OutAnnot.csv', index=False)
        else:
            print("Did not find recode_dir")


    def generateFSC2input(self, recode_dir, pops, output_name, bootstrap_block_size, bootstrap_reps, mem=16000, numcores=1, time='2-00:00'):
        '''Purpose:  Generate --multiSFS for fastsimcoal2 along with a given number of non-parametric block-bootstrapped replicates
           Notes: Must provide the block size for bootstrapping as well as number of bootstrap replicates
                  As of now, the necessary template files for FSC2 must come from elsewhere.  Beware of running this method with numerous populations'''

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        outdir = recode_dir + "FSC2input_" + output_name + "/"

        if os.path.exists(outdir) is False:
            os.mkdir(outdir)

        num_pops = len(pops)

        # num_inds = [self.samp_nums.get(x) for x in pops]
        if os.path.exists(recode_dir) is True:
            print("Concatenating input files")
            concat_file = open(recode_dir + output_name + '.repol.concat.txt', 'w')
            for file in os.listdir(recode_dir):
                if file.endswith(".repol.txt") and file.split(".")[0] in pops:
                    pops.remove(file.split(".")[0])
                    with open(recode_dir + file) as infile:
                        for line in infile:
                            concat_file.write(line)
            if len(pops) != 0:
                print("Did not find repolarized files for the following populations ", pops, ".  Aborting!!")
            else:
                print("Finished preparing input data")

                shfile4 = open(output_name + '.fsc2input.sh', 'w')

                shfile4.write('#!/bin/bash\n' +
                              '#SBATCH -J ' + output_name + '.fsc2input.sh' + '\n' +
                              '#SBATCH -e ' + self.oande + output_name + '.fsc2input.err' + '\n' +
                              '#SBATCH -o ' + self.oande + output_name + '.fsc2input.out' + '\n' +
                              '#SBATCH -p nbi-medium\n' +
                              '#SBATCH -n ' + str(numcores) + '\n' +
                              '#SBATCH -t ' + str(time) + '\n' +
                              '#SBATCH --mem=' + str(mem) + '\n' +
                              'source python-3.5.1\n' +
                              'source env/bin/activate\n' +
                              'python3 ' + self.code_dir + '/FSC2input.py -i ' + recode_dir + output_name + '.repol.concat.txt -o ' + outdir + ' -prefix ' + output_name + ' -ws ' + str(bootstrap_block_size) + ' -bs ' + str(bootstrap_reps) + ' -np ' + str(num_pops) + '\n')
                shfile4.close()

                cmd1 = ('sbatch ' + output_name + '.fsc2input.sh')
                p1 = subprocess.Popen(cmd1, shell=True)
                sts1 = os.waitpid(p1.pid, 0)[1]

                os.remove(output_name + '.fsc2input.sh')

        else:
            print("!!!Did not find recode_dir!!!!")
