#!/usr/bin/env python35
import os
import subprocess
import pandas
import math
import datetime
import time


class scantools:

    def __init__(self, WorkingDir, encoding="-9"):
        if WorkingDir.endswith("/") is False:
            WorkingDir += "/"
        if os.path.exists(WorkingDir) is False:
            os.mkdir(WorkingDir)
            os.mkdir(WorkingDir + "OandE/")
        if os.path.exists(WorkingDir + "OandE/") is False:
            os.mkdir(WorkingDir + "OandE/")
        if encoding != "-9":
            POP_file = pandas.read_csv(WorkingDir + "PopKey.csv", header=0, encoding=encoding)
        else:
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
        self.log_file = open(WorkingDir + "LogFile.txt", 'a+')
        time = 'New Instance at: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()) + "\n"
        self.log_file.write(time)
        self.split_dirs = []
        for path in os.listdir(self.dir):
            if path.split("_")[0] == "VCF":
                self.split_dirs.append(self.dir + path)
        print("Be sure that you did 'unset SBATCH_PARTITION' prior to using ScanTools!!!")


    def removePop(self, pops_to_be_removed):
        '''Purpose: remove population from all object and recalculate min_ind'''
        if isinstance(pops_to_be_removed, list) is False:
            pops = [pops_to_be_removed]
        for popname in pops:
            popname = str(popname)
            if popname in self.pops:
                self.pops.remove(popname)
                self.samps.pop(popname, None)
                self.samp_nums.pop(popname, None)
                self.log_file.write("Removed Pop: " + popname + "\n")
            else:
                print("Population does not exist")
        min_ind = min([self.samp_nums[pop] for pop in self.pops])
        self.min_ind = min_ind

    def removeInds(self, ind_list):
        if isinstance(ind_list, list) is False:
            ind_list = [ind_list]
        for pop in self.samps:
            for indname in ind_list:
                indname = str(indname)
                if indname in self.samps[pop]:
                    self.samps[pop].remove(indname)
                    self.samp_nums[pop] -= 1
                    self.log_file.write("Removed Ind: " + indname + "\n")
        min_ind = min([self.samp_nums[pop] for pop in self.pops])
        self.min_ind = min_ind



    def combinePops(self, pops, popname):
        new_samps = []
        for pop in pops:
            for samp in self.samps[pop]:
                new_samps.append(samp)
        self.pops.append(popname)
        self.samps[popname] = new_samps
        self.samp_nums[popname] = len(new_samps)
        self.log_file.write("Combined Pops: " + str(pops) + " as " + popname + "\n")


    def splitVCFs(self, vcf_dir, min_dp, mffg, ref_path="/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta", gatk_path="/nbi/software/testing/GATK/nightly.2016.09.26/x86_64/jars/GenomeAnalysisTK.jar", repolarization_key="-99", pops='all', mem=16000, time='0-04:00', numcores=1, print1=False, overwrite=False, partition="nbi-long", keep_intermediates=False, use_scratch=False, scratch_path="/nbi/scratch/monnahap/", partition2="nbi-medium", time2="0-12:00"):
        '''Purpose:  Find all vcfs in vcf_dir and split them by population according to samples associated with said population.
                    Then, take only biallelic snps and convert vcf to table containing scaff, pos, ac, an, dp, and genotype fields.
                    Finally, concatenate all per-scaffold tables to one giant table. Resulting files will be put into ~/Working_Dir/VCFs/
            Notes: mffg is maximum fraction of filtered genotypes.  Number of actual genotypes allowed will be rounded up.
                    If you want to print the batch scripts, you must set print1 and overwrite to True'''

        if vcf_dir.endswith("/") is False:
            vcf_dir += "/"
        vcf_dir_name = vcf_dir.split("/")[-2]
        if use_scratch is True:
            if scratch_path.endswith("/") is False:
                scratch_path += "/"
            outdir = scratch_path + "VCF_" + str(vcf_dir_name) + "_DP" + str(min_dp) + ".M" + str(mffg) + "/"
        else:
            outdir = self.dir + "VCF_" + str(vcf_dir_name) + "_DP" + str(min_dp) + ".M" + str(mffg) + "/"
        self.vcf_dir = vcf_dir
        if outdir not in self.split_dirs:
            self.split_dirs.append(outdir)

        mem1 = int(mem / 1000)

        if os.path.exists(outdir) is False:
            os.mkdir(outdir)
        elif overwrite is True:
            print("Overwriting files in existing VCF directory")

        if pops == 'all':
            pops = self.pops

        for pop in pops:

            if (os.path.exists(outdir + pop + '.table.repol.txt') is True or os.path.exists(outdir + pop + '.table.recode.txt') is True) and overwrite is not True:
                print("Found file for pop = " + pop + '.  Set overwrite = True to overwrite files.')
            else:
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
                    shfile1 = open(pop + vcf_dir_name + '.sh', 'w')
                    shfile1.write('#!/bin/bash\n' +
                                  '#SBATCH -J ' + pop + vcf_dir_name + '.sh' + '\n' +
                                  '#SBATCH -e ' + self.oande + pop + vcf + '.gatk.err' + '\n' +
                                  '#SBATCH -o ' + self.oande + pop + vcf + '.gatk.out' + '\n' +
                                  '#SBATCH -p ' + str(partition) + '\n' +
                                  '#SBATCH -n ' + str(numcores) + '\n' +
                                  '#SBATCH -t ' + str(time) + '\n' +
                                  '#SBATCH --mem=' + str(mem) + '\n' +
                                  'source GATK-nightly.2016.09.26\n' +
                                  'java -Xmx' + str(mem1) + 'g -jar ' + gatk_path + ' -T SelectVariants -R ' + ref_path + ' -V ' + vcf_dir + vcf + sample_string1 + ' -o ' + outdir + vcf_basenames[v] + '.' + pop + '.vcf\n' +
                                  'gunzip ' + outdir + vcf_basenames[v] + '.' + pop + '.vcf.gz\n' +
                                  'java -Xmx' + str(mem1) + 'g -jar ' + gatk_path + ' -T VariantFiltration -R ' + ref_path + ' -V ' + outdir + vcf_basenames[v] + '.' + pop + '.vcf --genotypeFilterExpression "DP < ' + str(min_dp) + '" --genotypeFilterName "DP" -o ' + outdir + vcf_basenames[v] + '.' + pop + '.dp' + str(min_dp) + '.1.vcf\n')
                    if keep_intermediates is False: shfile1.write('rm ' + outdir + vcf_basenames[v] + '.' + pop + '.vcf\nrm ' + outdir + vcf_basenames[v] + '.' + pop + '.vcf.idx\n')
                    shfile1.write('java -Xmx' + str(mem1) + 'g -jar ' + gatk_path + ' -T VariantFiltration -R ' + ref_path + ' -V ' + outdir + vcf_basenames[v] + '.' + pop + '.dp' + str(min_dp) + '.1.vcf --setFilteredGtToNocall -o ' + outdir + vcf_basenames[v] + '.' + pop + '.dp' + str(min_dp) + '.vcf\n')
                    if keep_intermediates is False: shfile1.write('rm ' + outdir + vcf_basenames[v] + '.' + pop + '.dp' + str(min_dp) + '.1.vcf\nrm ' + outdir + vcf_basenames[v] + '.' + pop + '.dp' + str(min_dp) + '.1.vcf.idx\n')
                    shfile1.write('java -Xmx' + str(mem1) + 'g -jar ' + gatk_path + ' -T SelectVariants -R ' + ref_path + ' -V ' + outdir + vcf_basenames[v] + '.' + pop + '.dp' + str(min_dp) + '.vcf --maxNOCALLnumber ' + str(mfg) + ' -o ' + outdir + vcf_basenames[v] + '.' + pop + '.m' + str(mffg) + '.dp' + str(min_dp) + '.bi.vcf\n')
                    if keep_intermediates is False: shfile1.write('rm ' + outdir + vcf_basenames[v] + '.' + pop + '.dp' + str(min_dp) + '.vcf\nrm ' + outdir + vcf_basenames[v] + '.' + pop + '.dp' + str(min_dp) + '.vcf.idx\n')
                    shfile1.write('java -Xmx' + str(mem1) + 'g -jar ' + gatk_path + ' -T VariantsToTable -R ' + ref_path + ' -V ' + outdir + vcf_basenames[v] + '.' + pop + '.m' + str(mffg) + '.dp' + str(min_dp) + '.bi.vcf -F CHROM -F POS -F REF -F AN -F DP -GF GT -o ' + outdir + vcf_basenames[v] + '.' + pop + '_raw.table\n')
                    if keep_intermediates is False: shfile1.write('rm ' + outdir + vcf_basenames[v] + '.' + pop + '.m' + str(mffg) + '.dp' + str(min_dp) + '.bi.vcf\nrm ' + outdir + vcf_basenames[v] + '.' + pop + '.m' + str(mffg) + '.dp' + str(min_dp) + '.bi.vcf.idx\n')
                    shfile1.close()

                    if print1 is False:  # send slurm job to NBI SLURM cluster
                        cmd1 = ('sbatch ' + pop + vcf_dir_name + '.sh')
                        p1 = subprocess.Popen(cmd1, shell=True)
                        sts1 = os.waitpid(p1.pid, 0)[1]
                        joblist.append(p1.pid)

                    else:
                        file1 = open(pop + vcf_dir_name + '.sh', 'r')
                        data1 = file1.read()
                        print(data1)

                    os.remove(pop + vcf_dir_name + '.sh')


                # combine all variants table for each scaffold within a population

                shfile3 = open(pop + vcf_dir_name + '.sh', 'w')

                shfile3.write('#!/bin/bash\n' +
                              '#SBATCH -J ' + pop + vcf_dir_name + '.sh' + '\n' +
                              '#SBATCH -e ' + self.oande + pop + vcf_dir_name + '.cat.err' + '\n' +
                              '#SBATCH -o ' + self.oande + pop + vcf_dir_name + '.cat.out' + '\n' +
                              '#SBATCH -p ' + partition2 + '\n' +
                              '#SBATCH -n ' + str(numcores) + '\n' +
                              '#SBATCH -t ' + time2 + '\n' +
                              '#SBATCH --mem=' + str(mem) + '\n' +
                              'source python-3.5.1\n' +
                              'cat ' + outdir + '*' + pop + '_raw.table | tail -n+2 > ' + outdir + pop + '.table\n')
                if keep_intermediates is False: shfile3.write('rm ' + outdir + '*' + pop + '_raw.table\n')
                shfile3.write('python3 ' + self.code_dir + '/recode012.py -i ' + outdir + pop + '.table -pop ' + pop + ' -o ' + outdir + '\n')
                if keep_intermediates is False: shfile3.write('rm ' + outdir + pop + '.table\n')

                if repolarization_key == "-99":
                    print("No repolarization key provided.  Repolarized input files will not be produced.  Must set 'use_repol' to False in subsequent steps")
                else:
                    shfile3.write('python3 ' + self.code_dir + '/repol.py -i ' + outdir + pop + '.table.recode.txt -o ' + outdir + pop + ' -r ' + repolarization_key + '\n')
                    if keep_intermediates is False: shfile3.write('rm ' + outdir + pop + '.table.recode.txt\n')

                shfile3.close()

                if print1 is False:
                    cmd3 = ('sbatch -d singleton ' + pop + vcf_dir_name + '.sh')
                    p3 = subprocess.Popen(cmd3, shell=True)
                    sts3 = os.waitpid(p3.pid, 0)[1]
                else:
                    file3 = open(pop + vcf_dir_name + '.sh', 'r')
                    data3 = file3.read()
                    print(data3)

                os.remove(pop + vcf_dir_name + '.sh')

            if print1 is False:
                self.log_file.write("###  Split VCFs  ###\n" +
                                    "VCF Directory: " + vcf_dir + "\n" +
                                    "Reference Path: " + ref_path + "\n" +
                                    "Repolarization Key: " + repolarization_key + "\n" +
                                    "Output Directory: " + outdir + "\n" +
                                    "Min Depth Per Individual: " + str(min_dp) + "\n" +
                                    "Max Fraction of Filtered Genotypes: " + str(mffg) + "\n" +
                                    "Populations: " + str(pops) + "\n")


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


    def repolarize(self, recode_dir, repolarization_key, pops="all", time='0-02:00', mem='8000', numcores=1, partition="medium"):
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
                              '#SBATCH -p nbi-' + str(partition) + '\n' +
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



    def getPloidies(self, recode_dir, use_repol=True):
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
                    if use_repol is True:
                        tmp = open(recode_dir + pop + '.table.repol.txt', 'r')
                    else:    
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


    def calcwpm(self, recode_dir, window_size, min_snps, pops="all", print1=False, mem=16000, numcores=1, sampind="-99", partition="nbi-medium", use_repol=True, time="0-02:00", overwrite=False):
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

        dir_name = recode_dir.split("/")[-2]

        if os.path.exists(recode_dir) is True and sind > 3:

            for pop in pops:

                prefix = pop + ".WS" + str(window_size / 1000) + "k_MS" + str(min_snps) + "_" + str(sind) + "ind"
                if os.path.exists(recode_dir + prefix + '_WPM.txt') is True and overwrite is not True:
                    print("Output file for pop " + pop + ' already exists.  Set overwrite = True to overwrite.  Aborting.')
                else:
                    if os.path.exists(recode_dir + pop + suffix) is True:
                        shfile3 = open(dir_name + "." + pop + '.sh', 'w')

                        shfile3.write('#!/bin/bash\n' +
                                      '#SBATCH -J ' + dir_name + "." + pop + '.sh' + '\n' +
                                      '#SBATCH -e ' + self.oande + dir_name + "." + prefix + '.wpm.err' + '\n' +
                                      '#SBATCH -o ' + self.oande + dir_name + "." + prefix + '.wpm.out' + '\n' +
                                      '#SBATCH -p ' + str(partition) + '\n' +
                                      '#SBATCH -n ' + str(numcores) + '\n' +
                                      '#SBATCH -t ' + str(time) + '\n' +
                                      '#SBATCH --mem=' + str(mem) + '\n' +
                                      'source python-3.5.1\n' +
                                      'python3 ' + self.code_dir + '/wpm.py -i ' + recode_dir + pop + suffix + ' -o ' + recode_dir + ' -p ' + prefix + ' -sampind ' + str(sind) + ' -ws ' + str(window_size) + ' -ms ' + str(min_snps) + '\n')
                        shfile3.close()

                        if print1 is False:
                            cmd3 = ('sbatch -d singleton ' + dir_name + "." + pop + '.sh')
                            p3 = subprocess.Popen(cmd3, shell=True)
                            sts3 = os.waitpid(p3.pid, 0)[1]

                        else:
                            file3 = open(dir_name + "." + pop + '.sh', 'r')
                            data3 = file3.read()
                            print(data3)

                        os.remove(dir_name + "." + pop + '.sh')
                    else:
                        print("Did not find input files for: ", pop)

            if print1 is False:
                self.log_file.write("###  Calculate Within-Population-Metrics  ###\n" +
                                    "Input Directory: " + recode_dir + "\n" +
                                    "Window Size: " + str(window_size) + "\n" +
                                    "Minimum SNPs in a window: " + str(min_snps) + "\n" +
                                    "Number of individuals to be downsampled to: " + str(sampind) + "\n" +
                                    "Use repolarized data: " + str(use_repol) + "\n" +
                                    "Populations: " + str(pops) + "\n")

        elif sind <= 3:
            print("Number of individuals to be used/downsampled to is <= 3.  Unable to calculate within-population-metrics on so few individuals.")
        else:
            print("Did not find recode_dir.  Must run splitVCFs followed by recode before able to calculate within population metrics")


    def concatWPM(self, recode_dir, suffix, outname, pops='all'):
        '''Purpose:  Concatenate _WPM.txt files corresponding to populations indicated in pops parameter.'''

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True:
            new = open(recode_dir + outname + "_WPM.txt", 'w')
            if pops == 'all':
                pops = self.pops
            head = False
            for i, pop in enumerate(pops):
                try:
                    with open(recode_dir + pop + suffix, 'r') as inf:
                        for j, line in enumerate(inf):
                            if j == 0 and head is False:
                                new.write(line)
                                head = True
                            elif j != 0:
                                new.write(line)
                except FileNotFoundError:
                    print("Did not find _WPM.txt file for population: ", pop)


    def concatBPM(self, recode_dir, suffix, outname, pops='all', strict=False):
        '''Purpose:  Concatenate _WPM.txt files corresponding to populations indicated in pops parameter.'''

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if os.path.exists(recode_dir) is True:
            new = open(recode_dir + outname + "_BPM.txt", 'w')
            if pops == 'all':
                pops = self.pops
            head = False
            for file in os.listdir(recode_dir):
                if file.endswith(suffix):
                    if strict is False:
                        if file.split("_")[0][:3] in pops or file.split("_")[0][3:] in pops:
                            with open(recode_dir + file, 'r') as inf:
                                for j, line in enumerate(inf):
                                    if j == 0 and head is False:
                                        new.write(line)
                                        head = True
                                    elif j != 0:
                                        new.write(line)
                    elif strict is True:
                        if file.split("_")[0][:3] in pops and file.split("_")[0][3:] in pops:
                            with open(recode_dir + file, 'r') as inf:
                                for j, line in enumerate(inf):
                                    if j == 0 and head is False:
                                        new.write(line)
                                        head = True
                                    elif j != 0:
                                        new.write(line)


    def calcbpm(self, recode_dir, pops, output_name, window_size, min_snps, print1=False, mem=16000, numcores=1, partition="nbi-medium", use_repol=True, keep_intermediates=False, time="0-12:00", use_scratch=False, scratch_path="/nbi/scratch/monnahap"):
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

        if use_scratch is True:
            if scratch_path.endswith("/") is False:
                scratch_path += "/"
            tmpdir = scratch_path + recode_dir.split("/")[-2] + "/"
            if os.path.exists(tmpdir) is False:
                os.mkdir(tmpdir)
        else:
            tmpdir = recode_dir

        output_name += "_WS" + str(window_size) + "_MS" + str(min_snps)
        if os.path.exists(recode_dir) is True and len(pops) > 1:
            pop_num = 0
            file_string = ""
            for pop in pops:
                try:
                    a = open(recode_dir + pop + suffix, 'r')
                    a.close()
                    file_string += recode_dir + pop + suffix + " "
                    pop_num += 1
                except IOError:
                    print("Did not find input file for pop ", pop)
            if len(pops) != pop_num:
                print("Did not find all input files!!  Aborting.")
                os.remove(recode_dir + output_name + '.concat.txt')
            else:
                shfile3 = open(recode_dir + output_name + '.bpm.sh', 'w')

                shfile3.write('#!/bin/bash\n' +
                              '#SBATCH -J ' + output_name + '.bpm.sh' + '\n' +
                              '#SBATCH -e ' + self.oande + output_name + '.bpm.err' + '\n' +
                              '#SBATCH -o ' + self.oande + output_name + '.bpm.out' + '\n' +
                              '#SBATCH -p ' + str(partition) + '\n' +
                              '#SBATCH -n ' + str(numcores) + '\n' +
                              '#SBATCH -t ' + str(time) + '\n' +
                              '#SBATCH --mem=' + str(mem) + '\n' +
                              'source python-3.5.1\n' +
                              'sort -k3,3 -k4,4n -m ' + file_string + '> ' + tmpdir + output_name + '.concat.txt\n' +
                              'python3 ' + self.code_dir + '/bpm.py -i ' + tmpdir + output_name + '.concat.txt' + ' -o ' + recode_dir + ' -prefix ' + output_name + ' -ws ' + str(window_size) + ' -ms ' + str(min_snps) + ' -np ' + str(pop_num) + '\n')
                if keep_intermediates is False:
                    shfile3.write('rm ' + recode_dir + output_name + '.concat.txt')
                shfile3.close()

                if print1 is False:
                    cmd3 = ('sbatch -d singleton ' + recode_dir + output_name + '.bpm.sh')
                    p3 = subprocess.Popen(cmd3, shell=True)
                    sts3 = os.waitpid(p3.pid, 0)[1]

                    self.log_file.write("###  Calculate Between-Population-Metrics  ###\n" +
                                        "Input Directory: " + recode_dir + "\n" +
                                        "Window Size: " + str(window_size) + "\n" +
                                        "Minimum SNPs in a window: " + str(min_snps) + "\n" +
                                        "Use repolarized data: " + str(use_repol) + "\n" +
                                        "Populations: " + str(pops) + "\n")

                else:
                    file3 = open(recode_dir + output_name + '.bpm.sh', 'r')
                    data3 = file3.read()
                    print(data3)

                os.remove(recode_dir + output_name + '.bpm.sh')
        elif len(pops) < 2:
            print("'pops' argument must be a list of strings specifiying two or more population names as they appear in input file prefixes.  len(pops) was < 2")
        else:
            print("Did not find recode_dir.  Must run splitVCFs followed by recode before able to calculate between population metrics")

    def calcPairwisebpm(self, recode_dir, pops, window_size, min_snps, print1=False, mem=16000, numcores=1, partition="nbi-medium", use_repol=True, keep_intermediates=False, time="0-01:00", overwrite=False, use_scratch=False, scratch_path="/nbi/scratch/monnahap/"):
        '''Purpose:  Calculate between population metrics including: Dxy, Fst (using Weir and Cockerham 1984), and Rho (Ronfort et al. 1998)
           Notes: User provides a list of populations to be included.  For pairwise estimates, simply provide two populations
                    Calculations are done for windows of a given bp size.  User also must specify the minimum number of snps (min_snps) in a window
                    for calculations to be made'''

        if use_repol is True:
            suffix = '.table.repol.txt'
        else:
            suffix = '.table.recode.txt'

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if use_scratch is True:
            if scratch_path.endswith("/") is False:
                scratch_path += "/"
            tmpdir = scratch_path + recode_dir.split("/")[-2] + "/"
            if os.path.exists(tmpdir) is False:
                os.mkdir(tmpdir)
        else:
            tmpdir = recode_dir

        if os.path.exists(recode_dir) is True and len(pops) > 1:
            # Concatenate input files and sort them
            for i, pop1 in enumerate(pops):  # Add data from all populations to single, huge listg
                for pop2 in pops[i + 1:]:
                    output_name = pop1 + pop2 + "_WS" + str(window_size) + "_MS" + str(min_snps)
                    skip = False
                    try:
                        a = open(recode_dir + pop1 + suffix, 'r')
                        a.close()
                    except IOError:
                        print("Did not find input file for pop ", pop1)
                        skip = True
                    try:
                        a = open(recode_dir + pop2 + suffix, 'r')
                        a.close()
                    except IOError:
                        print("Did not find input file for pop ", pop2)
                        skip = True
                    if skip is True:
                        print("Did not find all input files!!  Aborting pairwise bpm for contrast: ", output_name)
                    elif os.path.exists(recode_dir + output_name + '_BPM.txt') and overwrite is False:
                        print(recode_dir + output_name + '_BPM.txt already exists.  Set "overwrite" to True if you want to overwrite.')
                    else:
                        shfile3 = open(recode_dir + output_name + '.bpm.sh', 'w')

                        shfile3.write('#!/bin/bash\n' +
                                      '#SBATCH -J ' + output_name + '.bpm.sh' + '\n' +
                                      '#SBATCH -e ' + self.oande + output_name + '.bpm.err' + '\n' +
                                      '#SBATCH -o ' + self.oande + output_name + '.bpm.out' + '\n' +
                                      '#SBATCH -p ' + partition + '\n' +
                                      '#SBATCH -n ' + str(numcores) + '\n' +
                                      '#SBATCH -t ' + str(time) + '\n' +
                                      '#SBATCH --mem=' + str(mem) + '\n' +
                                      'source python-3.5.1\n' +
                                      'sort -k3,3 -k4,4n -m ' + recode_dir + pop1 + suffix + " " + recode_dir + pop2 + suffix + " > " + tmpdir + output_name + '.concat.txt\n' +
                                      'python3 ' + self.code_dir + '/bpm.py -i ' + tmpdir + output_name + '.concat.txt' + ' -o ' + recode_dir + ' -prefix ' + output_name + ' -ws ' + str(window_size) + ' -ms ' + str(min_snps) + ' -np 2\n')
                        if keep_intermediates is False:
                            shfile3.write('rm ' + recode_dir + output_name + '.concat.txt')
                        shfile3.close()

                        if print1 is False:
                            cmd3 = ('sbatch -d singleton ' + recode_dir + output_name + '.bpm.sh')
                            p3 = subprocess.Popen(cmd3, shell=True)
                            sts3 = os.waitpid(p3.pid, 0)[1]

                        else:
                            file3 = open(recode_dir + output_name + '.bpm.sh', 'r')
                            data3 = file3.read()
                            print(data3)

                        os.remove(recode_dir + output_name + '.bpm.sh')

            if print1 is False:
                self.log_file.write("###  Calculate PAIRWISE Between-Population-Metrics  ###\n" +
                                    "Input Directory: " + recode_dir + "\n" +
                                    "Window Size: " + str(window_size) + "\n" +
                                    "Minimum SNPs in a window: " + str(min_snps) + "\n" +
                                    "Use repolarized data: " + str(use_repol) + "\n" +
                                    "Populations: " + str(pops) + "\n")
        elif len(pops) < 2:
            print("'pops' argument must be a list of strings specifiying two or more population names as they appear in input file prefixes.  len(pops) was < 2")
        else:
            print("Did not find recode_dir.  Must run splitVCFs followed by recode before able to calculate between population metrics")


    def calcAFS(self, recode_dir, data_name, sampind="-99", pops="-99", time="0-00:30", mem="8000", use_repol=True, print1=False, allow_one_missing=True, partition="nbi-short"):
        """Purpose:  calculate allele frequency spectrum for each population in pops list and downsample to number specified in sampind.  If no downsampling is specified the number of indviduals in the population (minus one if allow_one_missing=True) will be used."""
        if recode_dir.endswith("/") is False:
            recode_dir += "/"
        if pops == "-99":
            pops = self.pops

        if use_repol is True:
            suffix = '.table.repol.txt'
        else:
            suffix = '.table.recode.txt'
        for pop in pops:
            if sampind == "-99":
                si = self.samp_nums[pop]
                if allow_one_missing is True:
                    si -= 1
            else:
                si = sampind
            infile = recode_dir + pop + suffix
            shfile3 = open(infile + '.afs.sh', 'w')
            prefix = pop + "_" + data_name
            shfile3.write('#!/bin/bash\n' +
                          '#SBATCH -J ' + pop + '.afs.sh' + '\n' +
                          '#SBATCH -e ' + self.oande + pop + '.afs.err' + '\n' +
                          '#SBATCH -o ' + self.oande + pop + '.afs.out' + '\n' +
                          '#SBATCH -p ' + partition + '\n'
                          '#SBATCH -n 1' + '\n' +
                          '#SBATCH -t ' + str(time) + '\n' +
                          '#SBATCH --mem=' + str(mem) + '\n' +
                          'source python-3.5.1\n' +
                          'python3 ' + self.code_dir + '/calcAFS.py -i ' + infile + ' -o ' + recode_dir + ' -p ' + prefix + ' -sampind ' + str(si) + '\n')
            shfile3.close()
            if print1 is False:
                cmd3 = ('sbatch ' + infile + '.afs.sh')
                p3 = subprocess.Popen(cmd3, shell=True)
                sts3 = os.waitpid(p3.pid, 0)[1]

            else:
                file3 = open(infile + '.afs.sh', 'r')
                data3 = file3.read()
                print(data3)
            os.remove(infile + '.afs.sh')


    def calcFreqs(self, recode_dir, outfile_name, sites_file, pops="-99", time="0-02:00", partition="nbi-short", mem="8000", use_repol=True, print1=False):
        """Calls calcFreqs_atSites.py.  Takes a list of sites in (sites_file, should be formatted so that each line simply has scaffold and position, with scaffold simply coded as an integer 0-8) and calculates the allele frequency in each population (list_of_populations) at each site"""
        if recode_dir.endswith("/") is False:
            recode_dir += "/"
        if pops == "-99":
            pops = self.pops

        if use_repol is True:
            suffix = '.table.repol.txt'
        else:
            suffix = '.table.recode.txt'

        pop_string = ""
        for pop in pops:
            if os.path.exists(recode_dir + pop + suffix) is True:
                pop_string += pop + ","
        pop_string.strip(",")

        shfile3 = open('CalcFreqs.sh', 'w')
        shfile3.write('#!/bin/bash\n' +
                      '#SBATCH -J CalcFreqs.sh' + '\n' +
                      '#SBATCH -e ' + self.oande + 'CalcFreqs.err' + '\n' +
                      '#SBATCH -o ' + self.oande + 'CalcFreqstest.out' + '\n' +
                      '#SBATCH -p ' + partition + '\n'
                      '#SBATCH -n 1' + '\n' +
                      '#SBATCH -t ' + str(time) + '\n' +
                      '#SBATCH --mem=' + str(mem) + '\n' +
                      'source python-3.5.1\n' +
                      'python3 ' + self.code_dir + '/calcFreqs_atSites.py -i ' + recode_dir + ' -of ' + outfile_name + ' -s ' + sites_file + ' -pops ' + pop_string + ' -suffix ' + suffix + '\n')
        shfile3.close()
        if print1 is False:
            cmd3 = ('sbatch CalcFreqs.sh')
            p3 = subprocess.Popen(cmd3, shell=True)
            sts3 = os.waitpid(p3.pid, 0)[1]

        else:
            file3 = open('CalcFreqs.sh', 'r')
            data3 = file3.read()
            print(data3)
        os.remove('CalcFreqs.sh')


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
            try:
                print("Finding outliers for " + recode_dir + in_file + "\n")
                data = pandas.read_table(recode_dir + in_file, header=0)
                metrics = []
                for i in column_index_list:
                    try:
                        metrics.append(list(data.columns.values)[i])
                    except IndexError:
                        print("IndexError in " + recode_dir + in_file + "\n")
                data.start = data.start.astype(int)
                data.end = data.end.astype(int)
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
                df_outlier.to_csv(recode_dir + in_file.replace(".txt", "") + '_' + str(percentile) + 'tile_OutOnly.bed', index=False, sep='\t', columns=["scaff", "start", "end"], header=False)
            except (UnicodeDecodeError, IndexError, KeyError):
                print('Error with file: ' + recode_dir + in_file + "\n")


    def annotateOutliers(self, recode_dir, in_file, annotation_file='/nbi/Research-Groups/JIC/Levi-Yant/GenomeScan/LyV2.gff', overlap_proportion=0.000001, print1=False, partition="nbi-short", time="0-01:00", mem=4000):
        '''Purpose: annotate bed file from findOutliers using information in annotation_file
           Notes: The output (suffix ol_genes.gff) only contains the window locations along with annotation info and does not contain
                    the original metric information used to determine outliers.  Use mergeAnnotation to merge original outlier file with annotation info'''

        if recode_dir.endswith("/") is False:
            recode_dir += "/"
        basename = in_file.strip(".bed")
        if os.path.exists(recode_dir) is True:
            shfile1 = open(recode_dir + in_file + 'bedtools_gff.sh', 'w')
            shfile1.write('#!/bin/bash\n' +
                          '#SBATCH -J ' + in_file + '.bedtools.sh' + '\n' +
                          '#SBATCH -e ' + self.oande + in_file + '.bedtools.err\n' +
                          '#SBATCH -o ' + self.oande + in_file + '.bedtools.out\n' +
                          '#SBATCH -p ' + partition + '\n' +
                          '#SBATCH -n 1\n' +
                          '#SBATCH -t ' + time + '\n' +
                          '#SBATCH --mem=' + str(mem) + '\n' +
                          'source bedtools-2.17.0\n' +
                          'bedtools intersect -a ' + recode_dir + in_file + ' -b ' + annotation_file + ' -f ' + str(overlap_proportion) + ' -wo | grep transcript | grep -v transcription | sort -u |' +
                          """awk '{print $1,$2,$3,$7,$8,$9,$10,$12}'""" +
                          '| ' +
                          """tr ' ' '\t' """
                          '> ' + recode_dir + basename + '_genes.gff')
            shfile1.close()

            if print1 is False:
                cmd1 = ('sbatch ' + recode_dir + in_file + 'bedtools_gff.sh')
                p1 = subprocess.Popen(cmd1, shell=True)
                sts1 = os.waitpid(p1.pid, 0)[1]

            else:
                file3 = open(recode_dir + in_file + 'bedtools_gff.sh', 'r')
                data3 = file3.read()
                print(data3)

            os.remove(recode_dir + in_file + 'bedtools_gff.sh')

        else:
            print("recode_dir not found")


    def mergeAnnotation(self, recode_dir, outlier_file):
        '''Purpose: Merge the annotation information with the original outlier file results from findOutliers.'''

        if recode_dir.endswith("/") is False:
            recode_dir += "/"
        annotated_outlier_file = outlier_file.strip('.csv') + '_genes.gff'
        if os.path.exists(recode_dir) is True:
            try:
                outliers = pandas.read_csv(recode_dir + outlier_file, header=0)
                annotation = pandas.read_table(recode_dir + annotated_outlier_file, names=["scaff", "start", "end", "gene_start", "gene_end", "overlap", "strand", "geneName"])
            except IOError:
                print("Did not find either original outlier file or the annotated outlier file")
            merged = outliers.merge(annotation, on=["scaff", "start", "end"],)
            merged.to_csv(recode_dir + outlier_file.replace("_OutOnly.csv", "_OutAnnot.csv"), index=False)
        else:
            print("Did not find recode_dir")

    def Outliers(self, recode_dir, in_file, column_index_list, percentile, tails='upper', annotation_file='/nbi/Research-Groups/JIC/Levi-Yant/GenomeScan/LyV2.gff', overlap_proportion=0.000001):
        """Purpose:  Wraps .findOutliers, .annotateOutliers, and .mergeAnnotation into a single function"""
        print("Be sure that no old versions of gff files are in this directory")
        self.findOutliers(recode_dir, in_file, column_index_list, percentile, tails)
        bed_file = in_file.replace(".txt", "") + '_' + str(percentile) + 'tile_OutOnly.bed'
        self.annotateOutliers(recode_dir, bed_file, annotation_file, overlap_proportion)
        outlier_file = in_file.replace(".txt", "") + '_' + str(percentile) + 'tile_OutOnly.csv'
        annotated_outlier_file = outlier_file.strip('.csv') + '_genes.gff'
        while not os.path.exists(recode_dir + annotated_outlier_file):
            time.sleep(1)

        if os.path.isfile(recode_dir + annotated_outlier_file):
            time.sleep(30)
            self.mergeAnnotation(recode_dir, outlier_file)
        else:
            raise ValueError("File error for mergeAnnotation: %s" % recode_dir + annotated_outlier_file)


    def generateFSC2input(self, recode_dir, pops, output_name, bootstrap_block_size=50000, bootstrap_reps=0, mem=16000, numcores=1, time='2-00:00', print1=False, use_repol=True, keep_intermediates=False, alphabetical_pop_order='false', use_scratch=True, partition="nbi-medium", scratch_path="/nbi/scratch/monnahap/"):
        '''Purpose:  Generate --multiSFS for fastsimcoal2 along with a given number of non-parametric block-bootstrapped replicates
           Notes: Must provide the block size for bootstrapping as well as number of bootstrap replicates
                  As of now, the necessary template files for FSC2 must come from elsewhere.  Beware of running this method with numerous populations'''

        if recode_dir.endswith("/") is False:
            recode_dir += "/"

        if use_scratch is True:
            if scratch_path.endswith("/") is False:
                scratch_path += "/"
            tmpdir = scratch_path
        else:
            tmpdir = recode_dir


        outdir = recode_dir + "FSC2input_" + output_name + "/"

        if os.path.exists(outdir) is False:
            os.mkdir(outdir)

        num_pops = len(pops)

        if os.path.exists(recode_dir) is True:
            if use_repol is True:
                suffix1 = '.table.repol.txt'
                suffix2 = 'repol.concat.txt'
            else:
                suffix1 = '.table.recode.txt'
                suffix2 = '.recode.concat.txt'
            concat_name = tmpdir + output_name + suffix2
            missing = []
            if os.path.exists(concat_name) is False:
                print("Concatenating input files")
                concat_file = open(tmpdir + output_name + suffix2, 'w')
                for pop in pops:
                    try:
                        with open(recode_dir + pop + suffix1) as infile:
                            for line in infile:
                                concat_file.write(line)
                    except FileNotFoundError:
                        missing.append(pop)
            if len(missing) != 0:
                print("Did not find input files for the following populations:", missing, ".  Aborting!!")
            else:
                print("Finished preparing input data")

                shfile4 = open(output_name + '.fsc2input.sh', 'w')

                shfile4.write('#!/bin/bash\n' +
                              '#SBATCH -J ' + output_name + '.fsc2input.sh' + '\n' +
                              '#SBATCH -e ' + self.oande + output_name + '.fsc2input.err' + '\n' +
                              '#SBATCH -o ' + self.oande + output_name + '.fsc2input.out' + '\n' +
                              '#SBATCH -p ' + str(partition) + '\n' +
                              '#SBATCH -n ' + str(numcores) + '\n' +
                              '#SBATCH -t ' + str(time) + '\n' +
                              '#SBATCH --mem=' + str(mem) + '\n' +
                              'source python-3.5.1\n' +
                              'source env/bin/activate\n' +
                              'python3 ' + self.code_dir + '/FSC2input.py -i ' + concat_name + ' -o ' + outdir + ' -prefix ' + output_name + ' -ws ' + str(bootstrap_block_size) + ' -bs ' + str(bootstrap_reps) + ' -np ' + str(num_pops) + ' -alpha ' + str(alphabetical_pop_order) + '\n')
                if keep_intermediates is False:
                    shfile4.write('rm ' + concat_name + "\n")
                shfile4.close()

                if print1 is False:
                    cmd1 = ('sbatch ' + output_name + '.fsc2input.sh')
                    p1 = subprocess.Popen(cmd1, shell=True)
                    sts1 = os.waitpid(p1.pid, 0)[1]

                    self.log_file.write("###  Generate FSC2 Input  ###\n" +
                                        "Input Directory: " + recode_dir + "\n" +
                                        "Bootstrap Block Size: " + str(bootstrap_block_size) + "\n" +
                                        "Bootstrap Replicate Data Sets: " + str(bootstrap_reps) + "\n" +
                                        "Output Name: " + output_name + "\n" +
                                        "Populations: " + str(pops) + "\n")

                else:
                    file3 = open(output_name + '.fsc2input.sh', 'r')
                    data3 = file3.read()
                    print(data3)

                os.remove(output_name + '.fsc2input.sh')

        else:
            print("!!!Did not find recode_dir!!!!")

    def FSC2(self, input_dir, num_jobs=400, sleep_time=0.5, num_reps=50, min_sims=10000, max_sims=100000, conv_crit=0.001, min_ecm=10, max_ecm=40, calc_CI=False, submit_partition="nbi-medium", job_partition="nbi-short", submit_time="2-00:00", job_time="0-02:00", numcores=1, time="0-02:00", mem="8000", print1=False, overwrite="None", fsc2_path="/nbi/Research-Groups/JIC/Levi-Yant/Patrick/fsc_linux64/fsc25221", cluster="JIC", verbose=False):
        """This method parallelises job submission of fastsimcoal2, but requires a very specific set up of input files.  The output of '.generateFSC2input' should be a folder that contains the multi-dimensional SFS.  Place this folder in a new folder that will be the FSC2_Data_Parent_Directory.  This directory should also contain one or more template (.tpl) and estimates (.est) files whose format can be found in the fastsimcoal2 documentation.  For each sub-directory containing input data, this method will re-format and rename the .tpl and .est files to reflect the necessary information in the sub-directory multi-dimensional SFS and then submit these jobs to the cluster.  I've tried to make the code as general as possible, but this is one method that will likely require the user to read and understand the code in order to get things working well for them.  Also, a major potential source of errors is in the correct formatting of the .tpl and .est files, so it is worthwhile to ensure that these are correct (by running FSC2 on a subset of your sub-directories) before launching full-scale"""
        # num_jobs is number of jobs to be queued.  the program will submit this number of jobs and then submit additional jobs as previous ones finish, but never exceeding the specified number of jobs at any given time.

        tpl_files = []
        est_files = []

        for path in os.listdir(input_dir):
            if os.path.isdir(input_dir + path) and path.startswith("FSC2input"):
                samp_name = path.split("_")[1]
                if samp_name + "_DSFS.obs" not in os.listdir(input_dir + path):
                    print("Did not find input data file for: ", samp_name)
                if calc_CI is True:
                    num_files = 0
                    for file in os.listdir(input_dir + path):
                        if file.endswith("_DSFS.obs") and file.split("_")[-2].split(".")[-1][0:3] == "rep" and file != samp_name + "_DSFS.obs":
                            num_files += 1
                    if len(num_files) < 1:
                        print("Did not find bootstrap replicates for: ", samp_name)
                    else:
                        print("Found ", num_files, " replicate dsfs files for CI calculation for ", samp_name)
            if path.endswith(".tpl"):
                tpl_files.append(path)
                est_files.append(path.split(".")[0])
        if len(tpl_files) == 0:
            print("Did not find any tpl files!! Aborting!!")
        else:
            if any(os.path.exists(input_dir + k.split(".tpl")[0] + ".est") is False for k in tpl_files):
                print("Did not find all est files.  Aborting!!")
            else:
                shfile5 = open("FSC2_submit.sh", 'w')
                shfile5.write('#!/bin/bash\n' +
                              '#SBATCH -J fsc2_submit.sh\n' +
                              '#SBATCH -e ' + self.oande + 'fsc2_submit.err\n' +
                              '#SBATCH -o ' + self.oande + 'fsc2_submit.out\n' +
                              '#SBATCH -p ' + submit_partition + '\n' +
                              '#SBATCH -n 1\n' +
                              '#SBATCH -t ' + str(submit_time) + '\n' +
                              '#SBATCH --mem=4000\n' +
                              'unset SBATCH_PARTITION\n' +
                              'source python-3.5.1\n' +
                              'source env/bin/activate\n' +
                              'python3 ' + self.code_dir + '/FSC2_submit.py -i ' + input_dir + ' -nj ' + str(num_jobs) + ' -reps ' + str(num_reps) + ' -minsims ' + str(min_sims) + ' -maxsims ' + str(max_sims) + ' -c ' + str(conv_crit) + ' -min_ecm ' + str(min_ecm) + ' -max_ecm ' + str(max_ecm) + ' -ci ' + str(calc_CI) + ' -p ' + job_partition + ' -nc ' + str(numcores) + ' -mem ' + str(mem) + ' -t ' + time + ' -print1 ' + str(print1) + ' -Ov ' + str(overwrite) + ' -fsc2path ' + fsc2_path + ' -clust ' + cluster + ' -oande ' + self.oande + ' -verbose ' + str(verbose) + '\n')

                shfile5.close()
                if print1 is False:
                    cmd1 = ("sbatch FSC2_submit.sh")
                    p1 = subprocess.Popen(cmd1, shell=True)
                    sts1 = os.waitpid(p1.pid, 0)[1]

                else:
                    file3 = open("FSC2_submit.sh", 'r')
                    data3 = file3.read()
                    print(data3)

                os.remove("FSC2_submit.sh")


    def gatherFSC2output(self, parent_dir):
        """This method collects all information from the '.bestlhood' output files of FSC2 that is buried in the sub-directories and outputs the information into one of two files:  Likelihoods file and parameters file."""
        if parent_dir.endswith("/") is False:
            parent_dir += "/"
        get_header = True
        dirname = parent_dir.strip("/").split("/")[-1]
        Lhoodfile = open(parent_dir + dirname + "_FSC2_Likelihoods.txt", 'w')
        paramsfile = open(parent_dir + dirname + "_FSC2_Params.txt", 'w')
        for root, dirs, files in os.walk(parent_dir):
            for file in files:
                if file.endswith(".bestlhoods"):
                    name = file.split(".best")[0]
                    samp_names = name.split("_")[0]
                    model = ".".join(name.split("_")[1:])
                    if get_header is True:
                        get_header = False
                        with open(os.path.join(root, file), 'r') as ff:
                            for i, line in enumerate(ff):
                                if i == 0:
                                    Lhoodfile.write("Model\tSampleNames\tnum_params\tLhood_est\tLhood_obs\tLhood_diff\tAIC\n")
                                else:
                                    info = line.strip("\n").split("\t")
                                    num_params = len(info) - 2
                                    Lhood_est = float(info[-2])
                                    AIC = (2 * num_params) - (2 * (Lhood_est * math.log(10)))
                                    Lhoodfile.write(model + "\t" + samp_names + "\t" + str(num_params) + "\t" + str(Lhood_est) + "\t" + info[-1] + "\t" + str(Lhood_est - float(info[-1])) + "\t" + str(AIC) + "\n")
                                    paramsfile.write(model + "\t" + samp_names + "\t" + line)
                    else:
                        with open(os.path.join(root, file), 'r') as ff:
                            for i, line in enumerate(ff):
                                if i > 0:
                                    info = line.strip("\n").split("\t")
                                    num_params = len(info) - 2
                                    Lhood_est = float(info[-2])
                                    AIC = (2 * num_params) - (2 * (Lhood_est * math.log(10)))
                                    Lhoodfile.write(model + "\t" + samp_names + "\t" + str(num_params) + "\t" + str(Lhood_est) + "\t" + info[-1] + "\t" + str(Lhood_est - float(info[-1])) + "\t" + str(AIC) + "\n")
                                    paramsfile.write(model + "\t" + samp_names + "\t" + line)
        Lhoodfile.close()
        paramsfile.close()


    def queryFSC2input(self, input_dir, index_list, outname):
        """PROVIDE SUMMARY HERE AND IN README"""
        outfile = open(input_dir + outname + ".txt", 'w')
        outfile.write("Outname\tPops\tProp\n")
        for path in os.listdir(input_dir):
            if os.path.isdir(input_dir + path) and path.startswith("FSC2input"):
                samp_name = path.split("_")[1]
                if samp_name + "_DSFS.obs" in os.listdir(input_dir + path):
                    with open(input_dir + path + '/' + samp_name + "_DSFS.obs") as fsc2input:
                        for i, line in enumerate(fsc2input):
                            if i == 2:
                                line = line.strip("\n").strip("\t").split("\t")
                                tot = sum([int(j) for j in line])
                                sp = 0  # Shared polymorphisms
                                for ix in index_list:
                                    sp += int(line[ix])
                    outfile.write(outname + '\t' + samp_name + '\t' + str(float(sp) / float(tot)) + '\n')

    # def calcMissing(self, vcf_dir, min_dp, window_size, ref_path="/nbi/Research-Groups/JIC/Levi-Yant/Lyrata_ref/alygenomes.fasta", gatk_path="/nbi/software/testing/GATK/nightly.2016.09.26/x86_64/jars/GenomeAnalysisTK.jar", pops='all', mem=16000, time='0-04:00', numcores=1, print1=False, partition="long", keep_intermediates=False):
    #     """NOT FUNCTIONAL"""
    #     if vcf_dir.endswith("/") is False:
    #         vcf_dir += "/"
    #     vcf_dir_name = vcf_dir.split("/")[-2]
    #     outdir = self.dir + "VCF_" + str(vcf_dir_name) + "_DP" + str(min_dp) + ".M" + str(mffg) + "/"
    #     self.vcf_dir = vcf_dir
    #     if outdir not in self.split_dirs:
    #         self.split_dirs.append(outdir)

    #     mem1 = int(mem / 1000)

    #     if os.path.exists(outdir) is False:
    #             os.mkdir(outdir)

    #     if pops == 'all':
    #         pops = self.pops

    #     for pop in pops:
    #         # Add samples to list for each population according to PF file
    #         sample_string1 = ""
    #         for samp in self.samps[pop]:
    #             sample_string1 += " -sn " + samp

    #         vcf_list = []
    #         vcf_basenames = []
    #         for file in os.listdir(vcf_dir):
    #             if file[-6:] == 'vcf.gz':
    #                 vcf_list.append(file)
    #                 vcf_basenames.append(file[:-7])
    #             elif file[-3:] == 'vcf':
    #                 vcf_list.append(file)
    #                 vcf_basenames.append(file[:-4])
    #         for v, vcf in enumerate(vcf_list):
    #             # Select single population and biallelic SNPs for each scaffold and convert to variants table
    #             shfile1 = open(vcf_dir_name + '.sh', 'w')
    #             shfile1.write('#!/bin/bash\n' +
    #                           '#SBATCH -J ' + vcf_dir_name + '.sh' + '\n' +
    #                           '#SBATCH -e ' + self.oande + pop + vcf + '.gatk.err' + '\n' +
    #                           '#SBATCH -o ' + self.oande + pop + vcf + '.gatk.out' + '\n' +
    #                           '#SBATCH -p nbi-' + str(partition) + '\n' +
    #                           '#SBATCH -n ' + str(numcores) + '\n' +
    #                           '#SBATCH -t ' + str(time) + '\n' +
    #                           '#SBATCH --mem=' + str(mem) + '\n' +
    #                           'source GATK-nightly.2016.09.26\n' +
    #                           'java -Xmx' + str(mem1) + 'g -jar ' + gatk_path + ' -T SelectVariants -R ' + ref_path + ' -V ' + vcf_dir + vcf + sample_string1 + ' -o ' + outdir + vcf_basenames[v] + '.' + pop + '.vcf\n' +
    #                           'gunzip ' + outdir + vcf_basenames[v] + '.' + pop + '.vcf.gz\n')
    #             shfile1.close()

    #             if print1 is False:  # send slurm job to NBI SLURM cluster
    #                 cmd1 = ('sbatch ' + pop + vcf_dir_name + '.sh')
    #                 p1 = subprocess.Popen(cmd1, shell=True)
    #                 sts1 = os.waitpid(p1.pid, 0)[1]

    #             else:
    #                 file1 = open(pop + vcf_dir_name + '.sh', 'r')
    #                 data1 = file1.read()
    #                 print(data1)

    #             os.remove(pop + vcf_dir_name + '.sh')


    #     # combine all variants table for each scaffold within a population

    #     shfile3 = open(vcf_dir_name + '.sh', 'w')

    #     shfile3.write('#!/bin/bash\n' +
    #                   '#SBATCH -J ' + vcf_dir_name + '.sh' + '\n' +
    #                   '#SBATCH -e ' + self.oande + vcf_dir_name + '.missing.err' + '\n' +
    #                   '#SBATCH -o ' + self.oande + vcf_dir_name + '.missing.out' + '\n' +
    #                   '#SBATCH -p nbi-medium\n' +
    #                   '#SBATCH -n ' + str(numcores) + '\n' +
    #                   '#SBATCH -t 2-00:00\n' +
    #                   '#SBATCH --mem=' + str(mem) + '\n' +
    #                   'source python-3.5.1\n' +
    #                   'python3 ' + self.code_dir + '/MissingData.py -v ' + outdir + ' -w ' + str(window_size) + ' -dp ' + str(min_dp) + ' -gz false -o ' + outdir + 'MissingData_PerPop.txt\n')

    #     if keep_intermediates is False:
    #         shfile3.write('rm ' + outdir + '*.' + pop + '.vcf\n')
    #         shfile3.write('rm ' + outdir + '*.' + pop + '.vcf.idx\n')
    #     shfile3.close()

    #     if print1 is False:
    #         cmd3 = ('sbatch -d singleton ' + pop + vcf_dir_name + '.sh')
    #         p3 = subprocess.Popen(cmd3, shell=True)
    #         sts3 = os.waitpid(p3.pid, 0)[1]
    #     else:
    #         file3 = open(pop + vcf_dir_name + '.sh', 'r')
    #         data3 = file3.read()
    #         print(data3)

    #     os.remove(vcf_dir_name + '.sh')
