import os
import argparse
import subprocess
import time

jobindex = 0
Processes = []


def StartNew(shlist, maxjobs, chk_file_list):
    """ Start a new subprocess if there is work to do """
    global jobindex
    global Processes

    if jobindex < len(shlist):
        proc = subprocess.Popen(['sbatch', shlist[jobindex]])
        print("Started to Process %s", shlist[jobindex])
        print("Job number ", jobindex + 1, " out of ", len(shlist), ".  ", (float(jobindex + 1) / float(len(shlist))) * 100.0, ' percent complete.')
        jobindex += 1
        Processes.append(chk_file_list[jobindex])


def CheckRunning(shlist, maxjobs, chk_file_list):
    """ Check any running processes and start new ones if there are spare slots."""
    global Processes
    global jobindex

    if len(Processes) > 0:
        for p in range(len(Processes) - 1, -1, -1):  # Check the processes in reverse order
            if os.path.exists(Processes[p]) is False:  # If the process hasn't finished will return None
                del Processes[p]  # Remove from list - this is why we needed reverse order

    while (len(Processes) < maxjobs) and (jobindex < len(shlist)):  # More to do and some spare slots
        StartNew(shlist, maxjobs, chk_file_list)


def FSC2(input_dir, num_reps=50, min_sims=10000, max_sims=100000, conv_crit=0.001, min_ecm=10, max_ecm=40, calc_CI=False, partition="nbi-short", numcores=1, time="0-02:00", mem="8000", print1=False, overwrite="None", fsc2_path="/nbi/Research-Groups/JIC/Levi-Yant/Patrick/fsc_linux64/fsc25221", cluster="JIC", oande_dir="-99", verbose=False):
    """This method parallelises job submission of fastsimcoal2, but requires a very specific set up of input files.  The output of '.generateFSC2input' should be a folder that contains the multi-dimensional SFS.  Place this folder in a new folder that will be the FSC2_Data_Parent_Directory.  This directory should also contain one or more template (.tpl) and estimates (.est) files whose format can be found in the fastsimcoal2 documentation.  For each sub-directory containing input data, this method will re-format and rename the .tpl and .est files to reflect the necessary information in the sub-directory multi-dimensional SFS and then submit these jobs to the cluster.  I've tried to make the code as general as possible, but this is one method that will likely require the user to read and understand the code in order to get things working well for them.  Also, a major potential source of errors is in the correct formatting of the .tpl and .est files, so it is worthwhile to ensure that these are correct (by running FSC2 on a subset of your sub-directories) before launching full-scale"""
    Data_Files = []
    tpl_files = []
    est_files = []
    CI_Data_Files = []
    shlist = []
    chk_file_list = []

    if input_dir.endswith("/") is False:
        input_dir += "/"

    if oande_dir == "-99":
        oande_dir = input_dir

    for path in os.listdir(input_dir):
        if os.path.isdir(input_dir + path) and path.startswith("FSC2input"):
            samp_name = path.split("_")[1]
            if samp_name + "_DSFS.obs" in os.listdir(input_dir + path):
                for i in range(0, num_reps):
                    new_file = open(input_dir + path + "/" + samp_name + str(i) + "_DSFS.obs", 'w')
                    with open(input_dir + path + "/" + samp_name + "_DSFS.obs") as data_file:
                        for line in data_file:
                            new_file.write(line)
                        new_file.close()
                    Data_Files.append(input_dir + path + "/" + samp_name + str(i) + "_DSFS.obs")
            else:
                print("Did not find input data file for: ", samp_name)
            if calc_CI == "True":
                num_files = 0
                for file in os.listdir(input_dir + path):
                    if file.endswith("_DSFS.obs") and file.split("_")[-2].split(".")[-1][0:3] == "rep" and file != samp_name + "_DSFS.obs":
                        for i in range(0, num_reps):
                            new_file = open(input_dir + path + "/" + samp_name + file.split("_")[-2].split(".")[-1].split("_")[0]+ "_" + str(i) + "_DSFS.obs", 'w')
                            with open(input_dir + path + "/" + file) as data_file:
                                for line in data_file:
                                    new_file.write(line)
                                new_file.close()
                            CI_Data_Files.append(input_dir + path + "/" + samp_name + file.split("_")[-2].split(".")[-1].split("_")[0]+ "_" + str(i) + "_DSFS.obs")
                            num_files += 1
                if len(CI_Data_Files) < 1:
                    print("Did not find bootstrap replicates for: ", samp_name)
                else:
                    print("Found ", num_files, " replicate dsfs files for CI calculation for ", samp_name)
        if path.endswith(".tpl"):
            tpl_files.append(path)
            est_files.append(path.split(".")[0])
    if len(tpl_files) == 0:
        print("Did not find any tpl files!! Aborting!!")
    else:
        if calc_CI == "True":
            Data_Files = CI_Data_Files
        for file in Data_Files:
            name = file.split("_DSFS")[0]
            samp_name = name.split("/")[-1]
            for tpl in tpl_files:
                tpl_name = tpl.split(".tpl")[0]
                if os.path.isdir(name + "_" + tpl_name) is False or overwrite == "hard":
                    new_tpl = open(name + "_" + tpl_name + ".tpl", 'w')
                    new_data = open(name + "_" + tpl_name + "_DSFS.obs", 'w')
                    chk_file = open(name + "_" + tpl_name + ".chk.txt", 'w')
                    chk_file.close()
                    chk_file_list.append(name + "_" + tpl_name + ".chk.txt")
                    with open(file, 'r') as data:
                        print(file)
                        for i, line in enumerate(data):
                            print(i, line)
                            if i == 1:
                                pop_info = line.strip("\n").strip("\t").split("\t")
                                pop_num = int(pop_info[0])
                                samp_nums = pop_info[-pop_num:]
                            new_data.write(line)
                    with open(input_dir + tpl, 'r') as template:
                        samp_num_lines = pop_num + 4
                        for i, line in enumerate(template):
                            if i < samp_num_lines:
                                new_tpl.write(line)
                            elif i == samp_num_lines:
                                for num in samp_nums:
                                    new_tpl.write(num + "\n")
                            elif i >= samp_num_lines + len(samp_nums):
                                new_tpl.write(line)
                    new_est = open(name + "_" + tpl_name + ".est", 'w')
                    try:
                        with open(input_dir + tpl_name + ".est") as est:
                            for line in est:
                                new_est.write(line)
                    except FileNotFoundError:
                        print("Did not find est file for: ", tpl)
                    if cluster == "JIC":
                        shname = name + "_" + tpl_name + ".sh"
                        shlist.append(shname)
                        shfile5 = open(name + "_" + tpl_name + ".sh", 'w')
                        shfile5.write('#!/bin/bash\n' +
                                      '#SBATCH -J ' + name.split("/")[-1] + "_" + tpl_name + ".fsc2.sh" + '\n' +
                                      '#SBATCH -e ' + oande_dir + name.split("/")[-1] + "_" + tpl_name + ".fsc2.err" + '\n' +
                                      '#SBATCH -o ' + oande_dir + name.split("/")[-1] + "_" + tpl_name + ".fsc2.out" + '\n' +
                                      '#SBATCH -p ' + str(partition) + '\n' +
                                      '#SBATCH -n ' + str(numcores) + '\n' +
                                      '#SBATCH -t ' + str(time) + '\n' +
                                      '#SBATCH --mem=' + str(mem) + '\n' +
                                      'cd ' + os.path.abspath(os.path.join(file, os.pardir)) + "\n" +
                                      fsc2_path + ' -t ' + samp_name + "_" + tpl_name + ".tpl" + ' -e ' + samp_name + "_" + tpl_name + '.est -n ' + str(min_sims) + ' -N ' + str(max_sims) + ' -u -d -q -l ' + str(min_ecm) + ' -L ' + str(max_ecm) + ' -M ' + str(conv_crit) + ' \n' + 
                                      'rm ' + name + "_" + tpl_name + ".chk.txt")
                        shfile5.close()

                    if print1 == "True":
                        file3 = open(name.split("/")[-1] + tpl_name + ".sh", 'r')
                        data3 = file3.read()
                        print(data3)

                elif os.path.exists(name + "_" + tpl_name + "/" + samp_name + "_" + tpl_name + ".bestlhoods") is False:  # Intended to catch instances where FSC2 had run previously (and therefore created the output directory), but did not converge (and therefore output directory does not contain .bestlhoods file)
                    new_tpl = open(name + "_" + tpl_name + ".tpl", 'w')
                    new_data = open(name + "_" + tpl_name + "_DSFS.obs", 'w')
                    chk_file = open(name + "_" + tpl_name + ".chk.txt", 'w')
                    chk_file.close()
                    chk_file_list.append(name + "_" + tpl_name + ".chk.txt")
                    with open(file) as data:
                        for i, line in enumerate(data):
                            if i == 1:
                                pop_info = line.strip("\n").strip("\t").split("\t")
                                pop_num = int(pop_info[0])
                                samp_nums = pop_info[-pop_num:]
                            new_data.write(line)
                    with open(input_dir + tpl) as template:
                        samp_num_lines = pop_num + 4
                        for i, line in enumerate(template):
                            if i < samp_num_lines:
                                new_tpl.write(line)
                            elif i == samp_num_lines:
                                for num in samp_nums:
                                    new_tpl.write(num + "\n")
                            elif i >= samp_num_lines + len(samp_nums):
                                new_tpl.write(line)
                    new_est = open(name + "_" + tpl_name + ".est", 'w')
                    try:
                        with open(input_dir + tpl_name + ".est") as est:
                            for line in est:
                                new_est.write(line)
                    except FileNotFoundError:
                        print("Did not find est file for: ", tpl)
                    shname = name + "_" + tpl_name + ".sh"
                    shlist.append(shname)
                    shfile5 = open(name + "_" + tpl_name + ".sh", 'w')
                    shfile5.write('#!/bin/bash\n' +
                                  '#SBATCH -J ' + name.split("/")[-1] + "_" + tpl_name + ".fsc2.sh" + '\n' +
                                  '#SBATCH -e ' + oande_dir + name.split("/")[-1] + "_" + tpl_name + ".fsc2.err" + '\n' +
                                  '#SBATCH -o ' + oande_dir + name.split("/")[-1] + "_" + tpl_name + ".fsc2.out" + '\n' +
                                  '#SBATCH -p ' + str(partition) + '\n' +
                                  '#SBATCH -n ' + str(numcores) + '\n' +
                                  '#SBATCH -t ' + str(time) + '\n' +
                                  '#SBATCH --mem=' + str(mem) + '\n' +
                                  'cd ' + os.path.abspath(os.path.join(file, os.pardir)) + "\n" +
                                  fsc2_path + ' -t ' + samp_name + "_" + tpl_name + ".tpl" + ' -e ' + samp_name + "_" + tpl_name + '.est -n ' + str(min_sims) + ' -N ' + str(max_sims) + ' -u -d -q -l ' + str(min_ecm) + ' -L ' + str(max_ecm) + ' -M ' + str(conv_crit) + ' \n' +
                                  'rm ' + name + "_" + tpl_name + ".chk.txt")
                    shfile5.close()
                    if print1 == "True":
                        file3 = open(name + "_" + tpl_name + ".sh", 'r')
                        data3 = file3.read()
                        print(data3)
                else:
                    if verbose == "True": print("Output for " + samp_name + "_" + tpl_name + " already exists.  Use hard_overwrite = True to overwrite.")
    return shlist, chk_file_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, metavar='input_directory', required=True, help='input file created with recode012.py')
    parser.add_argument('-reps', type=int, metavar='Number_Replicates', required=True, help='Number of replicates for each scenario')
    parser.add_argument('-minsims', type=int, metavar='minimum_simulations_for_esfs', required=True, help='')
    parser.add_argument('-maxsims', type=int, metavar='maximum_simulations_for_esfs', required=True, help='')
    parser.add_argument('-min_ecm', type=int, metavar='minimum_ECM_steps', required=True, help='')
    parser.add_argument('-max_ecm', type=int, metavar='number_populations', required=True, help='')
    parser.add_argument('-nj', type=int, metavar='Number_of_Jobs_to_run', required=True, help='number of jobs to be submitted')
    parser.add_argument('-s', type=float, metavar='sleep_time', required=True, help='float designating amount of time to check on status of running jobs')
    parser.add_argument('-p', type=str, metavar='partition', required=True, help='partition to use on the cluster')
    parser.add_argument('-mem', type=int, metavar='amount of memory', required=True, help='integer designating amount of memory to use for each job')
    parser.add_argument('-t', type=str, metavar='job_time_limit', required=True, help='')
    parser.add_argument('-c', type=float, metavar='convergence_criteria', required=True, help='deprecated...do not modify')
    parser.add_argument('-nc', type=int, metavar='number_of_cores', required=True, help='number of cores to request for each job')
    parser.add_argument('-ci', type=str, metavar='Calculate_Confidence_intervals', required=True, help='Must have requested bootstrap replicate data sets from generateFSC2input in order to calculate confidence intervals.')
    parser.add_argument('-Ov', type=str, metavar='Overwrite_type', required=True, help='set to hard if you want to overwrite .bestlhoods files')
    parser.add_argument('-fsc2path', type=str, metavar='absolute_path_to_fsc2_executable', required=True, help='')
    parser.add_argument('-print1', type=str, metavar='print_shell_scripts?', required=True, help='boolean designating whether shell scripts should be submitted or printed')
    parser.add_argument('-verbose', type=str, metavar='verbose', required=True, help='prints verbose output')
    parser.add_argument('-clust', type=str, metavar='Cluster', required=True, help='which cluster to submit to.  e.g. JIC')
    parser.add_argument('-oande', type=str, metavar='output_and_error_file_directory', required=True, help='')
    args = parser.parse_args()

    shell_script_list, Chk_files = FSC2(input_dir=args.i, num_reps=args.reps, min_sims=args.minsims, max_sims=args.maxsims, conv_crit=args.c, min_ecm=args.min_ecm, max_ecm=args.max_ecm, calc_CI=args.ci, partition=args.p, numcores=args.nc, mem=args.mem, time=args.t, print1=args.print1, overwrite=args.Ov, fsc2_path=args.fsc2path, cluster=args.clust, oande_dir=args.oande, verbose=args.verbose)
    print(shell_script_list)
    print(Chk_files)
    CheckRunning(shell_script_list, args.nj, Chk_files)  # This will start the max processes running
    while (len(Processes) > 0):  # Some thing still going on.
        time.sleep(float(args.s))  # You may wish to change the time for this
        CheckRunning(shell_script_list, args.nj, Chk_files)
