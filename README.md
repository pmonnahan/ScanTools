# ScanTools

----
## Introduction

  ScanTools is a collection of scripts for window-based genomic analyses and a wrapper that facilitates job submission on a slurm-based computing cluster.  The program begins with a set of VCF's to be analyzed as well as a population key (example template can be found on github repository) that assigns the individual names in the VCF's to populations.  VCF's are split according by a population and converted to a simplified format, which are used for all downstream analyses.  Currently, the major downstream analyses implemented are calculation of within population diversity metrics and neutrality tests, between population differentiation metrics, and demographic analyses using fastsimcoal2.

----
## Requirements

- python3.+ (program has been debugged for python-3.5.1)
- python modules: os, subprocess, pandas, math, datetime, time
- bedtools 2.17 or higher
- GATK 3.6 
- fastsimcoal2 version 2.5.2.21
- A directory containing a file named 'PopKey.csv'.

----
## Getting started

Download the program with:

    git clone https://github.com/pmonnahan/ScanTools.git


Prior to initiating the program for the first time, create a directory to store output for a project and, into this folder, place a file named 'PopKey.csv'.  This file associates individual names in the input VCF with population names.  See *PopKey_Template.csv* for an example of how to format this file.  

To initiate the program, enter the *ScanTools* directory and invoke an instance of python3:

    cd <path>/ScanTools
    python3
 
Import the program into python3 with:
    
    >>> from ScanTools import scantools 

Using the wrapper (ScanTools.py) is greatly aided by a rudimentary understanding of python, particularly regarding the use of 'methods' associated with 'objects'.  The first step is to create a 'ScanTools' object and all subsequent operations will be performed by calling the methods associated with this object. 

    >>> project = scantools("<path_to_directory_containing_PopKey.csv>", "encoding")

This creates a ScanTools object named 'project' and assigns several relevant pieces of information to the methods of 'project'.  If the *PopKey.csv* is encoded, you will have to specify the type of encoding as the second argument.  However, this should typically not be the case.  All output from any subsequent steps will be output to the directory specified in the initialization of the ScanTools object.  

Methods in python can both store information about the object as well as perform operations to the object.  For example, the list of populations and individuals are stored in the *.pops* and *.samps* methods, respectively.  This information can be retrieved by:

    >>> project.pops
    ['Pop1','Pop2','Pop3']
    >>> project.samps
    {'Pop1': ['Ind1', 'Ind2'], 'Pop2': ['Ind3', 'Ind4', 'Ind5'], 'Pop3': ['Ind6']}

The first step before doing any analysis is to split the input VCFs by population and convert them to the format used by the program. This can be done by calling:

    >>> project.splitVCFs(vcf_dir=<path_to_vcfs>, min_dp=<minimum_depth_for_genotype_call>, mffg=<maximum_fraction_filtered_genotypes>)

In this case, *splitVCFs* is a method of the *project* object that performs an operation based on the information associated with *project*.  Specifically, this method makes several calls to GATK to filter the VCFs and convert them to table format.  Subsequently, it calls either one or two custom python scripts to convert the genotype calls to numeric format.  All subsequent analyses can be carried out similarly by calling the relevant methods.  The full list of methods and a brief description is provided below, but look within the python scripts themselves for a more detailed description.  

## Methods:

Within python3, you can view all methods associated with an object by calling 

    dir(<object_name>)

Ignore methods surrounded by underscores...these are 'intrinsic' methods that are not directly relevant for our purposes.  One more important note on methods in python is that they are quite versatile.  For example, they can be functions that take arguments and carry out operations associated with the parent object.  Or, methods can objects that store information associated with the parent object.  

Several of the ScanTools methods will return information about the parent object (*project* in the above example) when called.  Below :

- **.pops**: a list of the populations identified in the PopKey.csv file.  These populations will be used as default for subsequent analyses.
- **.samps**: a dictionary of samples and their associated populations identified in the PopKey.csv file
- **.samp_nums**:  a dictionary of populations and their associated sample sizes
- **.min_ind**: the minimum observed number of individuals in a population
- **.dir**: the directory containing the PopKey and all output
- **.split_dirs**: list of directories containing reformatted input files, if any exist.


## Functions


**Important** 
Most of these methods include a *print1* argument that, if set to *True*, will print the shell scripts instead of submitting them to the cluster.  This is a useful check that you should always do before your first real execution of a method.  The argument names provided below are not the actual argument names required by each function.  Rather, the arguments below are descriptions of the actual arguments in the program.  Also, most of the methods have an option to change the partition, time, and memory requested for each job.  For methods that necessarily generate large intermediate files (.splitVCFS(), .calcbpm(), .generateFSC2input()), there is an option to use scratch (*use\_scratch*=True and *scratch\_path*=/path/to/scratch/directory).  However, this requires that the user to setup their directory on the scratch drive.

 - **.removePop**([pops\_to\_be\_removed]): takes a list of populations and removes them from '.pops'. E.g., 

        project.removePops(['Pop1','Pop2'])

  Even single populations should be specified in a list format.  '.min_ind' will be recalculated.

 - **.removeInds**([individuals\_to\_be\_removed]): takes a list of individuals and removes them from '.samps'.  Recalculates '.min_ind'

 - **.combinePops**(pop\_list, new\_pop\_name): combines two or more populations into new population with the name specified.  Original populations remain unchanged.

 - **.splitVCFs**(vcf\_directory, minimum\_individual\_depth, max\_fraction\_filtered\_genotypes):  Split VCF's by population, filter, and convert to input format for downstream analyses.

 - **.recode**(table\_directory): Should not be necessary for the most part.  This will typically be called during '.splitVCFs', but was left in code as standalone method in case there is need to convert a table (from GATK's VariantsToTable) to the reformatted input used in ScanTools.

 - **.repolarize**(recoded\_file\_directory):  Also should not be necessary for most part as it is typically called in '.splitVCFs' if a repolarization key is provided.  A repolarization key is simply a tab-delimited list of sites where each line contains the scaffold and position for a site.  This code will go through the recoded, numeric genotype files and flip the reference and alt alleles, modifying genotypes accordingly.

 - **.getPloidies**(recoded\_file\_directory):  Determines ploidy of each population from the recoded files and assigns them to a new method of the ScanTools parent object.

 - **.calcAFS**(recoded\_file\_direcory, output\_file\_suffix, number\_of\_individuals\_to\_downsample\_to, list\_of\_populations, allow\_one\_missing=True):  Calls calcAFS.py for each population in list_of_populations.  Calculates genome-wide allele frequency spectrum.

 - **.calcFreqs**(self, recode\_dir, outfile\_name, sites\_file, list\_of\_populations):  Calls *calcFreqs\_atSites.py*.  Takes a list of sites in (sites\_file, should be formatted so that each line simply has scaffold and position, with scaffold simply coded as an integer 0-8) and calculates the allele frequency in each population (list\_of\_populations) at each site.

 - **.calcwpm**(recoded\_file\_directory, window\_size, minimum\_snps):  Calls wpm.py, which calculates within-population diversity metrics and neutrality test statistics in windows along the genome as well as genome-wide for each population.  Windows are specified in terms of base pairs.  Windows with fewer than the specified minimum SNPs will not be reported.  All populations are downsampled by default to the '.min_ind' value unless specified otherwise using 'sampind' argument.

 - **.concatWPM**(directory\_containing\_wpm\_output, suffix, name):  Concatenates output of wpm.py across populations.  Suffix is a string that identifies the name of the output files to concatenate that follows the population name.  'name' is used to name the concatenated file.

 - **.calcBPM**(recoded\_file\_directory, pops, output\_name, window\_size, minimum\_snps):  Calls bpm.py, which calculates between-population differentiation metrics: Fst, Rho (for interploidy comparisons), dxy, fixed differences, and allele frequency difference in windows along the genome.  This can calculate differentiation for two or more populations, but will get seriously bogged down if the number of populations is large.  For this reason, you must explicitly provide a list of populations as an argument.

 - **.calcPairwisebpm**: same as above, but instead calculates metrics for every pair of populations passed in the 'pops' argument.

 - **.findOutliers**(directory, filename, column_index_list, percentile): Identifies outliers from one or more metrics in a provided output file.  A list of column indices (1-based; i.e. first column is column 1 not column 0) identifying the columns containing metrics to be used for determining outliers.

 - **.annotateOutliers**(directory,outlier\_file): annotates 'OutOnly.bed' file with annotation information from gff file using bedtools.

 - **.mergeAnnotation**(directory, filename):  filename should be the 'OutOnly.csv' file generated from '.findOutliers'.  This method merges the annotation information from '.annotateOutliers' with the metric info in OutOnly.csv.

 - **.Outliers**(directory, metric\_file, column\_index\_list, percentile):  wraps '.findOutliers', '.annotateOutliers', and '.mergeAnnotation' into one method

 - **.generateFSC2input**(recoded\_file\_directory, pops, output\_name, bootstrap\_block\_size, bootstrap\_reps): 
Generates the multi-dimensional site frequency spectrum input for fastsimcoal2 from the recoded, numeric genotype files (Output is labeled as output\_name\_DSFS.obs).  Has been tested for 5 populations, but may require much more memory as the number of populations are increased.  The bootstrap arguments are for performing confidence interval calculation as specified in the fastsimcoal manual.  Set bootstrap reps to 0 if you don't want this to happen.  Otherwise, you'll notice many additional files labelled with 'rep#' in the name, which are the bootstrapped replicate data sets.

 - **.FSC2**(FSC2\_Data\_Parent\_Directory):  
This method parallelises job submission of fastsimcoal2, but requires a very specific set up of input files.  The output of '.generateFSC2input' should be a folder that contains the multi-dimensional SFS.  Place this folder in a new folder that will be the FSC2_Data_Parent_Directory.  This directory should also contain one or more template (.tpl) and estimates (.est) files whose format can be found in the fastsimcoal2 documentation.  For each sub-directory containing input data, this method will re-format and rename the .tpl and .est files to reflect the necessary information in the sub-directory multi-dimensional SFS and then submit these jobs to the cluster.  I've tried to make the code as general as possible, but this is one method that will likely require the user to read and understand the code in order to get things working well for them.  Also, a major potential source of errors is in the correct formatting of the .tpl and .est files, so it is worthwhile to ensure that these are correct (by running FSC2 on a subset of your sub-directories) before launching full-scale

 - **.gatherFSC2output**(FSC2\_Data\_Parent\_Directory):  This method collects all information from the '.bestlhood' output files of FSC2 that is buried in the sub-directories and outputs the information into one of two files:  Likelihoods file and parameters file.


