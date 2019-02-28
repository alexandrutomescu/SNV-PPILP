# SNV-PPILP
SNV-PPILP: Refined SNV calling for tumor data using perfect phylogenies and ILP (version 1.2)

## Running SNV-PPILP

The command to run SNV-PPILP is:

    python SNV-PPILP.py

Usage: 

    SNV-PPILP.py [-h] -i VCF -o OUTPUTFILE [-ilp LPSOLVEPATH] [-f F]

Required arguments:
    
    -i VCF            GATK's multi-sample .vcf file
    -o OUTPUTFILE     output file

Optional arguments:

    -h, --help        Show this help message and exit
    -ilp LPSOLVEPATH  Path to the directory containing lp_solve (this argument
                      is not needed if lp_solve is in the PATH variable)
    -f F              An SNV non detected by GATK in one sample gets weight =
                      (the average quality score of that SNV over all samples) -
                      (sqrt(100/F)) * (the standard deviation of that SNV over
                      all samples) [default=50]
    -hc               Set this flag if the vcf file is the output of GATK's
                      HaplotypeCaller

## Input and output

SNV-PPILP takes as input .vcf file produced by GATK's
Unified Genotyper run in multi-sample mode.

The output is a unique file in .csv format containing, for each sample,
the list of CHROM names and POS ID's *ONLY* of the SNVs present in that 
sample. SNV-PPILP does not report indels.


## Requirements

SNV-PPILP uses the free ILP solver called 'lp_solve'. This is included 
with SNV-PPILP (for Linux, OS X and Windows), and alternatively can be 
downloaded from: http://sourceforge.net/projects/lpsolve/

We have tested SNV-PPILP on the latest version 5.5.0.*, but previous
versions may work without problems. If SNV-PPILP runs and produces an
output, then that output is correct and your version of lp_solve is
supported by SNV-PPILP (see below for an example usage and an example
output file).

Note: lp_solve may already be installed on your system. You can check
this by simply trying to run: lp_solve -h

## Example usage   

A directory sample/ containing GATK's Unified Genotyper .vcf file from
6 samples is included with this distribution. This is how to run the
tool on this sample:

    python SNV-PPILP.py -ilp ./lp_solve_5.5.2.0_exe_ux64/ -o ./sample/6samples.corrected.csv -i ./sample/6samples.vcf 

This command assumes that lp_solve is in the directory called
lp_solve_5.5.2.0_exe_ux64/ located in the current folder. (The 
folder lp_solve_5.5.2.0_exe_ux64/ included in this distribtuion 
contains the lp_solve binaries ofr linux). You may need to set
execution rights to the binaries of lp_solve.
If lp_solve is already installed on your system, then probably 
it is already in the PATH variable, in which case you can skip 
the parameter -ilp.

The output is the file ./sample/6samples.corrected.csv.

## HaplotypeCaller

If your .vcf file is the output of GATK's HaplotypeCaller, then set the flag **-hc** when running SNV-PPILP. For a test run, try running SNC-PPILP the sample file ./sample/haplotypecaller.vcf, as follows

    python SNV-PPILP.py -hc -ilp ./lp_solve_5.5.2.0_exe_osx32/ -o ./sample/haplotype.corrected.csv -i ./sample/haplotypecaller.vcf

## Acknowledgements

Thanks to Harald Detering (University of Vigo, Spain) for writing the support for HaplotypeCaller and providing the sample file samples/haplotypecaller.vcf.

## Citation and contact         

If you use SNV-PPILP, please cite

Karen E. van Rens, Veli Mäkinen, Alexandru I. Tomescu, SNV-PPILP: Refined SNV calling using perfect phylogenies and ILP, Bioinformatics, Volume 31, Issue 7, 1 April 2015, Pages 1133–1135, DOI: [10.1093/bioinformatics/btu755](https://doi.org/10.1093/bioinformatics/btu755)


https://www.cs.helsinki.fi/gsa/SNV-PPILP/
alexandru.tomescu@helsinki.fi
