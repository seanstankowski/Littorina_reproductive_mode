This file contains all of the code needed to infer ancestral recombination graphs using the program ARGweaver. 

ARG_inference/scripts
|-- analysis.Rmd                                    # R notebook for step-by-step analysis of ARGweaver output
|-- bin
|   |-- functions_from_argweaver.R                  # Functions for running analysis.Rmd (modified from original ARGweaver repo)
|   |-- functions_v2.R                              # Functions for running analysis.Rmd (version 2, written by Arka Pal)
|   `-- vcf2sites.R                                 # R function to convert .vcf to .sites format used by ARGweaver
|-- pipeline_Littorina.sh                           # Pipeline running ARGweaver on Littorina samples
`-- slurm                                           # Folder containing SLURM scripts for submitting HPC cluster jobs
    |-- job-resumeArgweaver_masked.slurm            # Resume ARGweaver run on sequneces masked with missing sites 
    |-- job-resumeArgweaver.slurm                   # Resume ARGweaver run
    |-- job-runArgweaver_masked.slurm               # Initiate ARGweaver run on sequneces masked with missing sites
    `-- job-runArgweaver.slurm                      # Initiate ARGweaver run