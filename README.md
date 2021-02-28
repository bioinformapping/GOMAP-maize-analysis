# The Repository for GOMAP maize analysis paper

This contains all the source code that is necessary to run the analysis that was conducted for the GOMAP paper that is was published in [bioRxiv](https://www.biorxiv.org/content/10.1101/809988v2)

The scripts are mainly written in [R](https://www.r-project.org/) with one script in [Python](https://www.python.org/)

The analysis is separated into 7 different categories in the order of the appearance of the figures and tables in the manuscript 

## Script Categories

1. SRA
	The code to download the SRA data and create the SRA dataset plot numbers
2. Walltime
	The code to create the table and figure related to the GOMAP performance 
3. Get Annotations
	The code process the GOMAP and Community annotations and prepare the expanded annotations by determining ancetral GO terms via the GO hierarchy
4. Analysis Metrics
	Code to process the fasta files and GO annotations to generate the analysis metrics for GOMAP and Community annotations
5. Evaluation Metrics
    The code to generate the CAFA metrics using the gold-standard standard datasets for the GOMAP and Community datasets
6. Get plant-specific annotations
    The code to clean-up and generate plant-specific annotations 
7. Others
	The scripts that were used to generate any other tables, figures, and citations

## Obtaining Plant-specific Annotations from All GOMAP Annotations

This functionality is provided via R scripts that have been containerized using singularity and the instructions are available at [GOMAP-PlantSpecific](https://github.com/bioinformapping/GOMAP-PlantSpecific).
