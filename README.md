
#Overview
CopywriteR is a tool used to extract CNAs from off-target reads in targeted sequencing data. This pipeline automates the process of running CopywriteR, processing its output, and visualizing the results. It handles both sample and control data, performs data cleaning, and produces segmented copy number profiles.

#Prerequisites
- R version 3.5 or higher
- CopywriteR package
- DNAcopy package
- Other R packages: Rcpp, data.table, reshape2
- Pipeline Workflow
  
#Setting up the environment:

- Specifying library paths and loading necessary R packages.
- Command line arguments for sample and normal data are parsed and used throughout the script.

#Data preparation:

- Paths to sample and control BAM files are constructed and combined into a dataframe.
- Unique entries are ensured to avoid redundant calculations.

#Running CopywriteR:

- The CopywriteR function is invoked with specified parameters, including paths to sample and reference data, and output directories.
- Post-run, the CNA profiles are plotted for visualization.

#Data post-processing:
- Segmented data files are renamed and organized.
- Additional segmentation adjustments are made using custom thresholds and conditions to refine CNA detection.

#Visualization and Output:

- Segmentation results are visualized in PDF format to assess the quality and distribution of CNAs.
- Segmentation data are saved in text format for further analysis.

#Advanced Analysis Options:
- Scripts include functions for re-estimating segment copy numbers using various statistical approaches.
- Functions are provided to collapse consecutive segments with the same copy number and to drop small segments.
Usage
Configuration:

#Ensure all file paths and parameters are correctly set in the script to match your computing environment and data locations.
Execution:

Run the script in an R environment capable of handling the specified dependencies.
Command line arguments must be provided for the script to run correctly, specifying the paths to sample and control data.

#Output:

- Check the specified output directory for results, including tables of CNA segments and visualizations in PDF format.
- Files and Directories
- Input files: BAM files for samples and controls.
- Output files:
- Segmented copy number data (*.seg.txt files).
- Visual plots of segmented CNAs (*.pdf files).
- Detailed logs and parameter estimates for the analysis (ParameterList.txt).

#Troubleshooting
- Ensure that all required libraries are installed and properly loaded.
- Validate file paths and permissions to read input files and write outputs.
- Check R console output for error messages that can indicate what might have gone wrong during execution.
