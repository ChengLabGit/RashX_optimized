## RashX R Package
The RashX package provides a function process_requests to classify cutaneous human immune cells taken from skin rashes as either psoriatiform or atopic dermatiform based on their single-cell gene expression profiles. This function takes an input RDS created from default Seurat pre-processing and performs the classification.  It utilizes dependendies Seurat, metap, MAST, ggplot2, cowplot, patchwork, stringr, RColorBrewer, ComplexHeatmap, circlize, tidyverse, concaveman, ggforce, monocle3, dplyr, Rmisc, ggrepel, raster, vegan, and plyr.  Output includes heatmaps of atopic dermatitis and psoriosis gene signature expression; figures depicting query-sample classification in a reduced gene expression space; and statistical results reporting the significance of the proximity of a query sample to reference atopic dermatitis and psoriasis samples in this reduced gene expression space.

## Installation
To install the Rashx package from GitHub, you can use the remotes package in R. Run the following code in R:

```
install.packages("remotes")  
remotes::install_github("username/repo")
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')

```

Replace "username/repo" with the actual repository name where the Rashx package is hosted.

## Usage
The main function in the Rashx package is process_requests. Here is a brief description of its parameters:

rds1: The path to the user's query RDS file containing Seurat-processed cutaneous immune cell scRNAseq for which classification will be performed.  The expression matrix must be log-normalized (LogNormalize c.f. SCT in Seurat).
exported_file_name: The name of the exported CSV file containing statistical results of the classification.
plot1_name: The name of the AD-specific genes heatmap plot.
plot2_name: The name of the PV-specific genes heatmap plot.
plot3_name: The name of the sample-level means plot without individual cells.
plot4_name: The name of the sample-level means plot with individual cells.
original_file_name: The name of the original file.
To use the process_requests function, provide the required parameters and call the function as shown below:

```
process_requests(rds1, exported_file_name, plot1_name, plot2_name, plot3_name, plot4_name, original_file_name)
```
Make sure to replace the parameter values with the actual file paths and names.

The function will perform various data processing and analysis tasks, generate plots, and export a CSV file with the statistical analysis results.

Note: The specific implementation and functionality of the process_requests function may vary depending on the data and analysis requirements.
