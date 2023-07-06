## Rashx R Package
The Rashx package provides a function process_requests for processing requests using various R packages. This function takes several parameters and performs data processing and analysis tasks related to gene expression data. It utilizes packages such as Seurat, metap, MAST, scDblFinder, ggplot2, cowplot, patchwork, stringr, RColorBrewer, ComplexHeatmap, circlize, tidyverse, concaveman, ggforce, monocle3, dplyr, Rmisc, ggrepel, raster, vegan, and plyr.

## Installation
To install the Rashx package from GitHub, you can use the remotes package in R. Run the following code in R:

```
install.packages("remotes")  
remotes::install_github("username/repo")  

```

Replace "username/repo" with the actual repository name where the Rashx package is hosted.

## Usage
The main function in the Rashx package is process_requests. Here is a brief description of its parameters:

rds1: The path to the RDS file.
exported_file_name: The name of the exported CSV file.
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
