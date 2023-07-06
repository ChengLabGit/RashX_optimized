\name{process_requests}
\alias{process_requests}

\title{Process Requests}

\description{
  This function processes requests using various R packages.
}

\usage{
  process_requests(rds1, exported_file_name, plot1_name, plot2_name, plot3_name, plot4_name, original_file_name)
}

\arguments{
  \item{rds1}{
    The path to the RDS file.
  }

  \item{exported_file_name}{
    The name of the exported CSV file.
  }

  \item{plot1_name}{
    The name of the AD-specific genes heatmap plot.
  }

  \item{plot2_name}{
    The name of the PV-specific genes heatmap plot.
  }

  \item{plot3_name}{
    The name of the sample-level means plot without individual cells.
  }

  \item{plot4_name}{
    The name of the sample-level means plot with individual cells.
  }

  \item{original_file_name}{
    The name of the original file.
  }
}

\details{
  The function requires the following packages:
  \itemize
  \item Seurat
  \item metap
  \item MAST
  \item scDblFinder
  \item ggplot2
  \item cowplot
  \item patchwork
  \item stringr
  \item RColorBrewer
  \item ComplexHeatmap
  \item circlize
  \item tidyverse
  \item concaveman
  \item ggforce
  \item monocle3
  \item dplyr
  \item Rmisc
  \item ggrepel
  \item raster
  \item vegan
  \item plyr
}

\value{
  A message indicating the processing is complete.
}
