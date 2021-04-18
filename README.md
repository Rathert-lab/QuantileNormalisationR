# QuantileNormalisationR
## Quantile Normalisation of .bedgraph files

Performs quantile normalizationf of ChIP-seq read data stored in bedgraph files using R. The script employs the quantile normalization function of the R package ‘preprocessCore’ for its main functionality. A version for use on Windows and Linux available.

- R based
- Windows and Linux Version available
- Quantile normalization of arbitrary number of bedgraph files
- Produces plots for quality control of normalization (Scatter plot and correlation matrix and Boxplot of all input files before and after normalization)


## Usage
The QuantileNormalisationR requires an installation of R and five R packages to run. All packages should be installed automatically when running the script if they are not installed already.
- ggplot2
- reshape2
- optparse
- GGally
- preprocessCore

ChIP-seq read data files used as input should be formatted in the bedgraph format with readcount values in the fourth coloumn. All files for normalization should contain the same number of datapoints (same bin size and same genome).
```
chromA  chromStartA  chromEndA  dataValueA
chromB  chromStartB  chromEndB  dataValueB
```

The script can be executed from a Windows/Linux command terminal or from a terminal in RStudio. Use the following command to run the script in Windows (for Linux exchange "\\" with "/" and use the linux version of the script 'quantile_normalization_bedgraph_linux.R'):

```sh
C:\Program Files\R\R-4.0.3\bin\Rscript.exe ~\quantile_normalization_bedgraph_windows.R -f ~\Control.bedgraph,~\Experiment1.bedgraph,~\Experiment2.bedgraph -o ~\Normalization_output
```
The QuantileNormalisationR script takes two parameters:

- -f : A comma separated list of any length of the bedgraph filepaths for normalization
- -o : The output path for the normalized bedgraph files and the quality control plots


## License

MIT
