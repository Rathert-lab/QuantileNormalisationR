#!/usr/bin/env Rscript

# author: "Jan Mueller (11.12.19)"
# Last modifictaion: 
# - Jan Mueller (12.12.19)
# - Jan Mueller/Philipp Rathert (13.12.19)
# - Jan Mueller (07.04.21)



# > sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 18363)

# Matrix products: default

# Random number generation:
#  RNG:     Mersenne-Twister 
#  Normal:  Inversion 
#  Sample:  Rounding 
 
# locale:
# [1] LC_COLLATE=German_Germany.1252  LC_CTYPE=German_Germany.1252    LC_MONETARY=German_Germany.1252
# [4] LC_NUMERIC=C                    LC_TIME=German_Germany.1252    

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] preprocessCore_1.52.1 GGally_2.1.1          optparse_1.6.6        reshape2_1.4.4       
# [5] ggplot2_3.3.2        

# loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.5         rstudioapi_0.13    magrittr_2.0.1     getopt_1.20.3      tidyselect_1.1.0  
#  [6] munsell_0.5.0      colorspace_2.0-0   R6_2.5.0           rlang_0.4.9        stringr_1.4.0     
# [11] plyr_1.8.6         dplyr_1.0.2        tools_4.0.3        grid_4.0.3         gtable_0.3.0      
# [16] withr_2.3.0        ellipsis_0.3.1     tibble_3.0.4       lifecycle_0.2.0    crayon_1.3.4      
# [21] RColorBrewer_1.1-2 purrr_0.3.4        vctrs_0.3.5        glue_1.4.2         stringi_1.5.3     
# [26] compiler_4.0.3     pillar_1.4.7       generics_0.1.0     scales_1.1.1       reshape_0.8.8     
# [31] pkgconfig_2.0.3



# Load required packages
## If a package is installed, it will be loaded. If any 
## are not, the missing package(s) will be installed 
## from CRAN and then loaded.
if (!require(pacman)) install.packages('pacman', repos='http://cran.us.r-project.org')
if (!require(BiocManager)) install.packages('BiocManager', repos='http://cran.us.r-project.org')

pacman::p_load(ggplot2, reshape2, optparse, GGally, preprocessCore)


##get arguments
option_list = list(
  make_option(c("-f", "--filelist"), type="character", default=NULL, 
              help="Identically sorted bedgraph files. Paths seperated by ',' (required)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="Output directory for normalized files and plots.", metavar="character")

); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
input_files = strsplit(opt$filelist,',')[[1]]

if (length(input_files)<2){
  print_help(opt_parser)
  stop("At least two bedgraph files must be supplied for normalization.", call.=FALSE)
}


##read input files
ldf = list()

for (k in 1:length(input_files)){
 ldf[[k]] = read.delim(input_files[k], header=FALSE)
}

##extract read counts from column 4
reads_raw = cbind(as.matrix(ldf[[1]])[,4])

for (k in 2:length(input_files)){
 reads_raw = cbind(reads_raw,as.matrix(ldf[[k]])[,4])
}

class(reads_raw) = 'numeric'

##Quantile normalization
reads_normalized = normalize.quantiles(reads_raw)

##Data export
file_name_list = c()
for (k in 1:length(input_files)){
 export = cbind.data.frame(ldf[[k]][,1:3],reads_normalized[,k])
 file_name = strsplit(tail(strsplit(input_files[k], '/')[[1]], 1), '.', fixed=TRUE)[[1]][1]
 file_name_list = c(file_name_list, file_name)
 write.table(export, paste(gsub('"', '',opt$out), paste(file_name,'.normalized.bedgraph', sep =''), sep ='/'), sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

###plotting needs to be corrected: factorization after melt + naming of boxes (saves file_name from upper loop to a list)

##QC plots
plot_raw = melt(reads_raw)
plot_raw$Var2 <- as.factor(plot_raw$Var2)
levels(plot_raw$Var2) = file_name_list

plot_normalized = melt(reads_normalized)
plot_normalized$Var2 <- as.factor(plot_normalized$Var2)
levels(plot_normalized$Var2) = file_name_list


reads_normalized_df = as.data.frame(reads_normalized)
names(reads_normalized_df) = file_name_list
reads_raw_df = as.data.frame(reads_raw)
names(reads_raw_df) = file_name_list

# compute lower and upper whiskers (of control) which are used for the scaling of the plot + converting the reads for matrix plots
ylim1 = boxplot.stats(reads_normalized[,1])$stats[c(1, 5)]
ylim2 = boxplot.stats(reads_raw[,1])$stats[c(1, 5)]


pdf(paste(gsub('"', '',opt$out),"QC_plots_quantile-norm.pdf" , sep ='/'))


# #Boxplots
ggplot(data = plot_raw, aes(x =Var2, y = value, group = Var2, fill = Var2)) + geom_boxplot(notch = TRUE, show.legend = FALSE) +
	coord_cartesian(ylim = ylim1*1.05) + 
	labs(title = 'Raw Data', x = ' ', y = 'Counts') + 
	theme_classic()	+ 
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
	scale_fill_brewer(palette="BuPu")

ggplot(data = plot_normalized, aes(x =Var2, y = value, group = Var2, fill = Var2)) + geom_boxplot(notch = TRUE, show.legend = FALSE) +
	coord_cartesian(ylim = ylim1*1.05) + 
	labs(title = 'Quantile normalized Data', x = ' ', y = 'Counts') + 
	theme_classic()	+ 
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
	scale_fill_brewer(palette="BuPu")

#Scatterplot Matrix
geom_bin_plot <- function(data, mapping, ..., low = "#132B43", high = "#56B1F7") {
  ggplot(data = data, mapping = mapping) +
    geom_bin2d(bins=500) +
    scale_fill_gradient(low = low, high = high) +
    scale_x_log10(limits = c(0.1,1e7)) +
    scale_y_log10(limits = c(0.1,1e7)) + 
    geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=0.5)
}

ggpairs(reads_raw_df, title = 'Raw Data', 
  upper = list(continuous = 'cor'),
  lower = list(continuous = geom_bin_plot),
  diag = 'blank'
)

ggpairs(reads_normalized_df, title = 'Quantile normalized Data', 
  upper = list(continuous = 'cor'),
  lower = list(continuous = geom_bin_plot),
  diag = 'blank'
)

dev.off()