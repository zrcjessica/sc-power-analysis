library(ggplot2)
library(dplyr)
library(magrittr)
library(optparse)

# load arguments
option_list <- list(
	make_option("--power_analysis", action = "store", 
		type = "character", 
		help = "path to file with power analysis results"),
	make_option("--sample", action = "store", type = "character",
		default = "Rat_Amygdala_787A_all_seq",
		help = "sample name"),
	make_option("--sample_variable", action = "store", type = "character",
		default = "replicates",
		help = "from {replicates, cells} - define sample size that is varied for power analysis"),
	make_option("--out", action = "store", 
		type = "character",
		help = "output filename"),
	make_option("--plot_width", action = "store", type = "integer",
		default = 5,
		help = "width of plot"),
	make_option("--plot_height", action = "store", type = "integer",
		default = 3,
		help = "height of plot")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!(opt$sample_variable %in% c('replicates','cells'))) {
	stop("sample_variable argument must be one of 'replicates' or 'cells'")
}

df <- read.table(opt$power_analysis, sep = ",", header = FALSE, 
	col.names = c('replicates','cells','effect','power'))

if (opt$sample_variable == "replicates") {
	df %<>% group_by(replicates, effect) %>% summarize(power = mean(power))

	df$effect %<>% as.factor

	p <- ggplot(data = df, aes(x = replicates, y = power)) + 
	geom_line(aes(color=effect)) + geom_point(aes(color=effect)) +
	theme_minimal()

	pdf(opt$out, width = opt$plot_width, height = opt$plot_height)
	print(p)
	dev.off()
} else {
	df %<>% group_by(cells, effect) %>% summarize(power = mean(power))

	df$effect %<>% as.factor

	p <- ggplot(data = df, aes(x = cells, y = power)) + 
	geom_line(aes(color=effect)) + geom_point(aes(color=effect)) +
	theme_minimal()

	pdf(opt$out, width = opt$plot_width, height = opt$plot_height)
	print(p)
	dev.off()
}

