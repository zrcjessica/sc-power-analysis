library(Matrix)
library(Seurat)
library(optparse)

option_list <- list(
	make_option("--count_mtx", action = "store", 
		default = '/iblm/netapp/home/jezhou/power_analysis/data/snRNA/Rat_Amygdala_787A_all_seq/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
		type = "character", 
		help = "path to cellranger filtered_feature_bc_matrix/matrix.mtx.gz file"),
	make_option("--params", action = "store",
		default = "/iblm/netapp/home/jezhou/power_analysis/outs/Rat_Amygdala_787A_all_seq_nb_params.csv",
		help = "path to file with NB params estimated from real rat data"),
	make_option("--out", action = "store",
		type = "character",
		help = "output file name"),
	make_option("--replicates", action = "store", type = "integer",
		help = "number of replicates to simulate for each condition"),
	make_option("--nGenes", action = "store", type = "integer",
		default = 10000,
		help = "number of genes to simulate"),
	make_option("--pctDiff", action = "store", type = "double",
		default = 0.05,
		help = "percentage of DEGs to simulate"),
	make_option("--nSims", action = "store", type = "integer",
		default = 15,
		help = "number of simulations"),
	make_option("--effectSize", action = "store", type = "double",
		help = "effect size to simulate"),
	make_option("--nCells", action = "store", type = "integer",
		default = 5000,
		help = "number of cells to simulate for each replicate"),
	make_option("--alpha", action = "store", type = "double",
		default = 0.05,
		help = "significance level for power analysis"),
	make_option("--test", action = "store", type = "character",
		default = "wilcox",
		help = "test to use for Seurat FindMarkers()")
)

simulateCounts <- function(reps, params, depths, nCells, degs, effect, nGenes) {
	ctrl.counts <- list()
	ctrl.bcs <- c()

	exp.counts <- list()
	exp.bcs <- c()

	for (i in c(1:reps)) {

		ctrl.counts[[i]] <- t(apply(params, 1, function(x) {
			rnbinom(nCells, mu = x[1]/depths, size = x[2])}))
		ctrl.bcs <- c(ctrl.bcs, paste("ctrl", c(1:nCells), i, sep = "_"))

		exp.de <- t(apply(params[c(1:degs),], 1, function(x) {
			rnbinom(nCells, mu = effect*x[1]/depths, size = x[2])}))
		exp.same <- t(apply(params[c((degs+1):nGenes),], 1, function(x) {
			rnbinom(nCells, mu = x[1]/depths, size = x[2])}))
		exp.counts[[i]] <- rbind(exp.de, exp.same)
		exp.bcs <- c(exp.bcs, paste("exp", c(1:nCells), i, sep = "_"))

	}

	ctrl.counts <- do.call(cbind, ctrl.counts)
	exp.counts <- do.call(cbind, exp.counts)

	return(list(ctrl_counts = ctrl.counts, 
		ctrl_bcs = ctrl.bcs,
		exp_counts = exp.counts,
		exp_bcs = exp.bcs))
}

opt <- parse_args(OptionParser(option_list = option_list))

# cat(sprintf('starting power analysis for %s with %d iterations\n', opt$sample, opt$nSims))
cat(sprintf('starting power analysis with %d iterations\n', opt$nSims))

# get sequencing depths from real data
counts <- readMM(opt$count_mtx)
seq.depths <- colSums(counts)

est.params <- read.table(opt$params, header = T)
# filter for detected genes
est.params <- est.params[est.params$mu!=0,]
cat(sprintf("sampling parameters from %d detected genes\n", dim(est.params)[1]))
# compute number of DEGs expected 
nDiff <- floor(opt$pctDiff * opt$nGenes)
# print simulation parameters
cat(sprintf("Simulation parameters:\n%d genes\n%d DEGs expected\n%d cells\nalpha = %.2f\neffect size = %.2f\n",
	opt$nGenes, nDiff, opt$nCells, opt$alpha, opt$effectSize))
cat(sprintf("testing for DEGs with %s test\n", opt$test))

results <- vector("numeric", length = opt$nSims)

for (i in c(1:opt$nSims)) {
	cat(sprintf("sim n=%d\n", i))
	# sample parameters for simulated data
	sim.params <- data.frame(mu = sample(est.params$mu, opt$nGenes, replace = TRUE), 
		size = sample(est.params$size, opt$nGenes, replace = TRUE))
	# sample sequencing depths 
	sim.depths <- sample(seq.depths, opt$nCells, replace = TRUE)

	cat("simulating counts\n")
	# simulate data from control case and generate cell barcodes
	sim.out <- simulateCounts(opt$replicates,
		sim.params,
		sim.depths,
		opt$nCells,
		nDiff,
		opt$effectSize,
		opt$nGenes)

	gene.names <- paste0("G", c(1:opt$nGenes))

	# load simulated data into Seurat and test for DE
	dimnames(sim.out$ctrl_counts) <- list(gene.names, sim.out$ctrl_bcs)
	dimnames(sim.out$exp_counts) <- list(gene.names, sim.out$exp_bcs)

	ctrl.seurat <- CreateSeuratObject(sim.out$ctrl_counts, 
		project = "ctrl", min.cells = 1, min.features = 1)
	exp.seurat <- CreateSeuratObject(sim.out$exp_counts, 
		project = "exp", min.cells = 1, min.features = 1)

	# merge 
	combined <- merge(ctrl.seurat, y = exp.seurat, 
		add.cell.ids = c("ctrl", "exp"), project = "power_analysis")

	cat('Normalizing merged data\n')
	combined <- NormalizeData(combined)

	cat('computing DEGs\n')
	# compute DEGs
	degs <- FindMarkers(combined, ident.1 = "ctrl", ident.2 = "exp", test.use = opt$test,
		logfc.threshold = 0, min.pct = 0, min.cells.feature = 0, min.cells.group=0)

	cat('computing power\n')
	# compute power
	tp <- sum(rownames(degs[degs$p_val_adj<opt$alpha,]) %in% gene.names[c(1:nDiff)])
	tpr <- tp/nDiff

	results[i] <- tpr
}

# output <- data.frame(replicates = opt$replicates, effect = opt$effectSize, power = results)
output <- data.frame(replicates = opt$replicates, cells = opt$nCells, effect = opt$effectSize, power = results)
write.table(output, file = opt$out, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = ',')
