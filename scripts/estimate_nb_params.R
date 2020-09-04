library(optparse)
library(Matrix)

optimize_NB_pars <- function(par, input.counts, scaling) {
    
    est.mean <- par[1]^2
    est.disp <- max(par[2]^2, 1)
    
    out.ll <- dnbinom(x = input.counts, mu = est.mean*scaling, size = est.disp, log = T)
    
    -sum(out.ll)
}

# load arguments
option_list <- list(
    make_option("--count_mtx", action = "store", 
        default = '/iblm/netapp/home/jezhou/power_analysis/data/snRNA/Rat_Amygdala_787A_all_seq/outs/filtered_feature_bc_matrix/matrix.mtx.gz',
        type = "character", 
        help = "path to cellranger filtered_feature_bc_matrix/matrix.mtx.gz file"),
    make_option("--out", action = "store", 
        type = "character",
        help = "output filename")
)

opt <- parse_args(OptionParser(option_list = option_list))

# cat(sprintf('estimating NB parameters for %s\n', opt$sample))

counts <- readMM(opt$count_mtx)

scaling.factor <- colSums(counts)/sum(colSums(counts))

res <- apply(counts, 1, function(x) {
    mean.counts <- mean(x)
    res <- optim(c(mean.counts, 5), optimize_NB_pars, input.counts = x, scaling = scaling.factor)
    
    params <- res$par^2
    return(params)
})

cat(sprintf('writing params to %s\n', opt$out))

write.table(t(res), file = opt$out, col.names = c('mu','size'), row.names = F, quote = F)
