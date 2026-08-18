library(spacexr)
DS <- "B19-25653_8um"
library(Matrix)
library(methods)
library(reticulate)
setwd("/data/projects/spatialTX/spacexr_deconvolution/")
use_python("/usr/bin/python", required = TRUE)
py_config() 
scipy_sparse = import("scipy.sparse")
njobs <- 5


##### ------------------------------------REF PREP
refdir <- sprintf("./Deconvolution_%s", DS)  # or your path with the exported files

counts <- read.csv(file.path(refdir, "reference_dge.csv"), check.names = FALSE)
rownames(counts) <- counts$gene
counts$gene <- NULL

meta_data <- read.csv(file.path(refdir, "reference_meta_data.csv"), stringsAsFactors = FALSE)

# named vectors as spacexr expects
cell_types <- setNames(as.factor(meta_data$cluster), meta_data$barcode)
nUMI       <- setNames(meta_data$nUMI,            meta_data$barcode)

# align counts columns to meta barcodes (very important)
counts <- counts[, names(cell_types), drop = FALSE]

# Create Reference
reference <- Reference(counts, cell_types, nUMI)


##### ------------------------------------SPATIAL PREP

counts = scipy_sparse$load_npz(sprintf("./Deconvolution_%s/%s_spots_counts.npz",DS, DS ))
genes <- read.table(sprintf("./Deconvolution_%s/%s_spots_var.tsv",DS, DS ), header=TRUE, sep=",")$X0
BCs <- read.table(sprintf("./Deconvolution_%s/%s_spots_obs.tsv",DS, DS ), header=TRUE, sep=",")$X0
colnames(counts) <- genes
rownames(counts) <- BCs

counts <- t(counts)

# Coords
coords <- read.table(sprintf("./Deconvolution_%s/%s_spots_coords.tsv",DS, DS ), header=TRUE, sep=",")
rownames(coords) <- coords$X
coords <- coords[ -c(1) ]



##### ------------------------------------PUCKDATA PREP
puck <- SpatialRNA(coords, counts)

## Examine SpatialRNA object (optional)

print(dim(puck@counts)) # observe Digital Gene Expression matrix
hist(log(puck@nUMI,2)) # histogram of log_2 nUMI

print(head(puck@coords)) # start of coordinate data.frame
barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 

# This list can be restricted if you want to crop the puck e.g. 
# puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
# on the plot:
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                     title ='plot of nUMI') 


##### ------------------------------------RCTD PREP AND DECONVOLUTION



myRCTD <- create.RCTD(puck, reference, max_cores = njobs,UMI_min=10)

system.time({
  myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
})

resultsdir <- sprintf('RCTD_res%s', DS)
dir.create(resultsdir)
save(myRCTD, file = sprintf("./RCTD_res/myRCTD_fit%s.RData",DS ))



##### ------------------------------------Collect results
results <- myRCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights <- normalize_weights(results$weights) 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA

##### ------------------------------------Plots
plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, results$results_df) 
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)


#obtain a dataframe of only doublets
doublets <- results$results_df[results$results_df$spot_class == "doublet_certain",] 
# Plots all doublets in space (saved as 
# 'results/all_doublets.pdf')
plot_doublets(spatialRNA, doublets, resultsdir, cell_type_names) 


plot_doublets_type(spatialRNA, doublets, resultsdir, cell_type_names)
# a table of frequency of doublet pairs 
doub_occur <- table(doublets$second_type, doublets$first_type) 
# Plots a stacked bar plot of doublet ocurrences (saved as 
# 'results/doublet_stacked_bar.pdf')

plot_doub_occur_stack(doub_occur, resultsdir, cell_type_names) 
##### ------------------------------------Save res

write.csv(results$results_df, sprintf("./%s/results_%s.csv", resultsdir, DS), row.names = TRUE)
