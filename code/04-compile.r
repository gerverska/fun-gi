# Apply LULU post-clustering curation and assign fungal taxonomy ####

# Accept an argument for the number of threads and check it for validity ####
threads <- commandArgs(T) |> as.integer()

if(is.na(threads) == T){
    stop('Argument was converted to NA')
}
if(length(threads) < 1){
    stop('Please specify the number of threads to launch')
}
if(length(threads) > 1){
    stop('Too many arguments have been provided')
}
if(is.numeric(threads) == F){
    stop('Only numeric arguments are accepted')
}
if(threads < 1){
    stop('At least one thread is needed')
} else {
    cat(threads, 'threads requested', '\n')
}

# Install the lulu github package ####
remotes::install_github('gerverska/lulu@v1.0.0', dependencies = F)

# Load packages ####
library(Biostrings)
library(dada2)
library(lulu) # Need to see how to download this through conda!!!

# Load functions ####
source(file.path('code', '00-functions.r'))

# Create output directories ####
out <- '04-compile'
unlink(out, recursive = T)
dir.create(out)
system(paste('touch', file.path(out, 'README.md')))
unlink('scratch', recursive = T)
dir.create('scratch')

# Load inputs ####
meta <- read.csv(file.path('data', 'meta.csv')) # Change csv to rds!!!
seq.tab <- readRDS(file.path('03-denoise', 'seq-tab.rds'))

# Define fungal and Gigantea samples ####
fun.samp <- meta$fun_id |> unique()
gi.samp <- meta$gi_id |> unique()

# Extract the fungal dataset from the DADA2 output ####
fun.tab <- seq.tab[rownames(seq.tab) %in% fun.samp, ]
fun.tab <- fun.tab[, colSums(fun.tab) > 0]

# Rename OTUs and sequences ####
fun.otus <- paste0('OTU.', 1:ncol(fun.tab))
fun.seq <- getSequences(fun.tab) |> DNAStringSet()
names(fun.seq) <- fun.otus
colnames(fun.tab) <- fun.otus

# Identify ASVs with complete fungal sequences ####
writeXStringSet(fun.seq, file = file.path('scratch', 'itsx.fa'))

itsx.flags <- paste('-i', file.path('scratch', 'itsx.fa'),
                    '-t "fungi,tracheophyta"',
                    '--preserve T',
                    '--cpu', threads, 
                    '-o', file.path('logs', '04-compile-itsx'),
                    '--only_full T')
system2('ITSx', args = itsx.flags)

# Remove ASVs with incomplete fungal sequences ####
itsx.otus <- readDNAStringSet(file.path('logs', '04-compile-itsx.ITS2.fasta')) |> names()
itsx.seq <- fun.seq[names(fun.seq) %in% itsx.otus]
itsx.tab <- fun.tab[, names(itsx.seq)]

# Perform LULU post-clustering on the fungal subset ####
fun.lulu <- lulu.clust(tab = itsx.tab, seq = itsx.seq, multi = threads, name = 'fun', min.match = 0.97)

# Predict fungal OTU taxonomy ####
set.seed(666)
taxa <- assignTaxonomy(fun.lulu$seq,
                       file.path('data', 'sh_general_release_dynamic_all_psme-noga_29.11.2022.fasta.gz'), # Potentially update with _s_ version!!!
                       minBoot = 80,
                       tryRC = T,
                       multithread = threads,
                       outputBootstraps = T)

# Set the sequence object as a dataframe ####
fun.lulu$seq <- fun.lulu$seq |> as.data.frame()
colnames(fun.lulu$seq) <- 'sequence'

# Add taxonomic information and metadata to the fungal output ####
tax <- taxa$tax |> as.data.frame()
boot <- taxa$boot |> as.data.frame()
trim.tax <- lapply(tax, gsub, pattern="^[[:alpha:]]__", replacement='') |> data.frame()

rownames(trim.tax) <- rownames(fun.lulu$seq)
rownames(boot) <- rownames(fun.lulu$seq)

fun.lulu$tax <- trim.tax
fun.lulu$boot <- boot
fun.lulu$meta <- meta

# Perform LULU post-clustering on the Gigantea subset ####
gi.tab <- seq.tab[rownames(seq.tab) %in% gi.samp, ]
gi.tab <- gi.tab[, colSums(gi.tab) > 0]

gi.otus <- paste0('OTU.', 1:ncol(gi.tab))
gi.seq <- getSequences(gi.tab) |> DNAStringSet()
names(gi.seq) <- gi.otus
colnames(gi.tab) <- gi.otus

gi.lulu <- lulu.clust(tab = gi.tab, seq = gi.seq, multi = threads, name = 'gi', min.match = 0.85)
gi.lulu$seq <- gi.lulu$seq |> as.data.frame()
colnames(gi.lulu$seq) <- 'sequence'

gi.lulu$meta <- meta

# Combine the fungal and Gigantea objects before output ####
fun.gi <- list(fun = fun.lulu,
               gi = gi.lulu)

file.path(out, 'fun-gi.rds') |> saveRDS(fun.gi, file = _)

# Remove temporary files ####
unlink('scratch', recursive = T)
