# Denoise sequences and identify amplicon sequence variants (ASVs) ####

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

# Load packages ####
library(ggplot2) # Should be imported by dada2, but I also thought this about BioStrings...
library(dada2)

# Create output directories ####
out <- '03-denoise'
logs <- file.path(out, 'logs')
unlink(out, recursive = T)
dir.create(logs, recursive = T)
system(paste('touch', file.path(out, 'README.md')))
unlink('scratch', recursive = T)
dir.create('scratch')

# Read in forward and reverse reads ####
in.fwd <- in.path |> list.files(pattern = '.R1.fq.gz', full.names = T) |> sort()
in.rev <- in.path |> list.files(pattern = '.R2.fq.gz', full.names = T) |> sort()

# Make paths for trimmed and filtered files ####
filt <- file.path(out, 'filter')
filt.fwd <- gsub(in.path, filt, in.fwd)
filt.rev <- gsub(in.path, filt, in.rev)

# Trim and quality filter ####
trim <- filterAndTrim(in.fwd, filt.fwd,
                      in.rev, filt.rev,
                      maxEE = c(2,2),
                      multithread = threads)

# # Update list of trimmed file paths to exclude samples with no reads passing filters ####
filt.fwd <- filt |> list.files(pattern = 'R1.fq.gz', full.names = T)
filt.rev <- filt |> list.files(pattern = 'R2.fq.gz', full.names = T)

# Make fastqc reports for each sample and read direction after filtering and trimming ####
fastqc.fwd <- paste(paste(filt.fwd, collapse = ' '),
                    '-t', threads,
                    '-o scratch')
fastqc.rev <- paste(paste(filt.rev, collapse = ' '),
                    '-t', threads,
                    '-o scratch')

system2('fastqc', args = fastqc.fwd)
system2('fastqc', args = fastqc.rev)

# Synthesize multiqc reports for each read direction ####
multiqc.fwd <- paste('-f -o', logs,
                     '-n filt-R1.html',
                     '-ip',
                     "-i 'Forward reads after quality filtering and trimming'",
                     'scratch/*R1_fastqc.zip')
multiqc.rev <- gsub('R1.html', 'R2.html', multiqc.fwd) |> 
		gsub('Forward', 'Reverse', x = _) |> 
		gsub('R1_fastqc.zip', 'R2_fastqc.zip', x = _)

system2('multiqc', args = multiqc.fwd)
system2('multiqc', args = multiqc.rev)

# Dereplicate sequences ####
derep.fwd <- derepFastq(filt.fwd, verbose = F)
derep.rev <- derepFastq(filt.rev, verbose = F)

# Trim names for each derep object ####
names(derep.fwd) <- gsub('-rc.R1.fq.gz', '', names(derep.fwd))
names(derep.rev) <- gsub('-rc.R2.fq.gz', '', names(derep.rev))

# Learn errors for each read direction. 1e09 bases will take a long time. ####
err.fwd <- learnErrors(filt.fwd, multithread = threads, randomize = T)  # What does the randomize flag do? Not in manual...
err.fwd.plot <- plotErrors(err.fwd, nominalQ = T) + ggtitle(paste('Forward read error model'))
file.path(logs, 'error-fwd.png') |> ggsave(err.fwd.plot, width = 12, height = 9)
err.rev <- learnErrors(filt.rev, multithread = threads, randomize = T)
err.rev.plot <- plotErrors(err.rev, nominalQ = T) + ggtitle(paste('Reverse read error model'))
file.path(logs, 'error-rev.png') |> ggsave(err.rev.plot, width = 12, height = 9)

# Denoise reads in both directions ####
dada.fwd <- dada(derep.fwd, err = err.fwd, multithread = threads, pool = 'pseudo')
dada.rev <- dada(derep.rev, err = err.rev, multithread = threads, pool = 'pseudo')

# Merge forward and reverse reads and make a sequence table ####
merged <- mergePairs(dada.fwd, derep.fwd,
                     dada.rev, derep.rev,
                     trimOverhang = T)
seq.tab <- merged |> makeSequenceTable()
file.path(out, 'seq-tab.rds') |> saveRDS(seq.tab, file = _)

# Make a summary log ####
get.n <- function(x){
	sum(getUniques(x))
}

trim.summary <- trim |> data.frame()
trim.summary$sample <- rownames(trim.summary)
trim.summary$sample <- gsub('-rc-R1.fq.gz', '', trim.summary$sample)

track <- cbind(sapply(dada.fwd, get.n),
               sapply(dada.rev, get.n),
               sapply(merged, get.n)) |>
		data.frame()
track$sample <- rownames(track)

log <- merge(trim.summary, track, by = 'sample', all.x = T)
colnames(log) <- c('sample', 'input', 'filtered', 'denoised.fwd', 'denoised.rev', 'merged')
file.path(logs, 'dada2.rds') |> saveRDS(log, file = _)

# Remove the scratch directory ####
unlink('scratch', recursive = T)
