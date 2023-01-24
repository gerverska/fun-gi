# Rarefy the data to reduce the influence of read depth on results ####

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

# Load functions ####
read.count <- function(x, otus, subsample) {
    # Count OTU reads ####
    data.frame(OTU = otus[x],
               reads = length(subsample[subsample %in% otus[x] == T]))
}
resample <- function(x, samp, depth){
    # Randomly sample reads for a sample ####
    subsample <- sample(samp$OTU, size = depth, replace = F)
    
    # Figure out which OTUs are present in the sample and how many there are ####
    otus <- unique(subsample)
    rich <- length(otus)
    
    # Obtain read counts for each OTU ####
    counts <- lapply(1:rich, read.count, otus = otus, subsample = subsample)
    
    # Combine count output into a single dataframe ####
    out <- do.call('rbind', counts)
    out$iteration <- x
    out
}
rarefy <- function(x, depth, subsamples, replicas){
    # Subset reads belonging to a given sample ####
    samp <- subset(replicas, combo == x)
    
    # Sub-sample these reads several times ####
    rare <- lapply(1:subsamples, resample, samp = samp, depth = depth)
    
    # Combine all the output ####
    out <- do.call('rbind', rare)
    out$combo <- x
    out
}

# Create output directories ####
out <- '05-rarefy'
unlink(out, recursive = T)
dir.create(out)
system(paste('touch', file.path(out, 'README.md')))

# Load and process inputs ####
fun.gi <- readRDS(file.path('04-compile', 'fun-gi.rds'))

fun.tab <- fun.gi$fun$tab
fun.tab$fun_id <- rownames(fun.tab)

gi.tab <- fun.gi$gi$tab
gi.tab$gi_id <- rownames(gi.tab)

# Combine fungal and Gigantea OTU tables and metadata into a single dataframe ####
full <- merge(fun.gi$fun$meta, fun.tab, by = 'fun_id', all = T)
full <- merge(full, gi.tab, by = 'gi_id', all = T)
full$combo <- paste(full$tags, full$fun_n, full$gi_n, sep = '-') # Ideally fix some of this upstream!!!
colnames(full) <- colnames(full) |> gsub('-', '', x = _) # Ideally fix this upstream!!!

# Convert NA values post merge ####
otus <- colnames(full)[grepl('OTU', colnames(full))]

full.tab <- full[, otus]
full.tab[is.na(full.tab) == T] <- 0
full[, otus] <- full.tab

# Remove samples beneath target sequencing depth prior to rarefaction ####
depth <- 750 # Need to show how we got at this number!!!

full$reads <- subset(full, select = otus) |> rowSums()
filt <- full |> subset(x = _, reads >= depth)
filt.tab <- filt[, otus]
filt.tab$combo <- filt$combo

# Convert the OTU table to long format ####
long <- reshape(filt.tab,
                direction = 'long', varying = otus,
                v.names = 'reads', times = otus,
                timevar = 'OTU', idvar = 'combo')

# Create a row for each read of each OTU in each sample ####
replicas <- lapply(long, rep, long$reads) |> as.data.frame() |>
    subset(x = _, select = -reads)

# Subsample each sample a number of times to a given depth ####
combo <- filt.tab$combo
subsamples <- 1000
rare <- lapply(combo, rarefy,
               depth = depth, subsamples = subsamples,
               replicas = replicas) |> do.call(rbind, args = _)

# Create a column for each OTU and remove NA values ####
wide <- reshape(rare,
                direction = 'wide', timevar = 'OTU',
                idvar = c('combo', 'iteration'),
                v.names = 'reads',
                sep = '_')
colnames(wide) <- colnames(wide) |> gsub('reads_', '', x = _)
wide[is.na(wide) == T] <- 0

# Subset the dataframe for OTU columns ####
wide.tab <- wide[, colnames(wide) %in% otus]

# Identify fungal OTUs ####
fun.otus <- colnames(wide.tab)[grepl('funOTU', colnames(wide.tab)) == T]

# Calculate statistics for each subsample ####
wide.tab$noga_ra <- wide.tab$funOTU.1 / depth
wide.tab$noga_load <- wide.tab$funOTU.1 / wide.tab$giOTU.1

wide.tab$fun_ra <- rowSums(wide.tab[, fun.otus]) / depth
wide.tab$fun_load <- rowSums(wide.tab[, fun.otus]) / wide.tab$giOTU.1

stat.tab <- wide.tab[!(colnames(wide.tab) %in% otus)]
stat.tab$combo <- wide$combo
stat.tab$iteration <- wide$iteration

stats <- colnames(stat.tab)[grepl("ra$|load$", colnames(stat.tab)) == T]

# Lengthen the dataframe to make summarization simpler ####
stat.long <- reshape(stat.tab,
                direction = 'long', varying = stats,
                v.names = 'value', times = stats,
                timevar = 'stat', idvar = c('combo', 'iteration'))

# Summarize each statistic for each sample using the subsamples ####
calc <- by(stat.long,
           INDICES = list(stat.long$combo, stat.long$stat),
           FUN = function(x){data.frame(combo = unique(x$combo),
                                        stat = unique(x$stat),
                                        mean = mean(x$value),
                                        sd = sd(x$value))
               }) |> do.call(rbind, args = _)

# Widen the dataframe to equalize the number of rows and number of samples ####
calc.wide <- reshape(calc,
                     direction = 'wide', timevar = 'stat',
                     idvar = 'combo',
                     sep = '_')

# Combine the rarefied output and the full metadata into a single dataframe ####
meta <- merge(full, calc.wide, by = 'combo', all = T)
fun.gi.rare <- list(meta = meta,
                    fun_seq = fun.gi$fun$seq,
                    tax = fun.gi$fun$tax,
                    boot = fun.gi$fun$boot,
                    gi_seq = fun.gi$gi$seq)

# Write the object out to an .rds file ####
file.path(out, 'fun-gi-rare.rds') |> saveRDS(fun.gi.rare, file = _)
