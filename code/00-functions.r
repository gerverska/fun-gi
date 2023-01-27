# Colorblind palette ####
color.pal <- c('#000000', '#E69F00', '#56B4E9',
               '#009E73', '#F0E442', '#0072B2',
               '#D55E00', '#CC79A7')

# 04-compile.r ####
lulu.clust <- function(tab, seq, multi = threads, name, min.match = 0.97){
    
    # Collapse identical ASVs ####
    colnames(tab) <- seq |> as.character() |> unname()
    unique.tab <- collapseNoMismatch(tab)
    
    # Remove chimeras ####
    nochim.tab <- removeBimeraDenovo(unique.tab,
                                     method = 'consensus',
                                     multithread = multi,
                                     verbose = T)
    
    # Reassign OTU names ####
    nochim.otus <- paste0(name, '_OTU.', 1:ncol(nochim.tab))
    nochim.seq <- getSequences(nochim.tab) |> DNAStringSet()
    names(nochim.seq) <- nochim.otus
    colnames(nochim.tab) <- nochim.otus
    
    # Prepare data for LULU clustering ####
    t.nochim.tab <- nochim.tab |> t() |> as.data.frame()
    
    # Make a matchlist by aligning all remaining sequences against each other
    writeXStringSet(nochim.seq,
                    file = file.path('scratch', paste0('vsearch-', name, '.fa'))
    )
    
    vsearch.flags <- paste('--usearch_global', file.path('scratch', paste0('vsearch-', name, '.fa')),
                           '--db', file.path('scratch', paste0('vsearch-', name, '.fa')),
                           '--self',
                           '--id', (min.match - 0.01),
                           '--strand plus',
                           '--iddef 1',
                           '--threads', multi,
                           '--userout', file.path('scratch', paste0('matchlist-', name, '.txt')),
                           '--log', file.path(logs, paste0('vsearch-', name, '.txt')),
                           '--userfields query+target+id',
                           '--maxaccepts 0',
                           '--query_cov 0.9',
                           '--maxhits 10')
    system2('vsearch', args = vsearch.flags)
    
    matchlist <- read.table(file.path('scratch', paste0('matchlist-', name, '.txt')),
                            header = F,
                            as.is = T,
                            stringsAsFactors = F)
    
    # Use LULU to make a curation object and update the OTU table and sequences ####
    curation <- lulu(t.nochim.tab, matchlist, minimum_match = (min.match * 100))
    lulu.tab <- curation$curated_table |> t() |> as.data.frame()
    lulu.seq <- nochim.seq[names(nochim.seq) %in% curation$curated_otus]
    
    # Move the LULU log from its original location to the log directory ####
    system(paste('mv lulu.log*', file.path(logs, paste0('lulu-', name, '.txt'))
    ))
    
    # Output the new objects ####
    return(list(tab = lulu.tab,
                seq = lulu.seq))
}

# 05-rarefy.r ####
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