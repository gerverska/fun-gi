# This script returns the following: ####

# Load packages ####
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(patchwork)

# Load trimming summary and experiment metadata ####
meta <- read.csv(file.path('logs', 'meta.csv')) |>
    pivot_longer(ends_with('id'), values_to = 'id', names_to = 'primers') |>
    select(template, dilution, id, primers, combo)

trim <- read.delim(file.path('logs', '02-trim3-reads.txt'), sep = ' ') |>
		select(id = trim3, reads) |>
		separate(id, c(NA, NA, NA, 'fwd', 'rev'), '-', F) |> 
		mutate(pair = paste0(fwd, '-', rev))

# Identify obvious Fun + Gi chimeras ####
fwd.chimeras <- trim |> filter(grepl('fun', fwd) & grepl('gi', rev)) |> select(reads) |> sum()
rev.chimeras <- trim |> filter(grepl('gi', fwd) & grepl('fun', rev)) |> select(reads) |> sum()
chimeras <- fwd.chimeras + rev.chimeras

# Track shifts in frameshift codes across the entire sequencing run ####
# 5.8S-Fun + ITS4-Fun ####
fun.shifts <- trim |>
		filter(grepl('fun', fwd) & grepl('fun', rev)) |>
		group_by(pair) |> 
		summarize(total = sum(reads)) |>
		arrange(pair)

fun.shifts$pair <- str_remove_all(fun.shifts$pair, "_fun_fwd|_fun_rev") |>
	str_replace_all('n', 'N')

# Calculate the proportion of correct frameshift demultiplexing ####
fun.total <- fun.shifts$total |> sum() + chimeras

fun.pairs <- c('N03-N03', 'N04-N04',
							 'N05-N05', 'N06-N06')
fun.correct <- fun.shifts |>
		filter(pair %in% fun.pairs) |>
		select(total) |>
		sum()

# Calculate rates at which single or double indels occurs in a read ####
fun.del.stop <- which(fun.shifts$pair == 'N03-N03') - 1 

fun.del <- fun.shifts |> slice(1:fun.del.stop)

fun.del.sum <- fun.del |>
		select(total) |>
		sum()

fun.in.start <- which(fun.shifts$pair == 'N06-N06') + 1

fun.in <- fun.shifts |> slice(fun.in.start:n())

fun.in.sum <- fun.in |>
		select(total) |>
		sum()

fun.text <- paste0(fun.correct, ' / ', fun.total, '\nITS2-trimmed reads\ncorrectly demultiplexed', '\n\n',
									 fun.del.sum, ' reads with deletions in N03-N03', '\n',
									 fun.in.sum, ' reads with insertions in N06-N06', '\n\n',
									 'Red and blue dashed lines indicate pairs\nused to calculate deletion and insertion\noccurrence rates, respectively')

# Prepare values for plotting ####
fun.y <- fun.shifts$total |> max() * 0.9

# Overall summary ####
fun.shift.reads <- ggplot(fun.shifts, aes(x = pair, y = total)) +
		geom_bar(stat = 'identity') +
		geom_vline(xintercept = 'N03-N03', color = 'red', linetype = 'dashed') +
		geom_vline(xintercept = 'N06-N06', color = 'blue', linetype = 'dashed') +
		scale_y_continuous(n.breaks = 8) +
		annotate('text', x = 'N02-N06', y = fun.y, label = fun.text) +
		xlab('\nITS2 frameshift pairs') +
		ylab('Total reads\n') +
		theme_cowplot() +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
					axis.title.x = element_text(face = 'bold'),
					axis.title.y = element_text(face = 'bold'))

# GiF + GiR ####
gi.shifts <- trim |>
		filter(grepl('gi', fwd) & grepl('gi', rev)) |>
		group_by(pair) |> 
		summarize(total = sum(reads)) |>
		arrange(pair)

gi.shifts$pair <- str_remove_all(gi.shifts$pair, "_gi_fwd|_gi_rev") |>
	str_replace_all('n', 'N')

# Calculate the proportion of correct frameshift demultiplexing ####
gi.total <- gi.shifts$total |> sum() + chimeras

gi.pairs <- c('N06-N06', 'N07-N07',
							'N08-N08', 'N09-N09')
gi.correct <- gi.shifts |>
	filter(pair %in% gi.pairs) |>
	select(total) |>
	sum()

# Calculate rates at which at least one insertion or deletion occurs in a read ####
gi.del.stop <- which(gi.shifts$pair == 'N06-N06') - 1 

gi.del <- gi.shifts |> slice(1:gi.del.stop)

gi.del.sum <- gi.del |>
	select(total) |>
	sum()

gi.in.start <- which(gi.shifts$pair == 'N09-N09') + 1

gi.in <- gi.shifts |> slice(gi.in.start:n())

gi.in.sum <- gi.in |>
	select(total) |>
	sum()

gi.text <- paste0(gi.correct, ' / ', gi.total, '\nGi-trimmed reads\ncorrectly demultiplexed', '\n\n',
									 gi.del.sum, ' reads with deletions in N06-N06', '\n',
									 gi.in.sum, ' reads with insertions in N09-N09')

# Calculate global chimera rate ####
chimera.text <- paste0('\n', chimeras, ' / ', (fun.total + gi.total - chimeras), ' reads\nas fungal-plant chimeras')
chimera.gi.text <- paste0(gi.text, '\n', chimera.text)

# Prepare values for plotting #### 
gi.y <- gi.shifts$total |> max() * 0.9

# Overall summary ####
gi.shift.reads <- ggplot(gi.shifts, aes(x = pair, y = total)) +
		geom_bar(stat = 'identity') +
		geom_vline(xintercept = 'N06-N06', color = 'red', linetype = 'dashed') +
		geom_vline(xintercept = 'N09-N09', color = 'blue', linetype = 'dashed') +
		annotate('text', x = 'N05-N07', y = gi.y, label = chimera.gi.text) +
		scale_y_continuous(n.breaks = 10) +
		xlab('\nGigantea frameshift pairs') +
		ylab('') +
		theme_cowplot() +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
					axis.title.x = element_text(face = 'bold'),
					axis.title.y = element_text(face = 'bold'))

# Combine both plots ####
shift.reads <- fun.shift.reads + gi.shift.reads

# Join the two datasets to each other ####
join <- left_join(meta, trim, by = 'id', keep = T) |>
    select(combo, template, primers, dilution, reads) |>
    pivot_wider(names_from = primers, values_from = reads) |> 
		mutate(reads = fun_id + gi_id,
					 ra = fun_id / reads,
					 load = fun_id / gi_id,
					 root = load^(1/4),
					 log = log(load)) |>
		separate(combo, c('fun_shift', 'gi_shift'), "[.]", F) |>
		separate(fun_shift, c(NA, NA, NA, 'fun_fwd', 'fun_rev'), '-', F) |> 
		separate(gi_shift, c(NA, NA, NA, 'gi_fwd', 'gi_rev'), '-', F) |> 
		mutate(fun_shift = paste0(fun_fwd, '-', fun_rev)) |> 
		mutate(gi_shift = paste0(gi_fwd, '-', gi_rev))

join$fun_shift <- str_remove_all(join$fun_shift, "_fun_fwd|_fun_rev") |>
	str_replace_all('n', 'N')
join$gi_shift <- str_remove_all(join$gi_shift, "_gi_fwd|_gi_rev") |>
	str_replace_all('n', 'N')

# Subset the dataset ####
standard <- filter(join, template == 'standard')
dfssmt <- filter(join, template == 'dfssmt')

# Examine whether sequencing depth could be correlated with hamPCR metrics ####
dilution.reads <- ggplot(standard, aes(x = log10(dilution), y = reads)) +
		geom_point(aes(color = as.factor(dilution*100)), size = 3) +
		scale_x_continuous(n.breaks = 10) +
		scale_y_continuous(n.breaks = 8) +
		xlab('\nLog-transformed fungal dilution') +
		ylab('Total reads\n') +
		labs(color = 'Dilution\n( prop. fungal gDNA )') +
		scale_color_colorblind() +
		theme_cowplot() +
		theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
					axis.title.x = element_text(face = 'bold'),
					axis.title.y = element_text(face = 'bold'))
	
reads.root <- ggplot(standard, aes(x = log10(reads), y = root)) +
	geom_point(aes(color = log2(dilution*100)), size = 3) +
	scale_x_continuous(n.breaks = 10) +
	scale_y_continuous(n.breaks = 8) +
	xlab('\nLog-transformed read depth') +
	ylab('Root-transformed fungal load\n') +
	labs(color = 'Log-transformed dilution\n( prop. fungal gDNA )\n') +
	scale_color_viridis_c() +
	theme_cowplot() +
	theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
				axis.title.x = element_text(face = 'bold'),
				axis.title.y = element_text(face = 'bold'))

# Combine both plots ####
dilution.reads.roots <- dilution.reads / reads.root

# Fit LOESS smooths to the data for ad hoc estimation of DFSSMT fungal:plant dilutions ####
# NOTE: The following plots are likely premature given that the previous plots
# show that sequencing depth may be correlated with both dilution treatment and hamPCR metrics!!!
ra.fit <- loess(dilution ~ ra, standard)
load.fit <- loess(dilution ~ load, standard)
root.fit <- loess(dilution ~ root, standard)
log.fit <- loess(dilution ~ log, standard)

dfssmt$ra_est <- predict(ra.fit, dfssmt$ra)
dfssmt$load_est <- predict(load.fit, dfssmt$load)
dfssmt$root_est <- predict(root.fit, dfssmt$root)
dfssmt$log_est <- predict(log.fit, dfssmt$log)

estimates <- dfssmt |> select(ends_with('_est'))
min.est <- estimates |> min() |> round(5)
mean.est <- estimates |> colMeans()
overall.mean <- mean.est |> mean() |> round(5)
max.est <- estimates |> max() |> round(5)

rel.text <- paste0('Minimum dilution esimate: ', min.est, '\n',
									 'Mean dilution estimate: ', overall.mean, '\n',
									 'Maximum dilution estimate: ', max.est, '\n',
									 'Squares represent DFSSMT samples')

# Plot the relationship between the dilution ratio ####
# Fungal relative abundance ####
ra <- ggplot(standard, aes(x = dilution, y = ra)) +
		geom_point(aes(color = as.factor(dilution)), size = 3) +
		geom_point(aes(x = ra_est, y = ra), dfssmt, size = 3, shape = 0, fill = 'white') +
		geom_smooth() +
		geom_vline(xintercept = mean.est[[1]], color = 'red', linetype = 'dashed') +
    ylab('Fungal relative abundance\n') +
    xlab('\n') +
    labs(color = 'Dilution ( prop. fungal gDNA )') +
		scale_x_continuous(n.breaks = 10) +
    theme_cowplot() +
    scale_color_colorblind() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
    			axis.title.x = element_text(face = 'bold'),
      		axis.title.y = element_text(face = 'bold'),
      		legend.title = element_text(face = 'bold'))

# Fungal load (fungal : plant) ####
load <- ggplot(standard, aes(x = dilution, y = load)) +
		geom_point(aes(color = as.factor(dilution)), size = 3) +
		geom_point(aes(x = load_est, y = load), dfssmt, size = 3, shape = 0, fill = 'white') +
		geom_smooth() +
		geom_vline(xintercept = mean.est[[2]], color = 'red', linetype = 'dashed') +
    ylab('Fungal load ( fungal : plant )\n') +
    xlab('\n') +
		labs(color = 'Dilution ( prop. fungal gDNA )') +
		scale_x_continuous(n.breaks = 10) +
    theme_cowplot() +
    scale_color_colorblind() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
    			axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'))

# Root-transformed fungal load ####
root <- ggplot(standard, aes(x = dilution, y = root)) +
		geom_point(aes(color = as.factor(dilution)), size = 3) +
		geom_point(aes(x = root_est, y = root), dfssmt, size = 3, shape = 0, fill = 'white') +
		geom_smooth() +
		geom_vline(xintercept = mean.est[[3]], color = 'red', linetype = 'dashed') +
    ylab('Root-transformed fungal load\n') +
    xlab('\nDilution') +
		labs(color = 'Dilution ( prop. fungal gDNA )') +
		scale_x_continuous(n.breaks = 10) +
		scale_y_continuous(n.breaks = 6) +
    theme_cowplot() +
    scale_color_colorblind() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
    			axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'))

# Log-transformed fungal load ####
log <- ggplot(standard, aes(x = dilution, y = log)) +
    geom_point(aes(color = as.factor(dilution)), size = 3) +
		geom_point(aes(x = log_est, y = log), dfssmt, size = 3, shape = 0, fill = 'white') +
		geom_smooth() +
		geom_vline(xintercept = mean.est[[4]], color = 'red', linetype = 'dashed') +
		annotate('text', x = 0.00625, y = 0, label = rel.text) +
    ylab('Log-transformed fungal load\n') +
    xlab('\nDilution') +
		labs(color = 'Dilution ( prop. fungal gDNA )') +
		scale_x_continuous(n.breaks = 10) +
		scale_y_continuous(n.breaks = 6) +
    theme_cowplot() +
    scale_color_colorblind() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
    			axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'))

# Combine all plots ####
rel <- (ra + load) / (root + log) + plot_layout(guides = 'collect')
