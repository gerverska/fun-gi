# Load packages ####
library(cowplot)
library(ggthemes)
library(ggplot2)

# Load functions ####

# Create output directories ####
out <- '06-analyze'
unlink(out, recursive = T)
dir.create(out)
system(paste('touch', file.path(out, 'README.md')))

# Load and process inputs ####
fun.gi <- readRDS(file.path('05-rarefy', 'fun-gi-rare.rds'))

meta <- fun.gi$meta
standard <- subset(meta, template == 'standard')
dfssmt <- subset(meta, template == 'dfssmt')

# Examine whether sequencing depth, dilution treatments, and/or frameshift codes are correlated ####
ggplot(standard, aes(x = log10(dilution), y = reads)) +
    geom_point(aes(color = fun_n, shape = gi_n), size = 3) +
    scale_y_continuous(n.breaks = 8, limits = c(0, 14000)) +
    scale_x_continuous(n.breaks = 12, limits = c(-4.25, -1.75)) +
    xlab('\nLog-transformed fungal dilution ( prop. fungal gDNA )') +
    ylab('Total reads\n') +
    scale_color_colorblind() +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'))

# Plot the relationship between the dilution ratio ####
# Root-transformed fungal load ####
ggplot(standard, aes(x = dilution, y = mean_noga_root)) +
    geom_point(aes(color = fun_n, shape = gi_n), size = 3) +
    stat_smooth(span = 1, n = 32) +
    # facet_wrap(~fun_n) +
    geom_hline(yintercept = min(dfssmt$mean_noga_root, na.rm = T), linetype = 'dotted') +
    geom_hline(yintercept = mean(dfssmt$mean_noga_root, na.rm = T), linetype = 'dashed') +
    geom_hline(yintercept = max(dfssmt$mean_noga_root, na.rm = T), linetype = 'solid') +
    ylab('Root-transformed NOGA load\n') +
    xlab('\n Dilution ( prop. fungal gDNA )') +
    scale_x_continuous(n.breaks = 10) +
    # ylim(0.5, 10) +
    scale_y_continuous(n.breaks = 10) +
    theme_cowplot() +
    scale_color_colorblind() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'))
