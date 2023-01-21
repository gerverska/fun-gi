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

# Examine whether sequencing depth could be correlated with hamPCR metrics ####
ggplot(standard, aes(x = log10(dilution), y = reads)) +
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

# Plot the relationship between the dilution ratio ####
# Fungal relative abundance ####
ggplot(standard, aes(x = dilution, y = mean_noga_ra)) +
    geom_point(aes(color = as.factor(dilution)), size = 3) +
    # geom_point(aes(x = ra_est, y = ra), dfssmt, size = 3, shape = 0, fill = 'white') +
    geom_smooth() +
    geom_hline(yintercept = c(min(dfssmt$mean_noga_ra, na.rm = T), max(dfssmt$mean_noga_ra, na.rm = T)), 
               color = 'red', linetype = 'dashed') +
    ylab('NOGA relative abundance\n') +
    xlab('\n Dilution ( prop. fungal gDNA )') +
    labs(color = 'Dilution') +
    scale_x_continuous(n.breaks = 10) +
    theme_cowplot() +
    scale_color_colorblind() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'))

# Fungal load (fungal : plant) ####
ggplot(standard, aes(x = dilution, y = mean_noga_load)) +
    geom_point(aes(color = as.factor(dilution)), size = 3) +
    geom_smooth() +
    geom_hline(yintercept = c(min(dfssmt$mean_noga_load, na.rm = T), max(dfssmt$mean_noga_load, na.rm = T)), 
               color = 'red', linetype = 'dashed') +
    ylab('NOGA load ( NOGA : Gigantea )\n') +
    xlab('\n Dilution ( prop. fungal gDNA )') +
    labs(color = 'Dilution') +
    scale_x_continuous(n.breaks = 10) +
    theme_cowplot() +
    scale_color_colorblind() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'))

# Root-transformed fungal load ####
ggplot(standard, aes(x = dilution, y = mean_noga_root)) +
    geom_point(aes(color = as.factor(dilution)), size = 3) +
    # geom_point(aes(x = root_est, y = root), dfssmt, size = 3, shape = 0, fill = 'white') +
    geom_smooth() +
    geom_hline(yintercept = c(min(dfssmt$mean_noga_root, na.rm = T), max(dfssmt$mean_noga_root, na.rm = T)), 
               color = 'red', linetype = 'dashed') +
    ylab('Root-transformed NOGA load\n') +
    xlab('\n Dilution ( prop. fungal gDNA )') +
    labs(color = 'Dilution') +
    scale_x_continuous(n.breaks = 10) +
    scale_y_continuous(n.breaks = 10) +
    theme_cowplot() +
    scale_color_colorblind() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'))

# Log-transformed fungal load ####
ggplot(standard, aes(x = dilution, y = mean_noga_log)) +
    geom_point(aes(color = as.factor(dilution)), size = 3) +
    # geom_point(aes(x = log_est, y = log), dfssmt, size = 3, shape = 0, fill = 'white') +
    geom_smooth() +
    geom_hline(yintercept = c(min(dfssmt$mean_noga_log, na.rm = T), max(dfssmt$mean_noga_log, na.rm = T)), 
               color = 'red', linetype = 'dashed') +
    ylab('Log-transformed NOGA load\n') +
    xlab('\n Dilution ( prop. fungal gDNA )') +
    labs(color = 'Dilution') +
    scale_x_continuous(n.breaks = 10) +
    scale_y_continuous(n.breaks = 6) +
    theme_cowplot() +
    scale_color_colorblind() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'))
