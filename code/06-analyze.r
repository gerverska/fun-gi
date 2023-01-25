# Load packages ####
library(cowplot)
library(ggplot2)
library(nlme)

# Set palette ####
color.pal <- c('#000000', '#E69F00', '#56B4E9',
               '#009E73', '#F0E442', '#0072B2',
               '#D55E00', '#CC79A7')

# Create output directories ####
out <- '06-analyze'
unlink(out, recursive = T)
dir.create(out)
system(paste('touch', file.path(out, 'README.md')))

# Load and process inputs ####
fun.gi <- readRDS(file.path('05-rarefy', 'fun-gi-rare.rds'))

meta <- fun.gi$meta
rownames(meta) <- meta$combo
meta$fun_n <- meta$fun_n |> gsub("_fun_fwd|_fun_rev", '', x = _) |> gsub('n', 'N', x = _)
meta$gi_n <- meta$gi_n |> gsub("_gi_fwd|_gi_rev", '', x = _) |> gsub('n', 'N', x = _)

standard <- subset(meta, template == 'standard')
dfssmt <- subset(meta, template == 'dfssmt')

# Examine whether sequencing depth, dilution treatments, and/or frameshift codes are correlated ####
depth <- ggplot(standard, aes(x = log10(dilution), y = log10(reads)
                           )) +
    geom_point(aes(color = fun_n, shape = gi_n), size = 3) +
    scale_x_continuous(n.breaks = 9) +
    scale_y_continuous(n.breaks = 8) +
    xlab('\nLog-transformed NOGA dilution ( prop. NOGA gDNA )') +
    ylab('Total reads\n') +
    scale_color_manual(values = color.pal) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'),
          legend.position = 'none')

# Design a full model, including frameshifts as random effects ####
full <- lme(log10(mean_noga_load) ~ log10(dilution),
            random = list(fun_n = ~ 1, gi_n = ~ 1),
            data = standard,
            method = 'ML',
            na.action = na.omit)

# Check for normally distributed residuals ####
full |> residuals() |> qqnorm()
full |> residuals() |> qqline()

# Check for normally distributed random intercepts ####
ranef(full)$fun_n$`(Intercept)` |> qqnorm()
ranef(full)$fun_n$`(Intercept)` |> qqline()
ranef(full)$gi_n$`(Intercept)` |> qqnorm()
ranef(full)$gi_n$`(Intercept)` |> qqline()

# Design a base model, only including dilution as an effect ####
base <- lm(log10(mean_noga_load) ~ log10(dilution), standard)

# Check for normally distributed residuals ####
base |> residuals() |> qqnorm()
base |> residuals() |> qqline()

# Test whether the full and base models are similar to each other ####
test <- anova(full, base) |> data.frame()

# Add model residuals to the standard dataframe ####
res <- data.frame(res = base$residuals)
res$combo <- rownames(res)
standard <- merge(standard, res, by = 'combo', all.x = T)

# Plot the relationship between load and the dilution ratio ####
load <- ggplot(standard, aes(x = log10(dilution), y = log10(mean_noga_load)
                           )) +
    geom_point(aes(color = fun_n, shape = gi_n), size = 3) +
    stat_smooth(method = 'lm') +
    geom_hline(yintercept = min(log10(dfssmt$mean_noga_load), na.rm = T), linetype = 'dotted') +
    geom_hline(yintercept = mean(log10(dfssmt$mean_noga_load), na.rm = T), linetype = 'dashed') +
    geom_hline(yintercept = max(log10(dfssmt$mean_noga_load), na.rm = T), linetype = 'solid') +
    scale_x_continuous(n.breaks = 9) +
    scale_y_continuous(n.breaks = 6) +
    xlab('\nLog-transformed NOGA dilution ( prop. NOGA gDNA )') +
    ylab('Log-transformed NOGA load\n') +
    labs(color = 'Fun frameshift', shape = 'Gi frameshift') +
    guides(color = guide_legend(order = 1), 
           shape = guide_legend(order = 2)) +
    theme_cowplot() +
    scale_color_manual(values = color.pal) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'))

# Plot the relationship between frameshift pairs and model residuals ####
shifts <- ggplot(standard, aes(x = fun_n, y = gi_n, color = res)) +
    geom_jitter(size = 3, width = 0.15, height = 0.1, alpha = 0.75) +
    xlab('\nFun frameshift') +
    ylab('Gi frameshift\n') +
    labs(color = 'Residuals') +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'))
