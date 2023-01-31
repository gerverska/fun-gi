# Determine the viability of hamPCR ####

# Load packages ####
library(ggplot2)
library(nlme)

# Load functions ####
source(file.path('code', '00-functions.r'))

# Create output directories ####
out <- '06-analyze'
logs <- file.path(out, 'logs')
unlink(out, recursive = T)
dir.create(logs, recursive = T)
system(paste('touch', file.path(out, 'README.md')))

# Load and process rarefied output ####
fun.gi <- readRDS(file.path('05-rarefy', 'fun-gi.rds'))

meta <- fun.gi$meta
rownames(meta) <- meta$combo
meta$fun_n <- meta$fun_n |> gsub("_fun_fwd|_fun_rev", '', x = _) |> gsub('n', 'N', x = _)
meta$gi_n <- meta$gi_n |> gsub("_gi_fwd|_gi_rev", '', x = _) |> gsub('n', 'N', x = _)

standard <- subset(meta, template == 'standard')
dfssmt <- subset(meta, template == 'dfssmt')
dfssmt$noga_load <- dfssmt$fun_OTU.1 / dfssmt$gi_OTU.1

# Design a full model, including frameshifts as random effects ####
full <- lme(log10(mean_noga_load) ~ log10(dilution),
            random = list(fun_n = ~ 1, gi_n = ~ 1),
            data = standard,
            method = 'ML',
            na.action = na.omit)

# Check for normally distributed residuals ####
file.path(logs, 'full-residuals.png') |> png()
full |> residuals() |> qqnorm()
full |> residuals() |> qqline()
dev.off()

# Check for normally distributed random intercepts ####
file.path(logs, 'fun-full-rand.png') |> png()
ranef(full)$fun_n$`(Intercept)` |> qqnorm()
ranef(full)$fun_n$`(Intercept)` |> qqline()
dev.off()
file.path(logs, 'gi-full-rand.png') |> png()
ranef(full)$gi_n$`(Intercept)` |> qqnorm()
ranef(full)$gi_n$`(Intercept)` |> qqline()
dev.off()

# Design a base model, only including dilution as an effect ####
base <- lm(log10(mean_noga_load) ~ log10(dilution), standard)

# Check for normally distributed residuals ####
file.path(logs, 'base-residuals.png') |> png()
base |> residuals() |> qqnorm()
base |> residuals() |> qqline()
dev.off()

# Test whether the full and base models are similar to each other ####
lrtp <- anova(full, base)$`p-value`[2] |> round(3)

# Extract model summary ####
summ <- base |> summary()
ftp <- summ$coefficients[8] |> round(3)
ftp <- ifelse(ftp < 0.001, '< 0.001', ftp)
r2 <- summ$r.squared |> round(digits = 3)
note <- paste0('Frameshift vs no-frameshift\nLikelihood ratio test: P-value = ', lrtp,
               '\n\nBase model F-test P-value = ', ftp,
               '\nBase model R-squared = ', r2) |> gsub('= <', '<', x = _)

# Add model residuals to the "standard" dataframe ####
res <- data.frame(res = base$residuals)
res$combo <- rownames(res)
standard <- merge(standard, res, by = 'combo', all.x = T)

# Calculate the minimum, mean, and maximum NOGA load obtained from the DFSSMT samples ####
dfssmt.stats <- data.frame(stat = c('minimum', 'mean', 'maximum'),
                           value = c(min(log10(dfssmt$noga_load), na.rm = T),
                                     mean(log10(dfssmt$noga_load), na.rm = T),
                                     max(log10(dfssmt$noga_load), na.rm = T)
                                     ))
                                      
# Plot the relationship between load and the dilution ratio ####
plot.y <- max(log10(standard$mean_noga_load), na.rm = T) * 0.75
plot.x <- min(log10(standard$dilution), na.rm = T) * 0.8

load <- ggplot(standard, aes(x = log10(dilution), y = log10(mean_noga_load)
                           )) +
    geom_point(aes(color = fun_n, shape = gi_n), size = 3) +
    stat_smooth(method = 'lm') +
    geom_hline(aes(yintercept = value, linetype = stat), dfssmt.stats) +
    scale_x_continuous(n.breaks = 9) +
    scale_y_continuous(n.breaks = 6) +
    xlab('\nLog-transformed NOGA dilution ( prop. NOGA gDNA )') +
    ylab('Log-transformed NOGA load\n') +
    labs(color = 'Fun frameshift',
         shape = 'Gi frameshift',
         linetype = 'DFSSMT\nNOGA load') +
    annotate('text', x = plot.x, y = plot.y, label = note, size = 5) +
    guides(color = guide_legend(order = 1), 
           shape = guide_legend(order = 2),
           linetype = guide_legend(order = 3)) +
    scale_color_manual(values = color.pal) +
    theme_classic() +
    theme(text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, vjust = 0.5),
          axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'))

file.path(out, 'dilution-load.png') |> ggsave(load, width = 12, height = 9)

# Plot the relationship between frameshift pairs and model residuals ####
set.seed(666)
shifts <- ggplot(standard, aes(x = fun_n, y = gi_n, color = res)) +
    geom_jitter(size = 3, width = 0.15, height = 0.1, alpha = 0.75) +
    xlab('\nFun frameshift') +
    ylab('Gi frameshift\n') +
    labs(color = 'Residuals') +
    scale_color_viridis_c() +
    theme_classic() +
    theme(text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, vjust = 0.5),
          axis.title.x = element_text(face = 'bold'),
          axis.title.y = element_text(face = 'bold'),
          legend.title = element_text(face = 'bold'))

file.path(out, 'residuals.png') |> ggsave(shifts, width = 9, height = 6)

# Comparing relative abundance and load metrics ####
fun.tab <- standard[, colnames(standard)[grepl('OTU', colnames(standard)) == T]]
fun.tab <- fun.tab[, colSums(fun.tab) > 0]

tax <- fun.gi$tax[rownames(fun.gi$tax) %in% colnames(fun.tab), ] |> t() |> data.frame()

fun.tab <- standard |> subset(select = colnames(standard)[grepl('OTU', colnames(standard)) == T])
fun.tab <- standard |> subset(select = colnames(standard)[grepl('OTU', colnames(standard)) == T])
