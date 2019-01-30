################################################################################
#
#   File name: calculate_significance.R
#
#   Authors: Stefano Pirro' ( s.pirro@qmul.ac.uk ) (*aka wynstep*)
#
#   Short Description:
#     Starting from a PRS score and a range (0, maxPRS), produces a gaussian distribution of values
#     and assess the significance of the score
#
################################################################################

#===============================================================================
#    Load libraries
#===============================================================================
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))

#===============================================================================
#    Catching the arguments
#===============================================================================

option_list = list(
  make_option(c("-p", "--prs"), action="store", default=NA, type='character',
              help="PRS score"),
  make_option(c("-m", "--max"), action="store", default=NA, type='character',
              help="max PRS score"),
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="Default directory")
)

opt = parse_args(OptionParser(option_list=option_list))
prs <- as.numeric(opt$prs)
max_prs <- as.numeric(opt$max)
outdir <- opt$dir

#===============================================================================
#   Generating a cutted, normal distribution
#===============================================================================
mn <- max_prs/2
# Generate numbers (mostly) from min to max
distribution <- rnorm(100000, mean = mn, sd = mn/3)
# Do something about the out-of-bounds generated values
distribution <- pmax(0, distribution)
distribution <- pmin(max_prs, distribution)

# Scale the prs according to the normal distribution (z-score)
prs_distribution_scaled <- data.frame("zscore" = scale(c(prs,distribution)))
scaled_prs <- round(prs_distribution_scaled$zscore[1], digits = 2)
breaks = c(scaled_prs,unique(round(prs_distribution_scaled$zscore), digits = 2))

report = paste0("Scaled PRS: ",scaled_prs, "\n", 
                "Original PRS: ", prs, "\n",
                "Max PRS: ",round(max_prs, digits = 2), "\n")

#===============================================================================
#   Plot values
#===============================================================================

p <- ggplot(data = prs_distribution_scaled, aes(x=zscore)) +
  geom_density(color="darkblue", fill="lightblue", alpha=0.5) +
  geom_vline(xintercept=scaled_prs, color="blue", linetype="dashed", size=1) +
  ggtitle("scaled PRS vs Normal Distribution (0, maxPRS)") +
  scale_x_continuous(breaks=breaks, labels = breaks) +
  theme(axis.text.x=element_text(color = "black", size=8, angle=45, vjust = 0.5)) +
  annotate("text", x=Inf, y = Inf, label = report, vjust=1, hjust=1)

plot_fn = paste(outdir,"scaled_prs_distribution.png", sep = "/")
ggsave(filename = plot_fn, plot = print(p))