library(ggplot2)
library(cowplot)
library(scales)
library(tidyqpcr)

# set default settings for plots

theme_set(
  theme_cowplot(font_size = 12, 
                font_family = "sans",
                rel_small = 10/12,
                rel_tiny = 9/12,
                rel_large = 12/12) %+replace% 
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          panel.border=element_rect(colour = "grey50",linetype = "solid",size=0.5),
          panel.grid.major = element_line(size = 0.15, linetype = 'solid',colour = "grey"),
          strip.background = element_blank(),
          plot.title=element_text(face = "bold", hjust = 0.5),
          strip.text=element_text(size = 10, face = "bold")
    )
)

# set default figure widths following cell press guidelines
fig_width_2column = 165
fig_dpi = 300

# Update geom defaults for consistent sizes

update_geom_defaults("point", list(size = 2))
update_geom_defaults("pointrange", list(size = 1, fatten = 1.5))
update_geom_defaults("linerange", list(size = 1))
update_geom_defaults("errorbar", list(size = 1))
update_geom_defaults("errorbarh", list(size = 1))
text_cor_size = 3.5

# default diagonal line
geom_diagline <- function(slope = 1, 
                          intercept = 0, 
                          linetype = "dashed", 
                          colour = "black",
                          size = 0.2,
                          ...) {
  geom_abline(slope = slope,
              intercept = intercept,
              linetype = linetype,
              colour = colour,
              size = size,
              ...)
}
# default vertical line with intercept 1
geom_vline1 <- geom_vline(xintercept=1, linetype="dashed", color = "black",size = 0.2)

# Set general variables for motif construct names and colour schemes

construct_to_label_dictionary_TSA1_RPS3 <-  
  tibble(mod = c("WT","modC","modE","modD","modA","modB","mod0"), 
         label = c("WT", "mod_NGG", "mod_HTH", "mod_HNH", "mod_NTN", "mod_NAA", "mod_NNN"))

construct_to_label_dictionary_PIR1 <- 
  tibble(construct = c("modG", "modF", "modE", "modD", "modC", "modA", "modB", "WT"), 
         label = c("mod_NTNNN", "mod_ANNNN", "mod_ATNNN", "mod_ATHNH", "mod_ATNHH", "mod_ANHHH", "mod_NTHHH", "WT"))

RPS3_TSA1_colour_scheme <- c( "WT"      = "#6c6c6c", 
                              "mod_NGG" = "#7FBF74",
                              "mod_HTH" = "#C288BB", 
                              "mod_HNH" = "#71519C",
                              "mod_NTN" = "#5978BA",
                              "mod_NAA" = "#C05558",
                              "mod_NNN" = "black")

PIR1_colour_scheme <- c("mod_NTNNN" = "#DD76A5",
                        "mod_ANNNN" = "#968BC2",
                        "mod_ATNNN" = "#71519C",
                        "mod_ATHNH" ="#D97F1D",
                        "mod_ATNHH" = "#4EB0B5", 
                        "mod_ANHHH" = "#5978BA",
                        "mod_NTHHH" = "#C15659",
                        "WT" = "#6c6c6c")

# Create plotting functions for  specific data sets

RNA_relative_abundance_figure_options <- list(
  geom_vline1,
  geom_point(aes(x=rel_abund_delta_deltacq,
                 y=mod_label,colour=mod_label),
             shape = 18, size = 2),
  scale_x_log2nice(omag = seq(-5,5),scilabels=FALSE),
  guides(colour=FALSE),
  theme(axis.text.x=element_text(angle=0,vjust=0.5),
        legend.position="bottom",
        legend.box.margin=margin(20,10,10,180),
        strip.text.x = element_text(vjust = 0.95)),
  stat_summary(aes(x=rel_abund_delta_deltacq,y=mod_label),
               fun="mean",colour="black",
               geom="crossbar", size=0.3, width=0.9)
  )

#protein_raw_abundance_figure_options <- list(
#  geom_vline1,
#  geom_point(aes(y=Terminator,x=fluo_per_OD_at_max_gr, colour = Terminator),
#             shape = 18, size = 2),
#  scale_colour_hue(h = c(0, 360)+20,l=60,c=60),
#  stat_summary(aes(y=Terminator,x=fluo_per_OD_at_max_gr),
#               fun="mean",colour="black",
#               geom="crossbar", size=0.3, width=0.9),
#  theme(legend.position = "none",
#        strip.text.x = element_text(vjust = 0.95)),
#  scale_x_continuous(oob=scales::squish(), limits = c(0,NA))
#)

#protein_relative_abundance_figure_options <- list(
#  scale_colour_hue(h = c(0, 360)+20,l=60,c=60),
#  geom_vline1,
#  stat_summary(aes(y=Terminator,x=fluo_per_OD_at_max_gr, colour = Terminator),
#               fun.data="mean_se",
#               geom="pointrange"),
#  theme(legend.position = "none", 
#        strip.text.x = element_text(vjust = 0.95)),
#  scale_x_continuous(oob=scales::squish(), limits = c(0,1.9))
#)
#
#protein_vs_RNA_figure_options <- list(
#    geom_diagline(),
#    geom_point(aes(y = mean_relative_abundance_protein, x = mean_relative_abundance_mrna, colour = label)),
#    geom_linerange(aes(ymax = mean_relative_abundance_protein + se_relative_abundance_protein,
#                      ymin = mean_relative_abundance_protein - se_relative_abundance_protein,
#                      x = mean_relative_abundance_mrna, 
#                      colour = label)),
#    geom_linerange(aes(xmax = mean_relative_abundance_mrna + se_relative_abundance_mrna,
#                       xmin = mean_relative_abundance_mrna - se_relative_abundance_mrna,
#                       y = mean_relative_abundance_protein, 
#                       colour = label)),
#    facet_wrap(~ promoter, ncol = 1),
#    theme(strip.text.x = element_text(vjust = 0.95)),
#    scale_x_log2nice(),
#    scale_y_log2nice()
#)

