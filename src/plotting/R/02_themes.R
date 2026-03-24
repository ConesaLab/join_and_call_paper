# 02_themes.R
# Color palettes, structural category labels, and reusable ggplot themes

xaxislevelsF1 <- c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "genic", "antisense", "fusion", "intergenic", "genic_intron")
xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic", "Antisense", "Fusion", "Intergenic", "Genic\nIntron")
cat.palette <- c("FSM" = "#6BAED6", "ISM" = "#FC8D59", "NIC" = "#78C679", "NNC" = "#EE6A50", "Genic\nGenomic" = "#969696", "Antisense" = "#66C2A4", "Fusion" = "goldenrod1", "Intergenic" = "darksalmon", "Genic\nIntron" = "#41B6C4")

theme_Publication <- function(base_size=14, base_family="sans") {
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5, margin = margin(0,0,0,0)),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(),
               axis.line.x = element_line(colour="black"),
               axis.line.y = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.box = "vertical",
               legend.key.size= unit(0.5, "cm"),
               #legend.margin = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(0,0,0,0),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold"),
       ))
}

bar_theme <- function(plot, title, ylims) {
  plot +
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    ) +
    coord_cartesian(ylim = ylims)
}

bar_theme_noaxis <- function(plot, title, ylims) {
  bar_theme(plot, title, ylims) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

compare_theme <- function(plot, title) {
  plot +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 12),
      legend.position = "none"
    ) +
    labs(x = "Strategy")
}

compare_theme_noaxis <- function(plot, title) {
  compare_theme(plot, title) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}
