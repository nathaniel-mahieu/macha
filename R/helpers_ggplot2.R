theme_nate <- function(base_size = 16)
{
  theme_grey(base_size = base_size) %+replace%
    theme(
      strip.background = element_rect(colour="grey95"),
      strip.text = element_text(size = base_size*0.5, colour="grey40"),

      panel.background = element_blank(),
      panel.grid.major = element_line(colour="grey95"),
      panel.grid.minor = element_line(colour="grey99"),
      plot.title = element_text(size = base_size * 0.75),

      axis.title.y = element_text(vjust=1, angle=90),
      axis.title = element_text(size = base_size*0.65, colour="grey40"),
      axis.line = element_line(colour="grey95"),
      axis.ticks = element_line(colour="grey95"),
      axis.text = element_text(size = base_size*0.5, colour = "grey40"),

      legend.position="top",
      legend.key = element_blank(),
      legend.text = element_text(colour="grey40"),
      legend.title  = element_text(colour="grey40")
    )
}

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)

  l <- gsub("0e\\+00","0",l)
  l <- gsub("e\\+","e",l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  #l <- gsub("e", "%*%10^", l)
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}


extrafont::loadfonts(device="pdf")
extrafont::loadfonts(device="postscript")
#extrafont::fonttable()
#extrafont::font_import("C:/Windows/Fonts/", pattern = "RobotoCondensed", prompt = F)
#C:/Program Files/R/R-3.3.1/library/extrafontdb/fontmap/
#Change lights to "Roboto Condensed Light"
#After ggsave(device="pdf") open and resave the file in Illustrator
library(gridExtra)
library(hrbrthemes)
library(ggplot2)
library(ggrepel)
theme_nate <- function(...)
{
  theme_ipsum_rc(...) %+replace%
    theme(
      panel.background = element_rect(fill="white", colour = "white"),
      panel.grid.major = element_line(colour="grey90"),
      panel.grid.minor = element_line(colour="grey95"),
      
      legend.position="top",
      legend.key = element_blank()
    )
}
theme_set(theme_nate())
