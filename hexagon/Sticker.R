#First generate Tikz Picture:

# \documentclass[tikz,convert={outfile=\jobname.svg}]{standalone}
# %\usepackage{tikz}
# \begin{document}
# \begin{tikzpicture}
# \node (bottomleft) [ draw, line width = 3pt, rectangle, minimum width=2cm, minimum height = 1.5cm] at (0,0cm) {};
# \node (topleft) [ draw, line width = 3pt, rectangle, minimum width=2cm, minimum height = 1.5cm] at (1.5,2.5cm) {};
# \node (topright) [ draw, line width = 3pt, rectangle, minimum width=2cm, minimum height = 1.5cm] at (5,1.5cm)  {};
# \node (bottomright) [ draw, line width = 3pt, rectangle, minimum width=2cm, minimum height = 1.5cm] at (3.5, -1cm) {    };
# \draw[->, line width = 3pt] (bottomleft) -- (topleft);
# \draw[->, line width = 3pt] (topleft) -- (topright);
# \draw[->, line width = 3pt] (bottomleft) -- (bottomright);
# \end{tikzpicture}
# \end{document}


library(showtext)
## Loading Google fonts (fonts.google.com)
#font_add_google("Source Code Pro", "myfont")
font_add_google("PT mono", "myfont")

library(hexSticker)

col = 'aquamarine2'
#col = 'greenyellow'
#col = 'deepskyblue4'
col = 'dodgerblue4'
asd <- getwd()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sticker("icmstate_hex.png", dpi = 1000,
        s_x=1, s_y=1, s_width = 0.7, s_height = 0.7,
        package = "mstate", p_size = 1, p_color = "black",
        p_family = 'myfont',
        p_x = 1, p_y = 1,
        h_fill = 'aquamarine2', h_color = "dodgerblue4",
        asp = 0.9, # default asp = 1 sformatta input figure!!!
        filename="icmstate_logo.png",
        white_around_sticker = FALSE)
setwd(asd)

asd <- getwd()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sticker("icmstate_hex2.png", dpi = 1000,
        s_x=1, s_y=1, s_width = 1.5, s_height = 1.5,
        package = "mstate", p_size = 1, p_color = "black",
        p_family = 'myfont',
        p_x = 1, p_y = 1,
        h_fill = 'aquamarine2', h_color = "dodgerblue4",
        asp = 0.9, # default asp = 1 sformatta input figure!!!
        filename="icmstate_logo2.png",
        white_around_sticker = TRUE)
setwd(asd)


empty_plot <- ggplot() + theme_void()


sticker(empty_plot, dpi = 1000,
        s_x=1, s_y=1, s_width = 1.5, s_height = 1.5,
        package = "mstate", p_size = 1, p_color = "black",
        p_family = 'myfont',
        p_x = 1, p_y = 1,
        h_fill = 'white', h_color = "black",
        asp = 0.9, # default asp = 1 sformatta input figure!!!
        filename="empty_hex.png",
        white_around_sticker = TRUE)




