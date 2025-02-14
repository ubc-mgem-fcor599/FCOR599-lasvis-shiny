library(shiny)
library(terra)
library(sf)
library(lidR)
library(tidyterra)
library(ggplot2)
library(ggspatial)
library(colourpicker)

# ---- Helper function for smoothing CHM (Gaussian kernel) ----
fgauss <- function(sigma, n = 5) {
  m <- matrix(ncol = n, nrow = n)
  col <- rep(1:n, n)
  row <- rep(1:n, each = n)
  x <- col - ceiling(n/2)
  y <- row - ceiling(n/2)
  m[cbind(row, col)] <- 1 / (2 * pi * sigma^2) * exp(-(x^2 + y^2) / (2 * sigma^2))
  m / sum(m)
}

# ---- Get paths to example CHM files ----
chm_paths <- list.files('data/chm_ex', pattern = '\\.tif$', full.names = TRUE)
# Create a named vector for the dropdown menu (display names are the file basenames)
example_choices <- c("Select an example" = "", setNames(chm_paths, basename(chm_paths)))

# ---- Plotting function with extended symbology options ----
plot_chm_ttops_tidy <- function(chm, ttops,
                                show_ttops   = TRUE,
                                ttops_color  = "red",
                                ttops_shape  = 16,
                                ttops_size   = 2,
                                ttops_alpha  = 1.0,
                                colormap     = "viridis")
{
  p <- ggplot() +
    geom_spatraster(data = chm) +
    labs(fill = "Canopy Height (m)") +
    scale_fill_viridis_c(option = colormap) +
    annotation_scale(
      height = unit(0.015, "npc"),
      width_hint = 0.25,
      pad_x = unit(0.07, "npc"),
      pad_y = unit(0.07, "npc"),
      text_cex = .8,
      text_col = 'white',
      line_col = 'white'
    ) +
    theme_void()

  if (show_ttops) {
    p <- p + geom_spatvector(
      data = vect(ttops),
      shape = ttops_shape,
      size = ttops_size,
      alpha = ttops_alpha,
      color = ttops_color
    )
  }
  return(p)
}

# Set max file size to 100MB
options(shiny.maxRequestSize = 100 * 1024^2)


# ---- Shiny UI ----
ui <- fluidPage(
  titlePanel("CHM Tree Top Detection"),

  fluidRow(
    column(width = 3,
           wellPanel(
             fileInput("chm_file", "Upload CHM TIF",
                       accept = c(".tif", ".tiff")),
             radioButtons("example_chm_radio", "Select Example CHM:",
                          choices = example_choices, selected = character(0)),
             br(),
             sliderInput("lmf_ws", "Window size (ws):",
                         min = 1, max = 10, value = 2),
             checkboxInput("variable_ws",
                           label = HTML("Variable Window Size (<a href='https://doi.org/10.14358/PERS.70.5.589' target='_blank'>Popescu and Wynne 2004</a>)"),
                           value = FALSE),
             numericInput("hmin", "Minimum height (hmin):",
                          value = 5, min = 0, max = 50, step = 1),
             checkboxInput("smooth_chm", "Smooth CHM before detection", value = FALSE),

             tabsetPanel(
               tabPanel("CHM Options",
                        selectInput("colormap", "Color Palette:",
                                    choices = c("viridis", "magma", "plasma", "inferno", "cividis"),
                                    selected = "viridis")
               ),
               tabPanel("Tree Top Options",
                        checkboxInput("show_ttops", "Show Tree Tops", value = TRUE),
                        colourInput("ttops_color", "Tree Top Color:", value = "red"),
                        selectInput("ttops_shape", "Tree Top Shape:",
                                    choices = c("circle" = 16, "square" = 15,
                                                "triangle" = 17, "diamond" = 18,
                                                "plus" = 3),
                                    selected = 16),
                        sliderInput("ttops_size", "Tree Top Size:", min = 0.5, max = 5,
                                    value = 2, step = 0.5),
                        sliderInput("ttops_alpha", "Tree Top Transparency:", min = 0, max = 1,
                                    value = 1, step = 0.1)
               )
             ),
             hr(),
             downloadButton("download_plot", "Download Plot (PNG)"),
             br(), br(),
             downloadButton("download_ttops", "Download Tree Tops (GPKG)"),
             hr(),
             actionButton("go_button", "Update Plot", class = "btn btn-success btn-lg")
           )
    ),
    column(width = 9,
           plotOutput("chm_plot", height = "700px")
    )
  ),

  absolutePanel(
    bottom = 5, right = 5,
    style = "background-color: rgba(255,255,255,0.7); padding: 5px; border-radius: 4px;",
    "Created by Liam Irwin 2025 for University of British Columbia MGEM Program"
  )
)
