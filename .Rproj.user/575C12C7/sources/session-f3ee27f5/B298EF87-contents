library(shiny)
library(terra)
library(sf)
library(lidR)
library(tidyterra)
library(ggplot2)
library(ggspatial)
library(colourpicker)
library(here)


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

# ---- Get paths to example CHM files (using relative path) ----
# 'here::here("data", "chm_ex")' will locate the 'data/chm_ex' folder relative to the script
chm_folder <- here::here("data", "chm_ex")
chm_paths  <- list.files(chm_folder, pattern = '\\.tif$', full.names = TRUE)
example_choices <- setNames(chm_paths, basename(chm_paths))

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

# Set max file size to 25MB
options(shiny.maxRequestSize = 25 * 1024^2)


# ---- Shiny Server ----
server <- function(input, output, session){

  rv <- reactiveValues(chm = NULL, ttops = NULL)

  # Show modal dialog when the app starts
  observe({
    showModal(modalDialog(
      title = "Welcome to CHM Tree Top Detection",
      "Please upload a CHM file or choose one of the example CHMs to begin.",
      easyClose = TRUE,
      footer = modalButton("OK")
    ))
  })

  # A helper function to update CHM and treetops
  updateChmAndTtops <- function(chmPath=NULL, chmFile=NULL,
                                smooth=FALSE, variable_ws=FALSE,
                                lmf_ws=2, hmin=5) {

    if (!is.null(chmFile) && chmFile != "") {
      chm <- rast(chmFile)
    } else if (!is.null(chmPath) && chmPath != "") {
      chm <- rast(chmPath)
    } else {
      showNotification("Please upload a CHM or select an example CHM", type = "error")
      return(NULL)
    }

    if (smooth) {
      chm <- terra::focal(chm, w = fgauss(sigma = 1, n = 5))
      names(chm) <- "Z"
    }

    if (variable_ws) {
      ttops <- locate_trees(
        chm,
        algorithm = lmf(ws = function(x) { x * 0.1 + 3 }, hmin = hmin)
      )
    } else {
      ttops <- locate_trees(
        chm,
        algorithm = lmf(ws = lmf_ws, hmin = hmin)
      )
    }

    list(chm = chm, ttops = ttops)
  }

  # Observer for "Update Plot" button
  observeEvent(input$go_button, {
    res <- updateChmAndTtops(
      chmPath   = input$example_chm_radio,
      chmFile   = if (!is.null(input$chm_file)) input$chm_file$datapath else NULL,
      smooth    = input$smooth_chm,
      variable_ws = input$variable_ws,
      lmf_ws    = input$lmf_ws,
      hmin      = input$hmin
    )
    if (!is.null(res)) {
      rv$chm   <- res$chm
      rv$ttops <- res$ttops
    }
  })

  # Automatically update plot when user picks an example CHM
  observeEvent(input$example_chm_radio, {
    if (input$example_chm_radio != "") {
      res <- updateChmAndTtops(
        chmPath   = input$example_chm_radio,
        chmFile   = NULL,
        smooth    = input$smooth_chm,
        variable_ws = input$variable_ws,
        lmf_ws    = input$lmf_ws,
        hmin      = input$hmin
      )
      if (!is.null(res)) {
        rv$chm   <- res$chm
        rv$ttops <- res$ttops
      }
    }
  })

  # Reactive expression for final ggplot
  plot_obj <- reactive({
    req(rv$chm, rv$ttops)
    plot_chm_ttops_tidy(
      chm         = rv$chm,
      ttops       = rv$ttops,
      show_ttops  = input$show_ttops,
      ttops_color = input$ttops_color,
      ttops_shape = as.numeric(input$ttops_shape),
      ttops_size  = input$ttops_size,
      ttops_alpha = input$ttops_alpha,
      colormap    = input$colormap
    )
  })

  output$chm_plot <- renderPlot({
    req(plot_obj())
    plot_obj()
  })

  # --- Download Plot (PNG) ---
  output$download_plot <- downloadHandler(
    filename = function() {
      req(rv$chm)
      res_vals <- res(rv$chm)
      res_str  <- paste0(round(res_vals[1], 2), "-", round(res_vals[2], 2))
      smoothed <- if (input$smooth_chm) "YES" else "NO"
      paste0("chm_plot_ws-",
             ifelse(input$variable_ws, "variable", input$lmf_ws),
             "_smooth-", smoothed,
             "_res-", res_str, ".png")
    },
    content = function(file) {
      ggsave(
        filename = file,
        plot = plot_obj(),
        device = "png",
        width = 10, height = 7, dpi = 300
      )
    }
  )

  # --- Download Treetops (GPKG) ---
  output$download_ttops <- downloadHandler(
    filename = function() {
      req(rv$chm)
      res_vals <- res(rv$chm)
      res_str  <- paste0(round(res_vals[1], 2), "-", round(res_vals[2], 2))
      smoothed <- if (input$smooth_chm) "YES" else "NO"
      paste0("treetops_ws-",
             ifelse(input$variable_ws, "variable", input$lmf_ws),
             "_smooth-", smoothed,
             "_res-", res_str, ".gpkg")
    },
    content = function(file) {
      req(rv$ttops)
      writeVector(vect(rv$ttops), file, filetype = "GPKG", overwrite = TRUE)
    }
  )
}
