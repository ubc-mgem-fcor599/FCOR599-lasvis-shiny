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


# ---- Shiny Server ----
server <- function(input, output, session){

  rv <- reactiveValues(chm = NULL, ttops = NULL)

  # Show modal dialog when the app starts
  observe({
    showModal(modalDialog(
      title = "Welcome to CHM Tree Top Detection",
      "Please upload a CHM file or select an example from the drop-down to begin.",
      easyClose = TRUE,
      footer = modalButton("OK")
    ))
  })

  # Observer for the Update Plot button
  observeEvent(input$go_button, {
    # Use the uploaded file if provided; otherwise use the selected example if available
    if (!is.null(input$chm_file) && input$chm_file$datapath != "") {
      chm <- rast(input$chm_file$datapath)
    } else if (input$example_chm != "") {
      chm <- rast(input$example_chm)
    } else {
      showNotification("Please upload a CHM or select an example CHM", type = "error")
      return(NULL)
    }

    # Optionally smooth the CHM before detection
    if (input$smooth_chm) {
      chm <- terra::focal(chm, w = fgauss(sigma = 1, n = 5))
      names(chm) <- "Z"
    }

    # Determine tree detection algorithm based on the Variable Window Size option
    if (input$variable_ws) {
      # Using a variable window size function: f(x) = 0.1*x + 3
      # Optionally, compute the maximum height (z_max) to illustrate the method:
      z_max <- terra::global(chm, "max", na.rm = TRUE)[1,1]
      # heights <- seq(0, z_max, by = 5)
      # ws_values <- (heights * 0.1 + 3)  # (For demonstration; the function below is passed to lmf)

      ttops <- locate_trees(chm, algorithm = lmf(ws = function(x) { x * 0.1 + 3 }, hmin = input$hmin))
    } else {
      ttops <- locate_trees(chm, algorithm = lmf(ws = input$lmf_ws, hmin = input$hmin))
    }

    rv$chm <- chm
    rv$ttops <- ttops
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
      paste0("chm_plot_ws-", ifelse(input$variable_ws, "variable", input$lmf_ws),
             "_smooth-", smoothed,
             "_res-", res_str,
             ".png")
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
      paste0("treetops_ws-", ifelse(input$variable_ws, "variable", input$lmf_ws),
             "_smooth-", smoothed,
             "_res-", res_str,
             ".gpkg")
    },
    content = function(file) {
      req(rv$ttops)
      writeVector(vect(rv$ttops), file, filetype = "GPKG", overwrite = TRUE)
    }
  )
}
