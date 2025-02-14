# FCOR 599 - Liam Irwin
# ShinyApps and Point Cloud Visualization
# How to make a 3D rotating GIF from a LAS/LAZ point cloud in R

# Function to create a rotating GIF of a point cloud
create_rotating_gif <- function(file_path, # path to LAS/LAZ file
                                col_by = "Z", # color by Z, RGB, Classification, etc.
                                decimate_to = NULL, # density to decimate to (NULL = no decimation)
                                bg_col = "white", # background color
                                photo_d = 4, # degrees to rotate per photo
                                fps = 10, # frames per second for output GIF
                                delay = NULL, # delay between frames in output GIF (overrides FPS)
                                rgl_width = 500, # width of rgl window (pixels)
                                rgl_height = 500, # height of rgl window (pixels)
                                out_dir = NULL, # output directory for GIF, defaults to file's directory
                                sleep_s = 0.001) # time to sleep between snapshots
{
  # Get name of las file for saving
  plot_name <- tools::file_path_sans_ext(basename(file_path))

  # Start timer
  tictoc::tic(msg = glue::glue("Generated gif for {plot_name} at {fps}fps"))

  # List of required packages
  req_pkgs <- c("lidR", "magick", "rgl", "glue", "tictoc", "stringr")

  # Identify any packages that are missing
  missing_pkgs <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]

  # Stop if any required packages are missing
  if (length(missing_pkgs) > 0) {
    stop(glue::glue("The following packages are required but not installed: {missing_pkgs}"))
  }

  # Load required packages
  lapply(req_pkgs, require, character.only = TRUE)

  # Set output directory to the file's directory if not provided
  if (is.null(out_dir)) {
    out_dir <- dirname(file_path)
  }

  # Create output directories if they don't exist
  snapshot_dir <- file.path(out_dir, "gif", "snapshot") # Where pngs are saved
  animation_dir <- file.path(out_dir, "gif", "out") # Where final GIF is saved

  if (!dir.exists(animation_dir)) {
    dir.create(animation_dir, recursive = TRUE)
  }

  # Remove any existing snapshot directory (delete intermediate files)
  if (dir.exists(snapshot_dir)) {
    unlink(snapshot_dir, recursive = TRUE)
    message(glue::glue("Deleted existing snapshot directory at {snapshot_dir}"))
  }

  # Create a new snapshot directory
  dir.create(snapshot_dir, recursive = TRUE)

  # Read LAS/LAZ file and optionally decimate points to chosen density
  las <- lidR::readLAS(file_path)
  if (!is.null(decimate_to)) {
    las <- lidR::decimate_points(las, lidR::random(decimate_to))
    message(glue::glue("Decimated point cloud to {decimate_to} points"))
  }

  # Set up rotation parameters
  num_frames <- round(360 / photo_d) # Number of frames to capture
  angle_list <- rep(photo_d * pi / 180, num_frames) # List of angles to rotate by

  # Plot the point cloud and set the window size
  lidR::plot(las, color = col_by, bg = bg_col)
  # Create an rgl device with the same window size
  par3d(windowRect = c(10, 5, rgl_width, rgl_height))

  # Capture snapshots for each rotation angle
  for (i in seq_along(angle_list)) {
    view3d(userMatrix = rotate3d(par3d("userMatrix"), angle_list[i], 0, 0, 1))
    Sys.sleep(sleep_s)
    snapshot_path <- file.path(snapshot_dir, sprintf("frame_%03d.png", i))
    rgl::rgl.snapshot(filename = snapshot_path)
  }

  # Close the rgl plotting device
  rgl::close3d()

  # Assemble snapshots into an animated GIF
  message(glue::glue("Finished capturing snapshots, assembling animated GIF... be patient"))
  # List all the snapshots RGL just took
  img_files <- list.files(snapshot_dir, full.names = TRUE)
  # Read them into R using magick
  img_list <- lapply(img_files, magick::image_read)

  img_joined <- magick::image_join(img_list)

  # Animate the series of images
  img_animated <- magick::image_animate(img_joined, fps = fps,
                                        delay = delay, optimize = TRUE)


  # Construct a file_name based on the input arguments
  file_name <- glue::glue("{plot_name}_{col_by}_{rgl_width}x{rgl_height}_{bg_col}_{photo_d}degree_{fps}fps.gif")

  # Save the animated GIF
  magick::image_write(image = img_animated, path = file.path(animation_dir, file_name))

  # Remove the old snapshot directory
  unlink(snapshot_dir, recursive = TRUE)

  tictoc::toc()  # End timer
  message(glue::glue("Saved animated GIF to {animation_dir}/{file_name}"))
}


# Change this path to where your LAS/LAZ file is located
las_path <- 'data/las_ex/1_pointcloud_norm.laz'

# Try the other example las point clouds (or with your own!)
# las_path <- 'data/las_ex/2_pointcloud_norm.laz'
# las_path <- 'data/las_ex/3_pointcloud_norm.laz'
# las_path <- 'data/las_ex/4_pointcloud_norm.laz'
# las_path <- 'data/las_ex/5_pointcloud_norm.laz'

# Tree colouring the point cloud by other attributes

col_by <- "Z"
# col_by <- "RGB"
# col_by <- "treeID"
# col_by <- "intensity"

bg_col <- "black"
# bg_col = 'white'
# Feel free to modify the function to add more custom parameters for plotting




# Run this to generate and save a rotating GIF!
create_rotating_gif(las_path, col_by = col_by, delay = 5,
                    bg_col = bg_col,
                    photo_d = 2,
                    rgl_width = 750, rgl_height = 750)
