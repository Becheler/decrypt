#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("At least two argument must be supplied (input file and max value)", call.=FALSE)
}
file=args[1]
max_N_value=strtoi(args[2])

print(file)

library("raster")
library("viridis")

make_movie <- function(history, ordered_times, max_N_value, working_folder){
  old_wd <- getwd()
  setwd(working_folder)
  for(i in ordered_times){
    png( paste0("N_", i, ".png"))
    print( plot(history[[i]], zlim=c(0,max_N_value), col = viridis(100), colNA="lightgrey", col) )
    dev.off() #only 129kb in size
  }
  system('n=0; ls -tr | while read i; do n=$((n+1)); mv -- "$i" "$(printf \'%04d\' "$n")"_"$i"; done')
  system("gm convert -delay 20 *.png animation.gif")
  setwd(old_wd)
}


message("Reading the demographic history geotiff file created by spatial_process")
history <- stack(file)

## MOVIE (may take some time like 3 minutes)
# Create a directory to store intermediary png files
dir.create("temp_movie")

# Precise the working directory to generate the demographic plots
working_folder <- paste0(getwd(),"/temp_movie")

# the time range to plot
ordered_times <- 1:nlayers(history)

# Standardize the plots legends with an expected maximal N value in the dataset
# like the maximal carrying capacity

message("Generating history movie")
# generate a MP4 file in the movie directory. Requires the ImageMagick package to be installed
make_movie(history, ordered_times, max_N_value, working_folder)

copied <- file.copy(from="temp_movie/animation.gif", to=getwd(), overwrite=TRUE)
message(paste("animation.gif generated at", getwd()))
removed <- file.remove(file.path("temp_movie", list.files("temp_movie")))
unlinked <- unlink(working_folder <- paste0(getwd(),"/temp_movie"), recursive=TRUE)
