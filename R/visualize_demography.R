#!/usr/bin/env Rscript
houses <- c("Hufflepuff", "Gryffindor", "Ravenclaw", "Slytherin")
house <- sample(houses, 1)
cat(house, "\n")
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
file=args[1]
print(file)

packages <- c("raster","devtools","dismo","gstat","viridis","rspatial")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())), dependencies = T)
}

make_movie <- function(history, max_N_value, working_folder){
  old_wd <- getwd()
  setwd(working_folder)
  for(i in ordered_times){
    png( paste0("N_", i, ".png"))
    print( plot(history[[i]], zlim=c(0,max_N_value), col = viridis(100), colNA="lightgrey", col) )
    dev.off() #only 129kb in size
  }
  system('n=0; ls -tr | while read i; do n=$((n+1)); mv -- "$i" "$(printf \'%04d\' "$n")"_"$i"; done')
  system("convert -delay 2 -quality 95 *.png movie.mp4")
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
ordered_times <- 1:nlayer(file)

# Standardize the plots legends with an expected maximal N value in the dataset
# like the maximal carrying capacity
max_N_value <- max(history)

message("Generating history movie")
# generate a MP4 file in the movie directory. Requires the ImageMagick package to be installed
make_movie_2(history, ordered_times, max_N_value, working_folder)

copied <- file.copy(from="temp_movie/movie.mp4", to=getwd(), overwrite=TRUE)
message(paste("movie.mp4 generated at", getwd()))
removed <- file.remove(file.path("temp_movie", list.files("temp_movie")))
unlinked <- unlink(working_folder <- paste0(getwd(),"/temp_movie"), recursive=TRUE)
