# Number of sampling scheme simulations under 1 demographic history"
n_sim_gen=5

# Number of loci
n_loci=5

# Introduction point latitude
lat_0=-20.0

# Introduction point longitude
lon_0=125.0

# Number of gene copies at introduction point
N_0=1000

# Number of generations to simulate
duration=500

# Fixed sampling point latitude
lat_1=-20.0

# Fixed sampling point longitude
lon_1=125.0

# Number of gene copies to sample around center 1
n_sample_1=30

# Sampling radius around center 1
radius_sample_1=30.0

# Number of gene copies to sample in population 2
n_sample_2=30

# Sampling radius around center 2
radius_sample_2=30.0

# Population size under which random sampling is not considered
sampling_threshold=30

# Environmental threshold delimiting suitable (above threshold) and unsuitable (under threshold) areas
suitability_threshold=26.4

# Carrying capacity in suitable areas
K_max=50

# Population persistance parameter in unsuitable areas
p=0.175

# Carrying capacity in unsuitable areas with probability 1-p
K_min_a=1

# Carrying capacity in unsuitable areas with probability p
K_min_b=20

# Constant growth rate
r=1

# Emigrant rate between the four neighboring cells
emigrant_rate=0.1

# Friction coefficient in suitable areas
friction_suitable=0.3

# Friction coefficient in unsuitable areas
friction_unsuitable=0.7

# File name for the simulated demography output
demography_out=output/N.tif

# File name for last demographic layer at sampling time
last_layer_out=output/last_N.tif

# File name for coordinates at which population size is positive at sampling time
distribution_area_out=output/distribution_area.shp

# File name for coordinates at which random sampling of individuals is operable
mask_out=output/mask.shp

# File name for the simulated sampling scheme output
sample_out=output/sample.shp

# Filename database storing the output
database=output/test.db
