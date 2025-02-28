using SpeciesDistributionToolkit
using CairoMakie

# Dispersal data from Schloss 2012 in km/year

v_raccoon = 3.48
v_skunk = 3.21

# TODO get current range, and expand by 20 years increments to get range limitation
raccoon = SDMLayer("rasters/Procyon_lotor_historical.tif"; bandnumber=2) .== 1
d_raccoon = pseudoabsencemask(DistanceToEvent, raccoon)

skunk = SDMLayer("rasters/Mephitis_mephitis_historical.tif"; bandnumber=2) .== 1
d_skunk = pseudoabsencemask(DistanceToEvent, skunk)

heatmap(d_skunk .<= (v_skunk * 100))
heatmap(d_raccoon .<= (v_raccoon * 100))