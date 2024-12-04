using SpeciesDistributionToolkit
using CairoMakie

rpath = joinpath(pwd(), "rasters")

# GADM data for QC municipalities
QC = SpeciesDistributionToolkit.gadm("CAN", 3)
MNAM = SpeciesDistributionToolkit.gadmlist("CAN", "Québec")

# Get a map
mpc = SDMLayers.read(joinpath(rpath, "Mephitis-mephitis-prediction-SSP245-2021-2040.tiff"))

# Zonal stat
