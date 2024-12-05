using SpeciesDistributionToolkit
using CairoMakie

rpath = joinpath(pwd(), "rasters")

# GADM data for QC municipalities
QC = SpeciesDistributionToolkit.gadm("CAN", "Québec")
MNAM = SpeciesDistributionToolkit.gadmlist("CAN", 3, "Québec")
MPOL = SpeciesDistributionToolkit.gadm("CAN", 3, "Québec")

# Get a map
mpc = SimpleSDMLayers._read_geotiff(joinpath(rpath, "Mephitis-mephitis-prediction-current.tiff"))
mpf = SimpleSDMLayers._read_geotiff(joinpath(rpath, "Mephitis-mephitis-prediction-SSP126-2021-2040.tiff"))

# Zonal stat
mmed = mosaic(median, mpc, MPOL)