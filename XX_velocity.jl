using SpeciesDistributionToolkit
using Statistics
using CairoMakie

provider = RasterData(WorldClim2, BioClim)
future = Projection(SSP126, CanESM5)

POL = SpeciesDistributionToolkit.openstreetmap("Montérégie")
bb = SpeciesDistributionToolkit.boundingbox(POL; padding=1.0)

cclim = [SDMLayer(provider; resolution=2.5, layer=i, bb...) for i in [1, 12]]
fclim = [SDMLayer(provider, future; resolution=2.5, layer=i, bb...) for i in [1, 12]]

heatmap(cclim[1])

# Get the pixels that are closest to the ref one

ref = fclim
tar = cclim
output = similar(ref[1], Float64)

for k in keys(output)
    
end