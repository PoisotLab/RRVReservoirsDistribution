using SpeciesDistributionToolkit
using CairoMakie
import DelimitedFiles
using PrettyTables
using Statistics
import CSV

include("S1_theme.jl")

QC = getpolygon(PolygonData(GADM, Countries), level=1, country="CAN")["Québec"]
provider = RasterData(WorldClim2, BioClim)

lnames = layers(provider)
descr = layerdescriptions(provider)

function _fromshap(taxname)
    v = Int64[DelimitedFiles.readdlm("data/$(replace(taxname, " " => "_")).params")...]
    S = [SDMLayer("rasters/02_predictions/$(replace(taxname, " " => "_"))_shapley.tif"; bandnumber=i) for i in eachindex(v)]
    mostimp = mosaic(argmax, map(x -> abs.(x), S))
    present = sort(unique(values(mostimp)))
    contrib = sum.(map(x -> abs.(x), S))
    contrib ./= sum(contrib)
    return (map=mostimp, vars=present, contrib=contrib, allvars=v)
end

function _abbreviator(x)
    x = replace(x, "Temperature" => "Temp.")
    x = replace(x, "Precipitation" => "Precip.")
    x = replace(x, "standard deviation" => "std. dev.")
    x = replace(x, "Coldest" => "Cold.")
    x = replace(x, "Warmest" => "Warm.")
    x = replace(x, "Wettest" => "Wett.")
    x = replace(x, "Quarter" => "Quart.")
    x = replace(x, " of " => " ")
    x = replace(x, r"\(.+\)" => "")
    x = replace(x, r"\s+" => " ")
    return x
end

Sp = _fromshap("Procyon lotor")
Sm = _fromshap("Mephitis mephitis")

_friendly_palette = [colorant"#0072B2", colorant"#56B4E9", colorant"#009E73", colorant"#F5C710", colorant"#E69F00", colorant"#D55E00"]

mIp = partialsortperm(Sp.contrib, 1:3, rev=true)
mIm = partialsortperm(Sm.contrib, 1:3, rev=true)

# Colors
Cp = _friendly_palette[1:2:5]
Cm = _friendly_palette[2:2:6]

# Get some figure action going
f = Figure(; size=(500, 400))

# Raccoon map
mp = Axis(f[1,1:2], aspect=DataAspect())
pal = fill(colorant"#cccccc", 12)
pal[mIp] .= Cp
heatmap!(mp, Sp.map, colormap=pal, colorrange=(1, 12))
lines!(mp, QC, color=:black)

# Skunk map
mm = Axis(f[2,1:2], aspect=DataAspect())
pal = fill(colorant"#cccccc", 12)
pal[mIm] .= Cm
heatmap!(mm, Sm.map, colormap=pal, colorrange=(1, 12))
lines!(mm, QC, color=:black)

# Axes
for ax in [mm, mp]
    hidedecorations!(ax)
    hidespines!(ax)
    tightlimits!(ax)
end

sp_uuid = Phylopic.imagesof("Procyon lotor"; items=1)
silhouetteplot!(mp, -61, 59, sp_uuid; markersize=42, label="P. lotor", color=:black)
sp_uuid = Phylopic.imagesof("Mephitis mephitis"; items=1)
silhouetteplot!(mm, -61, 59, sp_uuid; markersize=42, label="M. mephitis", color=:black)

# Variable importance

pstring = fill("WRONG VARIABLES", 3)#[descr[lnames[p]] for p in Sp.allvars[mIp]]
mstring = [descr[lnames[p]] for p in Sm.allvars[mIm]]
vstring = _abbreviator.(vcat(pstring, mstring))

Legend(
    f[1:2, 3],
    [PolyElement(; color=_friendly_palette[i]) for i in 1:length(_friendly_palette)],
    vstring;
    orientation=:horizontal,
    nbanks=6,
    framevisible=false,
    tellheight=false,
    tellwidth=true,
    vertical=false
)

# Show the figure
f
CairoMakie.save("figures/FINAL_shapmap.png", f)
