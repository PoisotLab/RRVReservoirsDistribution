using SpeciesDistributionToolkit
using CairoMakie
import DelimitedFiles
using PrettyTables
using Statistics
import CSV

include("S1_theme.jl")

taxname = "Mephitis mephitis"
if ~isempty(ARGS)
    taxname = join(ARGS[1:2], " ")
end
@info "Running for $(taxname)"

QC = SpeciesDistributionToolkit.gadm("CAN", "Québec")
provider = RasterData(WorldClim2, BioClim)

lnames = layers(provider)
descr = layerdescriptions(provider)

v = Int64[DelimitedFiles.readdlm("data/$(replace(taxname, " " => "_")).params")...]
S = [SDMLayer("rasters/$(replace(taxname, " " => "_"))_shapley.tif"; bandnumber=i) for i in eachindex(v)]

shapscale(S) = maximum(abs.(quantile(S, [0.05, 0.95]))) .* (-1, 1)

for (i, v) in enumerate(v)
    f = Figure(; size=(800, 700))
    ax = Axis(f[1, 1], aspect=DataAspect())
    hm = heatmap!(ax, S[2], colorrange=shapscale(S[i]), colormap=:RdYlBu)
    Colorbar(f[1, 1], hm, label=descr[lnames[v]], alignmode=Inside(), width=Relative(0.4), valign=:top, halign=:right, tellheight=false, tellwidth=false, vertical=false)
    lines!(QC, color=:black)
    hidedecorations!(ax)
    hidespines!(ax)
    current_figure()
    CairoMakie.save("paperfigures/$(replace(taxname, " " => "_"))_shapmap_BIO$(v).png", current_figure())
end

mostimp = mosaic(argmax, map(x -> abs.(x), S))
present = sort(unique(values(mostimp)))

dcols = CairoMakie.Colors.distinguishable_colors(length(present)+1)[2:end]

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

f = Figure(; size=(800, 700))
ax = Axis(f[1, 1], aspect=DataAspect())
heatmap!(ax, mostimp, colormap=dcols, colorrange=(1, length(present)))
lines!(QC, color=:black)
hidedecorations!(ax)
hidespines!(ax)
Legend(
    f[1, 1],
    [PolyElement(; color=dcols[i]) for i in 1:length(present)],
    _abbreviator.([descr[lnames[p]] for p in present]);
    orientation=:horizontal,
    nbanks=5,
    framevisible=false,
    alignmode=Inside(),
    width=Relative(0.4),
    valign=:top,
    halign=:right,
    tellheight=false,
    tellwidth=false,
    vertical=false
)
current_figure()
CairoMakie.save("paperfigures/$(replace(taxname, " " => "_"))_shapmap_mostimp.png", current_figure())