using SpeciesDistributionToolkit
using CairoMakie
import DelimitedFiles
using PrettyTables
import CSV

set_theme!()
CairoMakie.activate!(; type="png")
update_theme!(;
    backgroundcolor=:transparent,
    fontsize=12,
    Figure=(; backgroundcolor=:transparent),
    Axis=(
        backgroundcolor=:transparent,
    ),
    CairoMakie=(; px_per_unit=6),
)

QC = SpeciesDistributionToolkit.gadm("CAN", "Québec")

v = Int64[DelimitedFiles.readdlm("data/Procyon_lotor.params")...]
S = [SDMLayer("rasters/Procyon_lotor_shapley.tif"; bandnumber=i) for i in eachindex(v)]

function shapscale(S)
    return maximum(abs.(quantile(S, [0.05, 0.95]))) .* (-1, 1)
end

f = Figure(; size=(800, 700))
ax = Axis(f[1, 1], aspect=DataAspect())
hm = heatmap!(ax, S[2], colorrange=shapscale(S[2]), colormap=:RdYlBu)
Colorbar(f[1, 1], hm, label="TK VAR NAME", alignmode=Inside(), width=Relative(0.4), valign=:top, halign=:right, tellheight=false, tellwidth=false, vertical=false)
lines!(QC, color=:black)
hidedecorations!(ax)
hidespines!(ax)
current_figure()

mostimp = mosaic(argmax, map(x -> abs.(x), S))
present = sort(unique(values(mostimp)))

dcols = CairoMakie.Colors.distinguishable_colors(length(present)+1)[2:end]

f = Figure(; size=(800, 700))
ax = Axis(f[1, 1], aspect=DataAspect())
heatmap!(ax, mostimp, colormap=dcols, colorrange=(1, length(present)))
lines!(QC, color=:black)
hidedecorations!(ax)
hidespines!(ax)
Legend(
    f[1, 1],
    [PolyElement(; color=dcols[i]) for i in 1:length(present)],
    ["V$(i)" for i in present];
    orientation=:horizontal,
    nbanks=3,
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