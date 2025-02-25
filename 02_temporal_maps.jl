using CairoMakie
using Statistics
using SpeciesDistributionToolkit

set_theme!()
CairoMakie.activate!(; type = "png")
update_theme!(;
    backgroundcolor = :transparent,
    fontsize = 12,
    Figure = (; backgroundcolor = :transparent),
    Axis = (
        backgroundcolor = :transparent,
    ),
    CairoMakie = (; px_per_unit = 6),
)

# Load the current maps

function get_baseline()
    raccoon = SDMLayer("rasters/Procyon_lotor_historical.tif", bandnumber=2) .== 1.0
    skunk = SDMLayer("rasters/Mephitis_mephitis_historical.tif", bandnumber=2) .== 1.0
    return raccoon | skunk
end

function get_both_sp(scenario, timeframe)
    raccoon = SDMLayer("rasters/Procyon_lotor_$(scenario)_$(timeframe).tif", bandnumber=2) .== 1.0
    skunk = SDMLayer("rasters/Mephitis_mephitis_$(scenario)_$(timeframe).tif", bandnumber=2) .== 1.0
    return raccoon | skunk
end

scenario = "SSP585"
timeframes = ["2021-2040", "2041-2060", "2061-2080", "2081-2100"]

ranges = [get_both_sp(scenario, tf) for tf in timeframes]
pushfirst!(ranges, get_baseline())

# Now do the time expansion map
function _findfirst(x)
    out = findfirst(replace(x, nothing => false))
    return isnothing(out) ? 0 : out
end

earliest = nodata(mosaic(_findfirst, ranges), 0)

var_colors = collect(cgrad(:lapaz, length(ranges)+1, categorical=true))
deleteat!(var_colors, 2)

f = Figure(; size = (800, 700))
ax = Axis(f[1, 1], aspect=DataAspect())
heatmap!(ax, ranges[1], colormap=["#dfdfdf", "#dfdfdf"])
hm = heatmap!(ax, earliest; colormap = var_colors, colorrange=(1, length(ranges)))
hidedecorations!(ax)
hidespines!(ax)
Legend(
    f[2, 1],
    [PolyElement(; color = var_colors[i]) for i in 1:length(var_colors)],
    ["Current", "2030", "2050", "2070", "2090"];
    orientation = :horizontal,
    nbanks = 1,
)
current_figure()