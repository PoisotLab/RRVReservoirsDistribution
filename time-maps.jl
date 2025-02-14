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

sp = "Procyon lotor"
raster_path = joinpath(pwd(), "rasters", sp)
all_rasters = readdir(raster_path, join=true)
filter!(contains("range"), all_rasters)

function _read_to_bool(f)
    return convert(SDMLayer{Bool}, SDMLayer(f).>0)
end

current_range = _read_to_bool(all_rasters[findfirst(contains("current"), all_rasters)])
future_ranges = _read_to_bool.(all_rasters[findall(contains("SSP126"), all_rasters)])

pushfirst!(future_ranges, current_range)

function _findfirst(x)
    out = findfirst(replace(x, nothing => false))
    return isnothing(out) ? 0 : out
end

earliest = nodata(mosaic(_findfirst, future_ranges), 0)

var_colors = collect(cgrad(:lipari, length(future_ranges)+1, categorical=true))
deleteat!(var_colors, 2)

f = Figure(; size = (800, 700))
ax = Axis(f[1, 1])
heatmap!(ax, current_range, colormap=["#dfdfdf", "#dfdfdf"])
hm = heatmap!(ax, earliest; colormap = var_colors, colorrange=(1, length(future_ranges)))
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