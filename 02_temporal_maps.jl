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

QC = SpeciesDistributionToolkit.gadm("CAN", "Québec")

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

scenarios = ["SSP126", "SSP245", "SSP370", "SSP585"]
timeframes = ["2021-2040", "2041-2060", "2061-2080", "2081-2100"]

for scenario in scenarios
    ranges = [get_both_sp(scenario, tf) for tf in timeframes]
    pushfirst!(ranges, get_baseline())

    # Now do the time expansion map
    function _findfirst(x)
        out = findfirst(replace(x, nothing => false))
        return isnothing(out) ? 0 : out
    end

    earliest = nodata(mosaic(_findfirst, ranges), 0)

    #earliest = trim(interpolate(earliest; dest="+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs"))

    var_colors = collect(cgrad(:navia, length(ranges)+1, categorical=true))
    deleteat!(var_colors, 2)

    time_labels = ["Current", "2030", "2050", "2070", "2090"]

    f = Figure(; size = (800, 700))
    ax = Axis(f[1, 1], aspect=DataAspect())
    heatmap!(ax, ranges[1], colormap=["#dfdfdf", "#dfdfdf"])
    hm = heatmap!(ax, earliest; colormap = var_colors, colorrange=(1, length(ranges)))
    hidedecorations!(ax)
    hidespines!(ax)
    tightlimits!(ax)
    lines!(ax, QC, color=:black)
    ax_inset = Axis(f[1, 1],
        width=Relative(0.38),
        height=Relative(0.15),
        halign=1.0,
        valign=0.0,
        xticks = (1:5, time_labels)
    )
    hits = [sum(earliest.==v) for v in 1:length(ranges)]
    barplot!(ax_inset, hits, color = var_colors, strokewidth=1)
    hidespines!(ax_inset, :t, :r, :l)
    hideydecorations!(ax_inset)
    hidexdecorations!(ax_inset, ticklabels=false)
    tightlimits!(ax_inset)
    current_figure()

    CairoMakie.save("both_$(scenario).png", current_figure())
end