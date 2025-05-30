using CairoMakie
using Statistics
using SpeciesDistributionToolkit

include("S1_theme.jl")

QC = SpeciesDistributionToolkit.gadm("CAN", "Québec")

# Load the current maps

function get_baseline(; any=true)
    raccoon = SDMLayer("rasters/02_predictions/Procyon_lotor.tif", bandnumber=2)
    skunk = SDMLayer("rasters/02_predictions/Mephitis_mephitis.tif", bandnumber=2)
    return (raccoon + skunk) .>= (any ? 1 : 2)
end

function get_both_sp(scenario, timeframe; any=true)
    raccoon = SDMLayer("rasters/02_predictions/$(scenario)/$(timeframe)/Procyon_lotor.tif", bandnumber=2)
    skunk = SDMLayer("rasters/02_predictions/$(scenario)/$(timeframe)/Mephitis_mephitis.tif", bandnumber=2)
    return (raccoon + skunk) .>= (any ? 1 : 2)
end

scenarios = ["SSP126", "SSP245", "SSP370", "SSP585"]
timeframes = ["2021-2040", "2041-2060", "2061-2080", "2081-2100"]

for scenario in scenarios
    for any_sp in [true, false]
        ranges = [get_both_sp(scenario, tf; any=any_sp) for tf in timeframes]
        pushfirst!(ranges, get_baseline(; any=any_sp))

        # Now do the time expansion map
        function _findfirst(x)
            out = findfirst(replace(x, nothing => false))
            return isnothing(out) ? 0 : out
        end

        earliest = nodata(mosaic(_findfirst, ranges), 0)

        var_colors = collect(reverse(cgrad(:tempo, length(ranges) + 1, categorical=true)))
        deleteat!(var_colors, 2)

        time_labels = ["Current", "2030", "2050", "2070", "2090"]

        f = Figure(; size=(800, 700))
        ax = Axis(f[1, 1], aspect=DataAspect())
        heatmap!(ax, ranges[1], colormap=["#dfdfdf", "#dfdfdf"])
        hm = heatmap!(ax, earliest; colormap=var_colors, colorrange=(1, length(ranges)))
        hidedecorations!(ax)
        hidespines!(ax)
        tightlimits!(ax)
        lines!(ax, QC, color=:black)
        ax_inset = Axis(f[1, 1],
            width=Relative(0.38),
            height=Relative(0.15),
            halign=1.0,
            valign=0.0,
            xticks=(1:5, time_labels)
        )
        hits = [sum(earliest .== v) for v in 1:length(ranges)]
        barplot!(ax_inset, hits, color=var_colors, strokewidth=1)
        hidespines!(ax_inset, :t, :r, :l)
        hideydecorations!(ax_inset)
        hidexdecorations!(ax_inset, ticklabels=false)
        tightlimits!(ax_inset)
        current_figure()

        type = any_sp ? "either" : "both"

        CairoMakie.save("figures/03a_maps_temporal/arrival_date_$(type)_$(scenario).png", current_figure())
    end
end