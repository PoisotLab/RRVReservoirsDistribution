using CairoMakie
using Statistics
using SpeciesDistributionToolkit

include("S1_theme.jl")

QC = getpolygon(PolygonData(GADM, Countries); country="CAN", level=1)["Québec"]

shapscale(S) = maximum(abs.(quantile(S, [0.05, 0.95]))) .* (-1, 1)
shapscale(S) = (-0.5, 0.5) # This is useful for the final version of the figures (same range for all sub-plots)


function polarization(scenario, timeframe)
    raccoon = SDMLayer("rasters/02_predictions/$(scenario)/$(timeframe)/Procyon_lotor.tif", bandnumber=1)
    skunk = SDMLayer("rasters/02_predictions/$(scenario)/$(timeframe)/Mephitis_mephitis.tif", bandnumber=1)
    r_raccoon = SDMLayer("rasters/02_predictions/$(scenario)/$(timeframe)/Procyon_lotor.tif", bandnumber=2) .== 1
    r_skunk = SDMLayer("rasters/02_predictions/$(scenario)/$(timeframe)/Mephitis_mephitis.tif", bandnumber=2) .== 1

    in_range = union(keys(raccoon), keys(skunk))
    range_mask = similar(raccoon)
    for k in in_range
        range_mask[k] = true
        range_mask.indices[k] = true
    end

    # Mask for either species
    r_either = mosaic(any, [r_raccoon, r_skunk])
    nodata!(r_either, false)

    NDRI = (raccoon - skunk) / (raccoon + skunk)
    return raccoon, skunk, mask(NDRI, r_either), NDRI
end

scenarios = ["SSP126", "SSP245", "SSP370", "SSP585"]
timeframes = ["2021-2040", "2041-2060", "2061-2080", "2081-2100"]

f = Figure(; size=(1000, 1000))
axs = [Axis(f[i, j]; aspect=DataAspect(), title=scenarios[i], subtitle=timeframes[j]) for i in eachindex(scenarios), j in eachindex(timeframes)]
f

for (i, scenario) in enumerate(scenarios)
    for (j, timeframe) in enumerate(timeframes)
        raccoon, skunk, NDRI, fullval = polarization(scenario, timeframe)
        poly!(axs[i, j], QC, strokecolor=:black, strokewidth=1, color="#dfdfdf")
        hm = heatmap!(axs[i, j], NDRI, colormap=:diverging_bwg_20_95_c41_n256, colorrange=shapscale(NDRI))
        hidedecorations!(axs[i, j])
        hidespines!(axs[i, j])
        tightlimits!(axs[i, j])
        lines!(axs[i, j], QC, color=:black)
        Colorbar(f[i, j], hm, alignmode=Inside(), width=Relative(0.4), valign=:bottom, halign=:right, tellheight=false, tellwidth=false, vertical=false, flipaxis=false)
        ax_inset = Axis(f[i, j],
            aspect=1,
            width=Relative(0.23),
            height=Relative(0.23),
            halign=1.0,
            valign=1.0,
            xlabel="Raccoon",
            xlabelsize=12,
            ylabelsize=12,
            ylabel="Skunk",
        )
        lines!(ax_inset, [0.0, 1.0], [0.0, 1.0], color=:grey60, linestyle=:dash)
        scatter!(ax_inset, raccoon, skunk, color=values(fullval), markersize=1, alpha=0.15, colormap=:diverging_bwg_20_95_c41_n256, colorrange=shapscale(NDRI))
        hideydecorations!(ax_inset, label=false)
        hidexdecorations!(ax_inset, label=false)
        tightlimits!(ax_inset)
        current_figure()
    end
end

f
CairoMakie.save("figures/FINAL_polarization.png", f)