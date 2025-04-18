using CairoMakie
using Statistics
using SpeciesDistributionToolkit

include("S1_theme.jl")

QC = SpeciesDistributionToolkit.gadm("CAN", "Québec")

# Palette functions
include("S2_multivariate_palettes.jl")

shapscale(S) = maximum(abs.(quantile(S, [0.05, 0.95]))) .* (-1, 1)

scenarios = ["SSP126", "SSP245", "SSP370", "SSP585"]
timeframes = ["2021-2040", "2041-2060", "2061-2080", "2081-2100"]

for scenario in scenarios, timeframe in timeframes
    # Layers
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

    f = Figure(; size = (800, 700))
    ax = Axis(f[1, 1], aspect=DataAspect())
    poly!(ax, QC, strokecolor=:black, strokewidth=1, color="#dfdfdf")
    hm = heatmap!(ax, mask(NDRI, r_either), colormap=:balance, colorrange=shapscale(mask(NDRI, r_either)))
    hidedecorations!(ax)
    hidespines!(ax)
    tightlimits!(ax)
    Colorbar(f[1, 1], hm, label="NDRI", alignmode=Inside(), width=Relative(0.4), valign=:bottom, halign=:right, tellheight=false, tellwidth=false, vertical=false)
    lines!(ax, QC, color=:black)
    ax_inset = Axis(f[1, 1],
        aspect = 1,
        width=Relative(0.27),
        height=Relative(0.27),
        halign=1.0,
        valign=1.0,
        xlabel = "Raccoon",
        xlabelsize = 18,
        ylabelsize = 18,
        ylabel = "Skunk",
    )
    lines!(ax_inset, [0.0, 1.0], [0.0, 1.0], color=:grey60, linestyle=:dash)
    scatter!(ax_inset, raccoon, skunk, color=values(NDRI), markersize=1, alpha=0.15, colormap=:balance, colorrange=shapscale(mask(NDRI, r_either)))
    hideydecorations!(ax_inset, label=false)
    hidexdecorations!(ax_inset, label=false)
    tightlimits!(ax_inset)
    current_figure()
    
    CairoMakie.save("figures/03b_maps_bivariate/bivariate_suitability_$(scenario)_$(timeframe).png", current_figure())
end