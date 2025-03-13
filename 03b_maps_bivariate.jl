using CairoMakie
using Statistics
using SpeciesDistributionToolkit

include("S1_theme.jl")

QC = SpeciesDistributionToolkit.gadm("CAN", "Québec")

# Palette functions
include("S2_multivariate_palettes.jl")

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

    # Prepare the bivariate map
    nbins = 10

    rbiv = ((x) -> round(Int64, x * (nbins-1)) + 1).(raccoon)
    sbiv = ((x) -> round(Int64, x * (nbins-1)) + 1).(skunk)

    stevens1 = (p0 = colorant"#e8e8e8", p1 = colorant"#c85a5a", p2 = colorant"#64acbe")
    stevens2 = (p0 = colorant"#e8e8e8", p1 = colorant"#5ac8c8", p2 = colorant"#be64ac")
    stevens3 = (p0 = colorant"#e8e8e8", p1 = colorant"#6c83b5", p2 = colorant"#73ae80")
    stevens4 = (p0 = colorant"#e8e8e8", p1 = colorant"#c8b35a", p2 = colorant"#9972af")

    colpal = stevens3

    pal = _get_bivariate_colormap(; n_stops=nbins, colpal...)
    smpal = _get_bivariate_colormap(; n_stops=3, colpal...)

    lics = LinearIndices(pal)
    bivlayer = similar(range_mask, Integer)
    for k in keys(range_mask)
        try
            bivlayer[k] = lics[rbiv[k], sbiv[k]]
        catch
            bivlayer.indices[k] = false
        end
    end

    # Mask for either species
    r_either = mosaic(any, [r_raccoon, r_skunk])
    nodata!(r_either, false)

    f = Figure(; size = (800, 700))
    ax = Axis(f[1, 1], aspect=DataAspect())
    poly!(ax, QC, strokecolor=:black, color=colpal.p0, strokewidth=1)
    heatmap!(ax, mask(bivlayer, r_either), colormap=[pal[l] for l in vec(lics)])
    hidedecorations!(ax)
    hidespines!(ax)
    tightlimits!(ax)
    lines!(ax, QC, color=:black)
    ax_inset = Axis(f[1, 1],
        aspect = 1,
        width=Relative(0.25),
        height=Relative(0.25),
        halign=1.0,
        valign=1.0,
        xlabel = "Raccoon",
        ylabel = "Skunk",
    )
    heatmap!(ax_inset, 1:3, 1:3, smpal)
    scatter!(ax_inset, rescale(raccoon, ax_inset.xaxis.attributes.limits.val...), rescale(skunk, ax_inset.yaxis.attributes.limits.val...), color=:black, markersize=1, alpha=0.15)
    hideydecorations!(ax_inset, label=false)
    hidexdecorations!(ax_inset, label=false)
    tightlimits!(ax_inset)
    current_figure()
    
    CairoMakie.save("figures/03b_maps_bivariate/bivariate_suitability_$(scenario)_$(timeframe).png", current_figure())
end