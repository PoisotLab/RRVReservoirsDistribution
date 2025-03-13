using CairoMakie
using Statistics
using SpeciesDistributionToolkit

include("S1_theme.jl")

QC = SpeciesDistributionToolkit.gadm("CAN", "Québec")

# Palette functions
include("S2_multivariate_palettes.jl")

scenarios = ["SSP126"]#, "SSP245", "SSP370", "SSP585"]
timeframes = ["2021-2040", "2041-2060", "2061-2080", "2081-2100"]

for scenario in scenarios, timeframe in timeframes
    # Layers
    raccoon = SDMLayer("rasters/02_predictions/$(scenario)/$(timeframe)/Procyon_lotor.tif", bandnumber=1)
    skunk = SDMLayer("rasters/02_predictions/$(scenario)/$(timeframe)/Mephitis_mephitis.tif", bandnumber=1)
    r_raccoon = SDMLayer("rasters/02_predictions/$(scenario)/$(timeframe)/Procyon_lotor.tif", bandnumber=2) .== 1
    r_skunk = SDMLayer("rasters/02_predictions/$(scenario)/$(timeframe)/Mephitis_mephitis.tif", bandnumber=2) .== 1

    # Out of range!
    nodata!(r_raccoon, false)
    nodata!(r_skunk, false)

    in_range = union(keys(r_raccoon), keys(r_skunk))
    range_mask = similar(r_raccoon)
    for k in in_range
        range_mask[k] = true
        range_mask.indices[k] = true
    end

    mask!(skunk, range_mask)
    mask!(raccoon, range_mask)

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
    smpal = _get_bivariate_colormap(; n_stops=5, colpal...)

    lics = LinearIndices(pal)
    bivlayer = similar(range_mask, Integer)
    for k in keys(range_mask)
        try
            bivlayer[k] = lics[rbiv[k], sbiv[k]]
        catch
            bivlayer.indices[k] = false
        end
    end

    f = Figure(; size = (800, 700))
    ax = Axis(f[1, 1], aspect=DataAspect())
    heatmap!(ax, bivlayer, colormap=vec(pal))
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
    heatmap!(ax_inset, smpal)
    hideydecorations!(ax_inset, label=false)
    hidexdecorations!(ax_inset, label=false)
    tightlimits!(ax_inset)
    current_figure()
    
    CairoMakie.save("figures/03b_maps_bivariate/bivariate_suitability_$(scenario)_$(timeframe).png", current_figure())
end