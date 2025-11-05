using CairoMakie
using Statistics
using SpeciesDistributionToolkit

include("S1_theme.jl")

QC = getpolygon(PolygonData(GADM, Countries), level=1, country="CAN")["Québec"]

function _findfirst(x)
    out = findfirst(replace(x, nothing => false))
    return isnothing(out) ? 0 : out
end

# Load the current maps

function get_baseline(; any=true)
    raccoon = SDMLayer("rasters/02_predictions/Procyon_lotor.tif", bandnumber=2)
    skunk = SDMLayer("rasters/02_predictions/Mephitis_mephitis.tif", bandnumber=2)
    return (raccoon + skunk) .>= (any ? 1 : 2)
end

function get_one_baseline(; sp="Procyon_lotor")
    return SDMLayer("rasters/02_predictions/$(sp).tif", bandnumber=2) .> 0.0
end

function get_both_sp(scenario, timeframe; any=true)
    raccoon = SDMLayer("rasters/02_predictions/$(scenario)/$(timeframe)/Procyon_lotor.tif", bandnumber=2)
    skunk = SDMLayer("rasters/02_predictions/$(scenario)/$(timeframe)/Mephitis_mephitis.tif", bandnumber=2)
    return (raccoon + skunk) .>= (any ? 1 : 2)
end

function get_one_sp(scenario, timeframe; sp="Procyon_lotor")
    return SDMLayer("rasters/02_predictions/$(scenario)/$(timeframe)/$(sp).tif", bandnumber=2) .> 0.0
end

scenarios = ["SSP126", "SSP245", "SSP370", "SSP585"]
timeframes = ["2021-2040", "2041-2060", "2061-2080", "2081-2100"]

# Legend info
function create_palette(pal, n)
    var_colors = collect(reverse(cgrad(pal, n + 3, categorical=true)))
    deleteat!(var_colors, 2)
    deleteat!(var_colors, length(var_colors))
    return var_colors
end

pl_colors = create_palette(:Blues, length(timeframes))
mm_colors = create_palette(:Reds, length(timeframes))
bt_colors = create_palette(:Purples, length(timeframes))

time_labels = ["", "30s", "50s", "70s", "90s"]

# Function to change axis
function _ranges_in_axis!(ax, f, i, j, earliest; colormap=var_colors)
    poly!(ax, QC, color=:lightgrey)
    heatmap!(ax, earliest; colormap=colormap, colorrange=(1, 5))

    ax_inset = Axis(f[i, j],
        width=Relative(0.38),
        height=Relative(0.12),
        halign=1.0,
        valign=0.0,
        xticks=(1:5, time_labels),
        xlabelrotation = pi/3
    )
    hits = [sum(earliest .== v) for v in 1:5]
    barplot!(ax_inset, hits, color=colormap, strokewidth=1)
    hidespines!(ax_inset, :t, :r, :l)
    hideydecorations!(ax_inset)
    hidexdecorations!(ax_inset, ticklabels=false)
    tightlimits!(ax_inset)
    hidedecorations!(ax, label=false)
    hidespines!(ax)
    tightlimits!(ax)
    lines!(ax, QC, color=:black)
end

# The big ol figure

f = Figure(; size=(1000, 1000))

# Raccoon
ar_ssp1 = Axis(f[1, 1]; aspect=DataAspect(), ylabel="SSP126")
r = [get_one_sp("SSP126", tf; sp="Procyon_lotor") for tf in timeframes]
pushfirst!(r, get_one_baseline(; sp="Procyon_lotor"))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(ar_ssp1, f, 1, 1, earliest; colormap=pl_colors)

ar_ssp2 = Axis(f[2, 1]; aspect=DataAspect(), ylabel="SSP245")
r = [get_one_sp("SSP245", tf; sp="Procyon_lotor") for tf in timeframes]
pushfirst!(r, get_one_baseline(; sp="Procyon_lotor"))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(ar_ssp2, f, 2, 1, earliest; colormap=pl_colors)

ar_ssp3 = Axis(f[3, 1]; aspect=DataAspect(), ylabel="SSP370")
r = [get_one_sp("SSP370", tf; sp="Procyon_lotor") for tf in timeframes]
pushfirst!(r, get_one_baseline(; sp="Procyon_lotor"))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(ar_ssp3, f, 3, 1, earliest; colormap=pl_colors)

ar_ssp5 = Axis(f[4, 1]; aspect=DataAspect(), ylabel="SSP585")
r = [get_one_sp("SSP585", tf; sp="Procyon_lotor") for tf in timeframes]
pushfirst!(r, get_one_baseline(; sp="Procyon_lotor"))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(ar_ssp5, f, 4, 1, earliest; colormap=pl_colors)

# Skunk
am_ssp1 = Axis(f[1, 2]; aspect=DataAspect())
r = [get_one_sp("SSP126", tf; sp="Mephitis_mephitis") for tf in timeframes]
pushfirst!(r, get_one_baseline(; sp="Mephitis_mephitis"))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(am_ssp1, f, 1, 2, earliest; colormap=mm_colors)

am_ssp2 = Axis(f[2, 2]; aspect=DataAspect())
r = [get_one_sp("SSP245", tf; sp="Mephitis_mephitis") for tf in timeframes]
pushfirst!(r, get_one_baseline(; sp="Mephitis_mephitis"))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(am_ssp2, f, 2, 2, earliest; colormap=mm_colors)

am_ssp3 = Axis(f[3, 2]; aspect=DataAspect())
r = [get_one_sp("SSP370", tf; sp="Mephitis_mephitis") for tf in timeframes]
pushfirst!(r, get_one_baseline(; sp="Mephitis_mephitis"))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(am_ssp3, f, 3, 2, earliest; colormap=mm_colors)

am_ssp5 = Axis(f[4, 2]; aspect=DataAspect())
r = [get_one_sp("SSP585", tf; sp="Mephitis_mephitis") for tf in timeframes]
pushfirst!(r, get_one_baseline(; sp="Mephitis_mephitis"))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(am_ssp5, f, 4, 2, earliest; colormap=mm_colors)

# Either
ab_ssp1 = Axis(f[1, 3]; aspect=DataAspect())
r = [get_both_sp("SSP126", tf; any=true) for tf in timeframes]
pushfirst!(r, get_baseline(; any=true))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(ab_ssp1, f, 1, 3, earliest; colormap=bt_colors)

ab_ssp2 = Axis(f[2, 3]; aspect=DataAspect())
r = [get_both_sp("SSP245", tf; any=true) for tf in timeframes]
pushfirst!(r, get_baseline(; any=true))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(ab_ssp2, f, 2, 3, earliest; colormap=bt_colors)

ab_ssp3 = Axis(f[3, 3]; aspect=DataAspect())
r = [get_both_sp("SSP370", tf; any=true) for tf in timeframes]
pushfirst!(r, get_baseline(; any=true))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(ab_ssp3, f, 3, 3, earliest; colormap=bt_colors)

ab_ssp5 = Axis(f[4, 3]; aspect=DataAspect())
r = [get_both_sp("SSP585", tf; any=true) for tf in timeframes]
pushfirst!(r, get_baseline(; any=true))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(ab_ssp5, f, 4, 3, earliest; colormap=bt_colors)

sp_uuid = Phylopic.imagesof("Procyon lotor"; items=1)
silhouetteplot!(ar_ssp1, -62, 59, sp_uuid; markersize=40, label="P. lotor", color=:black)
sp_uuid = Phylopic.imagesof("Mephitis mephitis"; items=1)
silhouetteplot!(am_ssp1, -62, 59, sp_uuid; markersize=40, label="M. mephitis", color=:black)

sp_uuid = Phylopic.imagesof("Procyon lotor"; items=1)
silhouetteplot!(ab_ssp1, -62, 62, sp_uuid; markersize=40, label="P. lotor", color=:black)
sp_uuid = Phylopic.imagesof("Mephitis mephitis"; items=1)
silhouetteplot!(ab_ssp1, -60, 59, sp_uuid; markersize=40, label="M. mephitis", color=:black)

display(f)
CairoMakie.save("figures/FINAL_temporal.png", f)