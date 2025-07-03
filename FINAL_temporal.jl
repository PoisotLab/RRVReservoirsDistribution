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
var_colors = collect(reverse(cgrad(:tempo, length(timeframes) + 2, categorical=true)))
deleteat!(var_colors, 2)
time_labels = ["Current", "2030", "2050", "2070", "2090"]

# Function to change axis
function _ranges_in_axis!(ax, f, i, j, earliest)
    poly!(ax, QC, color=:lightgrey)
    heatmap!(ax, earliest; colormap=var_colors, colorrange=(1, 5))

    ax_inset = Axis(f[i, j],
        width=Relative(0.38),
        height=Relative(0.15),
        halign=1.0,
        valign=0.0,
        xticks=(1:5, time_labels),
        xlabelrotation = pi/2
    )
    hits = [sum(earliest .== v) for v in 1:5]
    barplot!(ax_inset, hits, color=var_colors, strokewidth=1)
    hidespines!(ax_inset, :t, :r, :l)
    hideydecorations!(ax_inset)
    hidexdecorations!(ax_inset, ticklabels=false)
    tightlimits!(ax_inset)
    hidedecorations!(ax)
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
_ranges_in_axis!(ar_ssp1, f, 1, 1, earliest)

ar_ssp2 = Axis(f[2, 1]; aspect=DataAspect(), ylabel="SSP245")
r = [get_one_sp("SSP245", tf; sp="Procyon_lotor") for tf in timeframes]
pushfirst!(r, get_one_baseline(; sp="Procyon_lotor"))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(ar_ssp2, f, 2, 1, earliest)

ar_ssp3 = Axis(f[3, 1]; aspect=DataAspect(), ylabel="SSP370")
r = [get_one_sp("SSP370", tf; sp="Procyon_lotor") for tf in timeframes]
pushfirst!(r, get_one_baseline(; sp="Procyon_lotor"))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(ar_ssp3, f, 3, 1, earliest)

ar_ssp5 = Axis(f[4, 1]; aspect=DataAspect(), ylabel="SSP585")
r = [get_one_sp("SSP585", tf; sp="Procyon_lotor") for tf in timeframes]
pushfirst!(r, get_one_baseline(; sp="Procyon_lotor"))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(ar_ssp5, f, 4, 1, earliest)

# Skunk
am_ssp1 = Axis(f[1, 2]; aspect=DataAspect())
r = [get_one_sp("SSP126", tf; sp="Mephitis_mephitis") for tf in timeframes]
pushfirst!(r, get_one_baseline(; sp="Mephitis_mephitis"))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(am_ssp1, f, 1, 2, earliest)

am_ssp2 = Axis(f[2, 2]; aspect=DataAspect())
r = [get_one_sp("SSP245", tf; sp="Mephitis_mephitis") for tf in timeframes]
pushfirst!(r, get_one_baseline(; sp="Mephitis_mephitis"))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(am_ssp2, f, 2, 2, earliest)

am_ssp3 = Axis(f[3, 2]; aspect=DataAspect())
r = [get_one_sp("SSP370", tf; sp="Mephitis_mephitis") for tf in timeframes]
pushfirst!(r, get_one_baseline(; sp="Mephitis_mephitis"))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(am_ssp3, f, 3, 2, earliest)

am_ssp5 = Axis(f[4, 2]; aspect=DataAspect())
r = [get_one_sp("SSP585", tf; sp="Mephitis_mephitis") for tf in timeframes]
pushfirst!(r, get_one_baseline(; sp="Mephitis_mephitis"))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(am_ssp5, f, 4, 2, earliest)

# Either
ab_ssp1 = Axis(f[1, 3]; aspect=DataAspect())
r = [get_both_sp("SSP126", tf; any=true) for tf in timeframes]
pushfirst!(r, get_baseline(; any=true))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(ab_ssp1, f, 1, 3, earliest)

ab_ssp2 = Axis(f[2, 3]; aspect=DataAspect())
r = [get_both_sp("SSP245", tf; any=true) for tf in timeframes]
pushfirst!(r, get_baseline(; any=true))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(ab_ssp2, f, 2, 3, earliest)

ab_ssp3 = Axis(f[3, 3]; aspect=DataAspect())
r = [get_both_sp("SSP370", tf; any=true) for tf in timeframes]
pushfirst!(r, get_baseline(; any=true))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(ab_ssp3, f, 3, 3, earliest)

ab_ssp5 = Axis(f[4, 3]; aspect=DataAspect())
r = [get_both_sp("SSP585", tf; any=true) for tf in timeframes]
pushfirst!(r, get_baseline(; any=true))
earliest = nodata(mosaic(_findfirst, r), 0)
_ranges_in_axis!(ab_ssp5, f, 4, 3, earliest)

sp_uuid = Phylopic.imagesof("Procyon lotor"; items=1)
silhouetteplot!(ar_ssp1, -62, 59, sp_uuid; markersize=40, label="P. lotor", color=:black)
sp_uuid = Phylopic.imagesof("Mephitis mephitis"; items=1)
silhouetteplot!(am_ssp1, -62, 59, sp_uuid; markersize=40, label="M. mephitis", color=:black)

f
CairoMakie.save("figures/FINAL_temporal.png", f)