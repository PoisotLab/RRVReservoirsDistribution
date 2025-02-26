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

scenario = "SSP585"
timeframe = "2081-2100"

# Layers
raccoon = SDMLayer("rasters/Procyon_lotor_$(scenario)_$(timeframe).tif", bandnumber=1)
skunk = SDMLayer("rasters/Mephitis_mephitis_$(scenario)_$(timeframe).tif", bandnumber=1)
r_raccoon = SDMLayer("rasters/Procyon_lotor_$(scenario)_$(timeframe).tif", bandnumber=2) .== 1
r_skunk = SDMLayer("rasters/Mephitis_mephitis_$(scenario)_$(timeframe).tif", bandnumber=2) .== 1

#mask!(raccoon, nodata(r_raccoon, false))
#mask!(skunk, nodata(r_skunk, false))

# Palette functions
include("multivariate-palettes.jl")

# Prepare the bivariate map
nbins = 3

rbiv = ((x) -> round(Int64, x * (nbins-1)) + 1).(raccoon)
sbiv = ((x) -> round(Int64, x * (nbins-1)) + 1).(skunk)

pal = _get_bivariate_colormap(n_stops=nbins, p1 = colorant"#73ae80", p2 = colorant"#6c83b5")
smpal = _get_bivariate_colormap(n_stops=3, p1 = colorant"#73ae80", p2 = colorant"#6c83b5")

lics = LinearIndices(pal)
bivlayer = similar(rbiv)
for k in keys(bivlayer)
    bivlayer[k] = lics[rbiv[k], sbiv[k]]
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
