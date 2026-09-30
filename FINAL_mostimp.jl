using SpeciesDistributionToolkit
using CairoMakie
import DelimitedFiles
using PrettyTables
using Statistics
import CSV

include("S1_theme.jl")

mm_sdm = SDeMo.loadsdm("models/Mephitis_mephitis.json")
mm_bar = mean(predict(mm_sdm; threshold=false))

pl_sdm = SDeMo.loadsdm("models/Procyon_lotor.json")
pl_bar = mean(predict(pl_sdm; threshold=false))

QC = getpolygon(PolygonData(GADM, Countries), level=1, country="CAN")["Québec"]
provider = RasterData(WorldClim2, BioClim)

lnames = layers(provider)
descr = layerdescriptions(provider)

function _fromshap(taxname)
  # This is the list of TRUE variables (i.e. BIOCLIM codes)
  v = Int64[DelimitedFiles.readdlm("data/$(replace(taxname, " " => "_")).params")...]
  # This is a map of INDICES, i.e. indices in the vector V, it goes from 1 to the number of variables
  S = [SDMLayer("rasters/02_predictions/$(replace(taxname, " " => "_"))_shapley.tif"; bandnumber=i) for i in eachindex(v)]
  mostimp = mosaic(argmax, map(x -> abs.(x), S))
  contrib = sum.(map(x -> abs.(x), S))
  contrib ./= sum(contrib)
  return (map=mostimp, contrib=contrib, variables=v)
end

function _abbreviator(x)
  x = replace(x, "Temperature" => "Temp.")
  x = replace(x, "Precipitation" => "Precip.")
  x = replace(x, "standard deviation" => "std. dev.")
  x = replace(x, "Coldest" => "Cold.")
  x = replace(x, "Warmest" => "Warm.")
  x = replace(x, "Wettest" => "Wett.")
  x = replace(x, "Quarter" => "Quart.")
  x = replace(x, " of " => " ")
  x = replace(x, r"\(.+\)" => "")
  x = replace(x, r"\s+" => " ")
  return x
end

Sp = _fromshap("Procyon lotor")
Sm = _fromshap("Mephitis mephitis")


# This gives the index in XX.variables for the most important variables
mIp = partialsortperm(Sp.contrib, 1:3, rev=true)
mIm = partialsortperm(Sm.contrib, 1:3, rev=true)

# Indices of the true most important variables
truvp = Sp.variables[mIp]
truvm = Sm.variables[mIm]
truv = unique(vcat(truvp, truvm))

# Colors
_friendly_palette = [
  colorant"#E69F00",
  colorant"#56B4E9",
  colorant"#009e73",
  colorant"#F0E442",
  colorant"#000000",
  colorant"#CC79A7"
]
cpal = vcat(_friendly_palette)

# Get some figure action going
f = Figure(; size=(500, 400))

# Raccoon map
mp = Axis(f[1:2, 1], aspect=DataAspect())
poly!(mp, QC, color=:lightgrey)
for v in mIp
  true_var = Sp.variables[v]
  idx = first(indexin(true_var, truv))
  heatmap!(mp, nodata(Sp.map, x -> x != v), colormap=[cpal[idx]])
end
lines!(mp, QC, color=:black, linewidth=0.5)
f

# Skunk map
mm = Axis(f[3:4, 1], aspect=DataAspect())
poly!(mm, QC, color=:lightgrey)
for v in mIm
  true_var = Sm.variables[v]
  idx = first(indexin(true_var, truv))
  heatmap!(mm, nodata(Sm.map, x -> x != v), colormap=[cpal[idx]])
end
lines!(mm, QC, color=:black, linewidth=0.5)
f

# Axes
for ax in [mm, mp]
  hidedecorations!(ax)
  hidespines!(ax)
  tightlimits!(ax)
end

sp_uuid = last(Phylopic.imagesof("Procyon lotor"; items=2))
silhouetteplot!(mp, -62, 59, sp_uuid; markersize=32, label="P. lotor", color=:black)
sp_uuid = last(Phylopic.imagesof("Mephitis mephitis"; items=2))
silhouetteplot!(mm, -62, 59, sp_uuid; markersize=32, label="M. mephitis", color=:black)

# Variable importance
pstring = [descr[lnames[p]] for p in truv]
vstring = _abbreviator.(pstring)

Legend(
  f[5, 1:3],
  [PolyElement(; color=c) for c in cpal[1:length(truv)]],
  vstring;
  orientation=:horizontal,
  nbanks=2,
  framevisible=false,
  tellheight=false,
  tellwidth=false,
  vertical=false
)
f

# Show the Shapley variables for one variable
vi = 5

label = descr[lnames[vi]]

pl_sh = Axis(f[1:2, 2:3], xlabel=label, ylabel="Effect on prediction")
scatter!(pl_sh, features(pl_sdm, vi), pl_bar .+ explain(pl_sdm, vi; threshold=false), color=:grey80)
lines!(pl_sh, partialresponse(pl_sdm, vi; threshold=false)..., color=:black, linewidth=2)
current_figure()

mm_sh = Axis(f[3:4, 2:3], xlabel=label, ylabel="Effect on prediction")
scatter!(mm_sh, features(mm_sdm, vi), mm_bar .+ explain(mm_sdm, vi; threshold=false), color=:grey80)
lines!(mm_sh, partialresponse(mm_sdm, vi; threshold=false)..., color=:black, linewidth=2)
current_figure()

# Show the figure
display(f)
CairoMakie.save("figures/FINAL_shapmap.png", f)

# Area
area = cellarea(Sp.map)

# Get the data for the table
for i in 1:19
  # i is the true bioclim code
  strings = [descr[lnames[i]]]
  if i in Sp.variables
    v_idx = findfirst(isequal(i), Sp.variables)
    push!(strings, string(round(Sp.contrib[v_idx]; digits=2)))
    v_area = sum(mask(area, nodata(Sp.map .== v_idx, false))) * 1e-3
    push!(strings, string(round(v_area; digits=2)))
  else
    append!(strings, ["", ""])
  end
  if i in Sm.variables
    v_idx = findfirst(isequal(i), Sm.variables)
    v_area = sum(mask(area, nodata(Sm.map .== v_idx, false))) * 1e-3
    push!(strings, string(round(Sm.contrib[v_idx]; digits=2)))
    push!(strings, string(round(v_area; digits=2)))
  else
    append!(strings, ["", ""])
  end
  @info strings
end
