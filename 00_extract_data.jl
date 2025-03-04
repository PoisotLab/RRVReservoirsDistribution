using SpeciesDistributionToolkit
using CairoMakie
import DelimitedFiles
import CSV

# Create the folders if they don't exist
for folder in ["models", "data", "rasters", "figures", "paperfigures"]
    if !ispath(folder)
        mkpath(folder)
    end
end

taxname = "Mephitis mephitis"
if ~isempty(ARGS)
    taxname = join(ARGS[1:2], " ")
end
@info "Running for $(taxname)"

@info "Loading the CSV file"
records = CSV.File("0022970-241126133413365.csv")
records = filter(r -> isequal(taxname)(r.species), records)

@info "Turning records in occurrences"
occ = Occurrences([Occurrence(r.species, true, (r.decimalLongitude, r.decimalLatitude), r.dateIdentified) for r in records])
bbox = SpeciesDistributionToolkit.boundingbox(occ; padding=2.0)

@info "Loading bioclim data for training"
provider = RasterData(WorldClim2, BioClim)
envirovars = [SDMLayer(provider; layer=i, bbox...) for i in eachindex(layers(provider))]

@info "Thinning the occurrences to the grid"
presencelayer = mask(envirovars[1], occ)

@info "Generating pseudo-absences mask"
event_dist = pseudoabsencemask(DistanceToEvent, presencelayer)
pa_mask = copy(event_dist)
nodata!(pa_mask, x -> x <= 20.0)
nodata!(pa_mask, x -> x >= 500.0)

@info "Sampling pseudo-absences"
absencelayer = backgroundpoints(pa_mask, 3sum(presencelayer))

@info "nodata for the points"
nodata!(absencelayer, false)
nodata!(presencelayer, false)

@info "Save the files"
function X(L, points)
    V = values.([mask(l, points) for l in L])
    return transpose(hcat(V...))
end

X₊ = X(envirovars, presencelayer)
X₋ = X(envirovars, absencelayer)

𝐲 = vcat(ones(Bool, sum(presencelayer)), zeros(Bool, sum(absencelayer)))
𝐗 = hcat(X₊, X₋)

fname = replace(taxname, " " => "_")

DelimitedFiles.writedlm("data/$(fname).X.dat", 𝐗)
DelimitedFiles.writedlm("data/$(fname).y.dat", 𝐲)