using SpeciesDistributionToolkit
using CairoMakie
import Dates

## Data

# Variables clim
bbox = (left = -78.288574, bottom = 43.572432, right=-68.653564,top=48.661943)

provider = RasterData(WorldClim2, BioClim)
envirovars = [SDMLayer(provider; layer=i, resolution=5.0, bbox...) for i in eachindex(layers(provider))]

heatmap(envirovars[1])

# Occurrences
sp = taxon("Mephitis mephitis")
occ = occurrences(sp, envirovars[1], "occurrenceStatus" => "PRESENT", "limit" => 300)
while length(occ) < 1400
    occurrences!(occ)
end

# Spatial thinning
presencelayer = mask(envirovars[1], occ)

# Pseudo-absences
pa_mask = pseudoabsencemask(DistanceToEvent, presencelayer)
nodata!(pa_mask, x -> x <= 10.0)
nodata!(pa_mask, x -> x >= 100.0)
absencelayer = backgroundpoints(pa_mask, 2sum(presencelayer))

# ???
nodata!(absencelayer, false)
nodata!(presencelayer, false)

heatmap(envirovars[1])
scatter!(presencelayer, color=:red)
scatter!(absencelayer, color=:grey)
current_figure()

## Modèle 
sdm = SDM(ZScore, DecisionTree, envirovars, presencelayer, absencelayer)

train!(sdm)
heatmap(predict(sdm, envirovars; threshold=false))
scatter!(presencelayer)
current_figure()

## Mesure performance
folds = kfold(sdm; k=10)
cv = crossvalidate(sdm, folds)

mcc(cv.validation)
mcc(cv.training)

balancedaccuracy(cv.validation)

## Optimisation
forwardselection!(sdm, folds, [1])

cv = crossvalidate(sdm, folds)

mcc(cv.validation)
mcc(cv.training)

## Importance des variables
variableimportance(sdm, folds)
variables(sdm)

## Résultats
heatmap(predict(sdm, envirovars; threshold=false), colormap=:Greys)
contour!(predict(sdm, envirovars), color=:red)
scatter!(presencelayer, color=:orange)
scatter!(absencelayer, color=:red, markersize=2)
res = current_figure()
save("mephitis_habitat.png", res)

## Bagging
ensemble = Bagging(sdm, 32)
bagfeatures!(ensemble)
train!(ensemble)

heatmap(predict(ensemble, envirovars; threshold=false), colormap=:Greys)
contour!(predict(ensemble, envirovars, consensus=majority), color=:red)
scatter!(presencelayer, color=:orange)
#scatter!(absencelayer, color=:red, markersize=2)
current_figure()

heatmap(predict(ensemble, envirovars; threshold=false, consensus=iqr), colormap=:Greys)
contour!(predict(ensemble, envirovars, consensus=majority), color=:red)
scatter!(presencelayer, color=:orange)
scatter!(absencelayer, color=:red, markersize=2)
current_figure()

heatmap(quantize(predict(ensemble, envirovars; threshold=false, consensus=iqr), 5), colormap=:Greys)
contour!(predict(ensemble, envirovars, consensus=majority), color=:red)
un = current_figure()
save("mephitis_uncertainty.png", un)

## Projection CC
futureclim = Projection(SSP370, CanESM5)
futurevars = [SDMLayer(provider, futureclim, timespan=Dates.Year(2041)=>Dates.Year(2060); layer=i, resolution=5.0, bbox...) for i in eachindex(layers(provider))]

predict(sdm, futurevars; threshold=false) |> heatmap
gl = heatmap(gainloss(predict(sdm, envirovars), predict(sdm, futurevars)), colormap=[:blue, :lightgrey, :red])
save("mephitis_gainloss_SSP370_CanESM5.png", gl)

heatmap(predict(ensemble, futurevars; threshold=false), colormap=:Greys)
contour!(predict(ensemble, futurevars, consensus=majority), color=:red)
proj = current_figure()
save("mephitis_proj_SSP370_CanESM5.png", proj)


## Species overlap

# Variables climatiques
bbox = (left = -78.288574, bottom = 43.572432, right = -68.653564, top = 48.661943)
provider = RasterData(WorldClim2, BioClim)
envirovars = [SDMLayer(provider; layer=i, resolution=5.0, bbox...) for i in eachindex(layers(provider))]

# Fonction pour créer un SDM pour une espèce
function create_sdm(taxon_name, envirovars)
    sp = taxon(taxon_name)
    occ = occurrences(sp, envirovars[1], "occurrenceStatus" => "PRESENT", "limit" => 300)
    while length(occ) < 1400
        occurrences!(occ)
    end
    presencelayer = mask(envirovars[1], occ)
    pa_mask = pseudoabsencemask(DistanceToEvent, presencelayer)
    nodata!(pa_mask, x -> x <= 10.0)
    nodata!(pa_mask, x -> x >= 100.0)
    absencelayer = backgroundpoints(pa_mask, 2sum(presencelayer))
    nodata!(absencelayer, false)
    nodata!(presencelayer, false)

    sdm = SDM(ZScore, DecisionTree, envirovars, presencelayer, absencelayer)
    train!(sdm)
    return sdm
end

# Création des SDM pour les deux espèces
sdm1 = create_sdm("Mephitis mephitis", envirovars)
sdm2 = create_sdm("Procyon lotor", envirovars)

# Projection des probabilités pour chaque espèce
proj1 = predict(sdm1, envirovars; threshold = false)
proj2 = predict(sdm2, envirovars; threshold = false)

# Fusion des projections pour montrer les zones de chevauchement
combined_proj = proj1 .+ proj2

# Visualisation
heatmap(combined_proj, colormap = [:blue, :lightgrey, :red])  # Bleu : Espèce 1, Rouge : Espèce 2, Gris clair : Chevauchement
#contour!(proj1 .> 0.5, color = :blue, linewidth = 1.5, label = "Mephitis mephitis")
#contour!(proj2 .> 0.5, color = :red, linewidth = 1.5, label = "Procyon lotor")

# Affichage
current_figure()

# Sauvegarde
save("species_overlap.png", current_figure())