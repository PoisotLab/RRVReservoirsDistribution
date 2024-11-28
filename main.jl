using SpeciesDistributionToolkit
using CairoMakie
CairoMakie.activate!(; type = "png", px_per_unit = 2) # Haute résolution pour les figures
import Dates
import Images
import Downloads

# Superficie pour l'entraînement du modèle
bbox = (left=-85.0,bottom=25.0,right=-55.0,top=56.0)

# Données environnementales pour l'entraînement
provider = RasterData(WorldClim2, BioClim)
envirovars = [SDMLayer(provider; layer=i, resolution=5.0, bbox...) for i in eachindex(layers(provider))]

# Température annuelle moyenne
heatmap(envirovars[1])

# Occurrences
sp = taxon("Mephitis mephitis")
occ = occurrences(sp, envirovars[1], "occurrenceStatus" => "PRESENT", "limit" => 300, "continent" => "NORTH_AMERICA")
while length(occ) < min(10_000, count(occ)) # Max. 10000, sinon toutes
    occurrences!(occ)
end

# Spatial thinning
presencelayer = mask(envirovars[1], occ)

# Pseudo-absences
pa_mask = pseudoabsencemask(DistanceToEvent, presencelayer)
nodata!(pa_mask, x -> x <= 10.0)
nodata!(pa_mask, x -> x >= 100.0)
absencelayer = backgroundpoints(pa_mask, 2sum(presencelayer))

# On garde uniquement les points qui sont soit une absence, soit une présence
nodata!(absencelayer, false)
nodata!(presencelayer, false)

heatmap(envirovars[1])
scatter!(presencelayer, color=:red, markersize=2)
scatter!(absencelayer, color=:grey, markersize=2)
current_figure()

## Modèle 
sdm = SDM(ZScore, DecisionTree, envirovars, presencelayer, absencelayer)

train!(sdm)
heatmap(predict(sdm, envirovars; threshold=false))
scatter!(presencelayer, markersize=2, color=:orange)
current_figure()

## Mesure performance
folds = kfold(sdm; k=10)
cv = crossvalidate(sdm, folds)

mcc(cv.validation)
mcc(cv.training)

balancedaccuracy(cv.validation)

## Optimisation
forwardselection!(sdm, folds, [1]; verbose=true)

cv = crossvalidate(sdm, folds)

mcc(cv.validation)
mcc(cv.training)

## Bagging (random forest)
ensemble = Bagging(sdm, 32)
bagfeatures!(ensemble)
train!(ensemble)

# Image pour les cartes (depuis Phylopic)
sp_uuid = Phylopic.imagesof(sp; items=1)
Phylopic.attribution(sp_uuid)
sp_thumbnail_url = Phylopic.thumbnail(sp_uuid)
sp_thumbnail_tmp = Downloads.download(sp_thumbnail_url)
sp_image = Images.load(sp_thumbnail_tmp)
sp_size = Vec2f(reverse(size(sp_image) ./ 2))

#### Prédiction seulement pour le QC
QC = SpeciesDistributionToolkit.gadm("CAN", "Québec")
qcbbox = SpeciesDistributionToolkit.boundingbox(QC; padding=1.5)

provider = RasterData(WorldClim2, BioClim)
futureclim = Projection(SSP245, CanESM5)
qccurrent = [SDMLayer(provider; layer=i, resolution=5.0, qcbbox...) for i in eachindex(layers(provider))]
qcfuture = [SDMLayer(provider, futureclim, timespan=Dates.Year(2041)=>Dates.Year(2060); layer=i, resolution=5.0, qcbbox...) for i in eachindex(layers(provider))]

bg = copy(qccurrent[1])
msk = mask!(copy(qccurrent[1]), QC)

# Masquage des données
qccurrent = [mask!(l, msk) for l in qccurrent]
qcfuture = [mask!(l, msk) for l in qcfuture]

# Prédictions
pr_qc_current = predict(ensemble, qccurrent; threshold=false)
rg_qc_current = predict(ensemble, qccurrent, consensus=majority; threshold=true)

pr_qc_future = predict(ensemble, qcfuture; threshold=false)
rg_qc_future = predict(ensemble, qcfuture, consensus=majority; threshold=true)

# Prédiction (climat historique)
fig_pred_qc = Figure(size=(800, 700))
ax = Axis(fig_pred_qc[1,1], aspect=DataAspect(), xlabel="Longitude", ylabel="Latitude")
heatmap!(ax, bg, colormap=[:lightgrey, :lightgrey])
hm = heatmap!(ax, pr_qc_current, colormap=:Oranges, colorrange=(0, 1))
lines!(ax, QC, color=:black, linewidth=1)
Colorbar(fig_pred_qc[1, 2], hm, height=Relative(0.7))
scatter!(ax, [-60.0], [60.0]; marker = sp_image, markersize = sp_size)
current_figure()

# Prédiction (climat futur)
fig_future_qc = Figure(size=(800, 700))
ax = Axis(fig_future_qc[1,1], aspect=DataAspect())
heatmap!(ax, bg, colormap=[:lightgrey, :lightgrey])
hm = heatmap!(ax, pr_qc_future, colormap=:Oranges, colorrange=(0, 1))
lines!(ax, QC, color=:black, linewidth=1)
Colorbar(fig_future_qc[1, 2], hm, height=Relative(0.7))
scatter!(ax, [-60.0], [60.0]; marker = sp_image, markersize = sp_size)
current_figure()

# Changement d'aire de distribution
fig_rangediff = Figure(size=(800, 700))
ax = Axis(fig_rangediff[1,1], aspect=DataAspect())
heatmap!(ax, bg, colormap=[:lightgrey, :lightgrey])
heatmap!(ax, gainloss(rg_qc_current, rg_qc_future), colormap=:diverging_bky_60_10_c30_n256, colorrange=(-1,1))
scatter!(ax, [-60.0], [60.0]; marker = sp_image, markersize = sp_size)
lines!(ax, QC, color=:black, linewidth=1)
#Colorbar(fig_pred_qc[1, 2], hm, height=Relative(0.7))
scatter!(ax, [-60.0], [60.0]; marker = sp_image, markersize = sp_size)
current_figure()

# Importance de la température (Shapley)
shap_temp_current = explain(ensemble, qccurrent, 1; threshold=false, samples=50)
shap_temp_future = explain(ensemble, qcfuture, 1; threshold=false, samples=50)

fig_shap_current = Figure(size=(800, 700))
ax = Axis(fig_shap_current[1,1], aspect=DataAspect(), xlabel="Longitude", ylabel="Latitude")
heatmap!(ax, bg, colormap=[:lightgrey, :lightgrey])
hm = heatmap!(ax, shap_temp_current, colormap=:diverging_bwr_20_95_c54_n256, colorrange=(-0.4, 0.4))
scatter!(ax, [-60.0], [60.0]; marker = sp_image, markersize = sp_size)
lines!(ax, QC, color=:black, linewidth=1)
Colorbar(fig_shap_current[1, 2], hm, height=Relative(0.7))
scatter!(ax, [-60.0], [60.0]; marker = sp_image, markersize = sp_size)
current_figure()

fig_shap_future = Figure(size=(800, 700))
ax = Axis(fig_shap_future[1,1], aspect=DataAspect(), xlabel="Longitude", ylabel="Latitude")
heatmap!(ax, bg, colormap=[:lightgrey, :lightgrey])
hm = heatmap!(ax, shap_temp_future, colormap=:diverging_bwr_20_95_c54_n256, colorrange=(-0.4, 0.4))
scatter!(ax, [-60.0], [60.0]; marker = sp_image, markersize = sp_size)
lines!(ax, QC, color=:black, linewidth=1)
Colorbar(fig_shap_future[1, 2], hm, height=Relative(0.7))
scatter!(ax, [-60.0], [60.0]; marker = sp_image, markersize = sp_size)
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