taxname = "Mephitis mephitis"
if ~isempty(ARGS)
    taxname = join(ARGS[1:2], " ")
end
@info "Running for $(taxname)"

using SpeciesDistributionToolkit
using CairoMakie
CairoMakie.activate!(; type="png", px_per_unit=2) # Haute résolution pour les figures
import Dates
import Images
import Downloads

# Superficie pour l'entraînement du modèle
bbox = (left=-90.0, bottom=35.0, right=-60.0, top=56.0)

# Données environnementales pour l'entraînement
provider = RasterData(WorldClim2, BioClim)
envirovars = [SDMLayer(provider; layer=i, resolution=5.0, bbox...) for i in eachindex(layers(provider))]

# Température annuelle moyenne
heatmap(envirovars[1])

# Fichier pour les figures
fpath = joinpath(pwd(), "figures")
if ~ispath(fpath)
    mkpath(fpath)
end

# Occurrences
sp = taxon(taxname)
fname = join(split(sp.name, " "), "-") # Pour sauvegarder les figures

occ = occurrences(sp, envirovars[1], "occurrenceStatus" => "PRESENT", "limit" => 300, "continent" => "NORTH_AMERICA")
while length(occ) < min(20_000, count(occ)) # Max. 10000, sinon toutes
    occurrences!(occ)
end

# Spatial thinning
presencelayer = mask(envirovars[1], occ)

# Pseudo-absences
event_dist = pseudoabsencemask(DistanceToEvent, presencelayer)
pa_mask = copy(event_dist)
nodata!(pa_mask, x -> x <= 8.0)
nodata!(pa_mask, x -> x >= 200.0)
absencelayer = backgroundpoints(pa_mask, 2sum(presencelayer))

# On garde uniquement les points qui sont soit une absence, soit une présence
nodata!(absencelayer, false)
nodata!(presencelayer, false)

fig_training = Figure()
ax = Axis(fig_training[1,1]; aspect=DataAspect())
heatmap(envirovars[1], colormap=[:lightgrey, :lightgrey])
scatter!(presencelayer, color=:red, markersize=3)
scatter!(absencelayer, color=:black, markersize=2)
current_figure()
save(joinpath(fpath, "$(fname)-data.png"), current_figure())

## Modèle 
sdm = SDM(ZScore, DecisionTree, envirovars, presencelayer, absencelayer)

train!(sdm)

heatmap(predict(sdm, envirovars; threshold=false))
scatter!(presencelayer, markersize=2, color=:orange)
current_figure()

## Optimisation
folds = kfold(sdm; k=10)
reset!(sdm)
forwardselection!(sdm, folds, [1]; verbose=true)

## Performance
cv = crossvalidate(sdm, folds)
mcc(cv.validation)
mcc(cv.training)

## Bagging (random forest)
ensemble = Bagging(sdm, 64)
bagfeatures!(ensemble)
train!(ensemble)

## OOB error rate
outofbag(ensemble) |> (x) -> 1 - accuracy(x)

# Image pour les cartes (depuis Phylopic)
sp_uuid = Phylopic.imagesof(sp; items=1)
Phylopic.attribution(sp_uuid)
sp_thumbnail_url = Phylopic.thumbnail(sp_uuid)
sp_thumbnail_tmp = Downloads.download(sp_thumbnail_url)
sp_image = Images.load(sp_thumbnail_tmp)
sp_size = Vec2f(reverse(size(sp_image) ./ 2))

#### Prédiction seulement pour le QC
QC = SpeciesDistributionToolkit.gadm("CAN", "Québec")
qcbbox = SpeciesDistributionToolkit.boundingbox(QC; padding=0.5)

provider = RasterData(WorldClim2, BioClim)
futureclim = Projection(SSP245, CanESM5)
qccurrent = [SDMLayer(provider; layer=i, resolution=5.0, qcbbox...) for i in eachindex(layers(provider))]
qcfuture = [SDMLayer(provider, futureclim, timespan=Dates.Year(2021) => Dates.Year(2040); layer=i, resolution=5.0, qcbbox...) for i in eachindex(layers(provider))]

bg = copy(qccurrent[1])
msk = mask!(copy(qccurrent[1]), QC)

# Masquage des données
qccurrent = [mask!(l, msk) for l in qccurrent]
# NE PAS CHANGER
for i in eachindex(qcfuture)
    qcfuture[i].x = qccurrent[i].x
    qcfuture[i].y = qccurrent[i].y
end
#
qcfuture = [mask!(l, msk) for l in qcfuture]

# Prédictions
pr_qc_current = predict(ensemble, qccurrent; threshold=false)
rg_qc_current = predict(ensemble, qccurrent, consensus=majority; threshold=true)

pr_qc_future = predict(ensemble, qcfuture; threshold=false)
rg_qc_future = predict(ensemble, qcfuture, consensus=majority; threshold=true)

# Prédiction (climat historique)
fig_pred_qc = Figure(size=(800, 700))
ax = Axis(fig_pred_qc[1, 1], aspect=DataAspect(), xlabel="Longitude", ylabel="Latitude")
heatmap!(ax, bg, colormap=[:lightgrey, :lightgrey])
hm = heatmap!(ax, pr_qc_current, colormap=:Oranges, colorrange=(0, 1))
lines!(ax, QC, color=:black, linewidth=1)
Colorbar(fig_pred_qc[1, 2], hm, height=Relative(0.7))
scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
current_figure()
save(joinpath(fpath, "$(fname)-pred-current.png"), current_figure())

# Prédiction (climat futur)
fig_future_qc = Figure(size=(800, 700))
ax = Axis(fig_future_qc[1, 1], aspect=DataAspect())
heatmap!(ax, bg, colormap=[:lightgrey, :lightgrey])
hm = heatmap!(ax, pr_qc_future, colormap=:Oranges, colorrange=(0, 1))
lines!(ax, QC, color=:black, linewidth=1)
Colorbar(fig_future_qc[1, 2], hm, height=Relative(0.7))
scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
current_figure()
save(joinpath(fpath, "$(fname)-pred-future.png"), current_figure())

# Changement d'aire de distribution
fig_rangediff = Figure(size=(800, 700))
ax = Axis(fig_rangediff[1, 1], aspect=DataAspect())
heatmap!(ax, bg, colormap=[:lightgrey, :lightgrey])
heatmap!(ax, gainloss(rg_qc_current, rg_qc_future), colormap=:diverging_bky_60_10_c30_n256, colorrange=(-1, 1))
scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
lines!(ax, QC, color=:black, linewidth=1)
#Colorbar(fig_pred_qc[1, 2], hm, height=Relative(0.7))
scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
current_figure()
save(joinpath(fpath, "$(fname)-range-shift.png"), current_figure())

# Importance de la température (Shapley)
shap_temp_current = explain(ensemble, qccurrent, 1; threshold=false, samples=50)
shap_temp_future = explain(ensemble, qcfuture, 1; threshold=false, samples=50)

fig_shap_current = Figure(size=(800, 700))
ax = Axis(fig_shap_current[1, 1], aspect=DataAspect(), xlabel="Longitude", ylabel="Latitude")
heatmap!(ax, bg, colormap=[:lightgrey, :lightgrey])
hm = heatmap!(ax, shap_temp_current, colormap=:diverging_bwr_20_95_c54_n256, colorrange=(-0.4, 0.4))
scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
lines!(ax, QC, color=:black, linewidth=1)
Colorbar(fig_shap_current[1, 2], hm, height=Relative(0.7))
scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
current_figure()
save(joinpath(fpath, "$(fname)-temperature-effect-current.png"), current_figure())

fig_shap_future = Figure(size=(800, 700))
ax = Axis(fig_shap_future[1, 1], aspect=DataAspect(), xlabel="Longitude", ylabel="Latitude")
heatmap!(ax, bg, colormap=[:lightgrey, :lightgrey])
hm = heatmap!(ax, shap_temp_future, colormap=:diverging_bwr_20_95_c54_n256, colorrange=(-0.4, 0.4))
lines!(ax, QC, color=:black, linewidth=1)
Colorbar(fig_shap_future[1, 2], hm, height=Relative(0.7))
scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
current_figure()
save(joinpath(fpath, "$(fname)-temperature-effect-future.png"), current_figure())

# Nouveauté climatique
using Statistics
μ = mean.(qccurrent)
σ = std.(qccurrent)

cr_current = (qccurrent .- μ) ./ σ;
cr_future = (qcfuture .- μ) ./ σ;

Δclim = similar(cr_current[1])
vals = values.(cr_future)
for position in keys(cr_current[1])
    dvar = [(cr_current[u][position] .- vals[u]) .^ 2.0 for u in variables(ensemble)]
    sm = vec(sum(hcat(dvar...); dims=2))
    sm[findall(isnan, sm)] .= 1000000.0
    Δclim[position] = minimum(sqrt.(sm))
end

# Nouveauté climatique
fig_novelty = Figure(size=(800, 700))
ax = Axis(fig_novelty[1, 1], aspect=DataAspect(), xlabel="Longitude", ylabel="Latitude")
heatmap!(ax, bg, colormap=[:lightgrey, :lightgrey])
hm = heatmap!(ax, Δclim, colormap=:linear_worb_100_25_c53_n256, colorrange=(0.0, floor(maximum(Δclim))))
lines!(ax, QC, color=:black, linewidth=1)
Colorbar(fig_novelty[1, 2], hm, height=Relative(0.7))
scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
current_figure()
save(joinpath(fpath, "$(fname)-climate-novelty.png"), current_figure())