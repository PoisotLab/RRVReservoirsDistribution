taxname = "Mephitis mephitis"
if ~isempty(ARGS)
    taxname = join(ARGS[1:2], " ")
end
@info "Running for $(taxname)"

using SpeciesDistributionToolkit
using CairoMakie
CairoMakie.activate!(; type="png", px_per_unit=3) # Haute résolution pour les figures
import Dates
import Images
import Downloads
using Statistics

# Superficie pour l'entraînement du modèle
bbox = (left=-90.0, bottom=32.5, right=-60.0, top=56.0)

# Données environnementales pour l'entraînement
provider = RasterData(WorldClim2, BioClim)
envirovars = [SDMLayer(provider; layer=i, resolution=10.0, bbox...) for i in eachindex(layers(provider))]

# Température annuelle moyenne
heatmap(envirovars[1])

# Fichier pour les figures
fpath = joinpath(pwd(), "figures")
if ~ispath(fpath)
    mkpath(fpath)
end

# Fichier pour les raster
rpath = joinpath(pwd(), "rasters")
if ~ispath(rpath)
    mkpath(rpath)
end

# Occurrences
sp = taxon(taxname)
fname = join(split(sp.name, " "), "-") # Pour sauvegarder les figures

occ = occurrences(sp, envirovars[1], "occurrenceStatus" => "PRESENT", "limit" => 300, "continent" => "NORTH_AMERICA")
while length(occ) < min(4_000, count(occ)) # Max. 10000, sinon toutes
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

# Image pour les cartes (depuis Phylopic)
sp_uuid = Phylopic.imagesof(sp; items=1)
Phylopic.attribution(sp_uuid)
sp_thumbnail_url = Phylopic.thumbnail(sp_uuid)
sp_thumbnail_tmp = Downloads.download(sp_thumbnail_url)
sp_image = Images.load(sp_thumbnail_tmp)
sp_size = Vec2f(reverse(size(sp_image) ./ 2))

# Carte avec les données
fig_training = Figure()
ax = Axis(fig_training[1, 1]; aspect=DataAspect())
heatmap!(ax, envirovars[1], colormap=[colorant"#efefef", colorant"#efefef"])
scatter!(ax, presencelayer, color=:orange, markersize=3)
scatter!(ax, absencelayer, color=:black, markersize=2)
scatter!(ax, [-65.0], [37.5]; marker=sp_image, markersize=sp_size)
current_figure()
save(joinpath(fpath, "$(fname)-data.png"), current_figure())

# On sauvegarde les raster
bg = !isnan.(envirovars[1])
SimpleSDMLayers.save(joinpath(rpath, "$(fname)-background.tiff"), convert(SDMLayer{UInt8}, bg))
SimpleSDMLayers.save(joinpath(rpath, "$(fname)-presence.tiff"), convert(SDMLayer{UInt8}, presencelayer))
SimpleSDMLayers.save(joinpath(rpath, "$(fname)-absence.tiff"), convert(SDMLayer{UInt8}, absencelayer))

## Modèle 
sdm = SDM(ZScore, DecisionTree, envirovars, presencelayer, absencelayer)

train!(sdm)

## Optimisation
folds = kfold(sdm; k=10)
reset!(sdm)
forwardselection!(sdm, folds, [1]; verbose=true)

## Performance
cv = crossvalidate(sdm, folds)
mcc(cv.validation)
mcc(cv.training)

## Bagging (random forest)
ensemble = Bagging(sdm, 32)
bagfeatures!(ensemble)
train!(ensemble)

## OOB error rate
outofbag(ensemble) |> (x) -> 1 - accuracy(x)

#### Prédiction seulement pour le QC
QC = SpeciesDistributionToolkit.gadm("CAN", "Québec")
qcbbox = SpeciesDistributionToolkit.boundingbox(QC; padding=1.5)

qccurrent = [SDMLayer(provider; layer=i, resolution=10.0, qcbbox...) for i in eachindex(layers(provider))]
bg = copy(qccurrent[1])
SimpleSDMLayers.save(joinpath(rpath, "$(fname)-qc.tiff"), convert(SDMLayer{UInt8}, !isnan.(bg)))
msk = mask!(copy(qccurrent[1]), QC)

# Masquage des données
qccurrent = [mask!(l, msk) for l in qccurrent]

# Prédictions
pr_qc_current = predict(ensemble, qccurrent; threshold=false)
SimpleSDMLayers.save(joinpath(rpath, "$(fname)-prediction-current.tiff"), convert(SDMLayer{Float32}, pr_qc_current))
rg_qc_current = predict(ensemble, qccurrent, consensus=majority; threshold=true)
SimpleSDMLayers.save(joinpath(rpath, "$(fname)-range-current.tiff"), convert(SDMLayer{Float32}, rg_qc_current))
S = explain(ensemble, qccurrent; threshold=false, samples=50)
minf = mosaic(x -> argmax(abs.(x)), S)

shaprange(v) = maximum(abs.(quantile(v, [0.05, 0.95]))) .* (-1, 1)

# Prédiction (climat historique)
fig_pred_qc = Figure(size=(800, 700))
ax = Axis(fig_pred_qc[1, 1], aspect=DataAspect())
heatmap!(ax, bg, colormap=[colorant"#efefef", colorant"#efefef"])
hm = heatmap!(ax, pr_qc_current, colormap=:Oranges, colorrange=(0, 1))
lines!(ax, QC, color=:black, linewidth=1)
Colorbar(fig_pred_qc[1, 2], hm, height=Relative(0.7))
scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
current_figure()
save(joinpath(fpath, "$(fname)-pred-current.png"), current_figure())

# Range (climat historique)
fig_pred_qc = Figure(size=(800, 700))
ax = Axis(fig_pred_qc[1, 1], aspect=DataAspect())
heatmap!(ax, bg, colormap=[colorant"#efefef", colorant"#efefef"])
hm = heatmap!(ax, rg_qc_current, colormap=:Oranges, colorrange=(0, 1))
lines!(ax, QC, color=:black, linewidth=1)
Colorbar(fig_pred_qc[1, 2], hm, height=Relative(0.7))
scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
current_figure()
save(joinpath(fpath, "$(fname)-range-current.png"), current_figure())

for (i, v) in enumerate(variables(ensemble))
    fig_S = Figure(size=(800, 700))
    ax = Axis(fig_S[1, 1], aspect=DataAspect())
    heatmap!(ax, bg, colormap=[colorant"#efefef", colorant"#efefef"])
    hm = heatmap!(ax, S[i], colormap=:diverging_bwg_20_95_c41_n256, colorrange=shaprange(S[i]))
    lines!(ax, QC, color=:black, linewidth=1)
    Colorbar(fig_S[1, 2], hm, height=Relative(0.7))
    scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
    current_figure()
    save(joinpath(fpath, "$(fname)-BIO$(v)-effect-current.png"), current_figure())
end

# Variable la plus importante
fig_maxinf = Figure(size=(800, 700))
ax = Axis(fig_maxinf[1, 1], aspect=DataAspect())
heatmap!(ax, bg, colormap=[colorant"#efefef", colorant"#efefef"])
var_colors = cgrad(:diverging_rainbow_bgymr_45_85_c67_n256, length(variables(ensemble)), categorical=true)
hm = heatmap!(ax, minf; colormap=var_colors, colorrange=(1, length(variables(ensemble))))
lines!(ax, QC, color=:black, linewidth=1)
scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
Legend(
    fig_maxinf[2, 1],
    [PolyElement(; color=var_colors[i]) for i in 1:length(variables(ensemble))],
    "BIO" .* string.(variables(ensemble));
    orientation=:horizontal,
    nbanks=1,
)
current_figure()
save(joinpath(fpath, "$(fname)-mostimportant-current.png"), current_figure())

for ssp in [SSP126, SSP245, SSP370, SSP585]
    futureclim = Projection(ssp, CanESM5)

    for tsp in SimpleSDMDatasets.timespans(provider, futureclim)
        range_begin = tsp.first.value
        range_end = tsp.second.value
        range_txt = "$(range_begin)-$(range_end)"


        qcfuture = [SDMLayer(provider, futureclim, timespan=tsp; layer=i, resolution=10.0, qcbbox...) for i in eachindex(layers(provider))]

        # NE PAS CHANGER
        for i in eachindex(qcfuture)
            qcfuture[i].x = qccurrent[i].x
            qcfuture[i].y = qccurrent[i].y
        end
        #
        qcfuture = [mask!(l, msk) for l in qcfuture]

        pr_qc_future = predict(ensemble, qcfuture; threshold=false)
        SimpleSDMLayers.save(joinpath(rpath, "$(fname)-prediction-$(ssp)-$(range_txt).tiff"), convert(SDMLayer{Float32}, pr_qc_future))
        rg_qc_future = predict(ensemble, qcfuture, consensus=majority; threshold=true)
        SimpleSDMLayers.save(joinpath(rpath, "$(fname)-range-$(ssp)-$(range_txt).tiff"), convert(SDMLayer{Float32}, rg_qc_future))

        # Prédiction (climat futur)
        fig_future_qc = Figure(size=(800, 700))
        ax = Axis(fig_future_qc[1, 1], aspect=DataAspect())
        heatmap!(ax, bg, colormap=[colorant"#efefef", colorant"#efefef"])
        hm = heatmap!(ax, pr_qc_future, colormap=:Oranges, colorrange=(0, 1))
        lines!(ax, QC, color=:black, linewidth=1)
        Colorbar(fig_future_qc[1, 2], hm, height=Relative(0.7))
        scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
        current_figure()
        save(joinpath(fpath, "$(fname)-pred-future-$(ssp)-$(range_txt).png"), current_figure())

        # Changement d'aire de distribution
        fig_rangediff = Figure(size=(800, 700))
        ax = Axis(fig_rangediff[1, 1], aspect=DataAspect())
        heatmap!(ax, bg, colormap=[colorant"#efefef", colorant"#efefef"])
        heatmap!(ax, gainloss(rg_qc_current, rg_qc_future), colormap=:diverging_isoluminant_cjm_75_c23_n256, colorrange=(-1, 1))
        scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
        lines!(ax, QC, color=:black, linewidth=1)
        #Colorbar(fig_pred_qc[1, 2], hm, height=Relative(0.7))
        scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
        current_figure()
        save(joinpath(fpath, "$(fname)-range-shift-$(ssp)-$(range_txt).png"), current_figure())

        # Importance de la température (Shapley)
        S = explain(ensemble, qcfuture; threshold=false, samples=50)
        minf = mosaic(x -> argmax(abs.(x)), S)


        for (i, v) in enumerate(variables(ensemble))
            fig_S = Figure(size=(800, 700))
            ax = Axis(fig_S[1, 1], aspect=DataAspect())
            heatmap!(ax, bg, colormap=[colorant"#efefef", colorant"#efefef"])
            hm = heatmap!(ax, S[i], colormap=:diverging_bwg_20_95_c41_n256, colorrange=shaprange(S[i]))
            lines!(ax, QC, color=:black, linewidth=1)
            Colorbar(fig_S[1, 2], hm, height=Relative(0.7))
            scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
            current_figure()
            save(joinpath(fpath, "$(fname)-BIO$(v)-effect-$(ssp)-$(range_txt).png"), current_figure())
        end

        # Variable la plus importante
        fig_maxinf = Figure(size=(800, 700))
        ax = Axis(fig_maxinf[1, 1], aspect=DataAspect())
        heatmap!(ax, bg, colormap=[colorant"#efefef", colorant"#efefef"])
        var_colors = cgrad(:diverging_rainbow_bgymr_45_85_c67_n256, length(variables(ensemble)), categorical=true)
        hm = heatmap!(ax, minf; colormap=var_colors, colorrange=(1, length(variables(ensemble))))
        lines!(ax, QC, color=:black, linewidth=1)
        scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
        Legend(
            fig_maxinf[2, 1],
            [PolyElement(; color=var_colors[i]) for i in 1:length(variables(ensemble))],
            "BIO" .* string.(variables(ensemble));
            orientation=:horizontal,
            nbanks=1,
        )
        current_figure()
        save(joinpath(fpath, "$(fname)-mostimportant-$(ssp)-$(range_txt).png"), current_figure())

        # Nouveauté climatique
        using Statistics
        μ = mean.(qccurrent)
        σ = std.(qccurrent)

        cr_current = (qccurrent .- μ) ./ σ
        cr_future = (qcfuture .- μ) ./ σ

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
        ax = Axis(fig_novelty[1, 1], aspect=DataAspect())
        heatmap!(ax, bg, colormap=[colorant"#efefef", colorant"#efefef"])
        hm = heatmap!(ax, Δclim, colormap=:linear_worb_100_25_c53_n256, colorrange=quantile(Δclim, [0.05, 0.95]))
        lines!(ax, QC, color=:black, linewidth=1)
        Colorbar(fig_novelty[1, 2], hm, height=Relative(0.7))
        scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
        current_figure()
        save(joinpath(fpath, "$(fname)-climate-novelty-$(ssp)-$(range_txt).png"), current_figure())
    end
end