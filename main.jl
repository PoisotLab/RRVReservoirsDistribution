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
import CSV
using Statistics

# Load CSV data
records = CSV.File("0022970-241126133413365.csv")
records = filter(r -> isequal(taxname)(r.species), records)

occ = Occurrences([Occurrence(r.species, true, (r.decimalLongitude, r.decimalLatitude), r.dateIdentified) for r in records])

bbox = SpeciesDistributionToolkit.boundingbox(occ; padding=2.0)

# Données environnementales pour l'entraînement
provider = RasterData(WorldClim2, BioClim)
envirovars = [SDMLayer(provider; layer=i, resolution=10.0, bbox...) for i in eachindex(layers(provider))]

# Image pour les cartes (depuis Phylopic)
sp_uuid = Phylopic.imagesof(taxname; items=1)
Phylopic.attribution(sp_uuid)
sp_thumbnail_url = Phylopic.thumbnail(sp_uuid)
sp_thumbnail_tmp = Downloads.download(sp_thumbnail_url)
sp_image = Images.load(sp_thumbnail_tmp)
sp_size = Vec2f(reverse(size(sp_image) ./ 2))

# Projection des données dans NAD 1983 Albers North America
envirovars = [interpolate(e; dest="+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs") for e in envirovars]

# Température annuelle moyenne
heatmap(envirovars[1])

# Fichier pour les figures
fpath = joinpath(pwd(), "figures", taxname)
if ~ispath(fpath)
    mkpath(fpath)
end

# Fichier pour les raster
rpath = joinpath(pwd(), "rasters", taxname)
if ~ispath(rpath)
    mkpath(rpath)
end

# Occurrences
fname = join(split(taxname, " "), "-") # Pour sauvegarder les figures

# Spatial thinning
presencelayer = mask(envirovars[1], occ)

# Pseudo-absences
event_dist = pseudoabsencemask(DistanceToEvent, presencelayer)
pa_mask = copy(event_dist)
nodata!(pa_mask, x -> x <= 10.0)
nodata!(pa_mask, x -> x >= 500.0)
absencelayer = backgroundpoints(pa_mask, 2sum(presencelayer))

# On garde uniquement les points qui sont soit une absence, soit une présence
nodata!(absencelayer, false)
nodata!(presencelayer, false)

# Carte avec les données
fig_training = Figure()
ax = Axis(fig_training[1, 1]; aspect=DataAspect())
heatmap!(ax, envirovars[1], colormap=[colorant"#efefef", colorant"#efefef"])
scatter!(ax, presencelayer, color=:orange, markersize=3)
scatter!(ax, absencelayer, color=:black, markersize=2)
#scatter!(ax, [-65.0], [37.5]; marker=sp_image, markersize=sp_size)
current_figure()
save(joinpath(fpath, "$(fname)-data.png"), current_figure())

# On sauvegarde les raster
bg = !isnan.(envirovars[1])
SimpleSDMLayers.save(joinpath(rpath, "$(fname)-background.tiff"), convert(SDMLayer{UInt8}, bg))
SimpleSDMLayers.save(joinpath(rpath, "$(fname)-presence.tiff"), convert(SDMLayer{UInt8}, presencelayer))
SimpleSDMLayers.save(joinpath(rpath, "$(fname)-absence.tiff"), convert(SDMLayer{UInt8}, absencelayer))

## Modèle 
sdm = SDM(RawData, DecisionTree, envirovars, presencelayer, absencelayer)
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
ensemble = Bagging(sdm, 16)
bagfeatures!(ensemble)
train!(ensemble)

## OOB error rate
outofbag(ensemble) |> (x) -> 1 - accuracy(x)

#### Prédiction seulement pour le QC
QC = SpeciesDistributionToolkit.gadm("CAN", "Québec")
qcbbox = SpeciesDistributionToolkit.boundingbox(QC; padding=0.0)

# Projection vers NAD 1983 Albers North America
qccurrent = [SDMLayer(provider; layer=i, resolution=5.0, qcbbox...) for i in eachindex(layers(provider))]
qccurrent = [interpolate(l; dest="+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs") for l in qccurrent]
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
ax = Axis(fig_pred_qc[1, 1])#, aspect=DataAspect())
heatmap!(ax, bg, colormap=[colorant"#efefef", colorant"#efefef"])
hm = heatmap!(ax, pr_qc_current, colormap=:Oranges, colorrange=(0, 1))
#lines!(ax, QC, color=:black, linewidth=1)
Colorbar(fig_pred_qc[1, 2], hm, height=Relative(0.7))
#scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
hidespines!(ax)
hidedecorations!(ax)
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

# Stemplot - importance des variables
s_varimp = mean.(map(l -> abs.(l), S))
s_varimp ./= sum(s_varimp)
o_varimp = sortperm(s_varimp; rev=true)

fig_varimp = Figure(size=(800, 700))
ax = Axis(fig_varimp[1, 1], ylabel="Cumulative variable importance")
scatter!(ax, cumsum(s_varimp[o_varimp]))
varlab = "BIO" .* string.(variables(ensemble)[o_varimp])
text!(ax, 1:length(variables(ensemble)), cumsum(s_varimp[o_varimp]); text=varlab, align=(:center, :top), offset=(0.0, -10.0))
hidexdecorations!(ax)
ylims!(ax, 0, 1)
current_figure()
save(joinpath(fpath, "$(fname)-varimp-current.png"), current_figure())

for ssp in [SSP126, SSP245, SSP370, SSP585]
    futureclim = Projection(ssp, CanESM5)

    for tsp in SimpleSDMDatasets.timespans(provider, futureclim)
        range_begin = tsp.first.value
        range_end = tsp.second.value
        range_txt = "$(range_begin)-$(range_end)"


        qcfuture = [SDMLayer(provider, futureclim, timespan=tsp; layer=i, resolution=2.5, qcbbox...) for i in eachindex(layers(provider))]

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

        # Range (climat futur)
        fig_pred_qc = Figure(size=(800, 700))
        ax = Axis(fig_pred_qc[1, 1], aspect=DataAspect())
        heatmap!(ax, bg, colormap=[colorant"#efefef", colorant"#efefef"])
        hm = heatmap!(ax, rg_qc_future, colormap=:Oranges, colorrange=(0, 1))
        lines!(ax, QC, color=:black, linewidth=1)
        Colorbar(fig_pred_qc[1, 2], hm, height=Relative(0.7))
        scatter!(ax, [-60.0], [60.0]; marker=sp_image, markersize=sp_size)
        current_figure()
        save(joinpath(fpath, "$(fname)-range-future-$(ssp)-$(range_txt).png"), current_figure())

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


        # Stemplot - importance des variables
        s_varimp = mean.(map(l -> abs.(l), S))
        s_varimp ./= sum(s_varimp)
        o_varimp = sortperm(s_varimp; rev=true)

        fig_varimp = Figure(size=(800, 700))
        ax = Axis(fig_varimp[1, 1], ylabel="Cumulative variable importance")
        scatter!(ax, cumsum(s_varimp[o_varimp]))
        varlab = "BIO" .* string.(variables(ensemble)[o_varimp])
        text!(ax, 1:length(variables(ensemble)), cumsum(s_varimp[o_varimp]); text=varlab, align=(:center, :top), offset=(0.0, -10.0))
        hidexdecorations!(ax)
        ylims!(ax, 0, 1)
        current_figure()
        save(joinpath(fpath, "$(fname)-varimp-$(ssp)-$(range_txt).png"), current_figure())
    end
end