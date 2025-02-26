using SpeciesDistributionToolkit
using CairoMakie
import DelimitedFiles
using PrettyTables
import CSV

taxname = "Mephitis mephitis"
if ~isempty(ARGS)
    taxname = join(ARGS[1:2], " ")
end
@info "Running for $(taxname)"

𝐗 = DelimitedFiles.readdlm("data/$(replace(taxname, " " => "_")).X.dat", Float32)
𝐲 = vec(DelimitedFiles.readdlm("data/$(replace(taxname, " " => "_")).y.dat", Bool))

@info "Train the SDM for all the known data"
sdm = SDM(ZScore, Logistic, 𝐗, 𝐲)
folds = kfold(sdm)

# Set some better training parameters
sdm.classifier.verbose = false
sdm.classifier.η = 1e-3
sdm.classifier.epochs = 5000

@info "Select variables"
forwardselection!(sdm, folds; verbose=true)

@info "Report variables"
DelimitedFiles.writedlm("data/$(replace(taxname, " " => "_")).params", variables(sdm))
DelimitedFiles.writedlm("data/$(replace(taxname, " " => "_")).threshold", threshold(sdm))

@info "Report on cross-validation"
cv = crossvalidate(sdm, folds)
@info "MCC val.", mcc(cv.validation)
@info "MCC trn.", mcc(cv.training)
@info "PPV val.", ppv(cv.validation)
@info "PPV trn.", ppv(cv.training)
@info "NPV val.", npv(cv.validation)
@info "NPV trn.", npv(cv.training)

@info "Loading bioclim data for prediction"
provider = RasterData(WorldClim2, BioClim)
QC = SpeciesDistributionToolkit.gadm("CAN", "Québec")
bbox = SpeciesDistributionToolkit.boundingbox(QC; padding=0.0)
envirovars = [SDMLayer(provider; layer=i, resolution=2.5, bbox...) for i in eachindex(layers(provider))]
mask!(envirovars, QC)

@info "Baseline prediction"
predict(sdm, envirovars; threshold=false) |> heatmap
predict(sdm, envirovars; threshold=true) |> heatmap

@info "Prepare to save stuff"
fname = replace(taxname, " " => "_")

@info "Plot the partial responses"
lnames = layers(provider)
ldescr = layerdescriptions(provider)

for v in variables(sdm)
    @info ldescr[lnames[v]]
    fpath = joinpath("figures", fname * "_" * "partialresponse_$(lnames[v]).png")
    f = Figure(; size=(600, 600))
    ax = Axis(f[1, 1], xlabel=ldescr[lnames[v]], ylabel="Score for $(taxname)")
    for i in 1:100
        lines!(ax, partialresponse(sdm, v; inflated=true, threshold=false)..., color=:grey, alpha=0.4)
    end
    lines!(ax, partialresponse(sdm, v; threshold=false)..., color=:black, linewidth=2)
    tightlimits!(ax)
    CairoMakie.save(fpath, f)
end

@info "Shapley values"
S = explain(sdm, envirovars; threshold=false)

@info "Save the Shapley values"
sname = joinpath("rasters", fname * "_" * "shapley.tif")
SimpleSDMLayers.save(sname, S)

@info "Save the range"
sname = joinpath("rasters", fname * "_" * "historical.tif")
proba = predict(sdm, envirovars; threshold=false)
range = predict(sdm, envirovars; threshold=true)
SimpleSDMLayers.save(sname, SDMLayer{Float64}[proba, range])

@info "Forecast"
for ssp in [SSP126, SSP245, SSP370, SSP585]
    futureclim = Projection(ssp, CanESM5)
    @info ssp
    for tsp in SimpleSDMDatasets.timespans(provider, futureclim)
        range_begin = tsp.first.value
        range_end = tsp.second.value
        range_txt = "$(range_begin)-$(range_end)"
        @info range_txt

        qcfuture = [SDMLayer(provider, futureclim, timespan=tsp; layer=i, resolution=2.5, bbox...) for i in eachindex(layers(provider))]
        qcfuture = [interpolate(qcf, first(envirovars)) for qcf in qcfuture] # Force compatibility
        mask!(qcfuture, QC)
        local proba = predict(sdm, qcfuture; threshold=false)
        local range = predict(sdm, qcfuture; threshold=true)
        sname = joinpath("rasters", fname * "_" * "$(ssp)_$(range_txt).tif")
        SimpleSDMLayers.save(sname, SDMLayer{Float64}[proba, range])
    end
end