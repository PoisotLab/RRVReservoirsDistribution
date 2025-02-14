using SpeciesDistributionToolkit
using CairoMakie
import DelimitedFiles
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

forwardselection!(sdm, folds; verbose=true)

cv = crossvalidate(sdm, folds)
mcc(cv.validation)
mcc(cv.training)

sdm.classifier.verbose = false
sdm.classifier.η = 1e-3
sdm.classifier.epochs = 12_000
train!(sdm)

@info "Loading bioclim data for prediction"
provider = RasterData(WorldClim2, BioClim)
QC = SpeciesDistributionToolkit.gadm("CAN", "Québec")
bbox = SpeciesDistributionToolkit.boundingbox(QC; padding=0.0)
envirovars = [SDMLayer(provider; layer=i, resolution=5.0, bbox...) for i in eachindex(layers(provider))]
mask!(envirovars, QC)

@info "Baseline prediction"
predict(sdm, envirovars; threshold=false) |> heatmap
predict(sdm, envirovars; threshold=true) |> heatmap

variableimportance(sdm, folds)

heatmap(explain(sdm, envirovars, 1; threshold=false))

lines(partialresponse(sdm, 7; threshold=false)...)
lines(partialresponse(sdm, 1; threshold=false)...)