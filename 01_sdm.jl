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

sdm = SDM(ZScore, Logistic, 𝐗, 𝐲)
folds = kfold(sdm)

train!(sdm)
cv = crossvalidate(sdm, folds)
