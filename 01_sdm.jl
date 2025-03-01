using SpeciesDistributionToolkit
using CairoMakie
import DelimitedFiles
using PrettyTables
import CSV

include("S1_theme.jl")

taxname = "Mephitis mephitis"
if ~isempty(ARGS)
    taxname = join(ARGS[1:2], " ")
end
@info "Running for $(taxname)"

𝐗 = DelimitedFiles.readdlm("data/$(replace(taxname, " " => "_")).X.dat", Float32)
𝐲 = vec(DelimitedFiles.readdlm("data/$(replace(taxname, " " => "_")).y.dat", Bool))

@info "Train the SDM for all the known data"
sdm = SDM(ZScore, NaiveBayes, 𝐗, 𝐲)
folds = kfold(sdm)

# Set some better training parameters
# classifier(sdm).verbose = false
# classifier(sdm).η = 1e-3
# classifier(sdm).epochs = 5000

@info "Select variables"
forwardselection!(sdm, folds; verbose=true)

@info "Report variables"
DelimitedFiles.writedlm("data/$(replace(taxname, " " => "_")).params", variables(sdm))
DelimitedFiles.writedlm("data/$(replace(taxname, " " => "_")).threshold", threshold(sdm))

@info "Save the trained model"
SDeMo.writesdm("models/$(replace(taxname, " " => "_")).json", sdm)

@info "Report on cross-validation"
cv = crossvalidate(sdm, folds)

measures = [mcc, ppv, npv, trueskill, markedness, plr, nlr]
M = permutedims([measure(c) for measure in measures, c in cv])
pt = pretty_table(M; header = measures)
open("data/$(replace(taxname, " " => "_")).crossvalidation", "w") do f
    pretty_table(f, M, header=measures)
end

@info "Learning curves"
τ = threshold(sdm)
thrs = LinRange(0.0, 1.0, 90)
C = zeros(ConfusionMatrix, length(thrs))
for (i, thr) in enumerate(thrs)
    threshold!(sdm, thr)
    C[i] =  ConfusionMatrix(predict(sdm), labels(sdm))
end
threshold!(sdm, τ)

f = Figure(; size=(900, 300))
ax1 = Axis(f[1, 1]; aspect = 1, xlabel = "Threshold", ylabel = "MCC")
lines!(ax1, thrs, mcc.(C))
ax2 = Axis(f[1, 2]; aspect = 1, xlabel = "Precision", ylabel = "Recall")
lines!(ax2, SDeMo.precision.(C), SDeMo.recall.(C))
ax3 = Axis(f[1, 3]; aspect = 1, xlabel = "TPR", ylabel = "FPR")
lines!(ax3, fpr.(C), tpr.(C))
for ax in [ax1, ax2, ax3]
    xlims!(ax, 0, 1)
    ylims!(ax, 0, 1)
    tightlimits!(ax)
end
current_figure()
CairoMakie.save("figures/01_sdm/$(replace(taxname, " " => "_"))-crossvalidation.png", current_figure())