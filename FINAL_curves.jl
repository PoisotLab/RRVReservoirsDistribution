using SpeciesDistributionToolkit
const SDT = SpeciesDistributionToolkit
using CairoMakie
import DelimitedFiles
using PrettyTables
using Statistics
import CSV

include("S1_theme.jl")

model_prolot = SDeMo.loadsdm("models/Procyon_lotor.json")
model_mepmep = SDeMo.loadsdm("models/Mephitis_mephitis.json")

# Plot the curve
function _forcurve(sdm; n=90)
    τ = threshold(sdm)
    thrs = LinRange(0.0, 1.0, n)
    C = zeros(ConfusionMatrix, length(thrs))
    for (i, thr) in enumerate(thrs)
        threshold!(sdm, thr)
        C[i] = ConfusionMatrix(predict(sdm), labels(sdm))
    end
    threshold!(sdm, τ)
    return thrs, C
end

# Confusion matrices
px, py = _forcurve(model_prolot)
mx, my = _forcurve(model_mepmep)

# Figures
f = Figure(; size=(500, 500))

# Axes
roc = Axis(f[1, 1]; xlabel="False Positive Rate", ylabel="True Positive Rate", aspect=DataAspect())
pr = Axis(f[1, 2]; xlabel="Precision", ylabel="Recall", aspect=DataAspect())
thr = Axis(f[2, 1], xlabel="Threshold", ylabel="MCC", aspect=DataAspect())

# Cleanup axes
for ax in [roc, pr, thr]
    tightlimits!(ax)
    xlims!(ax, 0, 1)
    ylims!(ax, 0, 1)
end

# Silhouettes - legend
sp_uuid = Phylopic.imagesof("Procyon lotor"; items=1)
silhouetteplot!(roc, -1.0, -1.0, sp_uuid; markersize=18, label="P. lotor", color=:black)
sp_uuid = Phylopic.imagesof("Mephitis mephitis"; items=1)
silhouetteplot!(roc, -1.0, -1.0, sp_uuid; markersize=18, label="M. mephitis", color=:orange)


# Legend
Legend(f[2, 2], roc, tellheight=false, tellwidth=false, framevisible=false)

# ROC curve
lines!(roc, fpr.(py), tpr.(py), color=:black)
lines!(roc, fpr.(my), tpr.(my), color=:orange)

# PR curve
lines!(pr, SDeMo.precision.(py), SDeMo.recall.(py), label="P. lotor", color=:black)
lines!(pr, SDeMo.precision.(my), SDeMo.recall.(my), label="M. mephitis", color=:orange)

# Learning curve
lines!(thr, px, mcc.(py), label="P. lotor", color=:black)
lines!(thr, mx, mcc.(my), label="M. mephitis", color=:orange)

# Show figure
f
CairoMakie.save("figures/SM_validation.png", f)
