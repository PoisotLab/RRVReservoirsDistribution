using SpeciesDistributionToolkit
using Statistics
using CairoMakie

include("S1_theme.jl")

taxname = "Mephitis mephitis"
if ~isempty(ARGS)
    taxname = join(ARGS[1:2], " ")
end
@info "Running for $(taxname)"

@info "Loading the SDM"
sdm = SDeMo.loadsdm("models/$(replace(taxname, " " => "_")).json")

@info "Loading bioclim data for prediction"
provider = RasterData(WorldClim2, BioClim)
QC = SpeciesDistributionToolkit.gadm("CAN", "Québec")
bbox = SpeciesDistributionToolkit.boundingbox(QC; padding=0.0)
envirovars = SDMLayer{Float32}[SDMLayer(provider; layer=i, resolution=2.5, bbox...) for i in eachindex(layers(provider))]
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
    fpath = joinpath("figures", "02_predictions", fname * "_" * "partialresponse_$(lnames[v]).png")
    f = Figure(; size=(600, 600))
    ax = Axis(f[1, 1], xlabel=ldescr[lnames[v]], ylabel="Score for $(taxname)")
    for i in 1:100
        lines!(ax, partialresponse(sdm, v; inflated=true, threshold=false)..., color=:grey, alpha=0.4)
    end
    lines!(ax, partialresponse(sdm, v; threshold=false)..., color=:black, linewidth=2)
    tightlimits!(ax)
    CairoMakie.save(fpath, f)
    fpath = joinpath("figures", "02_predictions", fname * "_" * "shapley_$(lnames[v]).png")
    f = Figure(; size=(600, 600))
    ax = Axis(f[1, 1], xlabel=ldescr[lnames[v]], ylabel="Score for $(taxname)")
    scatter!(ax, features(sdm, v), explain(sdm, v; threshold=false))
    tightlimits!(ax)
    CairoMakie.save(fpath, f)
end

@info "Plot the class assignment"
fpath = joinpath("figures", "02_predictions", fname * "_" * "outputs.png")
f = Figure(; size=(900, 450))
ax = Axis(f[1, 1], xlabel="Prediction", ylabel="Class")
boxplot!(ax, labels(sdm), predict(sdm; threshold=false), orientation=:horizontal, outliercolor=:black, color=:white, strokecolor=:black, strokewidth=2)
ylims!(ax, -0.5, 1.5)
xlims!(ax, 0, 1)
vlines!([threshold(sdm)], color=:red)
tightlimits!(ax)
current_figure()
CairoMakie.save(fpath, f)

@info "Shapley values"
S = explain(sdm, envirovars; threshold=false)

@info "Save the Shapley values"
sname = joinpath("rasters", fname * "_" * "shapley.tif")
SimpleSDMLayers.save(sname, S)

@info "Save the range"
sname = joinpath("rasters", fname * "_" * "historical.tif")
baseline_proba = predict(sdm, envirovars; threshold=false)
baseline_range = predict(sdm, envirovars; threshold=true)
SimpleSDMLayers.save(sname, SDMLayer{Float64}[baseline_proba, baseline_range])

SSPs = [SSP126, SSP245, SSP370, SSP585]
GCMs = [UKESM1_0_LL, MIROC6, ACCESS_CM2, CanESM5, EC_Earth3_Veg, MRI_ESM2_0, IPSL_CM6A_LR, CMCC_ESM2, CNRM_ESM2_1]
timespans = SimpleSDMDatasets.timespans(provider, Projection(SSPs[1], GCMs[1]))

function _predict(sdm, poly, template, provider, future, timespan; kwargs...)
    future = [SDMLayer(provider, future, timespan=timespan; layer=i, kwargs...) for i in eachindex(layers(provider))]
    future = [interpolate(fl, template) for fl in future]
    mask!(future, poly)
    return predict(sdm, future; threshold=false)
end

@info "Forecast - GCM averaging"
for ssp in SSPs
    for timespan in timespans
        range_begin = timespan.first.value
        range_end = timespan.second.value
        scores = SDMLayer{Float64}[]
        # Run for several GCMs
        for gcm in GCMs
            futureclim = Projection(ssp, gcm)
            @info ssp, "$(range_begin) → $(range_end)", gcm
            gcmscore = _predict(sdm, QC, baseline_proba, provider, Projection(ssp, gcm), timespan; resolution=2.5, bbox...)
            push!(scores, gcmscore)
        end
        # Average the predictions
        score = mosaic(mean, scores)
        range = score .>= threshold(sdm)
        sname = joinpath("rasters", fname * "_" * "$(ssp)_$(range_begin)-$(range_end).tif")
        SimpleSDMLayers.save(sname, SDMLayer{Float64}[score, range])
    end
end
