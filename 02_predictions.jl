using SpeciesDistributionToolkit
using Statistics
using CairoMakie

include("S1_theme.jl")

taxname = "Mephitis mephitis"
if ~isempty(ARGS)
    taxname = join(ARGS[1:2], " ")
end

@info "Running for $(taxname)"
taxcode = replace(taxname, " " => "_")
fslug = replace(last(splitpath(@__FILE__)), ".jl" => "")
figure_path = joinpath("figures", fslug, taxcode)
if !ispath(figure_path)
    mkpath(figure_path)
end

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

@info "Plot the partial responses"
lnames = layers(provider)
ldescr = layerdescriptions(provider)

for v in variables(sdm)
    @info ldescr[lnames[v]]

    fpath = joinpath(figure_path, taxcode * "_" * "partialresponse_$(lnames[v]).png")
    f = Figure(; size=(600, 600))
    ax = Axis(f[1, 1], xlabel=ldescr[lnames[v]], ylabel="Score for $(taxname)")
    for i in 1:100
        lines!(ax, partialresponse(sdm, v; inflated=true, threshold=false)..., color=:grey, alpha=0.4)
    end
    lines!(ax, partialresponse(sdm, v; threshold=false)..., color=:black, linewidth=2)
    tightlimits!(ax)
    CairoMakie.save(fpath, f)

    fpath = joinpath(figure_path, taxcode * "_" * "shapley_$(lnames[v]).png")
    f = Figure(; size=(600, 600))
    ax = Axis(f[1, 1], xlabel=ldescr[lnames[v]], ylabel="Score for $(taxname)")
    scatter!(ax, features(sdm, v), explain(sdm, v; threshold=false))
    tightlimits!(ax)
    CairoMakie.save(fpath, f)
end

@info "Plot the class assignment"
f = Figure(; size=(900, 450))
ax = Axis(f[1, 1], xlabel="Prediction", ylabel="Class", yticks=([0, 1], ["Absent", "Present"]))
col = Makie.wong_colors()[[i ? 1 : 2 for i in labels(sdm)]]
rainclouds!(ax, labels(sdm), predict(sdm; threshold=false), clouds=hist, orientation=:horizontal, color=col)
xlims!(ax, 0, 1)
vlines!([threshold(sdm)], color=:red)
tightlimits!(ax)
hidespines!(ax, :l, :r, :t)
current_figure()
CairoMakie.save(joinpath(figure_path, taxcode * "_" * "outputs.png"), f)

@info "Shapley values"
S = explain(sdm, envirovars; threshold=false)

@info "Save the Shapley values"
sname = joinpath("rasters", fslug, taxcode * "_" * "shapley.tif")
if !ispath(dirname(sname))
    mkpath(dirname(sname))
end
SimpleSDMLayers.save(sname, S)

@info "Save the range"
sname = joinpath("rasters", fslug, "$(taxcode).tif")
baseline_proba = predict(sdm, envirovars; threshold=false)
baseline_range = predict(sdm, envirovars; threshold=true)
SimpleSDMLayers.save(sname, SDMLayer{Float64}[baseline_proba, baseline_range])

SSPs = [SSP126, SSP245, SSP370, SSP585]
GCMs = [UKESM1_0_LL, MIROC6, ACCESS_CM2, CanESM5, EC_Earth3_Veg, MRI_ESM2_0, IPSL_CM6A_LR, CMCC_ESM2, CNRM_ESM2_1]
timespans = SimpleSDMDatasets.timespans(provider, Projection(SSPs[1], GCMs[1]))

function _predict(sdm, poly, template, provider, future, timespan; kwargs...)
    future = [SDMLayer(provider, future, timespan=timespan; resolution=2.5, layer=i, kwargs...) for i in eachindex(layers(provider))]
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
            future_fpath = joinpath("rasters", fslug, ssp, "$(range_begin)-$(range_end)")
            replace!(future_fpath, "figures" => "rasters")
            if !ispath(future_fpath)
                mkpath(future_fpath)
            end
            sname = joinpath(future_fpath, "$(taxcode)_$(gcm).tif")
            push!(scores, gcmscore)
        end
        # Save the projections for each GCM
        for i in eachindex(scores)

            SimpleSDMLayers.save(sname, scores[i])
        end
        # Average the predictions
        score = mosaic(mean, scores)
        range = score .>= threshold(sdm)
        future_fpath = joinpath("rasters", fslug, ssp, "$(range_begin)-$(range_end)")
        sname = joinpath(future_fpath, "$(taxcode).tif")
        SimpleSDMLayers.save(sname, SDMLayer{Float64}[score, range])
    end
end
