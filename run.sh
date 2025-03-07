species=("Procyon lotor" "Mephitis mephitis")

for sp in "${species[@]}"; do
    julia -t 6 --project 01_sdm.jl $sp
    julia -t 6 --project 02_predictions.jl $sp
done