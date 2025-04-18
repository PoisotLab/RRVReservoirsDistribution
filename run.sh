species=("Procyon lotor" "Mephitis mephitis")

for sp in "${species[@]}"; do
    julia -t auto --project 01_sdm.jl $sp
    julia -t auto --project 02_predictions.jl $sp
done