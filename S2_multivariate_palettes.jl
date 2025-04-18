using ColorBlendModes

function _vsup_grid(vbins, ubins, vpal; upal = colorant"#efefef55", s = 0.5, k = 1.0)
    pal = fill(upal, (vbins, ubins))
    for i in 1:ubins
        shrkfac = ((i - 1) / (ubins - 1))^k
        subst = 0.5 - shrkfac * s / 2
        pal[:, i] .= CairoMakie.cgrad(vpal)[LinRange(0.5 - subst, 0.5 + subst, vbins)]
        # Apply the mix to the uncertain color
        for j in 1:vbins
            pal[j, i] = ColorSchemes.weighted_color_mean(1 - shrkfac, pal[j, i], upal)
        end
    end
    return pal
end

function discretize(layer, n::Integer)
    categories = rescale(layer, 0.0, 1.0)
    n = n - 2
    map!(x -> round(x * (n + 1); digits = 0) / (n + 1), categories.grid, categories.grid)
    return n * categories
end


function _get_bivariate_colormap(;
    n_stops=3,
    p0 = colorant"#e8e8e8",
    p1 = colorant"#64acbe",
    p2 = colorant"#c85a5a",
    rev = false
)
    # Generate colormap
    cm1 = LinRange(p0, p1, n_stops)
    cm2 = LinRange(p0, p2, n_stops)
    if rev
        cmat = ColorBlendModes.BlendMultiply.(cm2, cm1')
    else
        cmat = ColorBlendModes.BlendMultiply.(cm1, cm2')
    end
    return cmat
end