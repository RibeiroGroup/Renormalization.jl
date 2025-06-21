export quncertainty

function quncertainty()

    fig = Figure(size=(800, 400))
    gd = GridLayout(fig[1,1])
    ax1 = Axis(gd[1,1], xlabel=L"$q$ ($\mu$m$^{-1}$)", ylabel="Energy (eV)", title="(a) Perovskite exciton-polaritons", xlabelsize=16, xticks=1:2:9)
    # ax2 = Axis(gd[1,2], title="(b) BODIPY exciton-BSW polaritons", xlabel=L"$q$ ($\mu$m$^{-1}$)", xlabelsize=16)

    # xlims!(ax1, 0, 10)
    # xlims!(ax2, 10, 12)

    # ylims!(ax1, 1.4, 3.5)
    # ylims!(ax2, 1.9, 2.25)

    ### Parameters from Nat. Commun. (Xu et al. 2023) ###

    # Bare exciton energy
    E0 = 2140.0 / 1000

    # Half of Rabi splitting
    Δ = 275.0 / 1000

    # Cavity fitting parameters
    EcD = 1570 / 1000
    LcD = 0.667 # μm
    # Cavity q values
    qD = range(start=0, stop=10, length=2000) |> collect
    # Bare cavity energies
    ecqD = EcD .* sqrt.(1 .+ ((qD .* LcD) / π).^2)

    # Ratios σ/ΩR = 2σ/Δ
    σratios = [0.1, 0.25, 0.5]

    quncertainty!(ax1, E0, Δ, ecqD, qD, σratios)

    axislegend(ax1)

    fig
end

function quncertainty!(ax, E0, Δ, Ecvals, qvals, σratios)

    # Create polariton object
    P = Polariton(E0, Δ, qvals, Ecvals)

    # Compute zeroth-order disperson, e.g. zero disorder
    LP0 = [lp0_energy(P, q) for q in P.qvals]
    UP0 = [up0_energy(P, q) for q in P.qvals]

    n = 1
    for σ in σratios
        lp = zeros(ComplexF64, length(P.qvals))
        up = zeros(ComplexF64, length(P.qvals))
    
        for i in eachindex(lp)
            lp[i] = find_renormalized_energy(P, P.qvals[i], 2*σ*Δ, LP0[i]) 
            up[i] = find_renormalized_energy(P, P.qvals[i], 2*σ*Δ, UP0[i]) 
        end

        dLP_dk = numerical_derivative_three_point(P.qvals, real.(lp))
        dUP_dk = numerical_derivative_three_point(P.qvals, real.(up))

        Δk_up = -imag.(up) ./ dUP_dk
        Δk_lp = -imag.(lp) ./ dLP_dk

        lines!(ax, P.qvals, Δk_lp ./ P.qvals, label=L"%$σ", color=Makie.wong_colors()[n])
        lines!(ax, P.qvals, Δk_up ./ P.qvals, color=Makie.wong_colors()[n], linestyle=:dash, label=L"%$σ")

        hlines!(ax, [1.0], linestyle=:dot, color=:gray)
        n += 1
    end
end

function find_boundaries(E0, Δ, Ecvals, qvals, σ; branch="UP", tol=1e-10)

    P = Polariton(E0, Δ, qvals, Ecvals)

    # Compute zeroth-order disperson, e.g. zero disorder
    LP0 = [lp0_energy(P, q) for q in P.qvals]
    UP0 = [up0_energy(P, q) for q in P.qvals]

    lp = [find_renormalized_energy(P, q, 2*σ*Δ, LP0[i]) for (i, q) in enumerate(P.qvals)]
    up = [find_renormalized_energy(P, q, 2*σ*Δ, UP0[i]) for (i, q) in enumerate(P.qvals)]

    dLP_dk = numerical_derivative_three_point(P.qvals, real.(lp))
    dUP_dk = numerical_derivative_three_point(P.qvals, real.(up))

    Δk_up = -imag.(up) ./ dUP_dk
    Δk_lp = -imag.(lp) ./ dLP_dk

    mask_up = (Δk_up ./ P.qvals) .< 1.0
    mask_lp = (Δk_lp ./ P.qvals) .< 1.0

    valid_q_up = P.qvals[mask_up]
    valid_q_lp = P.qvals[mask_lp]

    return [(valid_q_lp[1], valid_q_lp[end]), (valid_q_up[1], valid_q_up[end])]

end