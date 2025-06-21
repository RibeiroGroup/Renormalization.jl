using LaTeXStrings
using Makie

export renormalized_dispersion, exciton_contents, combined

function renormalized_dispersion(; σratios=[0.1, 0.25, 0.5], show_error=false, qfilter=false, verbose=false)

    fig = Figure(size=(800, 400))
    gd = GridLayout(fig[1,1])
    ax1 = Axis(gd[1,1], xlabel=L"$q$ ($\mu$m$^{-1}$)", ylabel="Energy (eV)", title="(a) Perovskite exciton-polaritons", xlabelsize=18, xticks=1:2:9)
    ax2 = Axis(gd[1,2], title="(b) BODIPY exciton-BSW polaritons", xlabel=L"$q$ ($\mu$m$^{-1}$)", xlabelsize=18)

    xlims!(ax1, 0, 10)
    xlims!(ax2, 10, 12)

    ylims!(ax1, 1.4, 3.5)
    ylims!(ax2, 1.9, 2.25)

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

    renormalized_dispersion!(ax1, E0, Δ, qD, ecqD, σratios; show_error=show_error, qfilter=qfilter, verbose=verbose)

    ### BSW Data from Tau et al. ###

    # Bare exciton energy
    E0 = 2130 / 1000

    # Half of Rabi splitting
    Δ = 142.02 / 2000

    # Load bare cavity data
    lines = readlines("../BSW.txt")
    qvals = zeros(length(lines))
    Ecvals = zeros(length(lines))

    for i in eachindex(qvals)
        a, b = split(lines[i])
        qvals[i] = parse(Float64, a)
        Ecvals[i] = parse(Float64, b)
    end

    plots = renormalized_dispersion!(ax2, E0, Δ, qvals, moving_average(Ecvals, 11), σratios; show_error=show_error, qfilter=qfilter, verbose=verbose)

    # Legend(gd[1,3], reverse([plots]), [0.0, σratios...], title=L"\sigma/\Omega_R", merge=true)

    fig
end

function renormalized_dispersion!(ax, E0, Δ, qvals, Ecvals, σratios; show_error=false, qfilter=false, verbose=false) 

    # Plot bare exciton and cavity energies
    hlines!(ax, [E0], linestyle=:dash, color=:gray)
    lines!(ax, qvals, Ecvals, linestyle=:dash, color=:gray)

    # Compute zeroth-order disperson, e.g. zero disorder
    LP0 = [lp0_energy(E0, Ec, Δ) for Ec in Ecvals]
    UP0 = [up0_energy(E0, Ec, Δ) for Ec in Ecvals]

    p0 = lines!(ax, qvals, LP0, label=L"0", linestyle=:dot, color=:black)
    lines!(ax, qvals, UP0, linestyle=:dot, color=:black)

    # Compute polariton objects
    Pols = reverse([Polariton(E0, Δ, qvals, Ecvals, σ*2*Δ) for σ in σratios])
    colors = reverse(Makie.wong_colors()[1:length(Pols)])

    # Plot bands for uncertainties
    if show_error
        for (i, P) in enumerate(Pols)

            # Create mask to remove large δq values
            m_up = ones(Bool, length(P.qvals))
            m_lp = ones(Bool, length(P.qvals))
            if qfilter
                m_up = (P.δq_up ./ P.qvals) .< 1.0
                m_lp = (P.δq_lp ./ P.qvals) .< 1.0
            end

            # Plot uncertainty energies as bands
            lp_lower_limit = real.(P.Elp) .- imag.(P.Elp)
            lp_upper_limit = real.(P.Elp) .+ imag.(P.Elp)
            up_lower_limit = real.(P.Eup) .- imag.(P.Eup)
            up_upper_limit = real.(P.Eup) .+ imag.(P.Eup)       

            band!(ax, P.qvals[m_lp], lp_lower_limit[m_lp], lp_upper_limit[m_lp], color=colors[i], alpha=0.4)
            band!(ax, P.qvals[m_up], up_lower_limit[m_up], up_upper_limit[m_up], color=colors[i], alpha=0.4)

            if verbose
                println("σ/ΩR = $(reverse(σratios)[i])")
                println("Maximum renormalization EUP - Ep(0) = $(maximum(real.(P.Eup) .- UP0))")
                println("Maximum renormalization Ep(0) - ELP = $(maximum(LP0 .- real.(P.Elp)))")
            end

        end
    end

    plots = []

    # Plot renormalized energies
    for (i, P) in enumerate(Pols)
        # Create mask to remove large δq values
        m_up = ones(Bool, length(P.qvals))
        m_lp = ones(Bool, length(P.qvals))
        if qfilter
            m_up = (P.δq_up ./ P.qvals) .< 1.0
            m_lp = (P.δq_lp ./ P.qvals) .< 1.0
        end

        p = lines!(ax, P.qvals[m_up], real.(P.Eup[m_up]), label=L"%$(reverse(σratios)[i])", linewidth=2, color=colors[i])
        lines!(ax, P.qvals[m_lp], real.(P.Elp[m_lp]), label=L"%$(reverse(σratios)[i])", linewidth=2, color=colors[i])

        push!(plots, p)
    end

    push!(plots, p0)

    return plots 
end

function exciton_contents()
    fig = Figure(size=(800, 400))
    gd = GridLayout(fig[1,1])
    ax1 = Axis(gd[1,1], xlabel=L"$q$ ($\mu$m$^{-1}$)", ylabel="Exciton Content (%)", title="(a) Perovskite exciton-polaritons", xlabelsize=16, xticks=1:2:9, yticks=0:25:100)
    ax2 = Axis(gd[1,2], title="(b) BODIPY exciton-BSW polaritons", xlabel=L"$q$ ($\mu$m$^{-1}$)", xlabelsize=16, yticks=0:25:100)

    xlims!(ax1, 0, 10)
    xlims!(ax2, 10, 12)

    linkyaxes!(ax1, ax2)
    ylims!(ax2, 0, 100)

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

    P = Polariton(E0, Δ, qD, ecqD, 0.0)
    lines!(ax1, P.qvals, 100P.up_exciton_content, label="UP", color=Makie.wong_colors()[4])
    lines!(ax1, P.qvals, 100P.lp_exciton_content, label="LP", color=Makie.wong_colors()[5], linestyle=:dash)


    ### BSW Data from Tau et al. ###

    # Bare exciton energy
    E0 = 2130 / 1000

    # Half of Rabi splitting
    Δ = 142.02 / 2000

    # Load bare cavity data
    lines = readlines("../BSW.txt")
    qvals = zeros(length(lines))
    Ecvals = zeros(length(lines))

    for i in eachindex(qvals)
        a, b = split(lines[i])
        qvals[i] = parse(Float64, a)
        Ecvals[i] = parse(Float64, b)
    end

    P = Polariton(E0, Δ, qvals, Ecvals, 0.0)
    lines!(ax2, P.qvals, 100P.up_exciton_content, label="UP", color=Makie.wong_colors()[4])
    lines!(ax2, P.qvals, 100P.lp_exciton_content, label="LP", color=Makie.wong_colors()[5], linestyle=:dash)

    colgap!(gd, 1, 30)

    Legend(gd[1,3], ax1)

    fig
end

function combined(; σratios=[0.1, 0.25, 0.5], show_error=false, qfilter=false)

    fig = Figure(size=(800, 700))
    gd = GridLayout(fig[1,1])
    ax1 = Axis(gd[1,1], xlabel=L"$q$ ($\mu$m$^{-1}$)", ylabel="Energy (eV)", title="Perovskite exciton-polaritons", xlabelsize=18, xticks=1:2:9)
    ax2 = Axis(gd[1,2], title="BODIPY exciton-BSW polaritons", xlabel=L"$q$ ($\mu$m$^{-1}$)", xlabelsize=18)
    ax3 = Axis(gd[2,1], xlabel=L"$q$ ($\mu$m$^{-1}$)", ylabel="Exciton Content (%)", xlabelsize=18, xticks=1:2:9, yticks=0:25:100)
    ax4 = Axis(gd[2,2], xlabel=L"$q$ ($\mu$m$^{-1}$)", xlabelsize=18, yticks=0:25:100)

    xlims!(ax1, 0, 10)
    xlims!(ax3, 0, 10)
    xlims!(ax2, 10, 12)
    xlims!(ax4, 10, 12)

    ylims!(ax1, 1.4, 3.5)
    ylims!(ax2, 1.9, 2.25)

    ylims!(ax3, 0, 105)
    ylims!(ax4, 0, 105)

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

    renormalized_dispersion!(ax1, E0, Δ, qD, ecqD, σratios; show_error=show_error, qfilter=qfilter)

    # Exciton contents
    P = Polariton(E0, Δ, qD, ecqD, 0.0)
    lines!(ax3, P.qvals, 100P.up_exciton_content, label="UP", color=Makie.wong_colors()[4])
    lines!(ax3, P.qvals, 100P.lp_exciton_content, label="LP", color=Makie.wong_colors()[5], linestyle=:dash)

    ### BSW Data from Tau et al. ###

    # Bare exciton energy
    E0 = 2130 / 1000

    # Half of Rabi splitting
    Δ = 142.02 / 2000

    # Load bare cavity data
    lines = readlines("../BSW.txt")
    qvals = zeros(length(lines))
    Ecvals = zeros(length(lines))

    for i in eachindex(qvals)
        a, b = split(lines[i])
        qvals[i] = parse(Float64, a)
        Ecvals[i] = parse(Float64, b)
    end

    plots = renormalized_dispersion!(ax2, E0, Δ, qvals, moving_average(Ecvals, 11), σratios; show_error=show_error, qfilter=qfilter)

    # Exciton contents
    P = Polariton(E0, Δ, qvals, Ecvals, 0.0)
    lines!(ax4, P.qvals, 100P.up_exciton_content, label="UP", color=Makie.wong_colors()[4])
    lines!(ax4, P.qvals, 100P.lp_exciton_content, label="LP", color=Makie.wong_colors()[5], linestyle=:dash)

    Legend(gd[1,3], reverse(plots), string.([0.0, σratios...]), L"\sigma/\Omega_R", merge=true)
    Legend(gd[2,3], ax3)

    for (ax, letter) in zip([ax1, ax2, ax3, ax4], ['a', 'b', 'c', 'd'])
        text!(ax, 0.02, 0.99, text="($letter)", color=:black, fontsize=18, align=(:left, :top), space=:relative)
    end

    fig
end
