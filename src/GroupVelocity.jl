using Statistics

export group_velocities, renorm_heatmap

function moving_average(data, window_size=5)
    n = length(data)
    smoothed = similar(data)
    smoothed .= data
    half_win = window_size ÷ 2

    for i in 5:(n-5)
        lo = max(1, i - half_win)
        hi = min(n, i + half_win)
        smoothed[i] = mean(data[lo:hi])
    end
    return smoothed
end

function group_velocities(;σratios=[0.1, 0.25, 0.5], qfilter=false)
    fig = Figure(size=(800,400))
    ax1 = Axis(fig[1,1], xlabel=L"$q$ ($\mu$m$^{-1}$)", ylabel=L"Group Velocity ($\mu$m$\cdot$ps$^{-1}$)", 
    title="Perovskite exciton-polaritons", xlabelsize=18, xticks=1:2:9, ylabelsize=18)
    ax2 = Axis(fig[1,2], title="BODIPY exciton-BSW polaritons", xlabel=L"$q$ ($\mu$m$^{-1}$)", xlabelsize=18)

    ### Parameters from Nat. Commun. (Xu et al. 2023) ###

    # Bare exciton energy
    E0 = 2140.0 / 1000

    # Half of Rabi splitting
    Δ = 275.0 / 1000

    # Cavity fitting parameters
    EcD = 1570 / 1000
    LcD = 0.667
    # Cavity q values (in cm⁻¹)
    qD = range(start=0, stop=10, length=2000) |> collect
    # Bare cavity energies
    ecqD = EcD .* sqrt.(1 .+ ((qD .* LcD) / π).^2)

    # Find resonant q
    _, idx = findmin(x->abs(x-E0), ecqD)

    vlines!(ax1, [qD[idx]], linestyle=:dash, color=:gray)

    group_velocities!(ax1, E0, Δ, qD, ecqD, σratios, qfilter=qfilter)

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

    group_velocities!(ax2, E0, Δ, qvals, moving_average(Ecvals, 11), σratios, qfilter=qfilter)

    # Find resonant q
    _, idx = findmin(x->abs(x-E0), Ecvals)

    vlines!(ax2, [qvals[idx]], linestyle=:dash, color=:gray)

    xlims!(ax1, 0.0, 10)
    xlims!(ax2, 10, 12)

    ylims!(ax1, 0, 210)
    ylims!(ax2, 25, 200)

    axislegend(ax2, L"\sigma/\Omega_R", position=:rt)

    for (ax, letter) in zip([ax1, ax2], ['a', 'b'])
        text!(ax, 0.02, 0.99, text="($letter)", color=:black, fontsize=18, align=(:left, :top), space=:relative)
    end

    fig
end

function group_velocities!(ax, E0, Δ, qvals, Ecvals, σratios; qfilter=false)

    # Given energy in eV and q in cm⁻¹, we need a conversion factor
    # eV -> s-1 using ħ = h/2π = 4.1356677e-15 eV⋅s / 2π
    ħ = 0.0041356677 / 2π

    # Create polariton object
    P0 = Polariton(E0, Δ, qvals, Ecvals, 0.0)

    # Plot zero disorder values
    lines!(ax, P0.qvals, P0.vg_lp, label=L"0", linestyle=:dot, color=:black)

    for σ in σratios
        P = Polariton(E0, Δ, qvals, Ecvals, 2*σ*Δ)

        m = ones(Bool, length(P.qvals))
        if qfilter
            m = (P.δq_lp ./ P.qvals) .< 1.0
        end
    
        lines!(ax, P.qvals[m], P.vg_lp[m], label=L"%$σ")
    end
end

function renormalized_vs_exc(; σratios=[0.1, 0.25, 0.5], qfilter=false)
    fig = Figure(size=(800,400))
    ax1 = Axis(fig[1,1], xlabel="Exciton Content (%)", ylabel="Renormalization (%)", title="(a) Perovskite exciton-polaritons", xlabelsize=16, xticks=0:25:100)
    ax2 = Axis(fig[1,2], title="(b) BODIPY exciton-BSW polaritons", xlabel="Exciton Content (%)", xlabelsize=16)

    ### Parameters from Nat. Commun. (Xu et al. 2023) ###

    # Bare exciton energy
    E0 = 2140.0 / 1000

    # Half of Rabi splitting
    Δ = 275.0 / 1000

    # Cavity fitting parameters
    EcD = 1570 / 1000
    LcD = 0.667
    # Cavity q values (in cm⁻¹)
    qD = range(start=0, stop=10, length=2000) |> collect
    # Bare cavity energies
    ecqD = EcD .* sqrt.(1 .+ ((qD .* LcD) / π).^2)

    # Create polariton object for the zero disorder
    P0 = Polariton(E0, Δ, qD, ecqD, 0.0)

    for σ in σratios

        P = Polariton(E0, Δ, qD, ecqD, 2*σ*Δ)

        m = ones(Bool, length(P.qvals))
        if qfilter
            m = (P.δq_lp ./ P.qvals) .< 1.0
        end

        lines!(ax1, 100*P.lp_exciton_content[m], 100 .* (P.vg_lp[m] ./ P0.vg_lp[m]), label=L"%$σ")
    end

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

    # Create polariton object for the zero disorder
    P0 = Polariton(E0, Δ, qvals, Ecvals, 0.0)

    for σ in σratios

        P = Polariton(E0, Δ, qvals, Ecvals, 2*σ*Δ)

        m = ones(Bool, length(P0.qvals)) 
        if qfilter
            m = (P.δq_lp ./ P.qvals) .< 1.0
        end

        lines!(ax2, 100*P.lp_exciton_content, 100 .* (P.vg_lp ./ P0.vg_lp), label=L"%$σ")
    end

    axislegend(ax2, L"\sigma/\Omega_R", position=:lb)

    fig
end

function renormalized_vs_σ(exc_vals=[10, 25, 50, 75], σratios=0.1:0.1:0.8)
    fig = Figure(size=(800,400))
    ax1 = Axis(fig[1,1], xlabel=L"Disorder Strength ($\sigma/\Omega_R$)", ylabel="Renormalization (%)", title="(a) Perovskite exciton-polaritons", xlabelsize=16, xticks=0:0.2:1)
    ax2 = Axis(fig[1,2], title="(b) BODIPY exciton-BSW polaritons", xlabel=L"Disorder Strength ($\sigma/\Omega_R$)", xlabelsize=16)

    ### Parameters from Nat. Commun. (Xu et al. 2023) ###

    # Bare exciton energy
    E0 = 2140.0 / 1000

    # Half of Rabi splitting
    Δ = 275.0 / 1000

    # Cavity fitting parameters
    EcD = 1570 / 1000
    LcD = 0.667
    # Cavity q values (in cm⁻¹)
    qD = range(start=0, stop=10, length=2000) |> collect
    # Bare cavity energies
    ecqD = EcD .* sqrt.(1 .+ ((qD .* LcD) / π).^2)

    renormalized_vs_σ!(ax1, E0, Δ, ecqD, qD, σratios, exc_vals)

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

    renormalized_vs_σ!(ax2, E0, Δ, Ecvals, qvals, σratios, exc_vals)

    axislegend(ax2, "Exciton Content (%)", position=:lb, merge=true)

    fig
end

function renormalized_vs_σ!(ax, E0, Δ, Ecvals, qvals, σratios, exc_vals)

    # Compute zeroth-order disperson, e.g. zero disorder
    P0 = Polariton(E0, Δ, qvals, Ecvals, 0.0)

    LP_vgs = zeros(length(σratios), length(exc_vals))

    for (n,σ) in enumerate(σratios)

        P = Polariton(E0, Δ, qvals, Ecvals, 2*σ*Δ)

        for (i, exc) in enumerate(exc_vals)

            # Find the exciton content for the given exciton energy
            _, idx = findmin(x->abs(100*x-exc), P.lp_exciton_content)

            # Compute the group velocity at this exciton content
            LP_vgs[n, i] = P.vg_lp[idx] #/ P0.vg_lp[idx]
        end
    end

    for (i, exc) in enumerate(exc_vals)
        scatter!(ax, σratios, 100 .* LP_vgs[:, i], label=L"%$exc")
        lines!(ax, σratios, 100 .* LP_vgs[:, i], label=L"%$exc")
    end
end

function renorm_heatmap()
    fig = Figure(size=(800,400))
    ax1 = Axis(fig[1,1], xlabel=L"Disorder Strength ($\sigma/\Omega_R$)", ylabel="Exciton Content (%)", title="(a) Perovskite exciton-polaritons", xlabelsize=16, xticks=0:0.2:1)
    ax2 = Axis(fig[1,3], title="(b) BODIPY exciton-BSW polaritons", xlabel=L"Disorder Strength ($\sigma/\Omega_R$)", xlabelsize=16)

    ### Parameters from Nat. Commun. (Xu et al. 2023) ###

    # Bare exciton energy
    E0 = 2140.0 / 1000

    # Half of Rabi splitting
    Δ = 275.0 / 1000

    # Cavity fitting parameters
    EcD = 1570 / 1000
    LcD = 0.667
    # Cavity q values (in cm⁻¹)
    qD = range(start=0, stop=10, length=2000) |> collect
    # Bare cavity energies
    ecqD = EcD .* sqrt.(1 .+ ((qD .* LcD) / π).^2)

    hm1 = renorm_heatmap!(ax1, E0, Δ, ecqD, qD)

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

    hm2 = renorm_heatmap!(ax2, E0, Δ, Ecvals, qvals)

     Colorbar(fig[1,2], hm1, label = "Renormalization")  # Add colorbar for first heatmap
    Colorbar(fig[1,4], hm2, label = "Renormalization")  # Add colorbar for second heatmap

    # # Find global min/max for colorrange
    # vmin = min(minimum(hm1[1][]), minimum(hm2[1][]))
    # vmax = max(maximum(hm1[1][]), maximum(hm2[1][]))

    # # Set colorrange for both heatmaps
    # hm1.colorrange[] = (vmin, vmax)
    # hm2.colorrange[] = (vmin, vmax)

    # # Add a single colorbar for both
    # Colorbar(fig[:,3], hm1, label = "Renormalization")

    fig
end

function renorm_heatmap!(ax, E0, Δ, Ecvals, qvals, σratios=0.1:0.01:0.70)

    # Compute zeroth-order disperson, e.g. zero disorder
    P0 = Polariton(E0, Δ, qvals, Ecvals, 0.0)

    LP_vgs = zeros(length(σratios), length(P0.qvals))

    boundary_ex = zeros(length(σratios))
    boundary_δq = zeros(length(σratios))

    lpguess = P0.Elp

    for (n,σ) in enumerate(σratios)
        P = Polariton(E0, Δ, qvals, Ecvals, 2*σ*Δ, lpguess=lpguess)
        # Compute the group velocity at this exciton content
        renorm = P.vg_lp ./ P0.vg_lp

        LP_vgs[n, :] = renorm

        idx = findfirst(y -> y < 0.9, renorm)
        if idx === nothing
            boundary_ex[n] = NaN  # No clear boundary found
        else
            boundary_ex[n] = 100*P.lp_exciton_content[idx]
        end

        idx2 = findlast(y -> y < 1.0, P.δq_lp ./ P.qvals)
        if idx2 === nothing
            boundary_δq[n] = NaN  # No clear boundary found
        else
            boundary_δq[n] = 100*P.lp_exciton_content[idx2]
        end

        lpguess = P.Elp  # Update guess for next iteration
    end

    # Prepare axes: x = σratios, y = exciton content (in %)
    x = collect(σratios)
    y = 100 .* P0.lp_exciton_content

    # Makie expects the matrix as (y, x), so transpose if needed
    hm = heatmap!(ax, x, y, LP_vgs; colormap=:viridis)

    lines!(ax, x, boundary_ex, color=:red, linestyle=:dash, label="90% Exciton Content Boundary")
    lines!(ax, x, boundary_δq, color=:blue, linestyle=:dash, label="1.0 δq/q Boundary")

    # ax.xlabel = "Disorder Strength (σ/Ωᵣ)"
    # ax.ylabel = "Exciton Content (%)"
    # ax.title = "Renormalization Heatmap"
    # ax.xticks = x
    # ax.yticks = range(round(minimum(y)), round(maximum(y)), length=5)

    return hm
end

function renorm_vs_product(; σratios=0.1:0.1:0.8, qfilter=false)

    fig = Figure(size=(800,400))
    ax1 = Axis(fig[1,1], xlabel=L"(Exciton Content) × σ/ΩR", ylabel="Renormalization (%)", title="(a) Perovskite exciton-polaritons", xlabelsize=16, xticks=0:25:100)
    ax2 = Axis(fig[1,2], title="(b) BODIPY exciton-BSW polaritons", xlabel=L"(Exciton Content) × σ/ΩR", xlabelsize=16)

    ### Parameters from Nat. Commun. (Xu et al. 2023) ###

    # Bare exciton energy
    E0 = 2140.0 / 1000

    # Half of Rabi splitting
    Δ = 275.0 / 1000

    # Cavity fitting parameters
    EcD = 1570 / 1000
    LcD = 0.667
    # Cavity q values (in cm⁻¹)
    qD = range(start=0, stop=10, length=2000) |> collect
    # Bare cavity energies
    ecqD = EcD .* sqrt.(1 .+ ((qD .* LcD) / π).^2)

    renorm_vs_product!(ax1, E0, Δ, ecqD, qD, σratios)

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

    renorm_vs_product!(ax2, E0, Δ, Ecvals, qvals, σratios)

    axislegend(ax2, L"\sigma/\Omega_R", position=:lb)

    fig
end

function renorm_vs_product!(ax, E0, Δ, Ecvals, qvals, σratios)

    # Compute zeroth-order disperson, e.g. zero disorder
    P0 = Polariton(E0, Δ, qvals, Ecvals, 0.0)

    lpguess = P0.Elp

    for (n,σ) in enumerate(σratios)

        P = Polariton(E0, Δ, qvals, Ecvals, 2*σ*Δ, lpguess=lpguess)

        # Compute the group velocity at this exciton content
        renorm = P.vg_lp ./ P0.vg_lp


        m = (P.δq_lp ./ P.qvals) .< 1.0

        scatter!(ax, σ .* P.lp_exciton_content[m], 100 .* renorm[m], label=L"%$σ")

        lpguess = P.Elp  # Update guess for next iteration
    end
end