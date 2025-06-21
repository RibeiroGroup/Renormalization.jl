module Renormalization
using SpecialFunctions
using DataInterpolations
using NLsolve

export Polariton
export cavity_energy, lp0_energy, up0_energy, up0_lp0_energy
export find_renormalized_energy

struct Polariton{T, C}
    E0::T
    Δ::T 
    σ::T
    ΩR::T
    qvals::Vector{T}
    Ecvals::Vector{T}
    Eup::Vector{C}
    Elp::Vector{C}
    up_exciton_content::Vector{T}
    lp_exciton_content::Vector{T}
    vg_up::Vector{T}
    vg_lp::Vector{T}
    δq_up::Vector{T}
    δq_lp::Vector{T}
end

function Polariton(E0, Δ, qvals, Ecvals, σ; ħ=0.0041356677/2π, lpguess=nothing, upguess=nothing)

    ΩR = 2*Δ

    Eup = zeros(ComplexF64, length(qvals))
    Elp = zeros(ComplexF64, length(qvals))
    lp_exciton_content = zeros(Float64, length(qvals))
    up_exciton_content = zeros(Float64, length(qvals))

    # Compute properties for each q mode
    for i in eachindex(qvals)

        # Compute renormalize UP and LP energies
        if lpguess !== nothing
            Elp[i] = find_renormalized_energy(E0, Ecvals[i], Δ, σ; branch="LP", guess=lpguess[i])
        else
            Elp[i] = find_renormalized_energy(E0, Ecvals[i], Δ, σ; branch="LP")
        end

        if upguess !== nothing
            Eup[i] = find_renormalized_energy(E0, Ecvals[i], Δ, σ; branch="UP", guess=upguess[i])
        else
            Eup[i] = find_renormalized_energy(E0, Ecvals[i], Δ, σ; branch="UP")
        end

        # Exciton content
        up_exciton_content[i] = exciton_content(E0, Ecvals[i], Δ; branch="UP")
        lp_exciton_content[i] = exciton_content(E0, Ecvals[i], Δ; branch="LP")
    end

    # Compute group_velocities
    vg_up = numerical_derivative_three_point(qvals, real.(Eup)) ./ ħ
    vg_lp = numerical_derivative_three_point(qvals, real.(Elp)) ./ ħ

    # Compute q uncertainties
    δq_up = -imag.(Eup) ./ (ħ * vg_up)
    δq_lp = -imag.(Elp) ./ (ħ * vg_lp)

    return Polariton(E0, Δ, σ, ΩR, qvals, Ecvals, Eup, Elp,
        up_exciton_content, lp_exciton_content,
        vg_up, vg_lp, δq_up, δq_lp)
end

function cavity_energy(Ecvals, qvals, q)

    idx = findfirst(isequal(q), qvals)

    if idx !== nothing
        return Ecvals[idx]
    end

    # If q is not in qvals, we need to interpolate
    itp = CubicSpline(Ecvals, qvals)

    return itp(q)
end

function lp0_energy(E0, Ec, Δ)
    return 0.5*(E0 + Ec) - sqrt(Δ^2 + (0.5*(E0 - Ec))^2)
end

function up0_energy(E0, Ec, Δ)
    return 0.5*(E0 + Ec) + sqrt(Δ^2 + (0.5*(E0 - Ec))^2)
end

function up0_lp0_energy(E0, Ec, Δ)
    A = 0.5*(E0 + Ec) 
    B = sqrt(Δ^2 + (0.5*(E0 - Ec))^2)

    return A + B, A - B
end

function photon_content(E0, Ec, Δ; branch="UP")

    ΩR = 2*Δ

    ω = 0.0

    if uppercase(string(branch)) == "UP"
        ω = up0_energy(E0, Ec, Δ)
    elseif uppercase(string(branch)) == "LP"
        ω = lp0_energy(E0, Ec, Δ)
    else
        error("Invalid branch option $up")
    end

    ωM = E0
    ωc = Ec

    Lq = ωM * ΩR^2 /(4ωc)

    return Lq / (Lq + (ω - ωc)^2)
end

function exciton_content(E0, Ec, Δ; branch="UP")
    return 1 - photon_content(E0, Ec, Δ; branch=branch)
end

### IMPLEMENT SOLVER T
function find_renormalized_energy(E0, Ec, Δ, σ; branch="UP", guess=nothing)

    if σ == 0.0
        # If σ is zero, return the zeroth-order energy
        if uppercase(string(branch)) == "UP"
            return up0_energy(E0, Ec, Δ)
        elseif uppercase(string(branch)) == "LP"
            return lp0_energy(E0, Ec, Δ)
        else
            error("Invalid branch option $branch")
        end
    end

    if guess === nothing
        if uppercase(string(branch)) == "UP"
            guess = up0_energy(E0, Ec, Δ)
        elseif uppercase(string(branch)) == "LP"
            guess = lp0_energy(E0, Ec, Δ)
        else
            error("Invalid branch option $branch")
        end
    end

    prefact = (Δ^2 * √π) / (σ)

    # Function for which the root must be found
    f(Ep) = im * erfcx(-im*((Ep - E0)/σ)) * prefact + Ep - Ec

    # This auxiliary function is necessary as a input to the `nlsolve` routine
    # which is used here to solve the non linear equation.
    function f!(F, xvec)
        Ep = complex(xvec[1], xvec[2])
        y = f(Ep)
        # println("y = $y")
        F[1] = real(y)
        F[2] = imag(y)
    end

    # Build a guess vector for the `nlsolve` routine.
    # Note that instead of using complex numbers
    # we use a two dimensional vector
    guess_vec = [real(guess), imag(guess)]

    sol = nlsolve(f!, guess_vec) 
    return complex(sol.zero[1], sol.zero[2])
end

function numerical_derivative_three_point(x, y)
    n = length(x)
    dy_dx = zeros(n)  # Initialize the derivative array

    for i in 2:n-1
        dy_dx[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
    end

    # Calculate the endpoints using a forward and a backward difference
    dy_dx[1] = (y[2] - y[1]) / (x[2] - x[1])
    dy_dx[n] = (y[n] - y[n-1]) / (x[n] - x[n-1])

    return dy_dx
end

include("Dispersion.jl")
include("GroupVelocity.jl")
include("Wavevector.jl")

end