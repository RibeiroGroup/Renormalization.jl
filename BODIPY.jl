using Plots
using Renormalization

function load_data(filename)

    lines = readlines(filename)
    x = zeros(length(lines))
    y = zeros(length(lines))

    for i in eachindex(x)
        a, b = split(lines[i])
        x[i] = parse(Float64, a)
        y[i] = parse(Float64, b)
    end

    return x, y
end

function dispersion0(E0, Δ)

    x, y = load_data("../BSW.txt")
    μm2cm = 1e4
    P = Polariton(E0, Δ, x .* μm2cm, y)

    LP = [lp0_energy(P, q) for q in P.qvals]
    UP = [up0_energy(P, q) for q in P.qvals]

    plot(P.qvals, LP, label="LP")
    plot!(P.qvals, UP, label="UP")
    hline!([P.E0], linestyle=:dash, color=:black, label="Bare Exciton")
    plot!(P.qvals, P.Ecvals, linestyle=:dash, color=:black, label="Bare cavity")

end

function dispersion(E0, Δ, σ, q)

    x, y = load_data("../BSW.txt")
    μm2cm = 1e4
    P = Polariton(E0, Δ, x .* μm2cm, y)

    m = P.qvals .> 10.0e4
    qvals = P.qvals[m]

    println(qvals[1])
    UP, LP = up0_lp0_energy(P, qvals[1])

    println(LP)

    return find_renormalized_energy(P, qvals[1], σ, complex(1933.80214654, 7.02069795e-055))

end