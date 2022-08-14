using ModelingToolkit
using DifferentialEquations
using Plots

# define variables and parameters of system
@variables t ρ[1:2, 1:2](t) δ(t) Ω(t)
@parameters Δ Ω0 ω ϕ δ0 icomplex
D = Differential(t)

# Hamiltonian
H = [
    0 Ω
    Ω Δ+δ
]

# bloch equations
eq = -icomplex * Symbolics.scalarize(H*ρ - ρ*H)
eqns = [D(ρ[idx, idy]) ~ eq[idx,idy] for idx in 1:2 for idy in 1:2]
# time varying detuning
append!(eqns, [δ ~ δ0*cos(ω*t+ϕ)])
append!(eqns, [Ω ~ Ω0*sin(ω*t+ϕ)])

# create symbolic system
@named bloch = ODESystem(eqns)
# structural_simplify simplifies a system of equations if possible by removing
# reduntancies etc.
bloch_simplified = structural_simplify(bloch)

# initial conditions
ρᵢ = zeros(ComplexF64,2,2)
ρᵢ[1,1] = 1
u0 = [ρ[idx,idy] => ρᵢ[idx,idy] for idx in 1:2 for idy in 1:2]

# initial parameters
p = [Ω0 => 1., Δ => 20., δ0 => 1., ω => 20., ϕ =>0, icomplex => 1im]

# setup and solve of ODE, t from 0 => π
prob = ODEProblem(bloch_simplified, u0, (0., π), p, jac = true)
sol = solve(prob, Tsit5(), dtmax = 1e-2)


plot(sol, vars = [abs(ρ[1,1]), abs(ρ[2,2])])

# scan ω
ωi = 0.0:0.01:30
function prob_func(prob, i, repeat)
    # change ω for each new iteration
    pnew = ModelingToolkit.varmap_to_vars([ω => ωi[i], Ω0 => 1., Δ=>20., δ0=>1., ϕ=>0., icomplex =>1im], parameters(bloch_simplified))
    remake(prob, p=pnew)
end

function output_func(sol, i)
    # only store last point for each solution in the ensemble
    (sol[end], false)
end

ens_prob = EnsembleProblem(prob, output_func = output_func, prob_func = prob_func)
ens_sol = solve(ens_prob, Tsit5(), EnsembleThreads(), trajectories = size(ωi)[1])
# extract gnd and exc
gnd = [real(ens_sol.u[i][1]) for i in 1:size(ωi)[1]]
exc = [real(ens_sol.u[i][4]) for i in 1:size(ωi)[1]]

pl = plot(ωi, gnd, label = "gnd")
plot!(pl, ωi, exc, label = "exc")

# scan ϕ
ϕi = 0.0:0.01:2π
function prob_func(prob, i, repeat)
    # change ω for each new iteration
    pnew = ModelingToolkit.varmap_to_vars([ω => 20, Ω0 => 1., Δ=>20., δ0=>1., ϕ=>ϕi[i], icomplex =>1im], parameters(bloch_simplified))
    remake(prob, p=pnew)
end

function output_func(sol, i)
    # only store last point for each solution in the ensemble
    (sol[end], false)
end

ens_prob = EnsembleProblem(prob, output_func = output_func, prob_func = prob_func)
ens_sol = solve(ens_prob, Tsit5(), EnsembleThreads(), trajectories = size(ϕi)[1])
# extract gnd and exc
gnd = [real(ens_sol.u[i][1]) for i in 1:size(ϕi)[1]]
exc = [real(ens_sol.u[i][4]) for i in 1:size(ϕi)[1]]

pl = plot(ϕi, gnd, label = "gnd")
plot!(pl, ϕi, exc, label = "exc")

# realistic parameters
p = [Ω0 => 2π*25, Δ => 2π*11e3, δ0 => 2π*100, ω => 2π*11e3, ϕ =>0, icomplex => 1im]

# setup and solve of ODE, t from 0 => π
prob = ODEProblem(bloch_simplified, u0, (0., 10e-3), p, jac = true)
sol = solve(prob, Tsit5(), dtmax = 1e-2, reltol=1e-8, abstol=1e-8)

plot(sol, vars = [abs(ρ[1,1]), abs(ρ[2,2])])

# scan ω
ωi = 10e3:10:12e3
function prob_func(prob, i, repeat)
    # change ω for each new iteration
    pnew = ModelingToolkit.varmap_to_vars([Ω0 => 2π*25, Δ => 2π*11e3, δ0 => 0*2π*100, ω => 2π*ωi[i], ϕ =>0, icomplex => 1im], parameters(bloch_simplified))
    remake(prob, p=pnew)
end

function output_func(sol, i)
    # only store last point for each solution in the ensemble
    (sol[end], false)
end

ens_prob = EnsembleProblem(prob, output_func = output_func, prob_func = prob_func)
ens_sol = solve(ens_prob, Tsit5(), EnsembleThreads(), trajectories = size(ωi)[1])
# extract gnd and exc
gnd = [real(ens_sol.u[i][1]) for i in 1:size(ωi)[1]]
exc = [real(ens_sol.u[i][4]) for i in 1:size(ωi)[1]]

pl = plot(ωi, gnd, label = "gnd")
plot!(pl, ωi, exc, label = "exc")

# scan ϕ
ϕi = 0.0:0.01:2π
function prob_func(prob, i, repeat)
    # change ω for each new iteration
    pnew = ModelingToolkit.varmap_to_vars([Ω0 => 2π*25, Δ => 2π*11e3, δ0 => 0*2π*100, ω => 2π*11e3, ϕ =>ϕi[i], icomplex => 1im], parameters(bloch_simplified))
    remake(prob, p=pnew)
end

function output_func(sol, i)
    # only store last point for each solution in the ensemble
    (sol[end], false)
end

ens_prob = EnsembleProblem(prob, output_func = output_func, prob_func = prob_func)
ens_sol = solve(ens_prob, Tsit5(), EnsembleThreads(), trajectories = size(ϕi)[1])
# extract gnd and exc
gnd = [real(ens_sol.u[i][1]) for i in 1:size(ϕi)[1]]
exc = [real(ens_sol.u[i][4]) for i in 1:size(ϕi)[1]]

pl = plot(ϕi, gnd, label = "gnd")
plot!(pl, ϕi, exc, label = "exc")