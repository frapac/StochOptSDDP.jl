# a brand new implementation of SDDP

# TODO
struct DynamicProgrammingModel
    umin::Vector{Float64}
    umax::Vector{Float64}
    xmin::Vector{Float64}
    xmax::Vector{Float64}
    cost::Function
    dynamics::Function

    pb::JuMP.Model
end

DynamicProgrammingModel(umin, umax, xmin, xmax, cost, dynamics) =
    DynamicProgrammingModel(umin, umax, xmin, xmax, cost, dynamics, JuMP.Model())

function optimize(dp::DynamicProgrammingModel)
    JuMP.optimize(dp.pb)
end


################################################################################
# Probability utilities
################################################################################
# we are able to implement different kind of transitions

struct IndependentProbability <: SOI.Probability
    laws::Scenarios.DiscreteLaw
    # we have an unique child!
    child::Int
end

# in this case, we just sample a realization of the marginal laws
sample(trans::IndependentProbability) = (trans.child, rand(trans.laws))


struct MarkovProbability <: SOI.Probability
    # transition from outgoing edges
    childproba::Scenarios.DiscreteLaw{Int}
    # independent realization
    laws::Scenarios.DiscreteLaw
end
# first attempt to build a constructor for MarkovTransition
# TODO: decide if it is convenient enough
MarkovProbability(μ::Scenarios.DiscreteLaw, childs, probaschild::Vector{Float64}) =
    MarkovTransition(Scenarios.DiscreteLaw(childs, probaschild), μ)

sample(trans::MarkovProbability) = (rand(trans.childproba), rand(trans.laws))


################################################################################
mutable struct NodeData
    timestep::Int
    pb::DynamicProgrammingModel
    cutstore # TODO
    cuts::Vector{AbstractCut}
    noises::NodeTransition
end


abstract type AbstractTransition end

struct Transition{T}
    id::Int # node id
    ξ::Vector{T} # noise
end
Scenario = Vector{<:AbstractTransition}

struct Path
    scenario::Scenario
    sol::Vector{Solution}
end



# cuts utilities
abstract type AbstractCut end

struct Cut <: AbstractCut
    β::Float64
    λ::Vector{Float64}
end

struct MultistageStochasticProgram <: SOI.AbstractStochasticProgram
    data::Vector{NodeData}
    num_stages::Int
end


struct Solution <: SOI.AbstractSolution
    status::Symbol
    objval::Float64
    ut::Vector{Float64} # control
    λt::Vector{Float64} # cut slope
    xf::Vector{Float64} # next optimal state
    θf::Vector{Float64} # value of θ
end

struct SDDP <: SOI.AbstractAlgorithm
    solvers::MOI.AbstractOptimizer
    options::Dict
end

################################################################################
# Cuts generator
################################################################################

abstract type AbstractCutGenerator end

####################
struct AverageCutGenerator <: AbstractCutGenerator end

# TODO: state `x` is missing

function gencuts(::AverageCutGenerator, sp, node, pos, pool)
    avgλ = zeros()
    avgβ = 0.

    ξ = probalaw(node)
    wghts = Scenarios.weights(ξ)

    # we have to consider every outgoing edges
    # and every independent noises
    for tr in SOI.get(sp, SOI.OutTransitions(), node)

        atomλ = zeros()
        atomβ = 0.

        for iξ in 1:length(ξ)

            # TODO
            sol = subsolve(node, tr, ξ.support[iξ, :])
            λ = sol.λt
            β = sol.objval - dot(λ, x)
            atomλ += wghts[iξ] * λ
            atomβ += wghts[iξ] * β
        end

        proba = SOI.get(sp, SOI.Probability(), tr)
        avgλ += proba * atomλ
        avgβ += proba * atomβ
    end

    return Cut(avgβ, avgλ)
end

####################
struct MultiCutGenerator <: AbstractCutGenerator end

function gencuts(::MultiCutGenerator, sp, node, pos, pool)

    multicuts = Cut[]

    ξ = probalaw(node)
    wghts = Scenarios.weights(ξ)

    # we have to consider every outgoing edges
    # and every independent noises
    for tr in SOI.get(sp, SOI.OutTransitions(), node)
        for iξ in 1:length(ξ)

            # TODO
            sol = subsolve(node, tr, ξ.support[iξ, :])
            λ = sol.λt
            β = sol.objval - dot(λ, x)
            push!(multicuts, Cut(β, λ))
        end
    end

    return multicuts
end


################################################################################

# for simulation
struct SolutionStore end

function forward_pass(sp::MultistageStochasticProgram, algo::SDDP)
    #TODO define
    # sample one noise scenarios
    scenario = sample(sp)
    #TODO one fwd pass only
    simulate(sp, scenario, algo)
end

function sample(sp::MultistageStochasticProgram)
    # TODO:
    nothing
end

function simulate(sp::MultistageStochasticProgram, scenario::Scenario, algo::SDDP)

    stats = SOI.SDDPStats()

    # TODO: define
    stats.solvertime += SOI.@_time mastersol = SOI.get(sp, SOI.Solution(), master)
    stats.nsolved += 1
    stats.niterations += 1
    infeasibility_detected = SOI.getstatus(mastersol) == :Infeasible

    #TODO
    current_node =
    current_state =

    num_stages = length(scenario)

    for t in 1:num_stages
        verbose >= 3 && println("Stage $t/$num_stages")

        current_noise = scenario[t].ξ
        sol = subsolve(node, current_state, current_noise)
        current_state = sol.xf
        current_node = scenario[t].id #TODO
        push!(path, sol)
    end

    # update stats
    stats.upperbound = z_UB
    stats.σ_UB = σ
    stats.npaths = algo.K
    stats.lowerbound = SOI.getobjectivevalue(mastersol)

    pathsd, mastersol, stats
end

function solve_one_stage_one_alea(model,
                                 param,
                                 m::JuMP.Model,
                                 t::Int64,
                                 xt::Vector{Float64},
                                 xi::Vector{Float64};
                                 relaxation=false::Bool,
                                 init=false::Bool,
                                 verbosity::Int64=0)

    # Get var defined in JuMP.model:
    x = getindex(m, :x)
    u = getindex(m, :u)
    w = getindex(m, :w)
    alpha = getindex(m, :alpha)

    # Update value of w:
    JuMP.fix.(w,xi)

    #update objective
    if isa(model.costFunctions, Function)
        @objective(m, Min, model.costFunctions(t, x, u, xi) + alpha)

    elseif isa(model.costFunctions, Vector{Function})
        cost = getindex(m, :cost)
        for i in 1:length(model.costFunctions)
            @constraint(m, cost >= model.costFunctions[i](t, x, u, xi))
        end
        @objective(m, Min, cost + alpha)
     end

    # Update constraint x == xt
    for i in 1:model.dimStates
        JuMP.setRHS(m.ext[:cons][i], xt[i])
    end

    if false
        println("One step one alea problem at time t=",t)
        println("for x =",xt)
        println("and w=",xi)
        print(m)
    end

    if model.IS_SMIP
        solved = relaxation ? solve_relaxed!(m, param,verbosity): solve_mip!(m, param,verbosity)
    else
        status = (verbosity>3) ? solve(m, suppress_warnings=false) : solve(m, suppress_warnings=false)
        solved = (status == :Optimal)
    end

    # get time taken by the solver:
    solvetime = try getsolvetime(m) catch 0 end

    if solved
        optimalControl = getvalue(u)
        # Return object storing results:
        result = NLDSSolution(
                          solved,
                          getobjectivevalue(m),
                          model.dynamics(t, xt, optimalControl, xi),
                          optimalControl,
                          getdual(m.ext[:cons]),
                          getvalue(alpha),
                          getcutsmultipliers(m))
    else
        println(m)
        println(status)
        error("Foo")
        # If no solution is found, then return nothing
        result = NLDSSolution()
    end

    return result, solvetime
end


function backward_pass!(sp::MultistageStochasticProgram, algo, trajectories)

    stats = SOI.SDDPStats()
    num_stages = SOI.get(sp, SOI.NumberOfStages())

    for t in (num_stages-1):-1:1
        # Make jobs
        endedpaths = SDDPPath{SOI.get(sp, SOI.TransitionType())}[]

        gencuts(pathsd[t], sp, stats)

        # OLD CODE
        #= jobsd = childjobs(sp, pathsd[t-1], algo.pathsampler, t, num_stages, endedpaths) # FIXME shouldn't need pathsampler here =#
        #= # Solve Jobs (parallelism possible here) =#
        #= solvejobs!(sp, jobsd, stats, algo.stopatinf) =#

        #= gencuts(pathsd[t-1], sp, stats, algo.ztol) =#
        #= # The cut is added after so that they are added by group and =#
        #= # duplicate detection in CutPruners works better =#
        #= applycuts(pathsd[t-1], sp) =#
    end
    stats
end
