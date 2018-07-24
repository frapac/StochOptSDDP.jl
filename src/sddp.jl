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

# cuts utilities
abstract type AbstractCut end

mutable struct Cut <: AbstractCut
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
    ut::Vector{Float64}
    λt::Vector{Float64}
    xf::Vector{Float64}
    θf::Vector{Float64}
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

# for forward pass
struct Path
    store::Vector{Solution}
end

function forwardpass!(sp::MultistageStochasticProgram, algo::SDDP)

    stats = SOI.SDDPStats()

    # TODO: define
    master = SOI.get(sp, SOI.MasterState())

    # TODO: define
    stats.solvertime += SOI.@_time mastersol = SOI.get(sp, SOI.Solution(), master)
    stats.nsolved += 1
    stats.niterations += 1
    infeasibility_detected = SOI.getstatus(mastersol) == :Infeasible

    # TODO: define
    num_stages = SOI.get(sp, SOI.NumberOfStages())

    for t in 1:num_stages
        verbose >= 3 && println("Stage $t/$num_stages")

        # Solve Jobs (parallelism possible here)
        infeasibility_detected |= solvejobs!(sp, jobsd, stats, algo.stopatinf)

        if algo[:forwardcuts]
            gencuts(pathsd[t-1], sp, stats, algo.ztol)
            # The cut is added after so that they are added by group and duplicate detection in CutPruners works better
            applycuts(pathsd[t-1], sp)
        end

        # Jobs -> Paths
        pathsd[t] = ???
    end

    if infeasibility_detected
        z_UB = Inf # FIXME assumes minimization
        σ = 0
    else
        for (_, paths) in pathsd[end]
            append!(endedpaths, paths)
        end
        z_UB, σ = meanstdpaths(endedpaths, algo.K)
    end

    # update stats
    stats.upperbound = z_UB
    stats.σ_UB = σ
    stats.npaths = algo.K
    stats.lowerbound = SOI.getobjectivevalue(mastersol)

    pathsd, mastersol, stats
end

function backwardpass!(sp::MultistageStochasticProgram, algo, trajectories)

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
