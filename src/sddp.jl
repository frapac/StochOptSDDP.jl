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

mutable struct NodeData
    timestep::Int
    pb::DynamicProgrammingModel
    cutstore # TODO
    cuts::Vector{AbstractCut}
    # TODO: where should we put proba vector?
    noises
end

mutable struct PathElement
    id::Int
    ξ::Vector{Float}
    sol::Solution
end

mutable struct Path
    path::Vector{PathElement}
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
# get number of stages inside a MultistageStochasticProgram
nstages(sp::MultistageStochasticProgram) = length(sp.data)

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

# TODO!
struct CutGenerator end


# TODO
function gencuts end


struct SolutionStore end

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
