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

mutable struct Transition
    id::Int # node id
    ξ::Vector{Float} # noise
end

mutable struct Path
    scenario::Vector{AbstractTransition}
    sol::Vector{Solution}
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
    ut::Vector{Float64} # control
    λt::Vector{Float64} # cut slope
    xf::Vector{Float64} # next optimal state
    θf::Vector{Float64} # value of θ
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

function forward_pass(sp::MultistageStochasticProgram, algo::SDDP)
    #TODO define
    scenario = sample_scenarios(sp::MultistageStochasticProgram,
    algo::SDDP,verbose::Int)
    #TODO one fwd pass only
    simulate_scenario(sp,scenario,SDDP)
end

function simulate_scenario(sp::MultistageStochasticProgram, algo::SDDP, scenario::Scenario)

    stats = SOI.SDDPStats()

    # TODO: define
    stats.solvertime += SOI.@_time mastersol = SOI.get(sp, SOI.Solution(),
    master)
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
        sol = solve(node, current_state, current_noise)
        current_state = sol.xf
        current_node = scenario[t].id #TODO
        push!(path,sol)
    end

    # update stats
    stats.upperbound = z_UB
    stats.σ_UB = σ
    stats.npaths = algo.K
    stats.lowerbound = SOI.getobjectivevalue(mastersol)

    pathsd, mastersol, stats
end

""" #TODO rewrite
Solve the Bellman equation at time t starting at state x under alea xi
with the current evaluation of Vt+1

# Description
The function solve
min_u current_cost(t,x,u,xi) + current_Bellman_Value_{t+1}(dynamic(t,x,u,xi))
and can return the optimal control and a subgradient of the value of the
problem with respect to the initial state x

# Arguments
* `model::SPmodel`:
    the stochastic problem we want to optimize
* `param::SDDPparameters`:
    the parameters of the SDDP algorithm
* `m::JuMP.Model`:
    The linear problem to solve, in order to approximate the
    current value functions
* `t::int`:
    time step at which the problem is solved
* `xt::Array{Float}`:
    current starting state
* `xi::Array{float}`:
    current noise value
* `relaxation::Bool`: default is false
    If problem is MILP, specify if it is needed to relax integer constraints.
* `init::Bool`:
    If specified, approximate future cost as 0

# Returns
* `solved::Bool`:
    True if the solution is feasible, false otherwise
* `NextStep`:
    Store solution of the problem
* `ts::Float64`:
   Solver's execution time
"""
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
