


struct SDDP <: SOI.AbstractAlgorithm
    solvers::MOI.AbstractOptimizer
    options::Dict
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
