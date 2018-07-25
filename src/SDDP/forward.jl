

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
