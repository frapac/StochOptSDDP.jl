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

function simulate_scenario(sp::MultistageStochasticProgram,
     scenario::Vector{<:AbstractTransition}, to::TimerOutput, verbose)

    current_node = get(sp, MasterNode())#FIXME called MasterState in SOI
    current_state = sp.initial_state

    num_stages = length(scenario)

    for t in 1:num_stages
        verbose >= 3 && println("Stage $t/$num_stages")

        current_noise = scenario[t].ξ
        sol = solve_one_stage_dp(node, current_state, current_noise)
        current_state = sol.xf
        current_node = scenario[t].child #TODO
        push!(path, sol)
    end

    # TODO update info

end
