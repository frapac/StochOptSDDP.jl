
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
