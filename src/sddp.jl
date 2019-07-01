# a brand new implementation of SDDP

# cuts utilities
abstract type AbstractCut end

struct Cut <: AbstractCut
    β::Float64
    λ::Vector{Float64}
end


function solve_one_stage_dpm(dp::DynamicProgrammingModel, xt,ξ)
    m = dp.pb
    dimstate = size(dp.xmax)
    dimcontrol = size(dp.umax)
    x = getindex(m, :x)
    u = getindex(m, :u)
    w = getindex(m, :w)
    alpha = getindex(m, :alpha)
    xf = getindex(m, :xf)

    # Update value of w:
    JuMP.fix.(w,ξ)

    #update objective
    if isa(dp.cost, Function)
        @objective(m, Min, dp.cost(t, x, u, ξ) + alpha)

    else
        error("`cost' is a$(typeof(dp.cost)) and not a Function")
    end

    # Update constraint x == xt
    if(size(xt)==dimstate)
        for i in 1:dimstate
            JuMP.setRHS(m.ext[:cons][i], xt[i])
        end
    else
        error("$(xt) is not of size $(dimstate)")
    end

    if false
        println("One step one alea problem at time t=",t)
        println("for x =",xt)
        println("and w=",ξ)
        print(m)
    end

    status = (verbosity>3) ? solve(m, suppress_warnings=false) : solve(m, suppress_warnings=false)
    solved = (status == :Optimal)

    # get time taken by the solver:
    solvetime = try getsolvetime(m) catch 0 end

    if solved
        # Return object storing results:
        result = Solution(
                          solved,
                          getobjectivevalue(m),
                          getvalue(u),
                          getdual(m.ext[:cons]),
                          getvalue(xf),
                          getvalue(alpha))
    else
        println(m)
        println(status)
        error("Foo")
        # If no solution is found, then return nothing
        result = Solution()
    end
    return Solution
end

struct SDDP <: SOI.AbstractAlgorithm
    solvers::MOI.AbstractOptimizer
    options::Dict
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
