################################################################################
mutable struct NodeData
    timestep::Int
    pb::DynamicProgrammingModel
    cutstore # TODO
    cuts::Vector{AbstractCut}
    noises::SOI.probability
end
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
    DynamicProgrammingModel(umin, umax, xmin, xmax,
                            cost, dynamics, JuMP.Model())

function optimize(dp::DynamicProgrammingModel)
    JuMP.optimize(dp.pb)
end


struct Solution <: SOI.AbstractSolution
    status::Symbol
    objval::Float64
    ut::Vector{Float64} # control
    λt::Vector{Float64} # cut slope
    xf::Vector{Float64} # next optimal state
    θf::Vector{Float64} # value of θ
end

function Node(timestep::Int64, cost::Function, dynamics::Function,
              umin::Vector{Float64}, umax::Vector{Float64},
              xmin::Vector{Float64}, xmax::Vector{Float64},
              noises::IndependentProbability)
    error("type of $(noises) is not supported")
end

function Node(timestep::Int64, cost::Function, dynamics::Function,
              umin::Vector{Float64}, umax::Vector{Float64},
              xmin::Vector{Float64}, xmax::Vector{Float64},
              noises::IndependentProbability)
    thisnode = NodeData()
    thisnode = build_node!(node, timestep, cost, dynamics,umin, umax, xmin, xmax,
                        noises)
    return = thisnode
end


function build_dynamic_problem(timestep::Int64, cost::Function,
                dynamics::Function, umin::Vector{Float64},
                umax::Vector{Float64}, xmin::Vector{Float64},
                xmax::Vector{Float64}, noises::IndependentProbability)
    Dimstate = size(xmin)
    Dimcontrol = size(umin)
    Dimalea = noises.laws.ndims
    # initialize w with dummy values
    w = Vector{Float64}(Dimalea)
    Modelpb = Model()
    @variable(Modelpb, xmin .<= x[1:Dimstate] .<= xmax)
    @variable(Modelpb, alpha[1:Dimstate])
    @variable(Modelpb, xf[1:Dimstate, 1:Dimalea])
    @variable(Modelpb, u[1:Dimcontrol, 1:Dimalea])

    # @objective(Modelpb, Min, cost(timestep, x, u, w))
    @constraint(Modelpb, dynamics(timestep, x, u, w) == 0)
    for i in [1:Dimalea]
        @constraint(Modelpb, xmin .<=xf[:,i] .<= xmax)
        @constraint(Modelpb, umin .<=u[:,i] .<= umax)
    end
    ##no constraint on alpha ????
    Nodepb = DynamicProgrammingModel(umin, umax, xmin, xmax, cost,
                                    dynamics, Modelpb )
    return Nodepb
end

function build_node!(node ::NodeData, timestep::Int64, cost::Function,
                dynamics::Function, umin::Vector{Float64},
                umax::Vector{Float64}, xmin::Vector{Float64},
                xmax::Vector{Float64}, noises::IndependentProbability)
    dp = build_dynamic_problem(timestep, cost, dynamics, umin, umax, xmin, xmax, noises)
    node.timestep = timestep
    node.pb = dp
    node.noises = noises
end

# Solves the DynamicProgrammingModel of a node, provided xt and ξ.
# Retruns a Solution
function solve_one_stage_dp(dp::DynamicProgrammingModel, xt,ξ)
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
