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


function build_dynamic_problem(timestep::Int64, cost::Function,
                dynamics::Function, umin::Vector{Float64},
                umax::Vector{Float64}, xmin::Vector{Float64},
                xmax::Vector{Float64}, noises::SOI.Probability,
                cuts::Vector{AbstractCut})
    Ncuts = length(cuts)
    Dimstate = size(xmin)
    Dimcontrol = size(umin)
    Dimalea = noises.laws.ndims
    # initialize w with dummy values
    w = Vector{Float64}(Dimalea)
    Modelpb = Model()
    @variable(Modelpb, xmin .<= x[1:Dimstate] .<= xmax)
    @variable(Modelpb, xmin .<= xf[1:Dimstate] .<= xmax)
    @variable(Modelpb, alpha)
    @variable(Modelpb, umin .<= u[1:Dimcontrol] .<= umax)

    @objective(Modelpb, Min, cost(timestep, x, u, w) + alpha)

    @constraint(Modelpb, dynamics(timestep, x, u, w) .== xf)
    for k in [1:Ncuts]
        @constraint(Modelpb, alpha >= dot(cuts[k].λ,xf) + cuts[k].β)
    end
    ##no constraint on alpha ????
    Nodepb = DynamicProgrammingModel(umin, umax, xmin, xmax, cost,
                                    dynamics, Modelpb )
    return Nodepb
end


function create_cost(c::Array{Float64,1}, d::Array{Float64,1},
                     f::Array{Float64,1})
    function cost(t, x, u, w)
        return dot(c, x) + dot(d, u) + dot(f, w)
    end
    return cost
end


function create_dynamics(A::Array{Float64,2}, B::Array{Float64,2},
                        C::Array{Float64,2})
    function dynamics(t, x, u, w)
        return A*x + B*u + C*w
    end
    return dynamics
end


# cost : min{ cx+du+fw + alpha}
# dynamics: Ax_{t} + Du_{t} + Ew_{t} = x_{t+1}
#cuts:      alpha >= <lambda, x_{t+1}> + beta
# bounds:   Fx_{t} <h
#           Gu_{t} <j
function build_dynamic_problem(timestep::Int64,
                c::Array{Float64,1}, d::Array{Float64,1}, f::Array{Float64,1},
                A::Array{Float64,2}, B::Array{Float64,2}, C::Array{Float64,2},
                F::Array{Float64,1}, G::Array{Float64,1}, h::Array{Float64,1},
                j::Array{Float64,1}, cuts::Vector{AbstractCut})
    Ncuts = length(cuts)
    Dimstate = length(c)
    Dimcontrol = length(d)
    Dimalea = length(f)

    cost = create_cost(c, d, f)
    dynamics = create_dynamics(A, B, C)
    w=Vector{Float64}(Dimalea)

    Modelpb = Model()
    @variable(Modelpb, xmin .<= x[1:Dimstate] .<= xmax)
    @variable(Modelpb, xmin .<= xf[1:Dimstate] .<= xmax)
    @variable(Modelpb, alpha)
    @variable(Modelpb, umin .<= u[1:Dimcontrol] .<= umax)
    @objective(Modelpb, Min, cost(t, x, u, w) + alpha)
    @constraint(Modelpb, dynamics(timestep, x, u, w) .== xf)
    for k in [1:Ncuts]
        @constraint(Modelpb, alpha >= dot(cuts[k].λ,xf) + cuts[k].β)
    end
    @contraint(Modelpb, F*x <= h)
    @contraint(Modelpb, G*x <= j)
    #extract bounds from h and j
    Nodepb = DynamicProgrammingModel(-1*j[Dimcontrol+1:2*Dimcontrol],
                                j[1:Dimcontrol], -1*h[Dimstate+1:2*Dimstate],
                                h[1:Dimstate], cost, dynamics, Modelpb )
    return Nodepb
end

# only fill in the pb field, the rest of the information is in the model
function build_dynamic_problem(timestep, model)
    Nodepb = DynamicProgrammingModel()
    Nodepb.pb = model
    return Nodepb
end


function build_node!(node ::NodeData, timestep::Int64, cost::Function,
                dynamics::Function, umin::Vector{Float64},
                umax::Vector{Float64}, xmin::Vector{Float64},
                xmax::Vector{Float64}, noises::SOI.Probability,
                cuts::Vector{AbstractCut})
    dp = build_dynamic_problem(timestep, cost, dynamics, umin, umax,
                                xmin, xmax, noises, cuts)
    node.cuts = cuts
    node.timestep = timestep
    node.pb = dp
    node.noises = noises
end


#build node according a set of matrixes such that:
# cost : min{ cx+du+fw + alpha}
# dynamics: Ax_{t} + Du_{t} + Ew_{t} = x_{t+1}
#cuts:      alpha >= <lambda, x_{t+1}> + beta
# bounds:   Fx_{t} <h
#           Gu_{t} <j
function build_node!(node ::NodeData, timestep::Int64,c::Array{Float64,1},
                d::Array{Float64,1}, f::Array{Float64,1}, A::Array{Float64,2},
                B::Array{Float64,2}, C::Array{Float64,2}, F::Array{Float64,1},
                G::Array{Float64,1}, h::Array{Float64,1}, j::Array{Float64,1},
                noises::SOI.Probability, cuts::Vector{AbstractCut})

    dp = build_dynamic_problem(timestep, c, d, f, A, B, C, F, G, h, j, cuts)
    node.cuts = cuts
    node.timestep = timestep
    node.pb = dp
    node.noises = noises
end

function build_node!(node ::NodeData, timestep::Int64, model::JuMP.Model,
                    noises::SOI.Probability, cuts::Vector{AbstractCut})
    dp = build_dynamic_problem(timestep, model)
    node.cuts = cuts
    node.timestep = timestep
    node.pb = dp
    node.noises = noises
end

function Node(timestep::Int64, cost::Function, dynamics::Function,
              umin::Vector{Float64}, umax::Vector{Float64},
              xmin::Vector{Float64}, xmax::Vector{Float64},
              noises::SOI.Probability, cuts::Vector{AbstractCut})
    if (isa(noises, IndependentProbability) ||isa(noises, MarkovProbability))
        thisnode = NodeData()
        build_node!(thisnode, timestep, cost, dynamics,umin, umax,
                            xmin, xmax, noises, cuts)
    else
        error("noise type is invalid")
    end
    return = thisnode
end


function Node(timestep::Int64,c::Array{Float64,1},
                d::Array{Float64,1}, f::Array{Float64,1}, A::Array{Float64,2},
                B::Array{Float64,2}, C::Array{Float64,2}, F::Array{Float64,1},
                G::Array{Float64,1}, h::Array{Float64,1}, j::Array{Float64,1},
                noises::SOI.Probability, cuts::Vector{AbstractCut})
    if (isa(noises, IndependentProbability) ||isa(noises, MarkovProbability))
        thisnode = NodeData()
        build_node!(thisnode, timestep, c, d, f, A, B, C, F, G, h, j,
                                    noises, cuts)
    else
        error("noise type is invalid")
    end
    return = thisnode
end

function Node(timestep::Int64, model::JuMP.Model, noises::SOI.Probability,
            cuts::Vector{AbstractCut})
    if (isa(noises, IndependentProbability) ||isa(noises, MarkovProbability))
        thisnode = NodeData()
        build_node!(thisnode, timestep, model, noises, cuts)
    else
        error("noise type is invalid")
    end
    return = thisnode
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

function solve_one_stage_node(node::NodeData, xt,ξ)
    solution = solve_one_stage_dp(node.pb, xt, ξ)
    return solution
end
