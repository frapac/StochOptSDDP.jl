################################################################################
abstract type AbstractNode end

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

optimize(dp::DynamicProgrammingModel) = JuMP.optimize(dp.pb)

mutable struct NodeData <: AbstractNode
    # time step number
    timestep::Int
    # parent nodes
    parents::Vector{Int}
    # optimization problem
    pb::DynamicProgrammingModel
    # cut pruners
    cutstore # TODO
    # cuts previously computed
    cuts::Vector{AbstractCut}
    # childs are stored inside noises
    noises::SOI.AbstractTransition
end

struct EmptyNode <: AbstractNode
    timesteps::Int
    parents::Vector{Int}
    noises::SOI.AbstractTransition
end




struct Solution <: SOI.AbstractSolution
    status::Symbol
    objval::Float64
    ut::Vector{Float64} # control
    λt::Vector{Float64} # cut slope
    xf::Vector{Float64} # next optimal state
    θf::Vector{Float64} # value of θ
end
