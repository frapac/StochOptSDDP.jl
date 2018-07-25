# a brand new implementation of SDDP







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



# cuts utilities
abstract type AbstractCut end

struct Cut <: AbstractCut
    β::Float64
    λ::Vector{Float64}
end

struct MultistageStochasticProgram <: SOI.AbstractStochasticProgram
    data::Vector{NodeData}
    num_stages::Int
end
