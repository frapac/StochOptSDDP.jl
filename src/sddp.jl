# a brand new implementation of SDDP

# cuts utilities
abstract type AbstractCut end

struct Cut <: AbstractCut
    β::Float64
    λ::Vector{Float64}
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
