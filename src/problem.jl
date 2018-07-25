

struct MultistageStochasticProgram <: SOI.AbstractStochasticProgram
    initial_state::Vector{Float}
    data::Vector{NodeData}
    num_stages::Int
end
