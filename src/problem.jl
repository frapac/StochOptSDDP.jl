

struct MultistageStochasticProgram <: SOI.AbstractStochasticProgram
    data::Vector{NodeData}
    num_stages::Int
end
