################################################################################

struct MultistageStochasticProgram <: SOI.AbstractStochasticProgram
    data::Vector{AbstractNode}
    num_stages::Int
end

# TODO: populate
