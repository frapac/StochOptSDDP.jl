################################################################################
# Probability utilities
################################################################################
# we are able to implement different kind of transitions

struct IndependentProbability <: SOI.Probability
    laws::Scenarios.DiscreteLaw
    # we have an unique child!
    child::Int
end

# in this case, we just sample a realization of the marginal laws
sample(trans::IndependentProbability) = (trans.child, rand(trans.laws))


struct MarkovProbability <: SOI.Probability
    # transition from outgoing edges
    childproba::Scenarios.DiscreteLaw{Int}
    # independent realization
    laws::Scenarios.DiscreteLaw
end
# first attempt to build a constructor for MarkovTransition
# TODO: decide if it is convenient enough
MarkovProbability(μ::Scenarios.DiscreteLaw, childs, probaschild::Vector{Float64}) =
    MarkovTransition(Scenarios.DiscreteLaw(childs, probaschild), μ)

sample(trans::MarkovProbability) = (rand(trans.childproba), rand(trans.laws))
