
push!(LOAD_PATH, "..")

using StochOptSDDP, Scenarios
const DP = StochOptSDDP


p1 = DP.IndependentProbability(Scenarios.DiscreteLaw([1.]), 2)
n1 = DP.EmptyNode(1, [-1], p1)

p2 = DP.IndependentProbability(Scenarios.DiscreteLaw([1.]), 3)
n2 = DP.EmptyNode(2, [1], p2)


sp = DP.MultistageStochasticProgram([n1, n2], 2)

println(DP.sample(n1.noises))
