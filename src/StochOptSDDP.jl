################################################################################

module StochOptSDDP

import JuMP
import MathOptInterface, StochOptInterface
using CutPruners, LightGraphs


const SOI = StochOptInterface
const MOI = MathOptInterface

include("sddp.jl")

end
