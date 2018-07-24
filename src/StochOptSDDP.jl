################################################################################

module StochOptSDDP

import JuMP
import MathOptInterface, StochOptInterface
using CutPruners, LightGraphs
import Scenarios


const SOI = StochOptInterface
const MOI = MathOptInterface

include("sddp.jl")

end
