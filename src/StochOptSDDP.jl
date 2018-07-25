################################################################################
module StochOptSDDP

import JuMP
import MathOptInterface, StochOptInterface
using CutPruners
import Scenarios


const SOI = StochOptInterface
const MOI = MathOptInterface

include("SDDP/sddp.jl")
include("Model/node.jl")
include("Model/proba.jl")
include("Model/problem.jl")

end
