################################################################################
module StochOptSDDP

import JuMP
import MathOptInterface, StochOptInterface
using CutPruners
import Scenarios


const SOI = StochOptInterface
const MOI = MathOptInterface

include("sddp.jl")
include("node.jl")
include("proba.jl")
include("problem.jl")

end
