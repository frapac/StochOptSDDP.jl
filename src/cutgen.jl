################################################################################
# Cuts generator
################################################################################

abstract type AbstractCutGenerator end

####################
struct AverageCutGenerator <: AbstractCutGenerator end

# TODO: state `x` is missing

function gencuts(::AverageCutGenerator, sp, node, pos, pool)
    avgλ = zeros()
    avgβ = 0.

    ξ = probalaw(node)
    wghts = Scenarios.weights(ξ)

    # we have to consider every outgoing edges
    # and every independent noises
    for tr in SOI.get(sp, SOI.OutTransitions(), node)

        atomλ = zeros()
        atomβ = 0.

        for iξ in 1:length(ξ)

            # TODO
            sol = subsolve(node, tr, ξ.support[iξ, :])
            λ = sol.λt
            β = sol.objval - dot(λ, x)
            atomλ += wghts[iξ] * λ
            atomβ += wghts[iξ] * β
        end

        proba = SOI.get(sp, SOI.Probability(), tr)
        avgλ += proba * atomλ
        avgβ += proba * atomβ
    end

    return Cut(avgβ, avgλ)
end

####################
struct MultiCutGenerator <: AbstractCutGenerator end

function gencuts(::MultiCutGenerator, sp, node, pos, pool)

    multicuts = Cut[]

    ξ = probalaw(node)
    wghts = Scenarios.weights(ξ)

    # we have to consider every outgoing edges
    # and every independent noises
    for tr in SOI.get(sp, SOI.OutTransitions(), node)
        for iξ in 1:length(ξ)

            # TODO
            sol = subsolve(node, tr, ξ.support[iξ, :])
            λ = sol.λt
            β = sol.objval - dot(λ, x)
            push!(multicuts, Cut(β, λ))
        end
    end

    return multicuts
end
