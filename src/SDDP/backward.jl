
function SOI.backward_pass!(sp::MultistageStochasticProgram, algo, paths, to)

    # TODO
    for node in reverse(paths)
        cut = gencuts(algo.options[:cutgenerator], sp, node, sol)
        addcut!(node, cut)
    end
    stats
end
