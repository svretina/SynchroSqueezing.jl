#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "experiments", "continuouswavelet_comparison.jl"))

const DEFAULT_LENGTHS = (2_048, 8_192, 16_384, 65_536)

function main()
    lengths = isempty(ARGS) ? DEFAULT_LENGTHS : Tuple(parse.(Int, ARGS))
    for n in lengths
        println("\n=== N = $n ===")
        compare_backends(n = n)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
