module betnet2

    using Distributions
    using BenchmarkTools
    using Base.Threads: @spawn, @threads
    using CSV
    using DataFrames
    using Random
    using StatsBase
    using FreqTables  
      
    const ERROR = 1e-5

    export
        #data and tree definition    
        TransNet,
        Tborg
    include("betnet2_code.jl")

end # module betnet2
