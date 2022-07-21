module DMSuite

export chebdif, chebdifft, chebint, cheb2bc, cheb4c, 
       fourdif, fourdifft, fourint, sincdif, sincdifft, 
       herroots, lagroots, legroots, poldif, polint, 
       herdif, lagdif
    
    include("chebdif.jl")
    include("chebdifft.jl")
    include("chebint.jl")
    include("cheb2bc.jl")
    include("cheb4c.jl")
    
    include("fourdif.jl")
    include("fourdifft.jl")
    include("fourint.jl")
    
    include("sincdif.jl")
    include("sincdifft.jl")
    
    include("herroots.jl")
    include("lagroots.jl")
    include("legroots.jl")

    include("poldif.jl")
    include("polint.jl")

    include("herdif.jl")
    include("lagdif.jl")

end