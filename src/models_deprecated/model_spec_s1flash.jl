
# Code from @carolsatye 20211112
# Minor modifications @PeterMaxwell 20211114 (compiles, but otherwise not tested correctness)


function mm_spec_s1flash_dimspec()::MPCCDimSpec
    dimspec = MPCCDimSpec(21, 2, 2, 18, 20, 1, 0)   
    return dimspec
end
# Dimension of spatial variables
# Columns in F (number of complementarities?)
# Number of F functions per complementarity variable (usually 2)
# Number of equality constraints
# Number of inequality constraints
# Number of continuous parameters for PMPCC
# Number of discrete parameters for PMPCC


function mm_spec_s1flash_defn(x, pr, ps)

    # PM20211114: some modifications are a work-around for a bug in Symbolics.jl - https://github.com/JuliaSymbolics/Symbolics.jl/issues/319

    # PM20211114: Constants. PM, modified to include prefix 'c' because F was conflicting with the F constraints
    cF = 1
    cA = [3.97786; 4.00139; 3.93002]
    cB = [1064.84; 1170.874; 1182.774]
    cC = [-41.136; -48.833; -52.532]
    cp = 5
    cz = [0.5; 0.3; 0.2]

    w = Symbolics.scalarize( x[1:3] )
    y = Symbolics.scalarize( x[4:6] )
    L = x[7]
    V = x[8]
    K = Symbolics.scalarize( x[9:11] )
    psat = Symbolics.scalarize( x[12:14] )          #pure component vapor pressure
    log_psat = Symbolics.scalarize( x[15:17] )
    a = x[18]                                       #corrected vapor fraction
    at = x[19]                                      #vapor fration calculated with the Rachford-rice equation
    sv = x[20]                                      #slack for vapor fraction
    sl = x[21]                                      #slack for liquid fraction

    f = 0.5*(a*cF - V)^2

    ce = Vector{Num}(undef, 18)
    ce[1:3] = Symbolics.scalarize( log_psat .- cA .+ cB./(pr[1] .+ cC) )        #Antoine's equation
    ce[4:6] = Symbolics.scalarize( psat .- 10 .^ log_psat )
    ce[7:9] = Symbolics.scalarize( K*cp .- psat )                               #Raoult's law
    ce[10:12] = Symbolics.scalarize( K.*w .- y )
    ce[13:15] = Symbolics.scalarize( L*w .+ V*y .- cF*cz )                      #Balance per component
    ce[16] = Symbolics.scalarize( L + V - cF )                                  #Total balance
    ce[17] = Symbolics.scalarize( sum((cz .* (K .- 1)) ./ (1 .+ at*(K .- 1))) ) #Rachford-Rice equation
    ce[18] = Symbolics.scalarize( a - sv + sl - at )

    ci = Vector{Num}(undef, 20)
    ci[1:3] = Symbolics.scalarize( w )
    ci[4:6] = Symbolics.scalarize( 1 .- w )
    ci[7:9] = Symbolics.scalarize( y )
    ci[10:12] = Symbolics.scalarize( 1 .- y )
    ci[13:15] = Symbolics.scalarize( K )
    ci[16:18] = Symbolics.scalarize( psat )
    ci[19] = Symbolics.scalarize( a )
    ci[20] = Symbolics.scalarize( 1 - a )

    F = Matrix{Num}(undef, 2, 2)
    F[1,1] = Symbolics.scalarize( L )
    F[2,1] = Symbolics.scalarize( sl )
    F[1,2] = Symbolics.scalarize( V )
    F[2,2] = Symbolics.scalarize( sv )

    return MPCCDefinition(f, ce, ci, F)
end



function mm_spec_s1flash_nzmask()::MPCCModelNZMask
    nczmask = MPCCModelNZMask(
        # ce
        [   Set{Int64}([15]),
            Set{Int64}([16]),
            Set{Int64}([17]),
            Set{Int64}([12, 15]),
            Set{Int64}([13, 16]),
            Set{Int64}([14, 17]),
            Set{Int64}([9, 12]),
            Set{Int64}([10, 13]),
            Set{Int64}([11, 14]),
            Set{Int64}([1, 4, 9]),
            Set{Int64}([2, 5, 10]),
            Set{Int64}([3, 6, 11]),
            Set{Int64}([1, 4, 7, 8, 9]),
            Set{Int64}([2, 5, 7, 8, 10]),
            Set{Int64}([3, 6, 7, 8, 111]),
            Set{Int64}([7, 8]),
            Set{Int64}([9, 10, 11, 19]),
            Set{Int64}([18, 19, 20, 21])
        ],
        # ci
        [   Set{Int64}([1]),
            Set{Int64}([2]),
            Set{Int64}([3]),
            Set{Int64}([1]),
            Set{Int64}([2]),
            Set{Int64}([3]),
            Set{Int64}([4]),
            Set{Int64}([5]),
            Set{Int64}([6]),
            Set{Int64}([4]),
            Set{Int64}([5]),
            Set{Int64}([6]),
            Set{Int64}([9]),
            Set{Int64}([10]),
            Set{Int64}([11]),
            Set{Int64}([12]),
            Set{Int64}([13]),
            Set{Int64}([14]),
            Set{Int64}([18]),
            Set{Int64}([18]), 
        ],
        # F
        [ 
            # F1
            [   Set{Int64}([7]), Set{Int64}([8]) ],
            # F2
            [   Set{Int64}([21]), Set{Int64}([20]) ] 
        ] )    
    return nczmask
end


function mm_spec_s1flash_testvectors()
    return Vector{MPCCModelTestVector}(undef, 0)
end


function mm_spec_s1flash_knownsols()

    return []
end



function mm_spec_s1flash_parameterisations(t)

    defn1 = MPCCParameterisationDefn(
            Vector{Num}([ t ]),
            (380.0, 400.0),               # tspan
            "Standard"    
        )

    return  Vector{MPCCParameterisationDefn}( [ defn1 ] )
end
