# for a given TripleAxis object (with its angles defining a point in Q-E space)
# calculate a number of possible spurious modes which could produce scattering
# at the same position

const _spurion_modes = ( [0 1;0 1], [1 0; 1 0], [2 0; 0 1], [0 2; 0 1], [2 0; 0 3], [3 0; 0 4], [0 1/2; 0 1/2], [3 0; 0 2], [1 0; 0 2], [1/2 0; 0 1] )
const _spurion_descr = ( "incoherent monochromator", "incoherent analyzer", 
                         "higher order 2ki", "incoherent higher order ki=2kf", 
                         "higher order 2ki,3kf", "higher order 3ki,4kf",
                         "incoherent monochromator, doubled analyzer d-spacing", 
                         "higher order 3ki,2kf", "higher order 2kf", "double analyzer d-spacing" )
function spurions(a::TripleAxis) # opitonal Symbol selection of which mode to calculate?
    for (n,descr) in enumerate(_spurion_descr)
        println(spurions(a,n),"\t",descr)
    end
end


function spurions(a::TripleAxis,n::Integer) 
    @assert 0<n<=length(_spurion_modes)
    skikf = _spurion_modes[n]*[getki(a),getkf(a)]
    Qspur=getQ_from_kikf(a,skikf[1],a.a[3],a.a[4],skikf[2])
    Espur=ħ²2mₙ*(skikf[1]^2-skikf[2]^2)
    return Lattice4Vector(Qspur,Espur)
end
