"""
The `TripleAxis` object contains all elements and angles for a single setting of
a triple-axis spectrometer instrument.

| Fieldname | Description |
|:---------:|:------------|
| a         | `Vector{Real}` of the instrument angles. |
|           | at least [θm,2θm,θs,2θs,θa,2θa] following extra angles are allowed |
| source    | `Source` object describing the physical extent of the incident beam |
| mono      | `Bragg` object describing the monochromator |
| sample    | `Sample` object describing the single-crystal sample |
| ana       | `Bragg` object describing the analyzer |
| detector  | `Detector` object describing the physical extend of the detector |
| arms      | `Vector{Real}` of distances between elements |
|           | [source-mono, mono-sample, sample-ana, ana-detector, mono-monitor]|
| horcol    | `Vector{Real}` of horizontal collimators FWHM in arc minutes |
| vercol    | `Vector{Real}` of vertical collimators FWHM in arc minutes |
| energy    | `Real` value of the fixed incident or final energy in meV |
| kforki    | `Bool` to indicate if kf (`false`) or ki (`true`) is fixed |

*Note:* `horcol` and `vercol` are defined in arc minutes to allow for specification of guide supermirror `m` values as negative numbers.
See `getαs` and `getβs` for the conversion from `m`/arc minutes to radian FWHM.
"""
type TripleAxis{A<:Real,R<:Real,H<:Real,V<:Real}
    a::Array{A,1} # array of angles; at least [θmono,2θmono,θsample,2θsample,θanalyzer,2θanalyzer] in order -- extra angles can follow
    source::Source
    mono::Bragg
    monitor::Detector
    sample::Sample
    # environment::SampleEnvironment
    ana::Bragg
    detector::Detector
    arms::Array{R,1}   # distance information
    horcol::Array{H,1} # horizontal collimation FWHM in arc minutes
    vercol::Array{V,1} #   vertical collimation FWHM in arc minutes
    energy::Real
    kforki::Bool # false==kf fixed, true==ki fixed
    extra::Nullable # for any extra information that the TripleAxis might need (pressure, field, temperature, ???)
    # initializer with vector length checks
    function TripleAxis{Ai,Ri,Hi,Vi}(angles::Array{Ai,1},src::Source,m::Bragg,moni::Detector,
                                     sam::Sample,n::Bragg,dect::Detector,distances::Array{Ri,1},
                                     hcol::Array{Hi,1},vcol::Array{Vi,1},energy::Real,kforki::Bool,ext::Nullable) where {Ai<:Real,Ri<:Real,Hi<:Real,Vi<:Real}
        @assert length(angles)>=6 "At least six angles are required, [θmono,2θmono,θsample,2θsample,θanalyzer,2θanalyzer] in order. Extra angles can follow."
        @assert length(distances)>=4 "At least four distances are required, [source-mono,mono-sample,sample-analyzer,analyzer-detector,(mono-monitor)]."
        @assert length(hcol)==4 "Four horizontal collimations are required."
        @assert length(vcol)==4 "Four vertical collimations are required."
        new(angles,src,m,moni,sam,n,dect,distances,hcol,vcol,energy,kforki,ext)
    end
end
function Base.copy(t::TripleAxis)
    a=copy(t.a)
    source=copy(t.source)
    mono=copy(t.mono)
    monitor=copy(t.monitor)
    sample=copy(t.sample)
    ana=copy(t.ana)
    detector=copy(t.detector)
    arms=copy(t.arms)
    horcol=copy(t.horcol)
    vercol=copy(t.vercol)
    energy=copy(t.energy)
    kforki=copy(t.kforki)
    extra=copy(t.extra)
    TripleAxis(a,source,mono,monitor,sample,ana,detector,arms,horcol,vercol,energy,kforki,extra)
end
function TripleAxis{T<:Real,R<:Real,S<:Real,V<:Real}(a::Array{T,1}=[0,0,0,0,0,0.],
                            src::Source=Source(),mono::Bragg=Bragg(),mn::Detector=Detector(),
                            smpl::Sample=Sample(),ana::Bragg=Bragg(),dt::Detector=Detector(),
                            ar::Array{R,1}=[Inf,Inf,Inf,Inf],
                            hc::Array{S,1}=[1,1,1,1]/60*pi,vc::Array{V,1}=[1,1,1,1]/60*pi,
                            en::Real=5.5,kb::Bool=false,ext::Nullable=Nullable())
    # use promote(...) to find the common type
#    (en,rest)=promote(en,a[1],ar[1],hc[1],vc[1])
#    t=typeof(en) # and map([common type],vector) to convert the vectors
#    TripleAxis{t}(map(t,a),src,mono,mn,smpl,ana,dt,map(t,ar),map(t,hc),map(t,vc),en,kb)
    TripleAxis{T,R,S,V}(a,src,mono,mn,smpl,ana,dt,ar,hc,vc,en,kb,ext)
end
TripleAxis(src::Source,o...)=TripleAxis([0,0,0,0,0,0.],src,o...)

fmtangle(angle)=@sprintf("%7.1f",angle)
#fmtvecpt(vecpt)=@sprintf("%+7.3f",vecpt)
function showTripleAxis(io::IO,t::TripleAxis,compact::Bool=false)
    ~compact && (Base.print(summary(t));Base.println(":"))
    ad=t.a/pi*180
    compact?(Base.print(io,"[");for ang in ad; Base.print(io,fmtangle(ang)); end;Base.print(io,"]°")):Base.println(io,"  angles: $(t.a/pi*180)°")
    ~compact && (Base.print(io,"  source: ");Base.show(io,t.source);Base.println(io,""))
    compact?(Base.showcompact(io,t.mono)):Base.println(io,"  mono: $(t.mono)")
    ~compact && (Base.print(io,"  monitor: ");Base.show(io,t.monitor);Base.println(io,""))
    compact?(Base.print(io,"{");Base.showcompact(io,t.sample);Base.print(io,"}")):(Base.print(io,"  samp: ");Base.showcompact(io,t.sample);Base.println(io,""))
    compact?(Base.showcompact(io,t.ana)):Base.println(io,"  ana.: $(t.ana)")
    ~compact && (Base.print(io,"  detector: ");Base.show(io,t.detector);Base.println(io,""))
    (~compact&&~all(map(x->isinf(x),t.arms))) && (Base.println(io,"  distances: $(t.arms)"))
    ~compact && (Base.print(io,"  collimation"))
    ~compact && (Base.print(io,"  horizontal: $(t.horcol)'"))
    ~compact && (Base.println(io,"  vertical: $(t.vercol)'"))
    compact?(Base.print(io," E"*(iskiFixed(t)?"ⁱ":"ᶠ")*"=");Base.showcompact(io,t.energy);Base.print(io,"meV")) : (Base.print(io,iskiFixed(t)?"  initial":"  final");Base.println(io," energy fixed at $(t.energy) meV"))
end
Base.show(io::IO,t::TripleAxis)=showTripleAxis(io,t,false)
Base.showcompact(io::IO,t::TripleAxis)=showTripleAxis(io,t,true)
sameLattice(a::TripleAxis,b::Union{Lattice,LatticeVector,Crystal,Sample})=sameLattice(a.sample,b)
sameLattice(b::Union{Lattice,LatticeVector,Crystal,Sample},a::TripleAxis)=sameLattice(b,a.sample)

function sameInstrument(a::TripleAxis,b::TripleAxis)
#    fields=[:energy,:source,:mono,:monitor,:sample,:ana,:detector,:arms,:horcol,:vercol,:kforki]
    fields=[:source,:mono,:monitor,:sample,:ana,:detector,:arms,:horcol,:vercol,:kforki] # FIXME add back in the energy check, and change it to be a Measured
    all([testfield(==,fld,a,b) for fld in fields])
end
(==)(a::TripleAxis,b::TripleAxis)=all(a.a.==b.a)&sameInstrument(a,b)
#
calck(τ::Real,twoθ::Real)=τ/sqrt(2.0-2.0*cos(twoθ))
# return the tuple (ki,kf) from specified (a2,a6) angles [or current values stored in the TripleAxis object t]
getki(t::TripleAxis,a2::Real=t.a[2])=calck(t.mono.tau,a2)
getkf(t::TripleAxis,a6::Real=t.a[6])=calck(t.ana.tau,a6)
# return the tuple (ei,ef) from specified E; relies upon fixed ki|kf information in t
calcEs(t::TripleAxis,E::Real)=iskiFixed(t)?(t.energy+0*E,t.energy-E):(t.energy+E,t.energy-0*E) # 0*E ensures both elements of tuple have the same promoted type
calcEi(t::TripleAxis,E::Real)=calcEs(t,E)[1]
calcEf(t::TripleAxis,E::Real)=calcEs(t,E)[2]
# The angle between ki and Q is given by tan(Φ)=[-kf sin(2θsample)]/[ki - kf cos(2θsample)]
getphi(ki::Real,kf::Real,a4::Real)=atan2(-kf*sin(a4),ki-kf*cos(a4))
#getphi(ki::Real,kf::Real,a4::Real)=atan2(kf*sin(a4),ki-kf*cos(a4))  #sin(a4)*scatteringsense(sample)?
getphi(t::TripleAxis)=getphi(getki(t),getkf(t),t.a[4])
function getphi{T<:TripleAxis}(ta::Array{T})
    phi1=getphi(ta[1])
    allphi=Array{typeof(phi1)}(size(ta)...)
    allphi[1]=phi1
    for i=2:length(ta); allphi[i]=getphi(ta[i]); end
    return allphi
end
#
function getTh2Th(s::Sample,Q::Lattice3Vector,ki::Real,kf::Real)
    (isDirect(Q)) && (warn("Passed direct- instead of reciprocal-lattice vector, reciprocating via star.");Q=star(Q))
    # the sample two theta can be calculated from the law of cosines, ki²+kf²=Q²+2 ki kf cos2θ, thanks to momentum conservation
    cos2θ = (ki^2+kf^2-norm(Q)^2)/2/ki/kf
    @assert abs(cos2θ)<=1 "The scattering triangle can not be closed for ki=$ki, kf=$kf, Q=$Q [cos2θ=$cos2θ]"
    tθ=acos(cos2θ)*getsense(s)
#    tθ=acos((ki^2+kf^2-norm(Q)^2)/2/ki/kf)*getsense(s)
    # sample2lab returns Q in the sample-table space,v=a*x+b*y+c*z, where c should be 0 for a triple-axis accessible Q
    x,y,z=(sample2lab(Q,s)...) # convert the Vector to a Tuple for asignment
    #@assert isa(z,Measured)?abs(z)<=Measured(0,2e-4):abs(z)<=2e-4 "Q has a component out of the scattering plane with magnitude $z Å⁻¹" ### Measured_asym.jl
    @assert isa(z,Measured)?abs(z)<=Measured(0,4e-8):abs(z)<=2e-4 "Q has a component out of the scattering plane with magnitude $z Å⁻¹"  ### Measured.jl
    # the sample angle depends on the angle between Q and ki (Φ), and the angle between Q and the sample orienting vector x
    θ=getphi(ki,kf,tθ)-atan2(y,x)
    (θ,tθ)
end
 #
iskiFixed(t::TripleAxis)=t.kforki
iskfFixed(t::TripleAxis)=~t.kforki
#
#const ħ²2mₙ=2.072142 # ħ²/(2mₙ) = 2.072142 meV Å² (ħ,² and ₙ used to avoid naming conflicts)
#
# return Q from specified (a2,a3,a4,a6) angles [or current values stored in the TripleAxis object t]
getQ(t::TripleAxis,a2::Real=t.a[2],a3::Real=t.a[3],a4::Real=t.a[4],a6::Real=t.a[6])=getQ_from_kikf(t,getki(t,a2),a3,a4,getkf(t,a6))
function getQ_from_kikf(t::TripleAxis,ki::Real,a3::Real,a4::Real,kf::Real)
    modQ=sqrt(ki^2+kf^2-2*ki*kf*cos(a4))
    Ψ=getphi(ki,kf,a4)-a3 # angle from first orienting vector to Q (phi is the angle from ki to Q)
#    Ψ-=pi/2 ## FIXME is this 90 degree shift right?
    lab2sample([cos(Ψ),sin(Ψ),0]*modQ,t.sample)
end
# getQi and getQf calculate Q positions for incoherent scattering with k=ki or k=kf, respectively
get_incoi_Q(t::TripleAxis)=getQ(t,t.a[2],t.a[3],t.a[4],t.a[2])
get_incof_Q(t::TripleAxis)=getQ(t,t.a[6],t.a[3],t.a[4],t.a[6])

# return E from specified (a2,a6) angles [or current values stored in the TripleAxis object t]
getEi(t::TripleAxis,a2::Real=t.a[2])=ħ²2mₙ*getki(t,a2)^2
getEf(t::TripleAxis,a6::Real=t.a[6])=ħ²2mₙ*getki(t,a6)^2
getE(t::TripleAxis,a2::Real=t.a[2],a6::Real=t.a[6])=ħ²2mₙ*(getki(t,a2)^2-getkf(t,a6)^2)
# Array versions of getQ,getEi,getEf,getE
for z in (:getEi,:getEf,:getE,:getQ,:get_incoi_Q,:get_incof_Q)
    @eval begin
        function $z{T<:TripleAxis}(v::Array{T})
            ot1=$z(v[1])
            out=Array{typeof(ot1)}(size(v)...)
            out[1]=ot1
            for j=2:length(v); out[j]=$z(v[j]); end
            return out
        end
    end
end

# return the tuple (Q,E) from specified (a2,a3,a4,a6) angles [or current values stored in the TripleAxis object t]
getQE(t::TripleAxis,a2::Real=t.a[2],a3::Real=t.a[3],a4::Real=t.a[4],a6::Real=t.a[6])=(getQ(t,a2,a3,a4,a6),getE(t,a2,a6))

getangle(t::TripleAxis,n::Integer)=t.a[n]
#function getangle{T<:TripleAxis}(v::Vector{T},n::Integer)
#    out=Array{eltype(v[1].a)}(length(v))
#    for i=1:length(v)
#        out[i]=v[i].a[n]
#    end
#    return out
#end
function getangle{T<:TripleAxis}(m::Array{T},n::Integer)
    out=Array{eltype(m[1].a)}(size(m)...)
    for i=1:length(m); out[i]=m[i].a[n]; end
    return out
end

getangles(t::TripleAxis)=(t.a[1:6]...)
function getangles{T<:TripleAxis}(v::Vector{T})
    out=Array{eltype(v[1].a)}(length(v),6)
    for i=1:length(v)
        out[i,:]=v[i].a[1:6]
    end
    return out
end
calcAngles(t::TripleAxis)=getangles(t)
function calcAngles(t::TripleAxis,Q::Lattice3Vector,E::Real=0.)
    # calcEs returns (Ei,Ef), collect converts it to [Ei,Ef],
    # ks=sqrt([Ei,Ef]/(ħ²/(2mₙ))), (ks...) splits the vector into a tuple for asignment to "ki,kf"
    Ei,Ef=calcEs(t,E)
    ki=sqrt(Ei/ħ²2mₙ)
    kf=sqrt(Ef/ħ²2mₙ)
    θm,tθm=getTh2Th(t.mono,ki)
    θs,tθs=getTh2Th(t.sample,Q,ki,kf)
    θa,tθa=getTh2Th(t.ana,kf)
    (θm,tθm,θs,tθs,θa,tθa)
end # this seems to work

#### in-place gotoQE!
function gotoQE!(t::TripleAxis,Q::Lattice3Vector,E::Real)
#    Ta=eltype(t.a);
    reta=collect(calcAngles(t,Q,E)) # collect converts the returned tuple into a vector
#    Tr=eltype(reta)
#    if Ta!=Tr && length(t.a)==length(reta)
#        t.a=reta
#    elseif Ta!=Tr && Tr<:Measured
#        # normally, t.a[1:6]=reta should suffice (with appropriate convert statements defined)
#        # but here reta might be a Measured array and I'm trying to avoid defining convert(Float64,Measured)
#        # because I want to avoid throwing away uncertainty information if at all possible
#        info("Discarding uncertainties in gotoQE! for angles $reta")
#        t.a[1:6]=value(reta)
#    else
        t.a[1:6]=reta
#    end
end
gotoQE!(t::TripleAxis,QE::Lattice4Vector)=gotoQE!(t,get3vector(QE),getE(QE))
## vector overloading: Vector{TripleAxis}, Vector{Lattice4Vector}
gotoQE!{T,L<:Lattice3Vector,R<:Real}(t::Array{T,1},Q::Array{L,1},E::Array{R,1})=(@assert compatible(t,Q,E);map((x,y,z)->gotoQE!(x,y,z),t,Q,E);)
gotoQE!{T,L<:Lattice4Vector}(t::Array{T,1},QE::Array{L,1})=(@assert compatible(t,QE);map((x,y)->gotoQE!(x,y),t,QE);)
#### copy gotoQE
gotoQE(t::TripleAxis,QE::Lattice4Vector)=(a=deepcopy(t);gotoQE!(a,QE);a) # copy and gotoQE
## single TripleAxis, multiple Lattice3Vectors&Real or Lattice4Vectors
#gotoQE{L<:Lattice4Vector}(t::TripleAxis,QE::Array{L,1})=map(x->gotoQE(t,x),QE)
gotoQE{L<:Lattice4Vector,N}(t,QE::Array{L,N})=map(x->gotoQE(t,x),QE)
#gotoQE{L<:Lattice3Vector,R<:Real}(t::TripleAxis,Q::Array{L,1},E::Array{R,1})=(@assert compatible(Q,E); map((x,y)->gotoQE(t,x,y),Q,E))
gotoQE{L<:Lattice3Vector,R<:Real,N}(t,Q::Array{L,N},E::Array{R,N})=(@assert compatible(Q,E); map((x,y)->gotoQE(t,x,y),Q,E))
## equal shaped TripleAxis and Lattice3Vector&Real or Lattice4Vector
gotoQE{T,L<:Lattice4Vector,N}(t::Array{T,N},QE::Array{L,N})=(@assert compatible(t,QE);map(gotoQE,t,QE))
gotoQE{T,L<:Lattice3Vector,R<:Real,N}(t::Array{T,N},Q::Array{L,N},E::Array{R,N})=(@assert compatible(t,Q,E);map(gotoQE,t,Q,E))
##
#### in-place gotoQ! (E=0)
gotoQ!(t::TripleAxis,Q::Lattice3Vector)=gotoQE!(t,Q,0.)
#gotoQ!{T<:TripleAxis,L<:Lattice3Vector}(t::Array{T,1},Q::Array{L,1})=(@assert compatible(t,Q); map(gotoQ!,t,Q))
gotoQ!{T,L<:Lattice3Vector,N}(t::Array{T,N},Q::Array{L,N})=(@assert compatible(t,Q); map(gotoQ!,t,Q))
#### copy gotoQ (E=0)
gotoQ(t,Q::Lattice3Vector)=(a=deepcopy(t);gotoQ!(a,Q);a)
gotoQ{L<:Lattice3Vector}(t,Q::Array{L})=map(x->gotoQ(t,x),Q)
#gotoQ{T<:TripleAxis,L<:Lattice3Vector}(t::Array{T,1},Q::Array{L,1})=(@assert compatible(t,Q); map(gotoQ,t,Q))
gotoQ{T,L<:Lattice3Vector,N}(t::Array{T,N},Q::Array{L,N})=(@assert compatible(t,Q); map(gotoQ,t,Q))

setangles!(t::TripleAxis,angles::Array{T,1}) where T<:Real =(@assert length(angles)==6; t.a[1:6]=angles[:])
setangles(t::TripleAxis,angles::Array{T,1}) where T<:Real =(a=deepcopy(t); setangles!(a,angles); a)
function setangles{T<:Real}(t::TripleAxis,angles::Array{T,2})
    @assert size(angles,1)==6
    out=[setangles(t,angles[:,1])]
    for i=2:size(angles,2); push!(out,setangles(t,angles[:,i])); end
    out::Array{typeof(t),1}
end
function setangles{T<:Real}(t::TripleAxis,angles::Array{T,3})
    @assert size(angles,1)==6
    out=cat(2,setangles(t,angles[:,:,1])) # make sure out is a Array{TripleAxis,2}
    for i=2:size(angles,3); out=cat(2,out,setangles(t,angles[:,:,i])); end
    out::Array{typeof(t),2}
end
function setangles{T<:Real}(t::TripleAxis,angles::Array{T,4})
    @assert size(angles,1)==6
    out=cat(3,setangles(t,angles[:,:,:,1])) # make sure out is a Array{TripleAxis,2}
    for i=2:size(angles,4); out=cat(3,out,setangles(t,angles[:,:,:,i])); end
    out::Array{typeof(t),3}
end

getmonoHR(t::TripleAxis)=getHR(t.mono,t.arms[1],t.arms[2],t.a[1])
getmonoVR(t::TripleAxis)=getVR(t.mono,t.arms[1],t.arms[2],t.a[1])
getanaHR(t::TripleAxis)=getHR(t.ana,t.arms[3],t.arms[4],t.a[5])
getanaVR(t::TripleAxis)=getVR(t.ana,t.arms[3],t.arms[4],t.a[5])

# TripleAxis object operator overloading: (probably) only used for interpolation between points in a scan
# Overload +,-,/,* for modifications to the current (Q,E) vector
+(a::TripleAxis,b::TripleAxis)=(@assert sameInstrument(a,b); gotoQE(a,get4vector(a)+get4vector(b)))
-(a::TripleAxis,b::TripleAxis)=(@assert sameInstrument(a,b); gotoQE(a,get4vector(a)-get4vector(b)))
+(a::TripleAxis,b::LatticeVector)=(@assert sameLattice(a,b); gotoQE(a,get4vector(a)+b))
+(a::LatticeVector,b::TripleAxis)=(@assert sameLattice(a,b); gotoQE(b,a+get4vector(b)))
-(a::TripleAxis,b::LatticeVector)=(@assert sameLattice(a,b); gotoQE(a,get4vector(a)-b))
-(a::LatticeVector,b::TripleAxis)=(@assert sameLattice(a,b); gotoQE(b,a-get4vector(b)))

# Scatterd {+,-,*,/} Number
for y in [:*,:/]
    @eval Base.$y(a::TripleAxis,b::Number)=gotoQE(a,Base.$y(get4vector(a),b))
end
# commutative operator
(*)(b::Number,a::TripleAxis)=(*)(a,b)

# overload mean to divide first instead of last to ensure all Q positions remain accessible (TripleAxis checks this!)
Base.mean{T<:TripleAxis}(a::T;k...)=Base.copy(a)
# This approach works just fine at zero energy transfer, but it's possible to
# have points at finite energy transfer for which, e.g., (Q,E)/2 is not an
# accessible point in reciprocal space.
#Base.mean{T<:TripleAxis}(v::Vector{T};k...)=Base.sum(v/Base.length(v);k...)
#Base.mean{T<:TripleAxis}(m::Array{T},d;k...)=Base.sum(m/Base.size(m,d);k...)
function Base.mean(v::Vector{T};k...) where T<:TripleAxis
    @assert all(sameInstrument(x,v[1]) for x in v)
    gotoQE(v[1],mean(get4vector.(v)))
end
# function Base.mean(m::Array{T},d;k...) where  T<:TripleAxis
#     @assert all(sameInstrument(x,m[1]) for x in m)
#     # reducing over a given dimension d requires more work
# end
function Base.mean(m::Array{T};k...) where  T<:TripleAxis
    @assert all(sameInstrument(x,m[1]) for x in m)
    gotoQE(m[1],mean(get4vector.(m)))
end


# Now it should be possible to write a function that inserts the convoluted function between measurement points in order to produce a smooth convoluted function line for plotting
# The simple case will only work for scans that are along defined directions in (Q,E) space.
# Interpolating rocking scans, which generally take a curved path through (Q,E) space, might give wacky results if the step size is large compared to \delta(Q,E)
# To help with such a special case, Abuse *,/ to overload +,- but acting on the angles for each TripleAxis
*(a::TripleAxis,b::TripleAxis)=(@assert sameInstrument(a,b); c=deepcopy(a); c.a=a.a+b.a; c)
/(a::TripleAxis,b::TripleAxis)=(@assert sameInstrument(a,b); c=deepcopy(a); c.a=a.a-b.a; c)


function combineRepeatedInstruments(v::AbstractArray{T},bno::AbstractArray{I},nb::Integer=maximum(bno);k...) where {T<:TripleAxis,I<:Integer}
    @assert size(v)==size(bno)
    [mean(v[bno.==i];k...) for i in 1:nb]
end



_FWHM2σ(c)=c./sqrt(8log(2))
_σ2FWHM(c)=c.*sqrt(8log(2))
""" Numbers expressed as FWHM can be related to a Gaussian variance, σ, by FWHM=sqrt(8ln2)σ. """
_FWHM2σ, _σ2FWHM
"""
Conversion from a guide `m` value to an angular standard deviation is wavelength dependent,
since the critical angle, `Θᶜ`, is wavelength dependent via
    Θᶜ = m Θ₀ λ
where `Θ₀` is the critical angle of natural Ni (0.1°/Å, or π/1800 radian/Å). Angular deviations
from -Θᶜ to +Θᶜ with equal probability gives a top-hat distribution, which has an angular standard
deviation of Θᶜ/sqrt(3).
"""
function _m2σ!{T<:Number}(t::TripleAxis,c::Array{T}) # convert guide m values to angular standard deviations (which are wavelength dependent)
    m=abs.(c)
    #Θ₀=pi/180/10 # the critical angle of natural Ni (in radian per angstrom)
    Θᶜ=(pi/1800)*m.*(2pi./[getki(t),getki(t),getkf(t),getkf(t)]) # the critical angle is mΘ₀λ
    σ=Θᶜ/sqrt(3) # and the angular standard deviation of a top-hat distribution between +/- Θᶜ is Θᶜ/sqrt(3)
    c[c.<0]=σ[c.<0]
end
"""
Given a `TripleAxis` object and an `Array` of numbers which are either FWHM angular deviations in
arc minutes or guide `m` values dependent on their sign (FWHM≥0, m<0), calculate the variance in
radian.
The FWHM angular deviations can be related to a Gaussian variance, σ, by FWHM=sqrt(8ln2)σ.
While conversion from a guide `m` value to an angular standard deviation is wavelength dependent,
since the critical angle, `Θᶜ`, is wavelength dependent via `Θᶜ = m Θ₀ λ`
where `Θ₀` is the critical angle of natural Ni (0.1°/Å, or π/1800 radian/Å). Angular deviations
from `-Θᶜ` to `+Θᶜ` with equal probability give a top-hat distribution, which has an angular standard
deviation of `Θᶜ/sqrt(3)`.
"""
function _arcmin_m_2σ{T<:Number}(t,c::Array{T})
    (c,crad)=promote(c,c/10800*pi) # _FWHM2σ will return a float, ensure c[c.>0]=... doesn't throw an InexactError()
    c[c.>0] = _FWHM2σ(crad[c.>0])
    any(c.<0) && (_m2σ!(t,c))
    c
end
"""
    getαs(t::TripleAxis)
Returns the horizontal guide/collimator-defined angular standard deviations in radian for a `TripleAxis` object.
See `_arcmin_m_2σ` for details.
"""
getαs(t)=_arcmin_m_2σ(t,deepcopy(t.horcol)) # t.horcol should be the horizontal collimation FWHM in arc minutes
"""
    getβs(t::TripleAxis)
Returns the vertical guide/collimator-defined angular standard deviations in radian for a `TripleAxis` object.
See `_arcmin_m_2σ` for details.
"""
getβs(t)=_arcmin_m_2σ(t,deepcopy(t.vercol)) # t.horcol should be the horizontal collimation FWHM in arc minutes

function geometricαs(t::TripleAxis)
    w=[getwidth(tasp.mono)+getwidth(tasp.source),getwidth(tasp.sample)+getwidth(tasp.mono),getwidth(tasp.ana)+getwidth(tasp.sample),getwidth(tasp.detector)+getwidth(tasp.ana)]/2
    atan2(w,tasp.arms[1:4])
end
function geometricβs(t::TripleAxis)
    w=[getheight(tasp.mono)+getheight(tasp.source),getheight(tasp.sample)+getheight(tasp.mono),getheight(tasp.ana)+getheight(tasp.sample),getheight(tasp.detector)+getheight(tasp.ana)]/2
    atan2(w,tasp.arms[1:4])
end
geteffictiveαs(t)=min(geometricαs(t),_arcmin_m_2σ(t,deepcopy(t.horcol))) # t.horcol should be the horizontal collimation FWHM in arc minutes
geteffectiveβs(t)=min(geometricβs(t),_arcmin_m_2σ(t,deepcopy(t.vercol))) # t.horcol should be the horizontal collimation FWHM in arc minutes
"""
Return effective angular standard deviations taking into account guides, collimators, and
physical dimensions of a `TripleAxis` object.
`geteffictiveαs` and `geteffectiveβs` are separate from `getαs` and `getβs` since it would be
inappropriate to use the effective angular standard deviations with the `Popovici` method as doing
so would double-count the influence of shape effects. Effective collimations *could* be used with
the `Cooper-Nathans` method, but this seems ill advised since the only advantage to that method
(a requisite knowledge of less parameters) is lost if one uses effective collimations.
"""
geteffictiveαs, geteffectiveβs

# generalized functions to access mosaics expressed as Gaussian variances
getηs(a::Union{Bragg,Sample})=_FWHM2σ(a.mosaic)
getηh(a::Union{Bragg,Sample})=_FWHM2σ(a.mosaic[1])
getηv(a::Union{Bragg,Sample})=_FWHM2σ(a.mosaic[2])

# Eventually each triple-axis Scatterd (TASPd, EIGERd, etc.) object will hold an Array of TripleAxis objects
# It will be advantagous to access at least some properties of the TripleAxis objects directly from the Array
getx(a::TripleAxis)=getx(a.sample);
getx{T<:TripleAxis}(v::Array{T,1})=(rx=[getx(v[1])];  for i=2:length(v);   push!(rx,getx(v[i]));   end; rx)
getx{T<:TripleAxis}(m::Array{T,2})=(rx= getx(m[:,1]); for i=2:size(m,2); rx=hcat(rx,getx(m[:,i])); end; rx)
gety(a::TripleAxis)=gety(a.sample);
gety{T<:TripleAxis}(v::Array{T,1})=(ry=[gety(v[1])];  for i=2:length(v);   push!(ry,gety(v[i]));   end; ry)
gety{T<:TripleAxis}(m::Array{T,2})=(ry= gety(m[:,1]); for i=2:size(m,2); ry=hcat(ry,gety(m[:,i])); end; ry)
getz(a::TripleAxis)=getz(a.sample);
getz{T<:TripleAxis}(v::Array{T,1})=(rz=[getz(v[1])];  for i=2:length(v);   push!(rz,getz(v[i]));   end; rz)
getz{T<:TripleAxis}(m::Array{T,2})=(rz= getz(m[:,1]); for i=2:size(m,2); rz=hcat(rz,getz(m[:,i])); end; rz)
#get3vector(t::TripleAxis)=getQ(t);
#get3vector{T<:TripleAxis}(v::Array{T,1})=(rQ=[get3vector(v[1])]; for i=2:length(v);   push!(rQ,get3vector(v[i]));  end; rQ)
#get3vector{T<:TripleAxis}(v::Array{T,2})=(rQ= get3vector(v[:,1]);for i=2:size(v,2); rQ=hcat(rQ,get3vector(v[:,i])); end; rQ)
get3vector{T<:TripleAxis}(t::Union{T,Array{T}})=getQ(t) # get3vector and getQ are identical.
#get4vector(t::TripleAxis)=((Q,E)=getQE(t); Lattice4Vector{typeof(geth(Q))}(Q,E))
get4vector(t::TripleAxis)=((Q,E)=getQE(t); Lattice4Vector(Q,E))
get4vector{T<:TripleAxis}(v::Array{T,1})=(rQ=[get4vector(v[1])]; for i=2:length(v);   push!(rQ,get4vector(v[i]));  end; rQ)
get4vector{T<:TripleAxis}(v::Array{T,2})=(rQ= get4vector(v[:,1]);for i=2:size(v,2); rQ=hcat(rQ,get4vector(v[:,i])); end; rQ)

geth{T<:TripleAxis}(t::Union{T,Array{T}})=geth(getQ(t))
getk{T<:TripleAxis}(t::Union{T,Array{T}})=getk(getQ(t))
getl{T<:TripleAxis}(t::Union{T,Array{T}})=getl(getQ(t))
#getE{T<:TripleAxis}(t::Union{T,Array{T}})=getE(getQ(t)) # already defined based on ki and kf

get_incoi_h{T<:TripleAxis}(t::Union{T,Array{T}})=geth(get_incoi_Q(t))
get_incoi_k{T<:TripleAxis}(t::Union{T,Array{T}})=getk(get_incoi_Q(t))
get_incoi_l{T<:TripleAxis}(t::Union{T,Array{T}})=getl(get_incoi_Q(t))
get_incof_h{T<:TripleAxis}(t::Union{T,Array{T}})=geth(get_incof_Q(t))
get_incof_k{T<:TripleAxis}(t::Union{T,Array{T}})=getk(get_incof_Q(t))
get_incof_l{T<:TripleAxis}(t::Union{T,Array{T}})=getl(get_incof_Q(t))

# a convenience function
lab2sample{T<:Real}(Q::Array{T},t::TripleAxis,o...)=lab2sample(Q,t.sample,o...)
lab2sample!{T<:Real}(Q::Array{T},t::TripleAxis,o...)=lab2sample!(Q,t.sample,o...)
# another convenience function (likely not used much)
sample2lab(t::TripleAxis)=sample2lab(getQ(t),t.sample)

# to help in the later construction of matching-size arrays in fuctions:
size(::TripleAxis)=() # size([any scalar])=()
ndims(::TripleAxis)=0 # ndims([any scaler])=0



# Common documentation
""" `getradius(obj)` returns the physical radius of a circular object. """
function getradius end
""" `getwidth(obj)` returns the physical width of an object. """
function getwidth end
""" `getheight(obj)` returns the physical height of an object. """
function getheight end
""" `setradius(obj,R)` takes the physical radius of a circular object and appropriately sets <r^2>. """
function setradius end
""" `setwidth(obj,W)` takes the physical width of an object and appropriately sets <y^2>. """
function setwidth end
""" `setheight(obj,H)` takes the physical height of an object and appropriately sets <z^2>. """
function setheight end
