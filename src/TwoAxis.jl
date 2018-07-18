"""
The `TwoAxis` object contains all elements and angles for a single setting of
a two-axis diffractometer instrument.

| Fieldname | Description |
|:---------:|:------------|
| a         | `Vector{Real}` of the instrument angles. |
|           | at least [θm,2θm,θs,2θs] following extra angles are allowed |
| source    | `Source` object describing the physical extent of the incident beam |
| mono      | `Bragg` object describing the monochromator |
| sample    | `Sample` object describing the single-crystal sample |
| detector  | `Detector` object describing the physical extend of the detector |
| arms      | `Vector{Real}` of distances between elements |
|           | [source-mono, mono-sample, sample-detector, mono-monitor]|
| horcol    | `Vector{Real}` of horizontal collimators FWHM in arc minutes |
| vercol    | `Vector{Real}` of vertical collimators FWHM in arc minutes |
| energy    | `Real` value of the incident energy in meV |

*Note:* `horcol` and `vercol` are defined in arc minutes to allow for specification of guide supermirror `m` values as negative numbers.
See `getαs` and `getβs` for the conversion from `m`/arc minutes to radian FWHM.
"""
type TwoAxis{A<:Real,R<:Real,H<:Real,V<:Real}
    a::Array{A,1} # array of angles; at least [θmono,2θmono,θsample,2θsample] in order -- extra angles can follow
    source::Source
    mono::Bragg
    monitor::Detector
    sample::Sample
    # environment::SampleEnvironment
    detector::Detector
    arms::Array{R,1}   # distance information
    horcol::Array{H,1} # horizontal collimation FWHM in arc minutes
    vercol::Array{V,1} #   vertical collimation FWHM in arc minutes
    energy::Real
    extra::Nullable # for any extra information that the TwoAxis might need (pressure, field, temperature, ???)
    # initializer with vector length checks
    function TwoAxis{Ai,Ri,Hi,Vi}(angles::Array{Ai,1},src::Source,m::Bragg,moni::Detector,sam::Sample,dect::Detector,
            distances::Array{Ri,1},hcol::Array{Hi,1},vcol::Array{Vi,1},energy::Real,ext::Nullable) where {Ai<:Real,Ri<:Real,Hi<:Real,Vi<:Real}
        @assert length(angles)>3 "At least four angles are required, [θmono,2θmono,θsample,2θsample] in order. Extra angles can follow."
        @assert length(distances)>2 "At least three distances are required, [source-mono,mono-sample,sample-detector,(mono-monitor)]."
        @assert length(hcol)==3 "Three horizontal collimations are required."
        @assert length(vcol)==3 "Three vertical collimations are required."
        new(angles,src,m,moni,sam,dect,distances,hcol,vcol,energy,ext)
    end
end
function Base.copy(t::TwoAxis)
    a=copy(t.a)
    source=copy(t.source)
    mono=copy(t.mono)
    monitor=copy(t.monitor)
    sample=copy(t.sample)
    detector=copy(t.detector)
    arms=copy(t.arms)
    horcol=copy(t.horcol)
    vercol=copy(t.vercol)
    energy=copy(t.energy)
    extra=copy(t.extra)
    TwoAxis(a,source,mono,monitor,sample,detector,arms,horcol,vercol,energy,extra)
end
function TwoAxis{T<:Real,R<:Real,S<:Real,V<:Real}(a::Array{T,1}=[0,0,0,0.],
                            src::Source=Source(),mono::Bragg=Bragg(),mn::Detector=Detector(),
                            smpl::Sample=Sample(),dt::Detector=Detector(),
                            ar::Array{R,1}=[Inf,Inf,Inf],
                            hc::Array{S,1}=[1,1,1]/60*pi,vc::Array{V,1}=[1,1,1]/60*pi,
                            en::Real=5.5,ext::Nullable=Nullable())
    TwoAxis{T,R,S,V}(a,src,mono,mn,smpl,dt,ar,hc,vc,en,ext)
end
TwoAxis(src::Source,o...)=TwoAxis([0,0,0,0.],src,o...)

function showTwoAxis(io::IO,t::TwoAxis,compact::Bool=false)
    ~compact && (Base.print(summary(t));Base.println(":"))
    ad=t.a/pi*180
    compact?(Base.print(io,"[");for ang in ad; Base.print(io,fmtangle(ang)); end;Base.print(io,"]°")):Base.println(io,"  angles: $(t.a/pi*180)°")
    ~compact && (Base.print(io,"  source: ");Base.show(io,t.source);Base.println(io,""))
    compact?(Base.showcompact(io,t.mono)):Base.println(io,"  mono: $(t.mono)")
    ~compact && (Base.print(io,"  monitor: ");Base.show(io,t.monitor);Base.println(io,""))
    compact?(Base.print(io,"{");Base.showcompact(io,t.sample);Base.print(io,"}")):(Base.print(io,"  samp: ");Base.showcompact(io,t.sample);Base.println(io,""))
    ~compact && (Base.print(io,"  detector: ");Base.show(io,t.detector);Base.println(io,""))
    (~compact&&~all(map(x->isinf(x),t.arms))) && (Base.println(io,"  distances: $(t.arms)"))
    ~compact && (Base.print(io,"  collimation"))
    ~compact && (Base.print(io,"  horizontal: $(t.horcol)'"))
    ~compact && (Base.println(io,"  vertical: $(t.vercol)'"))
    compact?(Base.print(io," E=");Base.showcompact(io,t.energy);Base.print(io,"meV")) : (Base.print(io,"  incident energy fixed at $(t.energy) meV"))
end
Base.show(io::IO,t::TwoAxis)=showTwoAxis(io,t,false)
Base.showcompact(io::IO,t::TwoAxis)=showTwoAxis(io,t,true)
sameLattice(a::TwoAxis,b::Union{Lattice,LatticeVector,Crystal,Sample})=sameLattice(a.sample,b)
sameLattice(b::Union{Lattice,LatticeVector,Crystal,Sample},a::TwoAxis)=sameLattice(b,a.sample)

function sameInstrument(a::TwoAxis,b::TwoAxis)
#    fields=[:energy,:source,:mono,:monitor,:sample,:detector,:arms,:horcol,:vercol]
    fields=[:source,:mono,:monitor,:sample,:detector,:arms,:horcol,:vercol] # FIXME add back in the energy check, and change it to be a Measured
    all([testfield(==,fld,a,b) for fld in fields])
end
(==)(a::TwoAxis,b::TwoAxis)=all(a.a.==b.a)&sameInstrument(a,b)

getki(t::TwoAxis,a2::Real=t.a[2])=calck(t.mono.tau,a2)
getkf(t::TwoAxis,a2::Real=t.a[2])=calck(t.mono.tau,a2)
# return the tuple (ei,ef) from specified E; relies upon fixed ki|kf information in t
calcEi(t::TwoAxis)=t.energy
calcEf(t::TwoAxis)=t.energy
# The angle between ki and Q is given by tan(Φ)=[-kf sin(2θsample)]/[ki - kf cos(2θsample)]
getphi(t::TwoAxis)=getphi(getki(t),getki(t),t.a[4]) # call the version defined in TripleAxis.jl
function getphi{T<:TwoAxis}(ta::Array{T})
    phi1=getphi(ta[1])
    allphi=Array{typeof(phi1)}(size(ta)...)
    allphi[1]=phi1
    for i=2:length(ta); allphi[i]=getphi(ta[i]); end
    return allphi
end

iskiFixed(t::TwoAxis)=true
iskfFixed(t::TwoAxis)=false
#
#const ħ²2mₙ=2.072142 # ħ²/(2mₙ) = 2.072142 meV Å² (ħ,² and ₙ used to avoid naming conflicts)
#
# return Q from specified (a2,a3,a4) angles [or current values stored in the TripleAxis object t]
function getQ(t::TwoAxis,a2::Real=t.a[2],a3::Real=t.a[3],a4::Real=t.a[4])
    k=getki(t,a2)
    modQ=k*sqrt(2-2*cos(a4))
    Ψ=getphi(k,k,a4)-a3 # angle from first orienting vector to Q (phi is the angle from ki to Q)
#    Ψ-=pi/2 ## FIXME is this 90 degree shift right?
    lab2sample([cos(Ψ),sin(Ψ),0]*modQ,t.sample)
end

# return E from specified (a2,a6) angles [or current values stored in the TripleAxis object t]
getEi(t::TwoAxis,a2::Real=t.a[2])=ħ²2mₙ*getki(t,a2)^2
getEf(t::TwoAxis,a2::Real=t.a[2])=ħ²2mₙ*getki(t,a2)^2
getE(t::TwoAxis,a2::Real=t.a[2])=ħ²2mₙ*(getki(t,a2)^2-getki(t,a2)^2) # always zero
# Array versions of getQ,getEi,getEf,getE
for z in (:getEi,:getEf,:getE,:getQ)
    @eval begin
        function $z{T<:TwoAxis}(v::Array{T})
            ot1=$z(v[1])
            out=Array{typeof(ot1)}(size(v)...)
            out[1]=ot1
            for j=2:length(v); out[j]=$z(v[j]); end
            return out
        end
    end
end

# return the tuple (Q,E) from specified (a2,a3,a4) angles [or current values stored in the TwoAxis object t]
getQE(t::TwoAxis,a2::Real=t.a[2],a3::Real=t.a[3],a4::Real=t.a[4])=(getQ(t,a2,a3,a4),getE(t,a2))

getangle(t::TwoAxis,n::Integer)=t.a[n]
function getangle{T<:TwoAxis}(m::Array{T},n::Integer)
    out=Array{eltype(m[1].a)}(size(m)...)
    for i=1:length(m); out[i]=m[i].a[n]; end
    return out
end

getangles(t::TwoAxis)=(t.a[1:4]...)
function getangles{T<:TwoAxis}(v::Vector{T})
    out=Array(eltype(v[1].a),length(v),4)
    for i=1:length(v)
        out[i,:]=v[i].a[1:4]
    end
    return out
end
calcAngles(t::TwoAxis)=getangles(t)
function calcAngles(t::TwoAxis,Q::Lattice3Vector,E::Real=0.)
    # calcEs returns (Ei,Ef), collect converts it to [Ei,Ef],
    # ks=sqrt([Ei,Ef]/(ħ²/(2mₙ))), (ks...) splits the vector into a tuple for asignment to "ki,kf"
    Ei=calcEi(t)
    ki=sqrt(Ei/ħ²2mₙ)
    θm,tθm=getTh2Th(t.mono,ki)
    θs,tθs=getTh2Th(t.sample,Q,ki,ki)
    (θm,tθm,θs,tθs)
end # this seems to work

#### in-place gotoQE!
function gotoQE!(t::TwoAxis,Q::Lattice3Vector,E::Real=0)
    reta=collect(calcAngles(t,Q,E)) # collect converts the returned tuple into a vector
    t.a[1:6]=reta
end
gotoQE!(t::TwoAxis,QE::Lattice4Vector)=gotoQE!(t,get3vector(QE),getE(QE))
#### in-place gotoQ! (E=0)
gotoQ!(t::TwoAxis,Q::Lattice3Vector)=gotoQE!(t,Q,0.)

setangles!{T<:Real}(t::TwoAxis,angles::Array{T,1})=(@assert length(angles)==4; t.a[1:4]=angles[:])
setangles(t::TwoAxis,o...)=(a=deepcopy(t); setangles!(a,o...); a)
function setangles{T<:Real}(t::TwoAxis,angles::Array{T,2})
    @assert size(angles,1)==4
    out=[setangles(t,angles[:,1])]
    for i=2:size(angles,2); push!(out,setangles(t,angles[:,i])); end
    out::Array{typeof(t),1}
end
function setangles{T<:Real}(t::TwoAxis,angles::Array{T,3})
    @assert size(angles,1)==4
    out=cat(2,setangles(t,angles[:,:,1])) # make sure out is a Array{TripleAxis,2}
    for i=2:size(angles,3); out=cat(2,out,setangles(t,angles[:,:,i])); end
    out::Array{typeof(t),2}
end

getmonoHR(t::TwoAxis)=getHR(t.mono,t.arms[1],t.arms[2],t.a[1])
getmonoVR(t::TwoAxis)=getVR(t.mono,t.arms[1],t.arms[2],t.a[1])

# TripleAxis object operator overloading: (probably) only used for interpolation between points in a scan
# Overload +,-,/,* for modifications to the current (Q,E) vector
+(a::TwoAxis,b::TwoAxis)=(@assert sameInstrument(a,b); gotoQE(a,get4vector(a)+get4vector(b)))
-(a::TwoAxis,b::TwoAxis)=(@assert sameInstrument(a,b); gotoQE(a,get4vector(a)-get4vector(b)))
+(a::TwoAxis,b::LatticeVector)=(@assert sameLattice(a,b); gotoQE(a,get4vector(a)+b))
+(a::LatticeVector,b::TwoAxis)=(@assert sameLattice(a,b); gotoQE(b,a+get4vector(b)))
-(a::TwoAxis,b::LatticeVector)=(@assert sameLattice(a,b); gotoQE(a,get4vector(a)-b))
-(a::LatticeVector,b::TwoAxis)=(@assert sameLattice(a,b); gotoQE(b,a-get4vector(b)))

# Scatterd {+,-,*,/} Number
for y in [:*,:/]
    @eval Base.$y(a::TwoAxis,b::Number)=gotoQE(a,Base.$y(get4vector(a),b))
    #@eval Base.$y{T<:TwoAxis}(a::Array{T},b::Number)=(out=Array{T}(size(a)); for i=1:length(a); out[i]=Base.$y(a[i],b); end; out)
end
# commutative operator
(*)(b::Number,a::TwoAxis)=(*)(a,b)

# overload mean to divide first instead of last to ensure all Q positions remain accessible (TripleAxis checks this!)
Base.mean{T<:TwoAxis}(a::T;k...)=Base.copy(a)
Base.mean{T<:TwoAxis}(v::Vector{T};k...)=Base.sum(v/Base.length(v);k...)
Base.mean{T<:TwoAxis}(m::Array{T},d;k...)=Base.sum(m/Base.size(m,d);k...)


# Now it should be possible to write a function that inserts the convoluted function between measurement points in order to produce a smooth convoluted function line for plotting
# The simple case will only work for scans that are along defined directions in (Q,E) space.
# Interpolating rocking scans, which generally take a curved path through (Q,E) space, might give wacky results if the step size is large compared to \delta(Q,E)
# To help with such a special case, Abuse *,/ to overload +,- but acting on the angles for each TripleAxis
*(a::TwoAxis,b::TwoAxis)=(@assert sameInstrument(a,b); c=deepcopy(a); c.a=a.a+b.a; c)
/(a::TwoAxis,b::TwoAxis)=(@assert sameInstrument(a,b); c=deepcopy(a); c.a=a.a-b.a; c)


function combineRepeatedInstruments(v::AbstractArray{T},bno::AbstractArray{I},nb::Integer=maximum(bno);k...) where {T<:TwoAxis,I<:Integer}
    @assert size(v)==size(bno)
    [mean(v[bno.==i];k...) for i in 1:nb]
end

function _m2σ!{T<:Number}(t::TwoAxis,c::Array{T}) # convert guide m values to angular standard deviations (which are wavelength dependent)
    m=abs(c)
    #Θ₀=pi/180/10 # the critical angle of natural Ni (in radian per angstrom)
    Θᶜ=(pi/1800)*m.*(2pi./[getki(t),getki(t)]) # the critical angle is mΘ₀λ
    σ=Θᶜ/sqrt(3) # and the angular standard deviation of a top-hat distribution between +/- Θᶜ is Θᶜ/sqrt(3)
    c[c.<0]=σ[c.<0]
end

function geometricαs(t::TwoAxis)
    w=[getwidth(tasp.mono)+getwidth(tasp.source),getwidth(tasp.sample)+getwidth(tasp.mono),getwidth(tasp.sample)+getwidth(tasp.detector)]/2
    atan2(w,tasp.arms[1:3])
end
function geometricβs(t::TwoAxis)
    w=[getheight(tasp.mono)+getheight(tasp.source),getheight(tasp.sample)+getheight(tasp.mono),getheight(tasp.sample)+getheight(tasp.detector)]/2
    atan2(w,tasp.arms[1:3])
end

# Eventually each triple-axis Scatterd (TASPd, EIGERd, etc.) object will hold an Array of TwoAxis objects
# It will be advantagous to access at least some properties of the TwoAxis objects directly from the Array
getx(a::TwoAxis)=getx(a.sample);
getx{T<:TwoAxis}(v::Array{T,1})=(rx=[getx(v[1])];  for i=2:length(v);   push!(rx,getx(v[i]));   end; rx)
getx{T<:TwoAxis}(m::Array{T,2})=(rx= getx(m[:,1]); for i=2:size(m,2); rx=hcat(rx,getx(m[:,i])); end; rx)
gety(a::TwoAxis)=gety(a.sample);
gety{T<:TwoAxis}(v::Array{T,1})=(ry=[gety(v[1])];  for i=2:length(v);   push!(ry,gety(v[i]));   end; ry)
gety{T<:TwoAxis}(m::Array{T,2})=(ry= gety(m[:,1]); for i=2:size(m,2); ry=hcat(ry,gety(m[:,i])); end; ry)
getz(a::TwoAxis)=getz(a.sample);
getz{T<:TwoAxis}(v::Array{T,1})=(rz=[getz(v[1])];  for i=2:length(v);   push!(rz,getz(v[i]));   end; rz)
getz{T<:TwoAxis}(m::Array{T,2})=(rz= getz(m[:,1]); for i=2:size(m,2); rz=hcat(rz,getz(m[:,i])); end; rz)
get3vector{T<:TwoAxis}(t::Union{T,Array{T}})=getQ(t) # get3vector and getQ are identical.
get4vector(t::TwoAxis)=((Q,E)=getQE(t); Lattice4Vector(Q,E))
get4vector{T<:TwoAxis}(v::Array{T,1})=(rQ=[get4vector(v[1])]; for i=2:length(v);   push!(rQ,get4vector(v[i]));  end; rQ)
get4vector{T<:TwoAxis}(v::Array{T,2})=(rQ= get4vector(v[:,1]);for i=2:size(v,2); rQ=hcat(rQ,get4vector(v[:,i])); end; rQ)

geth{T<:TwoAxis}(t::Union{T,Array{T}})=geth(getQ(t))
getk{T<:TwoAxis}(t::Union{T,Array{T}})=getk(getQ(t))
getl{T<:TwoAxis}(t::Union{T,Array{T}})=getl(getQ(t))

# a convenience function
lab2sample{T<:Real}(Q::Array{T},t::TwoAxis,o...)=lab2sample(Q,t.sample,o...)
lab2sample!{T<:Real}(Q::Array{T},t::TwoAxis,o...)=lab2sample!(Q,t.sample,o...)
# another convenience function (likely not used much)
sample2lab(t::TwoAxis)=sample2lab(getQ(t),t.sample)

# to help in the later construction of matching-size arrays in fuctions:
size(::TwoAxis)=() # size([any scalar])=()
ndims(::TwoAxis)=0 # ndims([any scaler])=0
