"""
The `TwoAxisBanana` object contains all elements and angles for a single setting of
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
type TwoAxisBanana{A<:Real,R<:Real,H<:Real,V<:Real}
    a1::A # θmono
    a2::A # 2θmono
    a3::A # θsample
    a4::Array{A,1} # array of 2θsample angles; [2θsample1,2θsample2,2θsample3,...] in order
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
    extra::Nullable # for any extra information that the TwoAxisBanana might need (pressure, field, temperature, ???)
    # initializer with vector length checks
    function TwoAxisBanana{Ai,Ri,Hi,Vi}(a1::Ai,a2::Ai,a3::Ai,a4::Array{Ai,1},src::Source,m::Bragg,moni::Detector,sam::Sample,dect::Detector,
            distances::Array{Ri,1},hcol::Array{Hi,1},vcol::Array{Vi,1},energy::Real,ext::Nullable) where {Ai<:Real,Ri<:Real,Hi<:Real,Vi<:Real}
        @assert length(a4)>0 "At least one 2θsample angle is required."
        @assert length(distances)>2 "At least three distances are required, [source-mono,mono-sample,sample-detector,(mono-monitor)]."
        @assert length(hcol)==3 "Three horizontal collimations are required."
        @assert length(vcol)==3 "Three vertical collimations are required."
        new(a1,a2,a3,a4,src,m,moni,sam,dect,distances,hcol,vcol,energy,ext)
    end
end
function Base.copy(t::TwoAxisBanana)
    a1=copy(t.a1)
    a2=copy(t.a2)
    a3=copy(t.a3)
    a4=copy(t.a4)
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
    TwoAxisBanana(a1,a2,a3,a4,source,mono,monitor,sample,detector,arms,horcol,vercol,energy,extra)
end
function TwoAxisBanana{T<:Real,R<:Real,S<:Real,V<:Real}(a1::T=zero(T),a2::T=zero(T),a3::T=zero(T),a4::Array{T,1}=[zero(T)],
                            src::Source=Source(),mono::Bragg=Bragg(),mn::Detector=Detector(),
                            smpl::Sample=Sample(),dt::Detector=Detector(),
                            ar::Array{R,1}=[Inf,Inf,Inf],
                            hc::Array{S,1}=[1,1,1]/60*pi,vc::Array{V,1}=[1,1,1]/60*pi,
                            en::Real=5.5,ext::Nullable=Nullable())
    TwoAxisBanana{T,R,S,V}(a1,a2,a3,a4,src,mono,mn,smpl,dt,ar,hc,vc,en,ext)
end
TwoAxisBanana(src::Source,o...)=TwoAxisBanana(0.,0.,0.,[0.],src,o...)

function showTwoAxis(io::IO,t::TwoAxisBanana,compact::Bool=false)
    ~compact && (Base.print(summary(t));Base.println(":"))
    ad=[t.a1;t.a2;t.a3;t.a4]/pi*180
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
Base.show(io::IO,t::TwoAxisBanana)=showTwoAxis(io,t,false)
Base.showcompact(io::IO,t::TwoAxisBanana)=showTwoAxis(io,t,true)
sameLattice(a::TwoAxisBanana,b::Union{Lattice,LatticeVector,Crystal,Sample})=sameLattice(a.sample,b)
sameLattice(b::Union{Lattice,LatticeVector,Crystal,Sample},a::TwoAxisBanana)=sameLattice(b,a.sample)

function sameInstrument(a::TwoAxisBanana,b::TwoAxisBanana)
#    fields=[:energy,:source,:mono,:monitor,:sample,:detector,:arms,:horcol,:vercol]
    fields=[:source,:mono,:monitor,:sample,:detector,:arms,:horcol,:vercol] # FIXME add back in the energy check, and change it to be a Measured
    all([testfield(==,fld,a,b) for fld in fields])
end
(==)(a::TwoAxisBanana,b::TwoAxisBanana)=a.a1==b.a1&&a.a2==b.a2&&a.a3==b.a3&&all(a.a4.==b.a4)&sameInstrument(a,b)

getki(t::TwoAxisBanana,a2::Real=t.a2)=calck(t.mono.tau,a2)
getkf(t::TwoAxisBanana,a2::Real=t.a2)=calck(t.mono.tau,a2)
# return the tuple (ei,ef) from specified E; relies upon fixed ki|kf information in t
calcEi(t::TwoAxisBanana)=t.energy
calcEf(t::TwoAxisBanana)=t.energy
# The angle between ki and Q is given by tan(Φ)=[-kf sin(2θsample)]/[ki - kf cos(2θsample)]
getphi(t::TwoAxisBanana)=getphi.(getki(t),getki(t),t.a4) # call the version defined in TripleAxis.jl
function getphi{T<:TwoAxisBanana}(ta::Array{T})
    phi1=getphi(ta[1])
    allphi=Array{typeof(phi1)}(size(ta)...)
    allphi[1]=phi1
    for i=2:length(ta); allphi[i]=getphi(ta[i]); end
    return allphi
end

iskiFixed(t::TwoAxisBanana)=true
iskfFixed(t::TwoAxisBanana)=false
#
#const ħ²2mₙ=2.072142 # ħ²/(2mₙ) = 2.072142 meV Å² (ħ,² and ₙ used to avoid naming conflicts)
#
# return Q from specified (a2,a3,a4) angles [or current values stored in the TripleAxis object t]
function getQ(t::TwoAxisBanana,a2::Real=t.a2,a3::Real=t.a3,a4::Array{T}=t.a4) where {T<:Real}
    k=getki(t,a2)
    modQ=k*sqrt.(2-2*cos.(a4))
    Ψ=getphi.(k,k,a4)-a3 # angle from first orienting vector to Q (phi is the angle from ki to Q)
#    Ψ-=pi/2 ## FIXME is this 90 degree shift right?
    lab2sample.([[cos(x),sin(x),0] for x in Ψ].*modQ,t.sample)
end

# return E from specified (a2,a6) angles [or current values stored in the TripleAxis object t]
getEi(t::TwoAxisBanana,a2::Real=t.a2)=ħ²2mₙ*getki(t,a2)^2
getEf(t::TwoAxisBanana,a2::Real=t.a2)=ħ²2mₙ*getki(t,a2)^2
getE(t::TwoAxisBanana,a2::Real=t.a2)=ħ²2mₙ*(getki(t,a2)^2-getki(t,a2)^2) # always zero
# Array versions of getQ,getEi,getEf,getE
for z in (:getEi,:getEf,:getE,:getQ)
    @eval begin
        function $z{T<:TwoAxisBanana}(v::Array{T})
            ot1=$z(v[1])
            out=Array{typeof(ot1)}(size(v)...)
            out[1]=ot1
            for j=2:length(v); out[j]=$z(v[j]); end
            return out
        end
    end
end

# return the tuple (Q,E) from specified (a2,a3,a4) angles [or current values stored in the TwoAxisBanana object t]
getQE(t::TwoAxisBanana,a2::Real=t.a2,a3::Real=t.a3,a4::Array{T}=t.a4) where {T<:Real}=(getQ(t,a2,a3,a4),getE(t,a2))

getangle(t::TwoAxisBanana,n::Integer)= 1==n ? t.a1 : 2==n ? t.a2 : 3==n ? t.a3 : 4==n ? t.a4 : throw(BoundsError(t,n))
function getangle(m::Array{T},n::Integer) where {T<:TwoAxisBanana}
    out=Array{eltype(m[1].a)}(size(m)...)
    for i=1:length(m); out[i]=m[i].a[n]; end
    return out
end

getangles(t::TwoAxisBanana)=(t.a1,t.a2,t.a3,t.a4...)
function getangles{T<:TwoAxisBanana}(v::Vector{T})
    out=Array(eltype(v[1].a4),length(v),3+length(v[1].a4))
    for i=1:length(v)
        out[i,:]=[v[i].a1;v[i].a2;v[i].a3;v[i].a4]
    end
    return out
end
calcAngles(t::TwoAxisBanana)=getangles(t)
function calcAngles(t::TwoAxisBanana,Q::Lattice3Vector,E::Real=0.)
    # calcEs returns (Ei,Ef), collect converts it to [Ei,Ef],
    # ks=sqrt([Ei,Ef]/(ħ²/(2mₙ))), (ks...) splits the vector into a tuple for asignment to "ki,kf"
    Ei=calcEi(t)
    ki=sqrt(Ei/ħ²2mₙ)
    θm,tθm=getTh2Th(t.mono,ki)
    θs,tθs=getTh2Th(t.sample,Q,ki,ki)
    (θm,tθm,θs,tθs)
end # this seems to work

#### in-place gotoQE!
function gotoQE!(t::TwoAxisBanana,Q::Lattice3Vector,E::Real=0)
    reta=collect(calcAngles(t,Q,E)) # collect converts the returned tuple into a vector
    t.a1,t.a2,t.a3,ta4=reta
    t.a4=(ta4-t.a4[1])+t.a4 # assume that we're specifying the a4 value of the first detector
end
gotoQE!(t::TwoAxisBanana,QE::Lattice4Vector)=gotoQE!(t,get3vector(QE),getE(QE))
#### in-place gotoQ! (E=0)
gotoQ!(t::TwoAxisBanana,Q::Lattice3Vector)=gotoQE!(t,Q,0.)

setangles!{T<:Real}(t::TwoAxisBanana,angles::Array{T,1})=(@assert length(angles)>=4; t.a1=angles[1];t.a2=angles[2];t.a3=angles[3];t.a4=angles[4:end])
setangles(t::TwoAxisBanana,o...)=(a=deepcopy(t); setangles!(a,o...); a)
function setangles{T<:Real}(t::TwoAxisBanana,angles::Array{T,2})
    @assert size(angles,1)>=4
    out=[setangles(t,angles[:,1])]
    for i=2:size(angles,2); push!(out,setangles(t,angles[:,i])); end
    out::Array{typeof(t),1}
end
function setangles{T<:Real}(t::TwoAxisBanana,angles::Array{T,3})
    @assert size(angles,1)>=4
    out=cat(2,setangles(t,angles[:,:,1])) # make sure out is a Array{TripleAxis,2}
    for i=2:size(angles,3); out=cat(2,out,setangles(t,angles[:,:,i])); end
    out::Array{typeof(t),2}
end

getmonoHR(t::TwoAxisBanana)=getHR(t.mono,t.arms[1],t.arms[2],t.a[1])
getmonoVR(t::TwoAxisBanana)=getVR(t.mono,t.arms[1],t.arms[2],t.a[1])

# TripleAxis object operator overloading: (probably) only used for interpolation between points in a scan
# Overload +,-,/,* for modifications to the current (Q,E) vector
+(a::TwoAxisBanana,b::TwoAxisBanana)=(@assert sameInstrument(a,b); gotoQE(a,get4vector(a)+get4vector(b)))
-(a::TwoAxisBanana,b::TwoAxisBanana)=(@assert sameInstrument(a,b); gotoQE(a,get4vector(a)-get4vector(b)))
+(a::TwoAxisBanana,b::LatticeVector)=(@assert sameLattice(a,b); gotoQE(a,get4vector(a)+b))
+(a::LatticeVector,b::TwoAxisBanana)=(@assert sameLattice(a,b); gotoQE(b,a+get4vector(b)))
-(a::TwoAxisBanana,b::LatticeVector)=(@assert sameLattice(a,b); gotoQE(a,get4vector(a)-b))
-(a::LatticeVector,b::TwoAxisBanana)=(@assert sameLattice(a,b); gotoQE(b,a-get4vector(b)))

# Scatterd {+,-,*,/} Number
for y in [:*,:/]
    @eval Base.$y(a::TwoAxisBanana,b::Number)=gotoQE(a,Base.$y(get4vector(a),b))
    #@eval Base.$y{T<:TwoAxisBanana}(a::Array{T},b::Number)=(out=Array{T}(size(a)); for i=1:length(a); out[i]=Base.$y(a[i],b); end; out)
end
# commutative operator
(*)(b::Number,a::TwoAxisBanana)=(*)(a,b)

# overload mean to divide first instead of last to ensure all Q positions remain accessible (TripleAxis checks this!)
Base.mean{T<:TwoAxisBanana}(a::T;k...)=Base.copy(a)
Base.mean{T<:TwoAxisBanana}(v::Vector{T};k...)=Base.sum(v/Base.length(v);k...)
Base.mean{T<:TwoAxisBanana}(m::Array{T},d;k...)=Base.sum(m/Base.size(m,d);k...)


# Now it should be possible to write a function that inserts the convoluted function between measurement points in order to produce a smooth convoluted function line for plotting
# The simple case will only work for scans that are along defined directions in (Q,E) space.
# Interpolating rocking scans, which generally take a curved path through (Q,E) space, might give wacky results if the step size is large compared to \delta(Q,E)
# To help with such a special case, Abuse *,/ to overload +,- but acting on the angles for each TripleAxis
*(a::TwoAxisBanana,b::TwoAxisBanana)=(@assert sameInstrument(a,b); c=deepcopy(a); c.a=a.a+b.a; c)
/(a::TwoAxisBanana,b::TwoAxisBanana)=(@assert sameInstrument(a,b); c=deepcopy(a); c.a=a.a-b.a; c)


function combineRepeatedInstruments{T<:TwoAxisBanana,I<:Integer}(v::AbstractArray{T},bno::AbstractArray{I},nb::Integer=maximum(bno);k...)
    @assert size(v)==size(bno)
    [mean(v[bno.==i];k...) for i in 1:nb]
end

function _m2σ!{T<:Number}(t::TwoAxisBanana,c::Array{T}) # convert guide m values to angular standard deviations (which are wavelength dependent)
    m=abs(c)
    #Θ₀=pi/180/10 # the critical angle of natural Ni (in radian per angstrom)
    Θᶜ=(pi/1800)*m.*(2pi./[getki(t),getki(t)]) # the critical angle is mΘ₀λ
    σ=Θᶜ/sqrt(3) # and the angular standard deviation of a top-hat distribution between +/- Θᶜ is Θᶜ/sqrt(3)
    c[c.<0]=σ[c.<0]
end

function geometricαs(t::TwoAxisBanana)
    w=[getwidth(tasp.mono)+getwidth(tasp.source),getwidth(tasp.sample)+getwidth(tasp.mono),getwidth(tasp.sample)+getwidth(tasp.detector)]/2
    atan2(w,tasp.arms[1:3])
end
function geometricβs(t::TwoAxisBanana)
    w=[getheight(tasp.mono)+getheight(tasp.source),getheight(tasp.sample)+getheight(tasp.mono),getheight(tasp.sample)+getheight(tasp.detector)]/2
    atan2(w,tasp.arms[1:3])
end

# Eventually each banana-detector Scatterd (WANDd, DMCd, etc.) object will hold an Array of TwoAxisBanana objects
# It will be advantagous to access at least some properties of the TwoAxisBanana objects directly from the Array
getx(a::TwoAxisBanana)=getx(a.sample);
getx{T<:TwoAxisBanana}(v::Array{T,1})=(rx=[getx(v[1])];  for i=2:length(v);   push!(rx,getx(v[i]));   end; rx)
getx{T<:TwoAxisBanana}(m::Array{T,2})=(rx= getx(m[:,1]); for i=2:size(m,2); rx=hcat(rx,getx(m[:,i])); end; rx)
gety(a::TwoAxisBanana)=gety(a.sample);
gety{T<:TwoAxisBanana}(v::Array{T,1})=(ry=[gety(v[1])];  for i=2:length(v);   push!(ry,gety(v[i]));   end; ry)
gety{T<:TwoAxisBanana}(m::Array{T,2})=(ry= gety(m[:,1]); for i=2:size(m,2); ry=hcat(ry,gety(m[:,i])); end; ry)
getz(a::TwoAxisBanana)=getz(a.sample);
getz{T<:TwoAxisBanana}(v::Array{T,1})=(rz=[getz(v[1])];  for i=2:length(v);   push!(rz,getz(v[i]));   end; rz)
getz{T<:TwoAxisBanana}(m::Array{T,2})=(rz= getz(m[:,1]); for i=2:size(m,2); rz=hcat(rz,getz(m[:,i])); end; rz)
get3vector{T<:TwoAxisBanana}(t::Union{T,Array{T}})=getQ(t) # get3vector and getQ are identical.
get4vector(t::TwoAxisBanana)=((Q,E)=getQE(t); Lattice4Vector(Q,E))
get4vector{T<:TwoAxisBanana}(v::Array{T,1})=(rQ=[get4vector(v[1])]; for i=2:length(v);   push!(rQ,get4vector(v[i]));  end; rQ)
get4vector{T<:TwoAxisBanana}(v::Array{T,2})=(rQ= get4vector(v[:,1]);for i=2:size(v,2); rQ=hcat(rQ,get4vector(v[:,i])); end; rQ)

geth{T<:TwoAxisBanana}(t::Union{T,Array{T}})=geth(getQ(t))
getk{T<:TwoAxisBanana}(t::Union{T,Array{T}})=getk(getQ(t))
getl{T<:TwoAxisBanana}(t::Union{T,Array{T}})=getl(getQ(t))

# a convenience function
lab2sample{T<:Real}(Q::Array{T},t::TwoAxisBanana,o...)=lab2sample(Q,t.sample,o...)
lab2sample!{T<:Real}(Q::Array{T},t::TwoAxisBanana,o...)=lab2sample!(Q,t.sample,o...)
# another convenience function (likely not used much)
sample2lab(t::TwoAxisBanana)=sample2lab(getQ(t),t.sample)

# to help in the later construction of matching-size arrays in fuctions:
size(::TwoAxisBanana)=() # size([any scalar])=()
ndims(::TwoAxisBanana)=0 # ndims([any scaler])=0
