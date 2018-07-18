"""
A `Bragg` object is a triple-axis spectrometer element that is *always* in Bragg scattering geometry.

| Fieldname | Description |
|:---------:|:------------|
| tau        | The wavevector magnitude of the reflection used to select a |
|            | wavelength from the `Bragg` object incident beam, τ, in Å⁻¹ |
| identifier | A string used to identify the `Bragg` object in a human-readable|
|            | format. A typical string is `"PG(002)"` indicates the `Bragg` |
|            | object is made of Pyrolytic Graphite and that its (002) |
|            | reflection is being used. |
| sense      | Scattering sense of the `Bragg` object, (only the sign matters)|
| kinref     | Kinematic reflectivity for reflectivity corrections. |
| mosaic     | A `Vector{Real}` of [horizontal, vertical] mosaic FWHM in radians |
| size       | A `Vector{SecondMoment}` of [width, height, depth] |
| radii      | A `Vector{Function}` to calculate the [horizontal, vertical] |
|            | radii of curvature using `getHR(::Bragg)` and `getVR(::Bragg)`|
"""
type Bragg{T<:Real,S<:SecondMoment,F<:Function}
    tau::T
    identifier::AbstractString
    sense::Int8 # scattering sense, +1 or -1
    kinref::T # kinematic reflectivity
    mosaic::Vector{T} # horizontal and vertical mosaic FWHM in radians
    size::Vector{S} # width,height,depth expresses as second moments
    radii::Vector{F} # horizontal and vertical radii of curvature functions
end
function Base.copy(a::Bragg)
    tau       =Base.copy(a.tau)
    identifier=Base.copy(a.identifier)
    sense     =Base.copy(a.sense)
    kinref    =Base.copy(a.kinref)
    mosaic    =Base.copy(a.mosaic)
    size      =Base.copy(a.size)
    radii     =Base.copy(a.radii)
    Bragg(tau,identifier,sense,kinref,mosaic,size,radii)
end
"""
    Bragg(τ,i,sn=+1,m=[30',30'],sz=[1,1,1],r=[flatBragg,flatBragg],kr=0)
Intializes a `Bragg` object given its wavevector magnitude and associated d spacing plus, optionally,
an identifier string, it's scattering sense, mosaic, physical extent, curvature functions,
and kinematic reflectivity.

    Bragg(d)
Calculates `τ` from `d=2π/τ` and initializes a `Bragg` option with all default values.

    Bragg(identifier::AbstractString="PG(002)")
Determines `τ` from the identifier string using `gettau(identifier)`, and then initializes a `Bragg` object.
"""
function Bragg{T<:Real,S<:Real,F<:Function}(τ::Real,i::AbstractString,sn::Real=1,
                               m::Vector{T}=[1,1]*pi/360,
                               sz::Vector{S}=[1,1,1],
                               r::Vector{F}=[flatBragg,flatBragg],kr::Real=0)
    (τ,kr)=promote(τ,kr)  # default to: kinematic reflectivity = 0
    sn=Int8(signbit(sn)?-1:1) #             sense = +1
    m=map(typeof(τ),m)        #             mosaic (horz./vert.) = 0
    #sz=map(typeof(τ),sz)      #             (w x h x d) = (1 x 1 x 1)
    sz=map(SecondMoment,sz.^2/12) # conver to second moments of a rectangular shape
    Bragg(τ,i,sn,kr,m,sz,r) #             (horz./vert.) radius = Inf <- i.e., flat
end
Bragg(d::Real)=Bragg(2*pi/d,"d=$(d)Å")
Bragg(s::AbstractString="PG(002)")=(τ=gettau(s); Bragg(τ,s))

Base.show(io::IO,m::Bragg)=showBragg(io,m,false)
Base.showcompact(io::IO,m::Bragg)=showBragg(io,m,true)

(==)(a::Bragg,b::Bragg)=(a.tau==b.tau)&(a.sense==b.sense)&(a.kinref==b.kinref)&all(a.mosaic.==b.mosaic)&all(a.size.==b.size)&all(a.radii.==b.radii)

flatBragg{T<:Real}(l1::T,l2::T,ang::T)=T(Inf)
verticallyCurvedBragg{T<:Real,N}(l1::AbstractArray{T,N},l2::AbstractArray{T,N},ang::AbstractArray{T,N})= 2 ./(1 ./l1 + 1 ./l2).*sin(ang)
horizontalyCurvedBragg{T<:Real,N}(l1::AbstractArray{T,N},l2::AbstractArray{T,N},ang::AbstractArray{T,N})=2 ./(1 ./l1 + 1 ./l2)./sin(ang)
rowlandCurvedBragg(l1::AbstractArray{T,N},l2::AbstractArray{T,N},ang::AbstractArray{T,N}) where {T<:Real,N} = (l1+l2)/2 ./sin(ang)
verticallyCurvedBragg{T<:Real}(l1::T,l2::T,ang::T)= 2/(1/l1+1/l2)*sin(ang)
horizontalyCurvedBragg{T<:Real}(l1::T,l2::T,ang::T)=2/(1/l1+1/l2)/sin(ang)
rowlandCurvedBragg(l1::T,l2::T,ang::T) where T<:Real = (l1+l2)/2sin(ang)

isCurved(x::Function)=~isinf(x(1.,1.,pi/2))
setHorzC!(b::Bragg,c::Function)=(b.radii[1]=(x,y,z)->1/c(x,y,z))
setVertC!(b::Bragg,c::Function)=(b.radii[2]=(x,y,z)->1/c(x,y,z))
setHorzR!(b::Bragg,c::Function)=(b.radii[1]=c)
setVertR!(b::Bragg,c::Function)=(b.radii[2]=c)
getHorzC(b::Bragg,l1::Real,l2::Real,ang::Real)=(1/b.radii[1](l1,l2,ang)*getsense(b))
getVertC(b::Bragg,l1::Real,l2::Real,ang::Real)=(1/b.radii[2](l1,l2,ang)*getsense(b))
getHR(b::Bragg,l1::Real,l2::Real,ang::Real)=(b.radii[1](l1,l2,ang)*getsense(b))
getVR(b::Bragg,l1::Real,l2::Real,ang::Real)=(b.radii[2](l1,l2,ang)*getsense(b))
function showBragg(io::IO,b::Bragg,compact::Bool)
    str=(signbit(b.sense)?"-":"+")
    if (isCurved(b.radii[1])||isCurved(b.radii[2]))
        if (isCurved(b.radii[1])&&isCurved(b.radii[2]))
            str*=compact?"d":"double "
        elseif isCurved(b.radii[1])
            str*=compact?"h":"horizontal "
        else
            str*=compact?"v":"vertical "
        end
        str*=compact?"f":"focusing "
    else
        str*=compact?"":"flat "
    end
    str*=b.identifier
    compact||(str*=" ")
    print(io,str)
end

function gettau(matref::AbstractString)
    # For a material,reflection specification of the form PG(002) or ge311, return τ=2π/d
    matref=lowercase(replace(matref,r"[(,)]","")) # allow for, e.g., PG(002) to match pg002
    knwMR=["pg002","pg004","ge111","ge220","ge311","be002","pg110","si111","he111"]
    τinvÅ=[1.87325,3.74650,1.92366,3.14131,3.68351,3.50702,5.49806,2.003936,1.8235248299941869]
    matched= matref .== knwMR
    any(matched) ? τinvÅ[matched][1] : error("Unknown material/reflection specification")
end
gettau(τ::Real)=τ
gettau(b::Bragg)=b.tau
resense!(b::Bragg,s::Number)=(b.sense=Int8(signbit(s)?-1:1))
setsense!(b::Bragg,s::Number)=resense!(b,s)
getsense(b::Bragg)=(signbit(b.sense)?-1:1)
# getTh2Th(b::Bragg,k::Real)=(th=asin(gettau(b)/2/k)*getsense(b);(th,2th))

function getTh2Th(b::Bragg,k::Real)
  arg=gettau(b)/2/k
  @assert 0<=arg<=1 "Passed wavevector k=$k Å¹ can not be selected by a $b Bragg device "
  th=asin(arg)
  return (th,2th)
end

function tauname(tau::Real,tol::Real=1e-3)
    knwMR=["PG(002)","PG(004)","Ge(111)","Ge(220)","Ge(311)","Be(002)","PG(110)","Si(111)","He(111)"]
    τinvÅ=[1.87325,3.74650,1.92366,3.14131,3.68351,3.50702,5.49806,2.003936,1.8235248299941869]
    matched= abs.(tau.-τinvÅ).<abs(tol)
    any(matched)? knwMR[matched][1] : "$(2pi/tau) Å"
end


# access widths and heights for Bragg objects
getdepth(b::Bragg)=sqrt(12*b.size[3])
getwidth(b::Bragg)=sqrt(12*b.size[1])
getheight(b::Bragg)=sqrt(12*b.size[2])
getshape(b::Bragg)=diagm(b.size[[3,1,2]]) # depth,width,height as Popovici's method expects, also second moments

setdepth(b::Bragg,d::SecondMoment)=(b.size[3]=d)
setdepth(b::Bragg,d::Real)=setwidth(b, SecondMoment(d^2/12))
setwidth(b::Bragg,w::SecondMoment)=(b.size[1]=w)
setwidth(b::Bragg,w::Real)=setwidth(b, SecondMoment(w^2/12))
setheight(b::Bragg,h::SecondMoment)=(b.size[2]=h)
setheight(b::Bragg,h::Real)=setwidth(b, SecondMoment(h^2/12))
