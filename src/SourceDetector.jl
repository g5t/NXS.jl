"""
`Source` is an abstract type to describe the source (or incident beam profile) of a
triple-axis scattering instrument. Three subtypes, `CircSource`, `ElliSource`, and `RectSource`
define circular, elliptical, and rectangular geometries.
"""
abstract type Source end
showSource(io::IO,s::Source,compact::Bool=false)=(compact?Base.print(io,"Source $(getwidth(s))w × $(getheight(s))h"):Base.print(io,"$(getwidth(s))×$(getheight(s))"))
Base.show(io::IO,s::Source)=showSource(io,s,false)
Base.showcompact(io::IO,s::Source)=showSource(io,s,true)
"""
    CircSource(r::T)
Define a circular `Source` object from its radius `r`.
If `T` is not a `SecondMoment` the radius is replaced by <x²>=r²/4.
"""
type CircSource <: Source
    r::SecondMoment # <x^2>=r, <y^2>=r
    CircSource{S<:Real}(r::SecondMoment{S})=new(r)
end
CircSource(r::Real)=CircSource(SecondMoment(r^2/4)) # r is the circle radius
"""
    ElliSource(w::T,h::T)
Define an eliptical `Source` object from `w` and `h` horizontal and vertical radii.
If `T` is not a `SecondMoment` the radii are replaced by <xᵢ²>=rᵢ²/4.
"""
type ElliSource <: Source
    w::SecondMoment # <x^2>=w
    h::SecondMoment # <y^2>=h
    ElliSource(w::SecondMoment,h::SecondMoment)=new(w,h)
end
ElliSource(a::Real,b::Real)=(t=typeof(promote(a^2,b^2)[1]);ElliSource(SecondMoment(t(a^2/4)),SecondMoment(t(b^2/4)))) # (a,b) are radii of the ellipse
"""
    RectSource(w::T,h::T)
Define a rectangular `Source` object from its total width `w` and height `h`.
If `T` is not a `SecondMoment` the dimensions are replaced by <xᵢ²>=dᵢ²/12.
"""
type RectSource <: Source
    w::SecondMoment # <x^2>=w
    h::SecondMoment # <y^2>=h
    RectSource(w::SecondMoment,h::SecondMoment)=new(w,h)
end
RectSource(w::Real,h::Real)=(t=typeof(promote(w^2,h^2)[1]);RectSource(SecondMoment(t(w^2/12)),SecondMoment(t(h^2/12)))) # (w,h) are total width and height
(==)(a::CircSource,b::CircSource)=(a.r==b.r)
(==)(a::Source,b::Source)=(a.w==b.w)&(a.h==b.h)

Base.copy(s::CircSource)=CircSource(Base.copy(s.r))
Base.copy(s::Source)=typeof(s)(Base.copy(s.w),Base.copy(s.h))
    

"""
`Detector` is an abstract type to describe the single-detector of a triple-axis scattering 
instrument. Three subtypes, `CircDetector`, `ElliDetecor`, and `RectDetector`
define circular, elliptical, and rectangular cross-section single-detectors.
"""
abstract type Detector end
showDetector(io::IO,s::Detector,compact::Bool=false)=(compact?Base.print(io,"Detector $(getwidth(s))w × $(getheight(s))h"):Base.print(io,"$(getwidth(s))×$(getheight(s))"))
Base.show(io::IO,s::Detector)=showDetector(io,s,false)
Base.showcompact(io::IO,s::Detector)=showDetector(io,s,true)
# separate Detector objects in case it's beneficial at some point to implement (monitor) efficiency corrections
"""
    CircDetector(r::T)
Define a circular `Detector` object from its radius `r`.
If `T` is not a `SecondMoment` the radius is replaced by <x²>=r²/4.
"""
type CircDetector <: Detector
    r::SecondMoment # <x^2>=r, <y^2>=r
    CircDetector{S<:Real}(r::SecondMoment{S})=new(r)
end
CircDetector(r::Real)=CircDetector(SecondMoment(r^2/4)) # r is the circle radius
"""
    ElliDetector(w::T,h::T)
Define an eliptical `Detector` object from `w` and `h` horizontal and vertical radii.
If `T` is not a `SecondMoment` the radii are replaced by <xᵢ²>=rᵢ²/4.
"""
type ElliDetector <: Detector
    w::SecondMoment # <x^2>=w
    h::SecondMoment # <y^2>=h
    ElliDetector(w::SecondMoment,h::SecondMoment)=new(w,h)
end
ElliDetector(a::Real,b::Real)=(t=typeof(promote(a^2,b^2)[1]);ElliDetector(SecondMoment(t(a^2/4)),SecondMoment(t(b^2/4)))) # (a,b) are radii of the ellipse
"""
    RectDetector(w::T,h::T)
Define a rectangular `Detector` object from its total width `w` and height `h`.
If `T` is not a `SecondMoment` the dimensions are replaced by <xᵢ²>=dᵢ²/12.
"""
type RectDetector <: Detector
    w::SecondMoment # <x^2>=w
    h::SecondMoment # <y^2>=h
    RectDetector(w::SecondMoment,h::SecondMoment)=new(w,h)
end
RectDetector(w::Real,h::Real)=(t=typeof(promote(w^2,h^2)[1]);RectDetector(SecondMoment(t(w^2/12)),SecondMoment(t(h^2/12)))) # (w,h) are total width and height
(==)(a::CircDetector,b::CircDetector)=(a.r==b.r)
(==)(a::Detector,b::Detector)=(a.w==b.w)&(a.h==b.h)

Base.copy(s::CircDetector)=CircDetector(Base.copy(s.r))
Base.copy(s::Detector)=typeof(s)(Base.copy(s.w),Base.copy(s.h))

getradius{T<:Union{CircSource,CircDetector}}(a::T)=sqrt(4*a.r)
getwidth{T<:Union{CircSource,CircDetector}}(a::T)=2*getradius(a)
getheight{T<:Union{CircSource,CircDetector}}(a::T)=2*getradius(a)
getwidth{T<:Union{ElliSource,ElliDetector}}(a::T)=2*sqrt(4*a.w)
getheight{T<:Union{ElliSource,ElliDetector}}(a::T)=2*sqrt(4*a.h)
getwidth{T<:Union{RectSource,RectDetector}}(a::T)=sqrt(12*a.w)
getheight{T<:Union{RectSource,RectDetector}}(a::T)=sqrt(12*a.h)

setradius{T<:Union{CircSource,CircDetector}}(a::T,r::Real)=T(r)
setwidth{T<:Union{CircSource,CircDetector}}(a::T,r::Real)=T(r/2)
setheight{T<:Union{CircSource,CircDetector}}(a::T,r::Real)=T(r/2)
setwidth{T<:Union{ElliSource,ElliDetector,RectSource,RectDetector}}(a::T,w::Real)=T(w/2,getheight(a)/2)
setheight{T<:Union{ElliSource,ElliDetector,RectSource,RectDetector}}(a::T,h::Real)=T(getwidth(a)/2,h/2)
setradius{T<:Union{CircSource,CircDetector}}(a::T,r::SecondMoment)=T(r)
setwidth{T<:Union{CircSource,CircDetector}}(a::T,r::SecondMoment)=T(r)
setheight{T<:Union{CircSource,CircDetector}}(a::T,r::SecondMoment)=T(r)
setwidth{T<:Union{ElliSource,ElliDetector,RectSource,RectDetector}}(a::T,w::SecondMoment)=T(w,a.h)
setheight{T<:Union{ElliSource,ElliDetector,RectSource,RectDetector}}(a::T,h::SecondMoment)=T(a.w,h)

# Default to rectangular sources, monitors, and detectors:
Source(w::Real=30.,h::Real=120.)=RectSource(w,h)
Detector(w::Real=30.,h::Real=160.)=RectDetector(w,h)

"""
    getshape(a::Union{Source,Detector})
returns a 2×2 diagonal matrix of the shape of `a` described in terms of its `SecondMoment` values.
"""
getshape(d::Union{CircSource,CircDetector})=diagm([d.r,d.r])
getshape(d::Union{Source,Detector})=diagm([d.w,d.h])  # uses second moment information
