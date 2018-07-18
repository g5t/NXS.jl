function _gaussian_variance(f,x,vx,A,vA,p,vp,σ,vσ)
  #                    FIXME +------- this 0 is due to binned x values where the bin width is *not* related to the variance!
  #                    FIXME V FIXME
  f.^2.*(((p-x)./σ.^2).^2.*vx*0 + vA./A.^2 + (x-p).^2.*vp./σ.^4 + ((x-p).^2./σ.^3-1./σ).^2.*vσ)
end
function _gaussian_nc_variance(f,x,vx,σ,vσ)
  f.^2.*((x./σ.^2).^2.*vx*0 + (x.^2./σ.^3-1./σ).^2.*vσ)
end
"""
    gaussian(area,position,fwhm,x)
A distribution of one variable, `x`, with `area`, central `position`,and full width at half maximum (`fwhm`).

    gaussian([area₁;position₁;fwhm₁;area₂;position₂;area₂;…],x)
The alternative syntax takes a vector `p` with length `3N` (or `3N+1`) which describes `N` gaussian
distributions that should be added.
If `length(p)%3≡1` the last element of the vector is taken to be a constant background.
"""
function gaussian{T<:Number}(area::MeasuredSymmetric,position::MeasuredSymmetric,fwhm::MeasuredSymmetric,x::AbstractVector{T})
    σ=abs(fwhm)/sqrt(8log(2)) # standard deviation
    A=value(area); vA=variance(area); p=value(position); vp=variance(position); w=value(σ); vw=variance(σ);
    f=gaussian(A,p,value(fwhm),value.(x)); vf=_gaussian_variance(f,value.(x),variance.(x),A,vA,p,vp,w,vw)
    MeasuredSymmetric(f,vf)
end
function gaussian{T<:Number}(area::Measured,position::Measured,fwhm::Measured,x::AbstractVector{T})
    σ=abs(fwhm)/sqrt(8log(2)) # standard deviation
    A=value(area);     vpA=positivevariance(area);     vmA=negativevariance(area);
    p=value(position); vpp=positivevariance(position); vmp=negativevariance(position);
    w=value(σ);        vpw=positivevariance(σ);        vmw=negativevariance(σ);
    f=gaussian(A,p,value(fwhm),value.(x));
    vpf=_gaussian_variance(f,value.(x),positivevariance.(x),A,vpA,p,vpp,w,vpw)
    vmf=_gaussian_variance(f,value.(x),negativevariance.(x),A,vmA,p,vmp,w,vmw)
    MeasuredAsymmetric(f,vpf,vmf)
end
function gaussian{R<:Real}(area::Real,position::Real,fwhm::Real,x::AbstractArray{R,1})
    σ=abs(fwhm)/sqrt(8log(2)) # standard deviation
    area/sqrt(2pi)/σ * exp.(-(x-position).^2/2/σ^2)
end
gaussian{T<:Measured}(area::Real,position::Real,fwhm::Real,x::AbstractVector{T})=gaussian(Measured(area,0),Measured(position,0),Measured(fwhm,0),x)
"""
    gaussian(area::Array,positon::Array,fwhm::Array,x::Array)
Calculates the Gaussian distribution where each point in x is allowed to have an independent `area`,
`position`, and `fwhm` (this is useful if the Gaussian distribution should be calculated for
multiple scans simultaneously, as when fitting simultaneously).
"""
function gaussian(area::Array,position::Array,fwhm::Array,x::Array)
    σ=abs.(fwhm)/sqrt(8log(2)) # standard deviation
    area./sqrt(2pi)./σ .* exp.(-(x-position).^2./2./σ.^2)
end
"""
    gaussian_nc(f,x)
A normalized and centered distribution of one variable, `x`, with either a single full-width at
half-maximum, `f`, or individual `fwhm` values for each value of `x`.
"""
gaussian_nc{R<:Real,T<:Real,N}(f::Union{R,Array{R,N}},x::Array{T,N})=sqrt.(4log(2)/pi./f.^2).*exp.(-x.^2.* 4log(2)./f.^2)
function gaussian_nc{R<:MeasuredSymmetric,N}(mw::Union{R,Array{R,N}},mx::Array{R,N})
    w=value.(mw); vw=variance.(mw)
    x=value.(mx); vx=variance.(mx)
    f=gaussian_nc(w,x)
    MeasuredSymmetric(f,_gaussian_nc_variance(f,x,vx,w,vw))
end
function gaussian_nc{R<:Measured,N}(mw::Union{R,Array{R,N}},mx::Array{R,N})
    w=value.(mw); pw=positivevariance.(mw); nw=negativevariance.(mw)
    x=value.(mx); px=positivevariance.(mx); nx=negativevariance.(mx)
    f=gaussian_nc(w,x)
    MeasuredSymmetric(f,_gaussian_nc_variance(f,x,px,w,pw),_gaussian_nc_variance(f,x,nx,w,nw))
end
gaussian_nc{R<:Measured,N,T}(mw::Union{R,Array{R,N}},x::Array{T,N})=gaussian_nc(mw,Measured(x,0))
gaussian_nc{R<:Measured,N,T}(w::Union{T,Array{T,N}},mx::Array{R,N})=gaussian_nc(Measured(w,0),mx)


function gaussian{T<:Number,R<:Number}(p::Vector{T},x::AbstractArray{R,1})
    (length(p)%3==1)?(b=p[end];p=p[1:end-1]):b=zero(eltype(p))
    out=zeros(x)+b
    for i=1:3:length(p)
        out+=gaussian(p[i:i+2]...,x)
    end
    out
end
"""
    absgaussian(area,position,fwhm,x::Array)
A variant of the Gaussian distribution `gaussian` with `abs(area)` used to restrict the distribution
to have a negative second derivative in `x` at `position`
"""
function absgaussian{R<:Real}(area::MeasuredSymmetric,position::MeasuredSymmetric,fwhm::MeasuredSymmetric,x::AbstractArray{R,1})
    σ=abs(fwhm)/sqrt(8log(2)) # standard deviation
    A=value(area); vA=variance(area);  p=value(position); vp=variance(position);  w=value(σ); vw=variance(σ);
    f=absgaussian(A,p,value(fwhm),x); vf=_gaussian_variance(f,x,0,A,vA,p,vp,w,vw)
    Measured(f,vf)
end
function absgaussian{R<:Real}(area::Measured,position::Measured,fwhm::Measured,x::AbstractArray{R,1})
    σ=abs(fwhm)/sqrt(8log(2)) # standard deviation
    A=value(area);     vpA=positivevariance(area);     vmA=negativevariance(area);
    p=value(position); vpp=positivevariance(position); vmp=negativevariance(position);
    w=value(σ);        vpw=positivevariance(σ);        vmw=negativevariance(σ);
    f=absgaussian(A,p,value(fwhm),x);
    vpf=_gaussian_variance(f,x,0,A,vpA,p,vpp,w,vpw)
    vmf=_gaussian_variance(f,x,0,A,vmA,p,vmp,w,vmw)
    MeasuredAsymmetric(f,vpf,vmf)
end
function absgaussian{R<:Real}(area::Real,position::Real,fwhm::Real,x::AbstractArray{R,1})
    σ=abs(fwhm)/sqrt(8log(2)) # standard deviation
    abs(area)/sqrt(2pi)/σ * exp.(-(x-position).^2/2/σ^2)
end
function absgaussian{T<:Number,R<:Real}(p::Vector{T},x::AbstractArray{R,1})
    (length(p)%3==1)?(b=p[end];p=p[1:end-1]):b=zero(eltype(p))
    out=zeros(x)+b
    for i=1:3:length(p)
        out+=absgaussian(p[i:i+2]...,x)
    end
    out
end

function _lorentzian_variance(x,f,A,vA,P,vP,G,vG)
  xmP2=(x-P).^2
  G2p4xmp22=(G.^2+4xmP2).^2
  f.^2.*(vA./A.^2+ 64*xmP2./G2p4xmp22.*vP + (G.^2-4xmP2).^2./G2p4xmp22.*vG./G.^2)
end
_lorentzian_nc_variance(x,f,G,vG)=f.^2.*( (G.^2-4x.^2).^2./(G.^2+4x.^2).^2.*vG./G.^2 )
"""
    lorentzian(area,position,Γ,x)
A Lorentzian (or Cauchy) distribution of one variable normalized to have integrated area `area`,
with its maximum at `position`, and full-width at half-maximum `Γ`, evaluated at the points in `x`

    lorentzian([area₁;position₁;Γ₁;area₂;position₂;Γ₂;…],x)
An alternative syntax evaluates the contributions from `N` Lorentzian distributions at the points
in `x`, where the properties of the `N` distributions are contained in a `3N`-long (or `3N+1` if the
opitonal last element is a constant added to the distributions) vector of parameters/
"""
function lorentzian{R<:MeasuredSymmetric,T<:Real,N}(a::Union{R,Array{R,N}},p::Union{R,Array{R,N}},γ::Union{R,Array{R,N}},x::Array{T,N})
    f=lorentzian(value.(a),value.(p),value.(γ),x)
    vf=_lorentzian_variance(x,f,value.(a),variance.(a),value.(p),variance.(p),value.(γ),variance.(γ))
    MeasuredSymmetric(v,vf)
end
function lorentzian{R<:Measured,T<:Real,N}(a::Union{R,Array{R,N}},p::Union{R,Array{R,N}},γ::Union{R,Array{R,N}},x::Array{T,N})
    f=lorentzian(value.(a),value.(p),value.(γ),x)
    vpf=_lorentzian_variance(x,f,value.(a),positivevariance.(a),value.(p),positivevariance.(p),value.(γ),positivevariance.(γ))
    vmf=_lorentzian_variance(x,f,value.(a),negativevariance.(a),value.(p),negativevariance.(p),value.(γ),negativevariance.(γ))
    MeasuredAsymmetric(v,vpf,vmf)
end
lorentzian{R<:Real,T<:Real,N}(a::Union{R,Array{R,N}},p::Union{R,Array{R,N}},γ::Union{R,Array{R,N}},x::Array{T,N})=2a./pi./γ./(1+(2(x-p)./γ).^2)

lorentzian_nc{R<:Real,T<:Real,N}(γ::Union{R,Array{R,N}},x::Array{T,N})=2./pi./γ./(1+(2x./γ).^2)
function lorentzian_nc{R<:MeasuredSymmetric,N}(mw::Union{R,Array{R,N}},mx::Array{R,N})
    w=value.(mw); vw=variance.(mw)
    x=value.(mx); vx=variance.(mx)
    f=lorentzian_nc(w,x)
    MeasuredSymmetric(f,_lorentzian_nc_variance(x,f,vx,w,vw))
end
function lorentzian_nc{R<:Measured,N}(mw::Union{R,Array{R,N}},mx::Array{R,N})
    w=value.(mw); pw=positivevariance.(mw); nw=negativevariance.(mw)
    x=value.(mx); px=positivevariance.(mx); nx=negativevariance.(mx)
    f=lorentzian_nc(w,x)
    MeasuredSymmetric(f,_lorentzian_nc_variance(x,f,px,w,pw),_lorentzian_nc_variance(x,f,nx,w,nw))
end
lorentzian_nc{R<:Measured,N,T}(mw::Union{R,Array{R,N}},x::Array{T,N})=lorentzian_nc(mw,Measured(x,0))
lorentzian_nc{R<:Measured,N,T}(w::Union{T,Array{T,N}},mx::Array{R,N})=lorentzian_nc(Measured(w,0),mx)

function lorentzian{T<:Number,R<:Real}(p::Vector{T},x::AbstractArray{R,1})
    (length(p)%3==1)?(b=p[end];p=p[1:end-1]):b=zero(eltype(p))
    out=zeros(x)+b
    for i=1:3:length(p)
        out+=lorentzian(p[i:i+2]...,x)
    end
    out
end





"""
A true Voigt function is the convolution of a Gaussian and a Lorentzian.
**See `gaussian` and `lorentzian`.**
This function returns an approximation to a true Voigtian via
[Humlićek's method](http://dx.doi.org/10.1016/0022-4073(82)90078-4)
with [Schreier's modifications](http://dx.doi.org/10.1016/0022-4073(92)90139-U)

    voigt(area,position,fwhm,Γ,x)
For scalars `area`, `position`, Gaussian `fwhm`, and Lorentzian `Γ` (fwhm), calculate an approximation
to the Voigtian at the points in x.
If all five inputs are arrays of the same size, treat the `i`th element of each as an independent
specification for the distribution to be evaluated: `voigt(area,position,fwhm,Γ,x)[i]==voigt(area[i],position[i],fwhm[i],Γ[i],x[i])`

"""
voigt{T<:Real}(gfwhm::Number,lfwhm::Number,x::AbstractArray{T,1})=voigt(one(T),zero(T),gfwhm,lfwhm,x)
function voigt{T<:Real}(area::Number,position::Number,gfwhm::Number,lfwhm::Number,x::AbstractArray{T,1})
    if lfwhm==0
        y=gaussian(area,position,gfwhm,x)
    elseif gfwhm==0
        y=lorentzian(area,position,lfwhm,x)
    else
        A=sqrt(log(2)/pi)*area*ones(size(x))
        γg=gfwhm/2*ones(size(x))
        hx=sqrt(log(2))*(x-position*ones(size(x)))./γg
        hy=sqrt(log(2))*ones(size(x))*abs(lfwhm/gfwhm)
        y=A./γg.*humlicek(hx,hy)
    end
    return y
end
voigt(g::Vector,l::Vector,x::Vector)=voigt(ones(x),zeros(x),g,l,x)
voigt(g::Matrix,l::Matrix,x::Matrix)=voigt(ones(x),zeros(x),g,l,x)
voigt(g::Matrix,l::Matrix,x::Vector)=voigt(ones(g),zeros(g),g,l,repeat(x,outer=[1,size(g,2)]))
voigt(g::Matrix,l::Vector,x::Matrix)=voigt(ones(x),zeros(x),g,repeat(l,outer=[1,size(x,2)]),x)
voigt(g::Vector,l::Matrix,x::Matrix)=voigt(ones(x),zeros(x),repeat(g,outer=[1,size(x,2)]),l,x)

function voigt(area::Array,position::Array,gfwhm::Array,lfwhm::Array,x::Array)
    y=zeros(x)
    voigt!(area,position,gfwhm,lfwhm,x,y)
    return y
end
function voigt!(area::Array,position::Array,gfwhm::Array,lfwhm::Array,x::Array,y::Array)
    @assert compatible(area,position,gfwhm,lfwhm,x,y)
    unsafevoigt!(area,position,gfwhm,lfwhm,x,y)
end
function unsafevoigt!(area::Array,position::Array,gfwhm::Array,lfwhm::Array,x::Array,y::Array)
    g=(lfwhm.==0); any(g)&&(y[g]=gaussian(area[g],position[g],gfwhm[g],x[g]))
    l=(gfwhm.==0); any(l)&&(y[l]=lorentzian(area[l],position[l],lfwhm[l],x[l]))
    v=.!(g.|l);    any(v)&&(y[v]=sqrt(log(2)/pi)./(gfwhm[v]/2).*area[v].*humlicek(sqrt(log(2))*(x[v]-position[v])./(gfwhm[v]/2),sqrt(log(2))*abs.(lfwhm[v]./gfwhm[v])))
    return y
end

function unsafevoigt_nc!{T<:MeasuredSymmetric}(𝒢::Array{T},ℒ::Array{T},x::Array{T},y::Array{T})
    g= ℒ.==0; l= 𝒢.==0; v=.!(g.|l)
    any(g) && (y[g]=gaussian_nc(𝒢[g],x[g]))
    any(l) && (y[l]=lorentzian_nc(ℒ[l],x[l]))
    any(g)&&status(:info,"some points have zero lorentzian width")
    any(l)&&status(:info,"some points have zero gaussian width")
    if any(v)
        val𝒢=value.(𝒢[v]); var𝒢=variance.(𝒢[v])
        valℒ=value.(ℒ[v]); varℒ=variance.(ℒ[v])
        valx=value.(x[v]); varx=variance.(x[v])
        z=sqrt(log(2))*(2valx./val𝒢  + im*abs.(valℒ./val𝒢))
        w=faddeeva(z)
        ∂V∂x=_∂humlicek∂x(z.*w)
        ∂V∂y=_∂humlicek∂y(z.*w)
        valV=2sqrt(log(2)/pi)./val𝒢.*real.(w)
        varV=4log(2)/pi*(4log(2)*∂V∂x.^2.*varx + log(2)*∂V∂y.^2.*varℒ + abs.(2sqrt(log(2)).*(valx.*∂V∂x+0.5*valℒ.*∂V∂y)./val𝒢+real.(w)).^2.*var𝒢)./val𝒢.^4
        y[v]=MeasuredSymmetric(valV,varV)
    end
    return y
end
function unsafevoigt_nc!{T<:Measured}(𝒢::Array{T},ℒ::Array{T},x::Array{T},y::Array{T})
    g= ℒ.==0; l= 𝒢.==0; v=.!(g.|l)
    any(g) && (y[g]=gaussian_nc(𝒢[g],x[g]))
    any(l) && (y[l]=lorentzian_nc(ℒ[l],x[l]))
    if any(v)
        val𝒢=value.(𝒢[v]); vap𝒢=positivevariance.(𝒢[v]); vam𝒢=negativevariance.(𝒢[v])
        valℒ=value.(ℒ[v]); vapℒ=positivevariance.(ℒ[v]); vamℒ=negativevariance.(ℒ[v])
        valx=value.(x[v]); vapx=positivevariance.(x[v]); vamx=negativevariance.(x[v])
        z=sqrt(log(2))*(2valx./val𝒢  + im*abs.(valℒ./val𝒢))
        w=faddeeva(z)
        ∂V∂x=_∂humlicek∂x(z.*w)
        ∂V∂y=_∂humlicek∂y(z.*w)
        valV=2sqrt(log(2)/pi)./val𝒢.*real.(w)
        vapV=4log(2)/pi*(4log(2)*∂V∂x.^2.*vapx + log(2)*∂V∂y.^2.*vapℒ + abs.(2sqrt(log(2)).*(valx.*∂V∂x+0.5*valℒ.*∂V∂y)./val𝒢+real.(w)).^2.*vap𝒢)./val𝒢.^4
        vamV=4log(2)/pi*(4log(2)*∂V∂x.^2.*vamx + log(2)*∂V∂y.^2.*vamℒ + abs.(2sqrt(log(2)).*(valx.*∂V∂x+0.5*valℒ.*∂V∂y)./val𝒢+real.(w)).^2.*vam𝒢)./val𝒢.^4
        y[v]=MeasuredAsymmetric(valV,vapV,vamV)
    end
    return y
end
unsafevoigt_nc!{T<:Measured}(𝒢::Array{T},ℒ::Array{T},x::Array,y::Array{T})=unsafevoigt_nc!(𝒢,ℒ,Measured(x,0),y)
unsafevoigt_nc!{T<:Measured}(𝒢::Array{T},ℒ::Array,x::Array{T},y::Array{T})=unsafevoigt_nc!(𝒢,Measured(ℒ,0),x,y)
unsafevoigt_nc!{T<:Measured}(𝒢::Array,ℒ::Array{T},x::Array{T},y::Array{T})=unsafevoigt_nc!(Measured(𝒢,0),ℒ,x,y)
unsafevoigt_nc!{T<:Measured}(𝒢::Array,ℒ::Array,x::Array{T},y::Array{T})=unsafevoigt_nc!(Measured(𝒢,0),Measured(ℒ,0),x,y)
unsafevoigt_nc!{T<:Measured}(𝒢::Array,ℒ::Array{T},x::Array,y::Array{T})=unsafevoigt_nc!(Measured(𝒢,0),ℒ,Measured(x,0),y)
unsafevoigt_nc!{T<:Measured}(𝒢::Array{T},ℒ::Array,x::Array,y::Array{T})=unsafevoigt_nc!(𝒢,Measured(ℒ,0),Measured(x,0),y)

function unsafevoigt_nc!(gfwhm::Array,lfwhm::Array,x::Array,y::Array) # normalized and centered
    g=(lfwhm.==0); any(g)&&(y[g]=gaussian_nc(gfwhm[g],x[g]))
    l=(gfwhm.==0); any(l)&&(y[l]=lorentzian_nc(lfwhm[l],x[l]))
    v=.!(g.|l);    any(v)&&(y[v]=sqrt(log(2)/pi)./(gfwhm[v]/2).*humlicek( sqrt(log(2))*(x[v])./(gfwhm[v]/2) , sqrt(log(2))*abs.(lfwhm[v]./gfwhm[v]) ) )
    return y
end
function voigt{T<:Number,R<:Real}(p::Vector{T},x::AbstractArray{R,1})
    (length(p)%4==1)?(b=p[end];p=p[1:end-1]):b=zero(eltype(p))
    out=zeros(x)+b
    for i=1:4:length(p)
        out+=voigt(p[i:i+3]...,x)
    end
    return out
end

# voigt takes full widths but for convoluting it's nicer to work in intrinsic widths:
       voigt_σγ(a,c,σ,γ,x,y)=       voigt(a,c,σ*sqrt(8log(2)),2γ,x)
      voigt_σγ!(a,c,σ,γ,x,y)=      voigt!(a,c,σ*sqrt(8log(2)),2γ,x,y)
unsafevoigt_σγ!(a,c,σ,γ,x,y)=unsafevoigt!(a,c,σ*sqrt(8log(2)),2γ,x,y)
unsafevoigt_σγ_nc!(σ,γ,x,y) =unsafevoigt_nc!(σ*sqrt(8log(2)),2γ,x,y)

"""
    humlicek(x,y)

The function `humlicek` calculates the real part of a Kramp function following Schreier's
method for calculating the convolution of a Gaussian and a Lorentzian lineshape.
It utilizes the C math library's `erfcx` to avoid under/overflow with small/large arguments.

For `z=x+im*y` the built-in Julia function `erfcx(-im*z)` computes the Faddeeva function `w(z)`,
which has the property `w(z)=V(x,y)+im*L(x,y)` where `V(x,y)` is the real Voigt function and
`L(x,y)` is the imaginary Voigt function.

This function is appropriately overloaded to calculate uncertainty if passed one or two `Measured`
values or arrays.

The name of this function reflects that it has supplanted an older version which achieved
the same end goal via Humlićek's approximations, which may have been slightly faster at the
expense of much more complicated source code and possible accuracy problems at large x.

See http://www.sciencedirect.com/science/article/pii/002240739290139U for Schreier's derivation,
and https://en.wikipedia.org/wiki/Faddeeva_function for information about the Faddeeva/Kramp function.
"""
function humlicek{T<:MeasuredSymmetric}(xm::AbstractArray{T},ym::AbstractArray{T})
    x=value.(xm); xv=variance.(xm); y=value.(ym); yv=variance.(ym);
    z=x+im*y
    w=faddeeva(z)
    return MeasuredSymmetric(real(w),_δhumlicek²(z.*w,xv,yv))
end
function humlicek{T<:Measured}(xm::AbstractArray{T},ym::AbstractArray{T})
    x=value.(xm); xpv=positivevariance.(xm); xnv=negativevariance.(xm)
    y=value.(ym); ypv=positivevariance.(ym); ynv=negativevariance.(ym)
    z=x+im*y
    w=faddeeva(z)
    return MeasuredAsymmetric(real(w),_δhumlicek²(z.*w,xpv,ypv),_δhumlicek²(z.*w,xnv,ynv))
end
humlicek{T<:Real}(x::AbstractArray{T},y::AbstractArray{T})=humlicek(x+im*y)
humlicek{T<:Complex}(z::AbstractArray{T})=real.(faddeeva(z))

humlicek{T<:Number}(x::AbstractArray{T},a::Number)=(t=promote_type(eltype(x),typeof(a));humlicek(map(t,x),a*ones(x)))
humlicek{T<:Number,R<:Number}(x::AbstractArray{T},y::AbstractArray{R})=(t=promote_type(T,R);humlicek(map(t,x),map(t,y)))
humlicek(x::Number,a::Number)=humlicek([x],[a])[1]

# """
#     _δhumlicek²(z*w(z),x,δx²,y,δy²)
#
# The partial derivatives of the Faddeeva function, `w(z)=V(x,y)+im*L(x,y)` are known:
#
#     ∂V(x,y)/∂x = -2*Re[z*w(z)]
#
#     ∂V(x,y)/∂y =  2*Im[x*w(z)] - 2/sqrt(pi)
#
# The variance of the real Voigt function is then
#
#     δV(x,y)² = 4*Re[z*w(z)]²*δx² + 4*(Im[z*w(z)]-1/sqrt(pi))²*δy²
# """
_δhumlicek²(zw,xv,yv)=_∂humlicek∂x(zw).^2.*xv+_∂humlicek∂y(zw).^2.*yv
_∂humlicek∂x(zw)=-2*real(zw)
_∂humlicek∂y(zw)=2*(imag.(zw)-1/sqrt(pi))


"""
    faddeeva(z)=w(z)

The Faddeeva distribution, w(z).
"""
faddeeva{T<:Complex}(z::AbstractArray{T})=erfcx.(-im*z)


function simpleincoherentline(p::Array,t)
    sz=size(t)
    ω=zeros(sz)
    S=ones(sz)*(length(p)>0?p[1]:1)
    γ=ones(sz)*(length(p)>1?p[2]:0)
    return (ω,S,γ)
end
#incoherentline=Dispersion(simpleincoherentline)
