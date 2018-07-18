#module M
immutable Measured <: Real
    val::Float64
    var::Float64
end
Measured(a::Real=0)=Measured(a,0.)
@vectorize_2arg Real Measured
@vectorize_1arg Real Measured
Measured(m::Measured)=m # no need to worry about copying because Measured values are immutable
@vectorize_1arg Measured Measured

zero(::Type{Measured})=Measured(0)
one(::Type{Measured})=Measured(1)

+(a::Measured,b::Measured)=Measured(a.val+b.val, a.var+b.var)
-(a::Measured,b::Measured)=Measured(a.val-b.val, a.var+b.var)
-(a::Measured)=Measured(-a.val,a.var)
/(x::Measured,y::Measured)=Measured(x.val./y.val, abs(1./y.val).^2.*x.var + abs(x.val./y.val.^2).^2.*y.var)
*(x::Measured,y::Measured)=Measured(x.val.*y.val, abs(y.val).^2.*x.var + abs(x.val).^2.*y.var)
.*(x::Measured,y::Measured)=x*y
./(x::Measured,y::Measured)=x/y

# comparisons, only including variance when comparing two Measured values
<(a::Measured,b::Measured)=isless(a,b)
>(a::Measured,b::Measured)=isless(b,a)
<(a::Measured,b::Real)=isless(a,b)
<(a::Real,b::Measured)=isless(a,b)
>(a::Measured,b::Real)=isless(b,a)
>(a::Real,b::Measured)=isless(b,a)
<=(a::Measured,b::Measured)=isless(a,b)|isapprox(a,b)
>=(a::Measured,b::Measured)=isless(b,a)|isapprox(a,b)
<=(a::Measured,b::Real)=a.val <= b
<=(a::Real,b::Measured)=a <= b.val
>=(a::Measured,b::Real)=a.val >= b
>=(a::Real,b::Measured)=a >= b.val
isless(a::Measured,b::Measured)=isless(a.val+sqrt(a.var),b.val-sqrt(b.var))
isapprox(a::Measured,b::Measured)=!(isless(a,b)|isless(b,a))
isless{T<:AbstractFloat}(a::Measured,b::T)=isless(a.val,b)
isless{T<:AbstractFloat}(a::T,b::Measured)=isless(a,b.val)
isless{T<:Real}(a::Measured,b::T)=isless(a.val,b)
isless{T<:Real}(a::T,b::Measured)=isless(a,b.val)
isapprox{T<:Real}(a::Measured,b::T)=isapprox(a.val,b)
isapprox{T<:Real}(a::T,b::Measured)=isapprox(a,b.val)

^(x::Measured,n::Integer)=Measured(x.val.^n, abs(n*x.val.^(n-1)).^2.*x.var)
^(x::Measured,n::Rational)=Measured(x.val.^n, abs(n*x.val.^(n-1)).^2.*x.var)
^(x::Measured,n::Real)=Measured(x.val.^n, abs(n*x.val.^(n-1)).^2.*x.var)
#sqrt(x::Measured)=x^(0.5)
cos(x::Measured)=Measured(cos(x.val), abs(sin(x.val)).^2.*x.var)
sin(x::Measured)=Measured(sin(x.val), abs(cos(x.val)).^2.*x.var)
tan(x::Measured)=sin(x)/cos(x)
acos(x::Measured)=Measured(acos(x.val),x.var/(1-x.val^2))
asin(x::Measured)=Measured(asin(x.val),x.var/(1-x.val^2))
atan(x::Measured)=Measured(atan(x.val),x.var/(1+x.val^2)^2)
atan2(y::Measured,x::Measured)=Measured(atan2(y.val,x.val), x.val^2*( (y.val/x.val)^2 * x.var + y.var)/(x.val^2+y.val^2)^2 )
sqrt(x::Measured)=Measured(sqrt(x.val),x.var/2/x.val)

log(x::Measured)=Measured(log(x.val), abs(1./x.val).^2 .* x.var)
log10(x::Measured)=log(x)./log(10.)
@vectorize_1arg Measured sqrt
@vectorize_1arg Measured cos
@vectorize_1arg Measured sin
@vectorize_1arg Measured tan
@vectorize_1arg Measured acos
@vectorize_1arg Measured asin
@vectorize_1arg Measured atan
@vectorize_2arg Measured atan2
@vectorize_1arg Measured log
@vectorize_1arg Measured log10

Base.abs(x::Measured)=Measured(Base.abs(x.val),x.var)
Base.isnan(x::Measured)=Base.isnan(x.val)
Base.signbit(a::Measured)=Base.signbit(a.val)

function discernablyUnique{T<:Measured}(x::AbstractArray{T,1})
    xval=val(x)
    xerr=err(x)
    # BitArrays are **not** initialized to zero! Use "falses" instead
    found=falses(size(x)...)
    du=falses(size(x)...)
    while sum(found)<length(x)
        i=findfirst(!found)
        withinerror = abs(xval-xval[i]) .<= xerr[i] # <= to ensure at least the current point is selected (if xerr[i]==0)
        found = found | withinerror
        du[i]=true
    end
    xval[du]
end
        


Base.size(x::Measured)=tuple() # like any single value, the size of a Measured is ()

# promote_type seems to be used by .* and ./ broadcasts
Base.promote_type(::Type{Measured},::Type{Measured})= Measured
Base.promote_type(::Type{Measured},::Type{Union{}})= Measured
Base.promote_type(::Type{Union{}},::Type{Measured})= Measured
Base.promote_type{T<:Real}(::Type{Measured},::Type{T})= Measured
Base.promote_type{T<:Real}(::Type{T},::Type{Measured})= Measured

# promote seems to be used by + and -
Base.promote(a::Measured,b::Measured)=(a,b)
Base.promote(a::Real,b::Measured)=(Measured(a),b)
Base.promote(a::Measured,b::Real)=(a,Measured(b))

# promote_rule appears to be used by .^
Base.promote_rule(::Type{Measured},::Type{Measured})= Measured
Base.promote_rule(::Type{Measured},::Type{Real})= Measured
Base.promote_rule(::Type{Real},::Type{Measured})= Measured

Base.convert(::Type{Bool},x::Measured)=convert(Bool,x.val)
Base.convert(::Type{Integer},x::Measured)=convert(Integer,x.val)
Base.convert(::Type{Rational{Base.GMP.BigInt}},x::Measured)=convert(Rational{GMP.BigInt},x.val)
Base.convert{T<:Integer}(::Type{Rational{T}},x::Measured)=convert(Rational{T},x.val)
Base.convert(::Type{Measured},x::Measured)=x
#Base.convert{S<:Real}(::Type{S},x::Measured)=convert(S,x.val)
Base.convert(::Type{Measured},x::Real)=Measured(convert(Float64,x)) # allows for, e.g., Measured[1,2,3,4]
Base.convert{T<:AbstractFloat}(::Type{T},x::Measured)=convert(T,x.val)


val(m::Measured)=m.val
var(m::Measured)=m.var
err(m::Measured)=sqrt(m.var)
@vectorize_1arg Measured val
@vectorize_1arg Measured var
@vectorize_1arg Measured err


pre(x)=(x<0?-1:1)*floor(Integer,abs(x))
zeropad(a::Real,b::Int)= ( n=a>0?round(Int,floor(log10(a)))+1:1; b-=n; (b>0?"$("0"^b)":"")*"$a" )
post(x,y)=zeropad(round(Integer,abs(x-pre(x))/10.0^y),abs(y))
function showMeasured(io::IO,m::Measured,compact::Bool=true)
    v,u = val(m),err(m)
    if u>0 # some uncertainty
#        if compact
            pv=(v!=0)?round(Integer,floor(log10(abs(v)))):0 # power of first digit in the value
            pu=round(Integer,floor(log10(u))) # power of first digit in the error
            dv=round(Integer,v/10.0^pu) # the value rounded to the precision of 1 digit in the error
            if dv==0 && v!=0 # the value isn't zero, but the rounded value is
                while dv==0; pu=pu-1; dv=round(Integer,v/10.0^pu); end # rounded value is zero but actual value is non-zero, increase the precision
            end
            du=round(Integer,u/10.0^pu) # the error rounded to the same precision
            # now actually do the formatting
            if (pv<0&&pu<0) # we need to try it both ways
                str1="$dv($du)ᴇ$pu"
                str2=((v<0)?"-":"")*"0."*"$(post(v,pu))($du)"
                str=(length(str1)<(length(str2)-3))?str1:str2 # pick str2 preferentially
            elseif pu>0 # two ways to try as well
                tdv=round(v/10.0^pv)
                tv=(v-tdv*10.0^pv)/10.0.^pv
                str1=((v<0)?"-":"")*"$(pre(tdv))."*"$(post(tv,pv-pu))($du)ᴇ$pv"
                str2="$(dv*10^pu)($(du*10^pu))"
                str=(length(str1)<(length(str2)-3))?str1:str2 # pick str2 preferentially
            elseif pu<0
                str="$(pre(v))"*"."*"$(post(v,pu))($du)"
            else # pu==0
                str="$dv($du)"
            end
            print(io,str)
#        else
#            show(io,v)
#            print(io," ± ")
#            showcompact(io,u)
#        end
    else
        compact ? showcompact(io,v) : show(io,v)
    end
end    
Base.show(io::IO,m::Measured)=showMeasured(io,m,false)
Base.showcompact(io::IO,m::Measured)=showMeasured(io,m,true)
function Base.alignment(x::Measured)
    m = match(r"^(.*?)((?:[(±].*)?)$", sprint(Base.showcompact_lim, x))
    m === nothing ? (length(sprint(Base.showcompact_lim, x)), 0) : (length(m.captures[1]), length(m.captures[2]))
end

## Plotting of Measured arrays
#function PyPlot.plot(x::AbstractArray{Measured,1},y::AbstractArray{Measured,1};marker="o",color="r",linestyle=" ",fillstyle="white")
#    (fillstyle == "white" || fillstyle == "w") ? (facecolor="white";fillstyle="full";):(facecolor=color);
#    h=PyPlot.errorbar(val(x),val(y),xerr=err(x),yerr=err(y))
#    PyPlot.setp(h[1],color=facecolor,marker=marker,markeredgecolor=color,markerfacecoloralt=color,fillstyle=fillstyle,linestyle=linestyle)
#    PyPlot.setp(h[2],marker=" ",color=color)
#    PyPlot.setp(h[3],color=color)
#    return h
#end
#PyPlot.plot{T<:Real}(x::AbstractArray{T,1},y::AbstractArray{T,1},o...)=PyPlot.plot(Measured(x),Measured(y),o...)
##plot{T<:Real}(x::AbstractArray{T,1},y::AbstractArray{Measured,1},o...)=plot(Measured(x),y,o...)
##plot{T<:Real}(x::AbstractArray{Measured,1},y::AbstractArray{T,1},o...)=plot(x,Measured(y),o...)

# Plotting of Measured arrays
function PyPlot.errorbar(x::AbstractArray{Measured,1},y::AbstractArray{Measured,1};marker="o",color="r",linestyle=" ",fillstyle="white")
    (fillstyle == "white" || fillstyle == "w") ? (facecolor="white";fillstyle="full";):(facecolor=color);
    h=PyPlot.errorbar(val(x),val(y),xerr=err(x),yerr=err(y))
    PyPlot.setp(h[1],color=facecolor,marker=marker,markeredgecolor=color,markerfacecoloralt=color,fillstyle=fillstyle,linestyle=linestyle)
    PyPlot.setp(h[2],marker=" ",color=color)
    PyPlot.setp(h[3],color=color)
    return h
end
#PyPlot.plot{T<:Real}(x::AbstractArray{T,1},y::AbstractArray{T,1},o...)=PyPlot.plot(Measured(x),Measured(y),o...)

function PyPlot.plot(x::AbstractArray{Measured,1},y::AbstractArray{Measured,1};marker="",color="b",linestyle="-",fillstyle="white",alpha=0.5,arealinestyle="",hatch="")
    (fillstyle=="white" || fillstyle=="w")? (facecolor="white";fillstyle="full"): (facecolor=color)
    hatch!=""? (fillfacecolor="none";edgecolor=color;filllinewidth=0.):(fillfacecolor=color;edgecolor="none";filllinewidth=0.)
    xlm=transpose(val(x).+hcat(-err(x),err(x)))[:]
    yll=transpose(val(y).-hcat(err(y),err(y)))[:]
    ymm=transpose(val(y).+hcat(err(y),err(y)))[:]
    hf=PyPlot.fill_between(xlm,yll,ymm,facecolor=fillfacecolor,alpha=alpha,hatch=hatch,edgecolor=edgecolor,linewidth=filllinewidth)
    h=PyPlot.plot(val(x),val(y),linestyle=linestyle,color=color)
#    PyPlot.setp(h[1],color=facecolor,marker=marker,markeredgecolor=color,markerfacecoloralt=color,fillstyle=fillstyle,linestyle=linestyle)
#    PyPlot.setp(h[2],marker=" ",color=color)
#    PyPlot.setp(h[3],color=color)
end

   


#end # module vv



