export Column, addunit, addunit!, getname, getunit
immutable Column
    name::AbstractString
    unit::AbstractString
    power::Int64
    value::Real
    norm::Nullable{Column}
    func::Nullable{AbstractString}
end
Column(a::Column;norm::Nullable{Column}=Nullable{Column}(),func::Nullable{AbstractString}=Nullable{AbstractString}())=Column(a.name,a.unit,a.power,a.value,norm,func)
function Column(name::AbstractString;unit::AbstractString="",power::Integer=0,value::Real=1,norm::Nullable{Column}=Nullable{Column}(),func::Nullable{AbstractString}=Nullable{AbstractString}())
    Column(name,unit,power,value,norm,func)
end
function showcolumn(io::IO,a::Column;compact::Bool=false,isnorm::Bool=false)
  hasfunc=~isnull(a.func)
  hasfunc && print(io,get(a.func),"(")
    if isnorm
        a.value != 1 && print(io,a.value)
        a.value != 1 && a.power != 0 && print(io,"×")
        a.power != 0 && print(io,"10",number2exponent(a.power)," ")
    end
    print(io,a.name)
    if isnorm
        print(io," ",a.unit)
    elseif !compact
        if (hasunit= a.unit != "")
            print(io," [")
            v=a.value; p=a.power
            if !haskey(powerprefix,p)
                pk=collect(keys(powerprefix))
                i=indmin(abs(pk.-p))
                newp=pk[i]
                v*=newp-p
                p=newp
            end
            v != 1 && print(io,v," ")
            print(io,powerprefix[p],a.unit)
        end
        if !isnull(a.norm)
            hasunit||(print(io," (1");hasunit=true)
            print(io," / ")
            showcolumn(io,get(a.norm);isnorm=true)
        end
        hasunit&&print(io,"]")
    end
    hasfunc && print(io,")")
end
function showcolumn_with_color(io::IO,a::Column;compact::Bool=false,isnorm::Bool=false,unitcolor::Symbol=:yellow)
  hasfunc=~isnull(a.func)
  hasfunc && print(io,get(a.func),"(")
    if isnorm
        a.value != 1 && print_with_color(unitcolor,io,"$(a.value)")
        a.value != 1 && a.power != 0 && print_with_color(unitcolor,io,"×")
        a.power != 0 && print_with_color(unitcolor,io,"10",number2exponent(a.power)," ")
        print_with_color(unitcolor,io,a.name)
        print_with_color(unitcolor,io," ",a.unit)
    else
        print(io,a.name)
    end
    if !isnorm && !compact
        hasunit=a.unit!=""
        hasnorm=!isnull(a.norm)
        hasunit&&!hasnorm&&print_with_color(unitcolor,io,"/")
        !hasunit&&hasnorm&&print_with_color(unitcolor,io," per ")
        hasunit&&hasnorm&&print_with_color(unitcolor,io," [")
        if hasunit
            v=a.value; p=a.power
            if !haskey(powerprefix,p)
                pk=collect(keys(powerprefix))
                i=indmin(abs(pk.-p))
                newp=pk[i]
                v*=newp-p
                p=newp
            end
            v != 1 && print_with_color(io,v," ")
            print_with_color(unitcolor,io,powerprefix[p],a.unit)
        end
        hasnorm&&hasunit&&print_with_color(unitcolor,io,"/")
        hasnorm&&showcolumn_with_color(io,get(a.norm);isnorm=true,unitcolor=:blue)
        hasunit&&hasnorm&&print_with_color(unitcolor,io,"]")
    end
    hasfunc && print(io,")")
end
Base.show(io::IO,a::Column)=Base.have_color?showcolumn_with_color(io,a;compact=false):showcolumn(io,a;compact=false)
Base.showcompact(io::IO,a::Column)=showcolumn(io,a;compact=true)
plainstring(a::Column;k...)=(io=IOBuffer();showcolumn(io,a;k...);String(take!(io)))

/(a::Column,b::Column)=Column(a,norm=Nullable{Column}(b))

Column{T<:AbstractString}(cv::Array{T,1})=Column[Column(x) for x in cv]
function addunit(a::Column,unit::AbstractString;power::Integer=a.power,value::Real=a.value,norm::Nullable{Column}=a.norm,func::Nullable{AbstractString}=a.func)
    Column(a.name,unit,power,value,norm,func)
end
function addunit{T<:Column,S<:AbstractString}(cv::Array{T,1},uv::Array{S,1};k...)
    @assert length(cv)==length(uv)
    Column[ addunit(c,u;k...) for (c,u) in zip(cv,uv)]
end
function addunit!{T<:Column,S<:AbstractString}(cv::AbstractArray{T},uv::AbstractArray{S};k...)
    @assert length(cv)==length(uv)
    for i=1:length(cv); cv[i]=addunit(cv[i],uv[i];k...); end
end
addunit!{T<:Column}(cv::AbstractArray{T},u::AbstractString;k...)=addunit!(cv,repeat([u],outer=[length(cv)]);k...)

getname(a::Column)=a.name
getunit(a::Column)=a.unit
for x in (:getname,:getunit)
    @eval $x{T<:Column}(av::Array{T})=[$x(a) for a in av]
end

const powerprefix= Dict(i=>k for (i,k) in zip([12,9,6,3,0,-2,-3,-6,-9,-12,-15],["T","G","M","k","","c","m","μ","n","p","f"]) )

number2exponent(a::Integer)=number2exponent("$a")
function number2exponent(a::AbstractString)
    n=["0","1","2","3","4","5","6","7","8","9","-"]
    p=["⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹","⁻"]
    for i=1:length(n);a=replace(a,n[i],p[i]);end
    return a
end

addfunc(c::Column,fn::Function)=addfunc(c,Symbol(fn))
function addfunc(c::Column,fn::Symbol)
  sf=[:log10] # add more special functions as necessary
  rt=["log₁₀"]
  i=findfirst(fn.==sf)
  name=i>0?rt[i]:String(fn)
  Column(c.name,c.unit,c.power,c.value,c.norm,Nullable{AbstractString}(name))
end
