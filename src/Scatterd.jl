# eventually one should add support for NLopt.jl style optimization!
# https://github.com/JuliaOpt/NLopt.jl
# which builds on top of the NLopt library out of MIT
# http://ab-initio.mit.edu/wiki/index.php/NLopt

# Overload some useful Base functions:
Base.size(a::Scatterd{0})=Base.size(a.data)
Base.size(a::Scatterd{0},n::Integer)=Base.size(a.data,n)
Base.length(a::Scatterd)=Base.size(a.data,1)

import Base: normalize,normalize!
# Normalization of (all) intensity column(s) by a named/numbered column
function normalize!(a::Scatterd{0},nc,r::Union{T,AbstractArray{T,1}}) where T<:Measured
    normidx=getIdx(a,nc)
    n=getDat(a,normidx) # Measured (value,variance) pairs of the normalization column
    vr=value( isa(r,Array) ? mean(r) : r )
    vrpow= round(Int,floor(log10(vr)))
    vrval= vr/10^vrpow
    normCol=Nullable( Column( name(getcolumns(a)[normidx]); unit=unit(getcolumns(a)[normidx]), value=vrval, power=vrpow ) )
    for c in [a.counters;a.timers] # We should normalize all counters and timers
        thisidx=getIdx(a,c)
        if normidx == thisidx # the normalization column is this counter
            setVal(a,thisidx,r)
            novalcol=Column( name(getcolumns(a)[thisidx]);unit=unit(getcolumns(a)[thisidx]))
            a.columns[thisidx]=Column(novalcol;norm=Nullable(novalcol))
        else
            setVal(a,thisidx, getDat(a,thisidx)./n.*r) # error propagation is handled correctly for the Measured type
            a.columns[thisidx]=Column(getcolumns(a)[getIdx(a,c)];norm=normCol)
        end
    end
end
function normalize!(a::Scatterd{N},nc,r::Union{T,AbstractArray{T,1}}) where {T<:Measured,N}
    colsofa=getcolumns(a)
    normidx=getIdx(a,nc)
    n=get_data_or_detector_column(a,nc) # might be (npts,) or (npts,ndetX,ndetY,...)
    vr=value( isa(r,Array) ? mean(r) : r )
    vrpow= round(Int,floor(log10(vr)))
    vrval= vr/10^vrpow
    normCol=Nullable( Column( name(colsofa[normidx]); unit=unit(colsofa[normidx]), value=vrval, power=vrpow ) )
    for c in [a.counters;a.timers] # We should normalize all counters and timers
        thisidx=getIdx(a,c)
        if normidx == thisidx # the normalization column is this counter)
            novalcol=Column( name(colsofa[thisidx]);unit=unit(colsofa[thisidx]))
            setDat(a,thisidx,r,Column(novalcol;norm=Nullable(novalcol)))
        else
            # each column may be (npts,) or (npts,ndetX,ndetY,...)
            # ./ *should* broadcast correctly by scaling-up as necessary
            # but we haven't done any checking to see if r is compatible FIXME (restrict r to be a scalar?)
            setDat(a,thisidx, get_data_or_detector_column(a,thisidx)./n.*r, Column(colsofa[getIdx(a,c)];norm=normCol))
        end
    end
end

function normalize!(a::Scatterd,nc,r::Union{T,AbstractArray{T,1}},dr::Union{T,AbstractArray{T,1}}=T(0)) where T<:Number
    normalize!(a,nc,Measured(r,dr.^2))
end
normalize(a::T,x...) where T<:Scatterd =(n=T(a);normalize!(n,x...);n)
# normalization of multiple Scatterd objects in an array:
#normalize!{T<:Scatterd{0}}(a::Array{T},x...)=map(y->normalize!(y,x...),a)
#normalize{T<:Scatterd{0}}(a::Array{T},x...)=map(y->normalize(y,x...),a)

"""
    mask!(a::Scatterd,range::AbstractVector{Number},whichcol=a.y)
    mask!(a::Scatterd,bound::Number,whichcol=a.y)
    mask!(a::Scatterd,badpoints::BitArray)

Mask the points in a scan using one of a `range`, an (typically upper) `bound`,
or a `BitArray` list of the points to be masked.

If a `range` is passed any value of the column specified by `whichcol` less than
the `minimum(range)` or more than `maximum(range)` has its mask flag set `true`.
If a `bound` is passed, a range of [0,`bound`] is created and then the same
procedure is performed.
In both cases `whichcol` can be anything that `getVal` accepts as a second input.
"""
# Mask the points in a scan outside of a specified range for some parameter
mask!(a::Scatterd{0},r::Number,wc=a.y)=mask!(a,[zero(r),r],wc)
function mask!(a::Scatterd{0},r::AbstractArray{T,1},wc=a.y) where T<:Number
    bad = (getVal(a,wc).<minimum(r)).|(getVal(a,wc).>maximum(r))
    a.mask[bad]=true; # bad points get masked
end
for i=0:2
    @eval mask!(a::Scatterd{$i},bad::BitArray{1})=(a.mask[bad]=true;)
end
mask{T<:Scatterd{0}}(a::T,x...)=(n=T(a);mask!(n,x...);n)

function Base.show(io::IO,p::Scatterd{D}) where D
    splitf=split(p.filename,"\n")
    splitc=split(p.command,"\n")
    names=["file","command","timer","counter","motor"]
    valus=[splitf,splitc,p.timers,p.counters,p.motors]
    lengt=map(length,valus)
    names[lengt.>1].*="s"
    lstrs=map(string,lengt)
    push!(names,"data")
    push!(lstrs,join(["$x" for x in size(p.data)],"×"))
    if D>0
        push!(names,"detector")
        push!(lstrs,join(["$x" for x in size(p.detector)],"×"))
        valus=vcat(valus,[getnames(p.data)],[getnames(p.detector)])
    else
        valus=vcat(valus,[p.columns]) # FIXME we can easily split "data" and "detector" columns once all instruments use NamedArrays
    end
    lstrs="(".*lstrs.*"): "
    lstrs[lstrs.=="(1): "]="   : "
    pre=map(length,names)
    post=map(length,lstrs)
    ex=maximum(pre+post)-pre-post+1
    pre=" ".*names.*map(x->" "^x,ex).*lstrs
    print(io,pre[1]);showpatharray(io,splitf);print(io,"\n")
    print(io,pre[2]);showscanarray(io,splitc);print(io,"\n")
    for i=3:length(pre); println(io,pre[i]*join(valus[i]," ")); end
end
#Base.show(io::IO,::MIME"text/plain",m::AbstractArray{T,1}) where {T<:Scatterd} = Base.show(io,"$(summary(m)):\n",m)
#Base.show(io::IO,::MIME"text/plain",m::AbstractArray{T,N}) where {T<:Scatterd,N<:Any} = Base.show(io,"$(summary(m)):\n",m)
Base.show(io::IO,::MIME"text/plain",s::Scatterd) = Base.print(io,summary(s),":\n",s)

function showpatharray{S<:AbstractString}(io::IO,sa::Array{S,1};color::Symbol=:blue)
    sat=copy(sa)
    l=length(sat)
    1==l && (Base.print(io,sat[1]); return)
    std(map(length,sat))==0 || error("I don't know how to handle different-length strings yet.")
    sat1s=split(sat[1],"")
    notunique=sat1s.==split(sat[2],"")
    for i=2:l; notunique.&=sat1s.==split(sat[i],""); end
    # notunique is true for character positions which are the same in all strings
    ctr=sum(abs,diff(notunique))
    if 3>ctr # one state change -> T...TF...F or F...FT...T
             # two state changes T...TF...FT...T or F...FT...TF...F
        preno=findfirst(.!notunique)
        postno=findlast(.!notunique)
        pre=sat[1][1:preno-1]
        post=sat[1][postno+1:end]
        satunique=Array{eltype(sat)}(l)
        for i=1:l; satunique[i]=sat[i][preno:postno]; end
        # satunique is *probably* a list of numbers
        Base.print(io,pre)
        try
            uniquenos=sort(map(parse,satunique)) # sort will throw an error if not all numbers
            eltype(uniquenos)<:Integer || error("Parsed result not integers")
            dun=diff(uniquenos)
            0==std(dun) || error("List does not have constant step size")
            Base.print_with_color(color,io,uniquenos[1],1==dun[1]?":":":$(dun[1]):",uniquenos[end])
        catch
            Base.print_with_color(color,io,"[",join(satunique,","),"]")
        end
        Base.print(io,post)
    else
        Base.print(io,"[",sat[1])
        Base.print_with_color(color,io," … ")
        Base.print(io,sat[end],"]")
    end
end

function showscanarray{S<:AbstractString}(io::IO,sa::Array{S,1};color::Symbol=:blue)
    l=length(sa)
    1==l && (Base.print(io,sa[1]);return)
    usa=unique(sa)
    ul=length(usa)
    if 1==ul # only one scan, just print that
        Base.print(io,usa[1])
    elseif l==ul # nothing to compact, just print the unique array
        Base.print(io,join(usa,", "))
    else # now the hard part
        uidx=Array{Int}(l)
        for i=1:l; uidx[i]=findfirst(sa[i].==usa); end
        Base.print_with_color(color,io,"$uidx")
        Base.print(io," in")
        for i=1:ul
            Base.print_with_color(color,io,"\n\t$i)")
            Base.print(io,usa[i])
        end
    end
end
