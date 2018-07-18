function uniqueBA(v)
    # create a false BitArray the same size as v.
    # size(v) is a tuple, ... unrolls it
    b=falses(size(v)...)
    for i=1:length(v)
        matching = v .== v[i]
        if sum(matching)>0
            # second matches only set the first value true
            b[findfirst(matching)]=true
        end
    end
    return b
end
"""
    matchBA(v1,v2)
Returns a BitArray, `b`, the same size and shape as `v1` with the elements of `b` indicating
the presence of the elements of `v1` in `v2`.

For example `matchBA([1,4,3,5,2],[3,1,2])` would return `[true,false,true,false,true]`

If `v1` is an `Array{AbstractString}` then a special form of `matchBA` exists that takes a `Regex` for `v2`.
This allows to form a `BitArray` indicating which strings in `v1` match the regular expression
(or expressions if an array) in `v2`.
"""
function matchBA{T<:AbstractString}(v1::Array{T},v2::Regex)
    b=BitArray{ndims(v1)}(size(v1)...)
    for i=1:length(v1); b[i]=!isa( match(v2,v1[i]), Void ); end
    return b
end
function matchBA(v1::Array,v2::Array)
    b=falses(size(v1)...) # make sure the BitArray is initialized
    for i in v2; b.|=matchBA(v1,i); end
    return b
end
matchBA(v1::Array,v2)= v1.==v2 # for the general case use a broadcast equivalence test

pad(a::Number,b::Int=6)= ( n=a>0?round(Int,floor(log10(a)))+1:1; b-=n; (b>0?"$("0"^b)":"")*"$a" )

function hasbeenfit(s::Scatterd) # utility to check if a scan has been fit
    (isdefined(s,:fitres)&&isdefined(s.fitres,:p))||(isdefined(s,:fits)&&length(s.fits)>0&&isdefined(s.fits[1],:p)) # FIXME needs to me modified again once only fitres or fits exists
end
hasbeenfit(s::FitResult)=isdefined(s,:p) # does the FitResult represent something which has been fit?
"""
    chi2(::FitResult)

The value minimized by `nlfit` for a function `f(p,x)` is

    Σ |( f(p,x) - y )/w|²

which is sometimes called the raw-χ², sum of the squared weighted error, or net residual.
Since it is usually possible to reduce the net residual by increasing the complexity
of `f` (and therefore the number of parameters that it depends on) one desires a means
by which to compare residuals between differently-complex functions that favors simplicity.
The normalized χ² is

    N⁻¹ Σ |( f(p,x) - y )/w|²

where `N` is the number of degrees of freedom, or the `length(x)-length(p)`.
This function computes and returns the normalized χ² from information stored within the
`FitResult` object.
"""
chi2(a::FitResult)= sum(abs2,a.residual)/a.dof
#chi2(a::Scatterd)=hasbeenfit(a)?chi2(a.fitres):NaN # push this onto the fit-result
function chi2(a::Scatterd)
    hasbeenfit(a)||(return NaN); f=a.fitres
    r=getErr(a,a.y) # the uncertainty in each point
    m=.!a.mask # points included in the fit
    length(f.residual)==sum(m) || (warn("Can not compute χ², was the mask modified after fitting?");return NaN)
    sum(abs2,f.residual./r[m])/f.dof
end

function newchi2(a::Scatterd)
    hasbeenfit(a)||(return NaN); f=a.fitres
    r=getErr(a,a.y); m=.!a.mask
    length(f.residual)==sum(m) || (status(:error,"can not compute χ²");status(:hint,"was the mas modified after fitting?");return NaN)
    sum(abs2,f.residual./r[m])/sum(m)
end
export newchi2


function createcolorlist(a,c=["r","o","y","g","b","i","v","k"],f=["full","left","right","bottom","top"],s=["o","d","^","<","v",">","s"],l=["none","--","-.",":"])
    la=(ndims(a)==0)?(a):(length(a))
    (c,f,s,l)=ndgridvecs(c,f,s,l)
    while length(c)<la; c=[c;c]; f=[f;f]; l=[l;l]; s=[s;s]; end
    c=c[1:la];f=f[1:la];l=l[1:la];s=s[1:la] # truncate output
    (c,f,s,l)
end

"""
Create vectors of color, fillstyle, marker, linestyle, and facecolor values that cycle through
all permutations of the first four with a separately cycling facecolor vector.
Optional keyword arguments `color`, `fillstyle`, `marker`, `linestyle` and `facecolor` accepted,
these can be single strings or arrays of strings, with arrays treated as vectors.
"""
function createcolorlistkeyed(a;color=["r","o","y","g","b","i","v","k"],
                                fillstyle=["full","left","right","bottom","top"],
                                marker=["o","d","^","<","v",">","s"],
                                linestyle=["none","--","-.",":"],
                                facecolor=["light","white"], ignoredkeys...)
    # ensure we have at Arrays
    isa(color,AbstractString)     && (color=[color])
    isa(fillstyle,AbstractString) && (fillstyle=[fillstyle])
    isa(marker,AbstractString)    && (marker=[marker])
    isa(linestyle,AbstractString) && (linestyle=[linestyle])
    isa(facecolor,AbstractString) && (facecolor=[facecolor])
    # and that those arrays are vectors
    1<ndims(color)     && (color=vec(color))
    1<ndims(fillstyle) && (fillstyle=vec(fillstyle))
    1<ndims(marker)    && (marker=vec(marker))
    1<ndims(linestyle) && (linestyle=vec(linestyle))
    1<ndims(facecolor) && (facecolor=vec(facecolor))
    # determine the number of output elements we need in each vector
    la= 0==ndims(a) ? a : length(a)
    (c,f,s,l)=ndgridvecs(color,fillstyle,marker,linestyle) # facecolor is handled separately
    lm = length(c) # the total number of elements in each c,f,s,l
    if lm<la
        rep=cld(la,lm)
        c=repeat(c,outer=rep)
        f=repeat(f,outer=rep)
        s=repeat(s,outer=rep)
        l=repeat(l,outer=rep)
    end
    x=vec(repeat(facecolor,outer=[1,cld(la,length(facecolor))]))
    c=c[1:la]; f=f[1:la]; s=s[1:la]; l=l[1:la]; x=x[1:la] # truncate vectors
    (c,f,s,l,x)
end


## getBytes modified from http://stackoverflow.com/questions/28402989/check-size-in-bytes-of-variable-using-julia
#getBytes(x::DataType) = sizeof(x);
#function getBytes(x)
#   total = 0;
#   fieldNames = fieldnames(typeof(x));
#   if fieldNames == []
#      return sizeof(x);
#   else
#     for fieldName in fieldNames
#        if !isa(getfield(x,fieldName),Function)
#            total += getBytes(getfield(x,fieldName));
#        end
#     end
#     return total;
#   end
#end

getMultiple(a::Dict,v::AbstractArray,default)=map(x->Base.get(a,x,default),v)
compatible(a,b)=(ndims(a)==ndims(b))&&all(collect(size(a)).==collect(size(b)))
compatible(a,b,c...)=compatible(a,b)&&compatible(b,c...)

# ndgrid code modified from the Julia examples, plus new ndgridvecs
ndgrid(v::AbstractVector) = copy(v)
function ndgrid(v1::AbstractVector, v2::AbstractVector)
    m, n = length(v1), length(v2)
    v1 = reshape(v1, m, 1)
    v2 = reshape(v2, 1, n)
    (repmat(v1, 1, n), repmat(v2, m, 1))
end
function ndgrid_fill(a, v, s, snext)
    for j = 1:length(a)
        a[j] = v[div(rem(j-1, snext), s)+1]
    end
end
function ndgrid(vs::AbstractVector...)
    n = length(vs)
    sz = map(length, vs)
    out = ntuple(i->Array{eltype(vs[i])}(sz),n)
    s = 1
    for i=1:n
        a = out[i]::Array
        v = vs[i]
        snext = s*size(a,i)
        ndgrid_fill(a, v, s, snext)
        s = snext
    end
    out
end
ndgridvecs(vs::AbstractVector...)=map(vec,ndgrid(vs...))

# these rotation matricies are rotations that change the direction of a vector within the x-y plane
for (f,i) in zip((:vectorrotation3x,:vectorrotation3y,:vectorrotation3z),(:_vr3x,:_vr3y,:_vr3z))
    @eval $f(θ::T) where T<:AbstractFloat = $i(one(T),zero(T),cos(θ),sin(θ))
    @eval $f(θ::Number) = $f(float(θ))
end
for (f,i) in zip((:vectorrotation3xd,:vectorrotation3yd,:vectorrotation3zd),(:_vr3x,:_vr3y,:_vr3z))
    @eval $f(θ::T) where T<:AbstractFloat = $i(one(T),zero(T),cosd(θ),sind(θ))
    @eval $f(θ::Number) = $f(float(θ))
end
@deprecate vectorrotation3xy(θ::Number) vectorrotaion3z(θ::Number)

_vr3x(o::T,z::T,c::T,s::T) where T =  T[ o  z  z;  z  c -s;  z  s  c]
_vr3y(o::T,z::T,c::T,s::T) where T =  T[ c  z  s;  z  o  z; -s  z  c]
_vr3z(o::T,z::T,c::T,s::T) where T =  T[ c -s  z;  s  c  z;  z  z  o]

vectorrotation4xy(θ::Number)=typeof(θ)[cos(θ) -sin(θ) 0 0; sin(θ) cos(θ) 0 0; 0 0 1 0; 0 0 0 1]
# these rotation matricies are basis-changes that keep the direction of the vector unchanged while the coordinate system is rotated
basisrotation3xy(θ::Number)=typeof(θ)[cos(θ) sin(θ) 0; -sin(θ) cos(θ) 0; 0 0 1]
basisrotation4xy(θ::Number)=typeof(θ)[cos(θ) sin(θ) 0 0; -sin(θ) cos(θ) 0 0; 0 0 1 0; 0 0 0 1]

isInvertable{T}(x::AbstractMatrix{T})=((det(x)==0)|isnan(det(x))?false:true)
