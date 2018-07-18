import Base: ones,zeros,+,-,*,/,get,find
export ScatteringFunction,FittingFunction,Bare,Convoluted,Full,Dispersion,Delta
export has,add,set,replace! # get,find

null=Nullable{Union{}}() # this definition won't be necessary in the future, supposedly
scale_one(p,x)=ones(size(x)...) # use size(x)... in case x is a vector of non-numbers
background_zero(p,x)=zeros(size(x)...)
"""
A specialized composite function type to hold internal functions and information 
relevant for fitting/evaluating the composite function. The type is function-input agnostic but has
been written with the intention that it will be passed an array of function parameter values and an
array of objects which are a subtype of `ScatteringInstrument` (with initial support for only `TripleAxis`).

The internal functions of the composite `ScatteringFunction` are 

  + `S` -- the core function, intended to be S(Q,E)
  
  + `P` -- a prefactor function
  
  + `B` -- a background function
  
with [`ScatteringFunction`]`(p,x) = P(p,x)*[∝S(p,x)]+B(p,x)` where `∝S(p,x)` is an indication that
the `S` function is not necessarily evaluated directly.

The other fields are `Nullable` or possibly `#undef` and should be accessed with caution.
The `Nullable` fields are

  + `useMonteCarlo` -- a `Nullable{Bool}` flag indication of whether the evaluation of `S` should utilize a Monte Carlo integration
  
  + `multiplicity` -- a `Nullable{Integer}` indicating the number of outputs of `S` for a scalar input x
  
  + `returneltype` -- a `Nullable{DataType}` indicating the element type of the (Array) output of `S`

The fields `multiplicity` and `returneltype` are used to construct the output of the composite 
function when the evaluation-type of `S` is such that some level of internal looping is required.

The last field `keywords` is `Array{Tuple{Symbol,Any},1}` and is intended to store combined keyword 
value pairs combined (or as if combined) by the slurp operator (`...`). The intended use for these
keyword value pairs is in the evaluation of `S`, but extended usage is not precluded.

Most `FittingFunction` subtypes include a `ScatteringFunction` composite function and it is the 
type of the `FittingFunction` object that ultimately determines how `S` is evaluated.
"""
type ScatteringFunction
    S::Function
    P::Function
    B::Function
    useMonteCarlo::Nullable{Bool}
    multiplicity::Nullable{Integer}
    returneltype::Nullable{DataType}
    keywords::Array{Tuple{Symbol,Any},1}
    ScatteringFunction(s::Function,p::Function,b::Function,u=null,m=null,r=null;k...)=new(s,p,b,u,m,r,k)
end
"""
    multiplicity(a::ScatteringFunction)
returns the value of the `multiplicity` field of `a`, which is a `Nullable{Integer}`. 
If a.multiplicity has not been set, this will (intentionally) throw a NullException.
"""
multiplicity(a::ScatteringFunction)=get(a.multiplicity)
"""
    returneltype(a::ScatteringFunction)
returns the value of the `returneltype` field of `a`, which is a `Nullable{DataType}`. 
If a.returneltype has not been set, this will (intentionally) throw a NullException.
"""
returneltype(a::ScatteringFunction)=get(a.returneltype)
"""
    multiplicity(a::ScatteringFunction,p,x)
returns the number of ouputs of `a.S(p,x)` (forcing `x` to be scalar, if necessary).
*This function ignores the value of a.multiplicity but a rewrite of the `S`-evaluation code would 
enable a more logically consistent behavior that returns the value of a.multiplicity or
calculates the value (storing it in a.multiplicity) if a.multiplicity is Null.*
"""
multiplicity(a::ScatteringFunction,p,x::Array)=multiplicity(a,p,x[1])
multiplicity(a::ScatteringFunction,p,x)=(r=a.S(p,x);isa(r,Tuple)?length(r[2]):length(r))
"""
    returneltype(a::ScatteringFunction,p,x)
returns the element type of the ouput of `a.S(p,x)` (forcing `x` to be scalar, if necessary).
*This function ignores the value of `a.returneltype` but a rewrite of the `S`-evaluation code would 
enable a more logically consistent behavior that returns the value of `a.returneltype` or
calculates the value (storing it in `a.returneltype`) if `a.returneltype` is Null.*
"""
returneltype(a::ScatteringFunction,p,x::Array)=returneltype(a,p,x[1])
returneltype(a::ScatteringFunction,p,x)=(r=a.S(p,x);isa(r,Tuple)?eltype(r[2]):eltype(r))
"""
    useMonteCarlo(a::ScatteringFunction)
returns the value of the `useMonteCarlo` field of `a`, or `true` if the field is `Null`.
"""
useMonteCarlo(a::ScatteringFunction)=get(a.useMonteCarlo,true)

"""
    has(a::ScatteringFunction,key::Symbol)
returns true if the symbol `key` is found in the field `keywords` of `a`
"""
has(a::ScatteringFunction,key::Symbol)         =find(a,key)>0
"""
    get(a::ScatteringFunction,key::Symbol,def)
Similar to the `get` function acting on a `Dict` or `Nullable`,
returns `def` or the last associated value of `key` if `key` is present in the `keywords` field of `a`.
"""
get(a::ScatteringFunction,key::Symbol,def)     =(pos=find(a,key);pos>0?a.keywords[pos][2]:def)
"""
    add(a::ScatteringFunction,key::Symbol,val)
Adds a key-value tuple to the `keywords` array field of `a`
"""
add(a::ScatteringFunction,key::Symbol,val)     =push!(a.keywords,(key,val))
"""
    set(a::ScatteringFunction; [key=value]xN)
replaces the `keywords` field of `a` with the new `Array{Tuple{Symbol,Any},1}` created by slurping
the `N` key-value pairs specified by `[key=value]xN`.
"""
set(a::ScatteringFunction;kwds...)             =(a.keywords=kwds)
"""
    find(a::ScatteringFunction,key::Symbol)
searches the `keywords` field of `a` for the presence of the symbol `key`.

returns 0 or the position of `key` in `a.keywords`
"""
find(a::ScatteringFunction,key::Symbol)        =findlast([k[1] for k in a.keywords].==key)
"""
    replace!(a::ScatteringFunction,key::Symbol,val)
If `key` is a keyword in the `keywords` field of `a`, its value is replaced by `val`, otherwise
an error is thrown.
"""
replace!(a::ScatteringFunction,key::Symbol,val)=(p=find(a,key);p>0?(a.keywords[p]=(a.keywords[p][1],val)):error("keyword $key not found in $(a.keywords)"))

# Implement ScatteringFunction types that should be evaluated in different ways:
#   Bare       = no convolution [defined for all (Q,E)]
#   Full       = 4D convolution [defined for all (Q,E)]
#   Dispersion = 3D convolution, with the energy integrated analytically
#   Delta      = 0D convolution, with all dimensions integrated analytically

#|     Type     |      `eltype(x)`           | `SF.S` Output | Convolution | P/B |
#|:------------:|:--------------------------:|:-------------------:|:--:|:---:|
#|   `Single`   | `Lattice4Vector` or `Real` |                     | no | no  |
#|    `Bare`    | `ScatteringInstrument`     |        S(Q,E)       | no | yes | 
#|    `Full`    | `ScatteringInstrument`     |        S(Q,E)       | 4D | yes |
#| `Dispersion` | `ScatteringInstrument`     | (ω(Q), S(Q), γ(Q))  | 3D | yes |
#|   `Delta`    | `ScatteringInstrument`     | (Q-Qₒ, S(Qₒ), γ(Qₒ)) | 0D | yes | # This line appears to be wrong
"""
`FittingFunction` is an abstract type to cover all functions used in fitting data.
Every fitting function has the fields:

| Fieldname | Description |
|:---------:|:------------|
|`SF`     | a `Function` or `ScatteringFunction` which is involved in fitting    |
|`scale`  | an overall scaling parameter to enable rescaling of the dependent    |
|         | parameter without knowledge of the parameters of `SF`                |
|`offset` | a parameter to enable adjusting the zero of the dependent parameter  |
|         | (e.g., a plotting offset) without knowledge of the parameters of `SF`|

The fields `scale` and `offset` are ignored when fitting `SF` to a dataset.

The concrete subtypes of `FittingFunction` are

|     Type     |      `eltype(x)`           | `SF.S` Output | Convolution | P/B |
|:------------:|:--------------------------:|:-------------------:|:--:|:---:|
|   `Single`   | `Lattice4Vector` or `Real` |                     | no | no  |
|    `Bare`    | `ScatteringInstrument`     |        S(Q,E)       | no | yes | 
|    `Full`    | `ScatteringInstrument`     |        S(Q,E)       | 4D | yes |
| `Dispersion` | `ScatteringInstrument`     | (ω(Q), S(Q), γ(Q))  | 3D | yes |
|   `Delta`    | `ScatteringInstrument`     |  (Qₒ, S(Qₒ), γ(Qₒ))  | 0D | yes |

Convolution of `SF.S` with an estimate of the instrumental resolution requires that it is defined 
for all (Q,E), but a `ScatteringInstrument` may not be able to access all (Q,E) due to kinematic
limitations. As a result, `SF.S` must take one or more `Lattice4Vector` as input while `SF.P` and 
`SF.B` can take the raw `ScatteringInstrument` objects as they are not convoluted and represent 
purely instrumental effects.

The convolution is performed in 4D (momentum, energy) space.
Reduced convolution dimensionality indicates that one or more dimensions are convoluted analytically.
In the case of `Dispersion` the analytical dimension is E and `SF.S` returns a tuple of arrays of
the characteristic energy, intensity, and energy width at the Q points specified in its input.
In the case of `Delta` all dimensions are handled analytically so `SF.S` must return a tuple of 
arrays of the position of the (nearest) δ-function, the integrated intensity of the δ-function, and
the intrinsic Lorentzian width of that δ-function (if it is not a true δ-function). To allow for
efficient handling of the 4D positions and widths, `Qₒ` must be returned as an `Array{Lattice4Vector}` 
and `γ(Qₒ)` as either an `Array{Real}` or `Array{Lattice4Tensor}` where the former is the special
case of a scalar times the identity rank-2 4×4 tensor.

If `x` is a `ScatteringInstrument`, `SF.S` is allowed to return Array(s) with one more dimension 
than `x`, where the extra dimension encodes information from multiple modes/branches/δ-functions and 
is summed over to produce a single value for each element of `x`.
"""
abstract type FittingFunction end
type Single <:FittingFunction; SF::Function; scale::Number; offset::Number; end
Single(S::Function)=Single(S,1,0)
type Bare  <: FittingFunction; SF::ScatteringFunction; scale::Number; offset::Number; end
abstract type Convoluted <: FittingFunction end
type Full       <: Convoluted; SF::ScatteringFunction; scale::Number; offset::Number; end
type Dispersion <: Convoluted; SF::ScatteringFunction; scale::Number; offset::Number; end
type Delta      <: Convoluted; SF::ScatteringFunction; scale::Number; offset::Number; end
# 
type MetaConvoluted <: FittingFunction; C::Convoluted; M::Function; end # TODO finish implementing meta fitting functions to allow for fitting convolution parameters

# Use metaprograming and delegation to redirect all (appropriate) functions to their inner ScatteringFunction
for x in (:Bare,:Full,:Dispersion,:Delta)
#    # create constructors that insert P(p,x)=ones(size(x)...) and B(p,x)=zeros(size(x)...)
#    @eval $x(S::Function,P::Function=scale_one,B::Function=background_zero,s=1.,o=0.;kwds...)=$x(ScatteringFunction(S,P,B;kwds...),s,o)
# FIXME we can't have it both ways in v0.5 :(
    # and similar constructors that take everything as keyworkd arguments
    @eval $x(S::Function;P::Function=scale_one,B::Function=background_zero,s=1.,o=0.,kwds...)=$x(ScatteringFunction(S,P,B;kwds...),s,o)
    # dispatch non-type-preserving functions that take any number of additional inputs and any number of keywords
    for y in (:useMonteCarlo,:multiplicity,:returneltype,:has,:get,:add,:set,:find,:replace!)
        @eval $y(a::$x,o...;k...)=$y(a.SF,o...;k...)
    end
end
# Introduce select conversion functions:
Bare(a::Full)=Bare(a.SF,a.scale,a.offset)
Full(a::Bare)=Full(a.SF,a.scale,a.offset)
# Allow for simplified calling: evaluate is overloaded for each type of FittingFunction
(f::Single)(p,x;kwds...)=f.scale*f.SF(p,x;kwds...)+f.offset
for z in (:Bare,:Full,:Dispersion,:Delta)
    @eval (f::$z)(p,x;kwds...)=f.scale*evaluate(f,p,x;kwds...)+f.offset
end
# basic arithmetic operations overloaded to modify scale and offset
(+){T<:FittingFunction}(a::T,b::Number)=(c=deepcopy(a);c.offset+=b;c)
(+){T<:FittingFunction}(a::Number,b::T)=(c=deepcopy(b);c.offset+=a;c)
(-){T<:FittingFunction}(a::T,b::Number)=(c=deepcopy(a);c.offset-=b;c)
(-){T<:FittingFunction}(a::T)=(c=deepcopy(a);c.offset*=-1;c.scale*=-1;c)
(-){T<:FittingFunction}(a::Number,b::T)=a+(-b)
(*){T<:FittingFunction}(a::T,b::Number)=(c=deepcopy(a);c.scale*=b;c)
(*){T<:FittingFunction}(a::Number,b::T)=(c=deepcopy(b);c.scale*=a;c)
(/){T<:FittingFunction}(a::T,b::Number)=(c=deepcopy(a);c.scale/=b;c)

# The MetaConvoluted relies on its wrapped function
(f::MetaConvoluted)(p,x;k...)=f.C.scale*evaluate(f,p,x;k...)+f.C.offset 
for z in (:+,:-,:*,:/)
    @eval $z(a::MetaConvoluted,b::Number)=MetaConvoluted($z(a.C,b),a.M)
    @eval $z(a::Number,b::MetaConvoluted)=MetaConvoluted($z(a,b.C),b.M)
end
