_integrateProject(r::ResFloat,m::Matrix{ResFloat},i::Integer)=(sqrt(2pi/m[i,i])*r,_integrateProject(m,i))
function _integrateProject(m::Matrix{ResFloat},i::Integer)
    # takes a matrix M, removes the ith row and column, and returns Mp
    # with the property |M|/M[i,i] = |Mp|
    (a,b)=size(m); @assert a==b "m is not square"
    (a==1) && (return m)
    k=1:a.!=i
    mp=m[k,k]
    b=(m[k,i]+m[i,k])/2
    mp-=(b*b')/m[i,i]
end
_integrateProject{T<:Matrix{ResFloat}}(mv::Array{T,1},i::Integer)=map(x->_integrateProject(x,i),mv)
_integrateProject{T<:Matrix{ResFloat},N}(mv::Array{T,N},i::Integer)=map(x->_integrateProject(x,i),mv)
_integrateProject{T<:ResRM}(rm::Array{T,1},i::Integer)= map(x->_integrateProject(_getM(x),i),rm)

function _getMij(a::Matrix{ResFloat},i::Integer,j::Integer)
    @assert (i<=size(a,1))&(i>0)&(j<=size(a,2))&(j>0)
    a[i,j]
end
function _getMij{T<:Matrix{ResFloat}}(a::Array{T,1},i::Integer,j::Integer)
    @assert (i<=size(a[1],1))&(i>0)&(j<=size(a[1],2))&(j>0)
    la=length(a)
    v=zeros(ResFloat,la)
    for k=1:la; v[k]=a[k][i,j]; end
    v
end
function _getMij{T<:Matrix{ResFloat},N}(a::Array{T,N},i::Integer,j::Integer)
    @assert (i<=size(a[1],1))&(i>0)&(j<=size(a[1],2))&(j>0)
    sa=size(a)
    v=zeros(ResFloat,sa...)
    for k=1:prod(sa...); v[k]=a[k][i,j]; end
    v
end

function verifymultiplicityandreturneltype!(p,t,f::FittingFunction)
    isnull(f.SF.multiplicity) && (f.SF.multiplicity=multiplicity(f.SF,p,get4vector(t)))
    isnull(f.SF.returneltype) && (f.SF.returneltype=returneltype(f.SF,p,get4vector(t)))
end

"""
    applyPB(ff::FittingFunction,p,t::TripleAxis,conv)
A `FittingFunction` is a composite type mainly comprised of three `Function`s, `S`, `P`, and `B`.
Different subtypes of `FittingFunction` exist to allow for disparate evaluation methods of the
`S` function, but all `FittingFunction` objects need to evaluate `P` and `B` in the same manner.
This function takes the `FittingFunction` object, `ff`, a fitting parameter (`Array`, likely), `p`,
a `TripleAxis` object, `t`, and the result of evaluating `S(p,t)`→`conv` and returns
`P(p,t)*conv+B(p,t)` while ensuring that the result has the same shape as `t`.
"""
function applyPB{T<:TripleAxis,N}(ff::FittingFunction,p,t::Array{T,N},conv;
        method::Function=get(ff,:method,popoviciRmM),
        RM::Array{ResRM,N}=get(ff,:RM,method.(t)::Array{ResRM,N}),kwds...)
    @assert ndims(conv)<3 "conv should only have size (N,) or (N,M)"
    conv.*=ff.SF.P(p,t[:]) # Multiply by the slow prefactor (the .* broadcasting multiplication scales-up PQω as necessary)
    (ndims(conv)==2)&&(conv=sum(conv,2)) # if necessary, sum over all modes; conv should only have size (N,) or (N,modes)
    conv=reshape(conv,size(t)) # reshape conv to be the same shape as the input TripleAxis Array
    conv.*=_getR.(RM)    # rescale conv by R0 (or Rm)
    conv+=ff.SF.B(p,t)      # add on the slow background
end
["""
The cross section expected at a point `x⃗ₒ=(Q⃗ₒ,Eₒ)` from a scattering function `S(x⃗)` is proportional
to `∫∫∫∫dx⃗ S(x⃗) R(x⃗-x⃗ₒ)` where `R(x⃗)=Rₒ sqrt(|M|)/(2π)²×exp(-½ x⃗'Mx⃗)` is the resolution function
and `M` is the resolution matrix. Both `Rₒ` and `M` depend on the instrument configuration at `x⃗ₒ`.
The functions `convolute` compute the 4-dimensional convolution for `R(x⃗)` and `S(x⃗)` for a list of
given `x⃗ₒ` points according to the method specified by the type and options of the `FittingFunction`
object which contains `S(x⃗)`.
"""]

"""
The functions `evaluate` use multiple dispatch to
"""
evaluate(S,p,s::Scatterd{0};k...)  =evaluate(S,p,s.instrument;k...)
evaluate(S,p,s::Scatterd;k...)=error("Please extend evaluation to Scatterd{>0}")
evaluate(S,p,t::TripleAxis;k...)=evaluate(S,p,[t];k...)[1]
function evaluate{T<:TripleAxis,N}(ff::Bare,p,t::Array{T,N};kwds...)
    conv=ff.SF.S(p,get4vector(t[:]))
    conv=applyPB(ff,p,t,conv;kwds...)
end
function evaluate{T<:TripleAxis,N}(ff::Convoluted,p,t::Array{T,N};kwds...)
    conv=convolute(ff,p,t;kwds...)
    conv=applyPB(ff,p,t,conv;kwds...)
end
"""
`MetaConvoluted` functions modify their resolution parameters via a meta function `mff.M`
and then convolute the `Convoluted` function `mff.C` with newly calculated resolution matricies.
Passing a precalculated resolution matrix array via the keyword `RM` will prevent the recalculation!!
"""
function evaluate(mff::MetaConvoluted,p,t::Array;kwds...)
    mt=mff.M(p,t)
    conv=convolute(mff.C,p,mt;kwds...)
    conv=applyPB(mff.C,p,mt,conv;kwds...)
end
function convolute{T<:TripleAxis,N}(ff::Union{Full,Dispersion},p,t::Array{T,N};uMC::Bool=useMonteCarlo(ff),kwds...)
    verifymultiplicityandreturneltype!(p,t,ff) # ensure that ff.modes and ff.returntype are set
    conv=uMC?convolute_MC(ff,p,t;kwds...):convolute_GH(ff,p,t;kwds...) # TODO put in a test for _GH case to prevent large n causing out-of-memory issues
end
function convolute_MC{T<:TripleAxis,N}(ff::Full,p,t::Array{T,N};
        seed::Integer=get(ff,:seed,-1),
        rng::AbstractRNG=get(ff,:rng,has(ff,:rng)?nothing:signbit(seed)?MersenneTwister(rand(UInt32)):MersenneTwister(seed)),
        points::Integer=get(ff,:points,4*151),
        rv::Array{ResFloat,2}=get(ff,:rv,has(ff,:rv)?nothing:zeros(ResFloat,4,points)),
        xv::Array{ResFloat,2}=get(ff,:xv,has(ff,:xv)?nothing:zeros(ResFloat,4,points)),
        qv::Array{Lattice4Vector{ResFloat},1}=get(ff,:qv,has(ff,:qv)?nothing:Array{Lattice4Vector{ResFloat}}(points)),
        method::Function=get(ff,:method,popoviciRmM),
        RM::Array{ResRM,N}=get(ff,:RM, has(ff,:RM) ? nothing : method.(t)::Array{ResRM,N}),
        kwds...)
    # the number of dimensions and size of RM should be the same as the size of t
    @assert compatible(t,RM)
    net=length(t) # number of elements t contains
    # The convolution output is stored in conv, which will later be reshaped to the size of t
    modes=multiplicity(ff) # this will throw a NullException if ff.SF.multiplicity isn't set!
    rtype=returneltype(ff) # this will throw a NullException if ff.SF.returneltype isn't set!
    conv=(modes>1)?zeros(rtype,net,modes):zeros(rtype,net)
    for i=1:net # loop over the number of elements in t
        # get random points distributed following a normal distribution with σ=1, centered on zero
        randn!(rng,rv);
        # Rescale the random values by eigen-Gaussian-widths and move from the eigen-coordinate
        # system to Popovici's (x,y,z,E) frame, this is x⃗-x⃗ₒ.
        # The eigen values of M are  Mᵢ=1/(σ²ᵢ) so σᵢ=1/sqrt(Mᵢ). Multiplying the random vectors
        # by [σ₁,σ₂,σ₃,σ₄] (elementwise) changes the standard deviation of the distribution
        # appropriately. The eigenbasis of the matrix M is defined by M=VSV⁻¹; where S is a diagonal
        # matrix with elements Mᵢ, and V is a square matrix comprised of column vectors of
        # the eigenvectors of M, vᵢ.
        # A vector, x⃗, can be transformed from the basis of M (the Popovici basis) to the eigenbasis
        # by  v⃗=V⁻¹x⃗; and so a vector in the eigenbasis, v⃗, can be transformed to
        # the basis of M by x⃗=Vv⃗
        rv./=sqrt.(eigvals(RM[i])) # uses broadcast to scale-up singleton dimensions
        _multMv!(eigvecs(RM[i]),rv,xv) # _multMv!(M,v,Mv) performs Mv=M*v for the n column vectors in v
        # perform the basis change from Popvici's basis to the Sample reciprocal lattice basis,
        lab2sample!(rv,t[i],qv) # stores the resultant Lattice4Vectors in qv
        Sx⃗ = ff.SF.S(p,get4vector(t[i])+qv)
        # a general MonteCarlo integral is the sum of the function [here S(x⃗)R(x⃗-x⃗ₒ)] evaluated
        # at N randomly chosen points divided by N. The normalization by N shall come later
        (modes>1)?conv[i,:]+=sum(Sx⃗,1):conv[i]+=sum(Sx⃗)
    end
    conv./=points  # rescale by the overall convolution normalization (1/#points for MonteCarlo)
end
function convolute_GH{T<:TripleAxis,N}(ff::Full,p,t::Array{T,N};
        n::Integer=get(ff,:n,6),
        xv::Array{ResFloat,2}=get(ff,:xv,has(ff,:xv)?nothing:zeros(ResFloat,4,n^4)),
        qv::Array{Lattice4Vector{ResFloat},1}=get(ff,:qv,has(ff,:qv)?nothing:Array{Lattice4Vector{ResFloat}}(n^4)),
        method::Function=get(ff,:method,popoviciRmM),
        RM::Array{ResRM,N}=get(ff,:RM, has(ff,:RM) ? nothing : method.(t)::Array{ResRM,N}),
        kwds...)
# function convolute_GH{T<:TripleAxis,N}(ff::Full,p,t::Array{T,N};
#         n::Integer=get(ff,:n,151),
#         xv::Array{ResFloat,2}=get(ff,:xv,zeros(ResFloat,4,4n)),
#         qv::Array{Lattice4Vector{ResFloat},1}=get(ff,:qv,Array{Lattice4Vector{ResFloat}}(4n)),
#         method::Function=get(ff,:method,popoviciRmM),
#         RM::Array{ResRM,N}=get(ff,:RM,method.(t)::Array{ResRM,N}),kwds...)
    DEBUG && status(:debug,"convolute_GH called with n=$n")
    @assert compatible(t,RM)
    net=length(t) # number of elements t contains
    # The convolution output is stored in conv, which will later be reshaped to the size of t
    modes=multiplicity(ff) # this will throw a NullException if ff.SF.multiplicity isn't set!
    rtype=returneltype(ff) # this will throw a NullException if ff.SF.returneltype isn't set!
    conv=(modes>1)?zeros(rtype,net,modes):zeros(rtype,net)
    # Use the Gauss-Hermite quadrature method to evaluate the convolution. For this we need roots
    # and weights for the 4-directions. Since we use the eigenvectors of M, there's no easy way
    # to single out a direction as being less important -- all are integrated with the same number
    # of points. r and w are either passed in or calculated in the function call and describe
    # the roots and weights for each direction;
    # we need to evaluate Sum(w₁w₂w₃w₄S(σ₁r₁,σ₂r₂,σ₃r₃,σ₄r₄)) over the 4D grid for each
    # element of the Array{TripleAxis}
    (r,w)=gausshermite(n)
    rv=hcat(ndgridvecs(r,r,r,r)...)' # 4xlength(r) Array{ResFloat,2}
    wv=prod(hcat(ndgridvecs(w,w,w,w)...),2) # length(r) Array{ResFloat,1}.
    # there's no need to scale up wv to size(SQω()), since .* later will do the equivalent for us
    #
    # # alternatively, save time and just evaluate the function along the eigen-directions
    # rv=hcat([r;0r;0r;0r],[0r;r;0r;0r],[0r;0r;r;0r],[0r;0r;0r;r])'
    # wv=prod(hcat([w-1;0w;0w;0w],[0w;w-1;0w;0w],[0w;0w;w-1;0w],[0w;0w;0w;w-1])+1,2)
    for i=1:net
        # The function SQω needs to be evaluated at the rescaled Hermitian-root positions
        # along the eigen-vectors of M. The rescaling comes from the convolution being
        # ∫∫∫∫dx⃗ S(x⃗) R(x⃗-x⃗ₒ) where R(x⃗)=Rₒ sqrt(|M|)/(2π)²×exp(-½ x⃗'Mx⃗)
        # ie, ∝∫∫∫∫dx⃗ S(x⃗)×exp(-½ [x⃗-x⃗ₒ]'M[x⃗-x⃗ₒ]), while Gauss-Hermite quadrature is an
        # approximation to an integral of the form ∫∫∫∫dy⃗ S(y⃗)×exp(-½y²). If we replace
        # [x⃗-x⃗ₒ] by the eigen-vectors of M, e⃗, then [x⃗-x⃗ₒ]'M[x⃗-x⃗ₒ] becomes
        # e⃗'Me⃗= e₁M₁₁e₁ + e₂M₂₂e₂ + e₃M₃₃e₃ + e₄M₄₄e₄, and
        # exp(-½ [x⃗-x⃗ₒ]'M[x⃗-x⃗ₒ]) → exp(-½M₁₁e₁²)exp(-½M₂₂e₂²)exp(-½M₃₃e₃²)exp(-½M₄₄e₄²).
        # Each exponential is then independent and we can define yᵢ²=Mᵢᵢeᵢ² to enable use
        # of Gauss-Hermite quadradture. yᵢ=sqrt(Mᵢᵢ)eᵢ → dyᵢ=sqrt(Mᵢᵢ)deᵢ; → eᵢ=yᵢ/sqrt(Mᵢᵢ)
        # Gauss-Hermite quadrature states that we only need to evaluate S(y⃗) at the roots
        # of the Hermite polynomial (of rank n) and take an appropriately weighted sum
        # to estimate ∫∫∫∫dy⃗ S(y⃗)×exp(-½y²). The rᵢ calculated above are the roots along
        # each yᵢ-axis; first rescale them to be along ei and multiply by the eigen-vectors
        # of M to go from eᵢ to [x⃗-x⃗ₒ]
        trv = rv./sqrt.(eigvals(RM[i])) # broadcasts to scale-up singleton dimensions to the size of rv
        _multMv!(eigvecs(RM[i]),trv,xv)
        # perform the basis change from Popvici's basis to the Sample reciprocal lattice basis,
        lab2sample!(xv,t[i],qv) # stores the resultant Lattice4Vectors in qv
        # adding the converted x⃗-x⃗ₒ to Q⃗ₒ and evaluating S(Q⃗)
        Sx⃗ = wv.*ff.SF.S(p,get4vector(t[i])+qv)
        (modes>1)?conv[i,:]=sum(Sx⃗,1):conv[i]=sum(Sx⃗) # perform the Gauss-Hermite sum
    end
    conv/=pi^2  # rescale by the overall convolution normalization
end

function convolute_MC{T<:TripleAxis,N}(ff::Dispersion,p,t::Array{T,N};
        seed::Integer=get(ff,:seed,-1),
        points::Integer=get(ff,:points,100),
        rng::AbstractRNG=get(ff,:rng,has(ff,:rng)?nothing:signbit(seed)?MersenneTwister(rand(UInt32)):MersenneTwister(seed)),
        rv::Array{ResFloat,2}=get(ff,:rv,has(ff,:rv)?nothing:zeros(ResFloat,3,points)),
        xv::Array{ResFloat,2}=get(ff,:xv,has(ff,:xv)?nothing:zeros(ResFloat,3,points)),
        qv::Array{Lattice3Vector{ResFloat},1}=get(ff,:qv,has(ff,:qv)?nothing:Array{Lattice3Vector{ResFloat}}(points)),
        method::Function=get(ff,:method,popoviciRmM),
        RM::Array{ResRM,N}=get(ff,:RM, has(ff,:RM) ? nothing : method.(t)::Array{ResRM,N}),
        kwds...)
    @assert compatible(t,RM)
    net=length(t) # number of elements t contains
    modes=multiplicity(ff) # this will throw a NullException if ff.SF.multiplicity isn't set!
    rtype=returneltype(ff) # this will throw a NullException if ff.SF.returneltype isn't set!
    # The convolution output is stored in conv, which will later be reshaped to the size of t
    conv=(modes>1)?zeros(rtype,net,modes):zeros(rtype,net)
    vy=(modes>1)?zeros(rtype,points,modes):zeros(rtype,points)
    RM3=_integrateProject(RM,4)
    σω=1./sqrt.(_getMij(RM,4,4))
    hklmat=BitArray([1 1 0 1; 1 1 0 1; 0 0 0 0; 1 1 0 1])

    for i=1:net # loop over the number of elements in t
        # get random points distributed following a normal distribution with σ=1, centered on zero
        randn!(rng,rv);
        # Rescale the random values by eigen-Gaussian-widths and move from the eigen-coordinate
        # system to Popovici's (x,y,z) frame, this is x⃗-x⃗ₒ.
        rv./=sqrt.(hkleigvals(RM[i])) # uses broadcast to scale-up singleton dimensions
        _multMv!(hkleigvecs(RM[i]),rv,xv) # _multMv!(M,v,Mv) performs Mv=M*v for the n column vectors in v
        # perform the basis change from Popvici's basis to the Sample reciprocal lattice basis,
        # stores the resultant Lattice3Vectors via promotion to Lattice4Vectors in qv
        lab2sample!(rv,t[i],qv)
        QE=get4vector(t[i])+qv # (Q,E) points to evaluate the model at
        (rE0,rS,rγ)=ff.SF.S(p,QE) # returns E0(Q), I(Q), gamma(Q)
        unsafevoigt_σγ_nc!(σω[i]*ones(rS),rγ,getE(QE)-rE0,vy);
        (modes>1)?conv[i,:]+=sum(rS.*vy,1):conv[i]+=sum(rS.*vy)
    end
    conv./=points       # rescale by the overall convolution normalization (1/#points for MonteCarlo)
end
# function convolute_GH{T<:TripleAxis,N}(ff::Dispersion,p,t::Array{T,N};
#         n::Integer=get(ff,:n,151),
#         xv::Array{ResFloat,2}=get(ff,:xv,has(ff,:xv) ? nothing : zeros(ResFloat,3,3n)),
#         qv::Array{Lattice3Vector{ResFloat},1}=get(ff,:qv,has(ff,:qv)?nothing:Array{Lattice3Vector{ResFloat}}(3n)),
#         method::Function=get(ff,:method,popoviciRmM),
#         RM::Array{ResRM,N}=get(ff,:RM, has(ff,:RM) ? nothing : method.(t)::Array{ResRM,N}),
#         kwds...)
function convolute_GH{T<:TripleAxis,N}(ff::Dispersion,p,t::Array{T,N};
        n::Integer=get(ff,:n,6),
        xv::Array{ResFloat,2}=get(ff,:xv,has(ff,:xv) ? nothing : zeros(ResFloat,3,n^3)),
        qv::Array{Lattice3Vector{ResFloat},1}=get(ff,:qv,has(ff,:qv)?nothing:Array{Lattice3Vector{ResFloat}}(n^3)),
        method::Function=get(ff,:method,popoviciRmM),
        RM::Array{ResRM,N}=get(ff,:RM, has(ff,:RM) ? nothing : method.(t)::Array{ResRM,N}),
        kwds...)
    @assert compatible(t,RM)
    net=length(t) # number of elements t contains
    modes=multiplicity(ff) # this will throw a NullException if ff.SF.multiplicity isn't set!
    rtype=returneltype(ff) # this will throw a NullException if ff.SF.returneltype isn't set!
    # The convolution output is stored in conv, which will later be reshaped to the size of t
    conv=(modes>1)?zeros(rtype,net,modes):zeros(rtype,net)
    # Use the Gauss-Hermite quadrature method to evaluate the convolution. For this we need roots
    # and weights for the 4-directions. Since we use the eigenvectors of M, there's no easy way
    # to single out a direction as being less important -- all are integrated with the same number
    # of points. r and w are either passed in or calculated in the function call and describe
    # the roots and weights for each direction;
    (r,w)=gausshermite(n)
    # we need to evaluate Sum(w₁w₂w₃w₄S(σ₁r₁,σ₂r₂,σ₃r₃,σ₄r₄)) over the 4D grid for each
    # element of the Array{TripleAxis}
    rv=hcat(ndgridvecs(r,r,r)...)' # 4xlength(r) Array{ResFloat,2}
    wv=prod(hcat(ndgridvecs(w,w,w)...),2) # length(r) Array{ResFloat,1}.
    # there's no need to scale up wv to size(SQω()), since .* later will do the equivalent for us
    #
    # alternatively, save time and just evaluate the function along the eigen-directions
    # rv=hcat([r;0r;0r],[0r;r;0r],[0r;0r;r])'
    # wv=prod(hcat([w-1;0w;0w],[0w;w-1;0w],[0w;0w;w-1])+1,2)
    vy=(modes>1)?zeros(returntype,length(rv),modes):zeros(returntype,length(rv))
    RM3=_integrateProject(RM,4)
    σω=1./sqrt.(_getMij(RM,4,4))
    for i=1:net
        # The function SQω needs to be evaluated at the rescaled Hermitian-root positions
        # along the eigen-vectors of M. The rescaling comes from the convolution being
        # ∫∫∫∫dx⃗ S(x⃗) R(x⃗-x⃗ₒ) where R(x⃗)=Rₒ sqrt(|M|)/(2π)²×exp(-½ x⃗'Mx⃗)
        # ie, ∝∫∫∫∫dx⃗ S(x⃗)×exp(-½ [x⃗-x⃗ₒ]'M[x⃗-x⃗ₒ]), while Gauss-Hermite quadrature is an
        # approximation to an integral of the form ∫∫∫∫dy⃗ S(y⃗)×exp(-½y²). If we replace
        # [x⃗-x⃗ₒ] by the eigen-vectors of M, e⃗, then [x⃗-x⃗ₒ]'M[x⃗-x⃗ₒ] becomes
        # e⃗'Me⃗= e₁M₁₁e₁ + e₂M₂₂e₂ + e₃M₃₃e₃ + e₄M₄₄e₄, and
        # exp(-½ [x⃗-x⃗ₒ]'M[x⃗-x⃗ₒ]) → exp(-½M₁₁e₁²)exp(-½M₂₂e₂²)exp(-½M₃₃e₃²)exp(-½M₄₄e₄²).
        # Each exponential is then independent and we can define yᵢ²=Mᵢᵢeᵢ² to enable use
        # of Gauss-Hermite quadradture. yᵢ=sqrt(Mᵢᵢ)eᵢ → dyᵢ=sqrt(Mᵢᵢ)deᵢ; → eᵢ=yᵢ/sqrt(Mᵢᵢ)
        # Gauss-Hermite quadrature states that we only need to evaluate S(y⃗) at the roots
        # of the Hermite polynomial (of rank n) and take an appropriately weighted sum
        # to estimate ∫∫∫∫dy⃗ S(y⃗)×exp(-½y²). The rᵢ calculated above are the roots along
        # each yᵢ-axis; first rescale them to be along ei and multiply by the eigen-vectors
        # of M to go from eᵢ to [x⃗-x⃗ₒ]
        trv=rv./sqrt.(hkleigvals(RM[i])) # broadcasts to scale-up the singleton dimensions to the size of rv
        _multMv!(hkleigvecs(RM[i]),trv,xv)
        # perform the basis change from Popvici's basis to the Sample reciprocal lattice basis,
        # stores the resultant Lattice3Vectors via promotion to E=0 Lattice4Vectors in qv
        lab2sample!(xv,t[i],qv)
        # adding the converted x⃗-x⃗ₒ to Q⃗ₒ and evaluating S(Q⃗)
        QE=get4vector(t[i])+qv
        (rE0,rS,rγ)=ff.SF.S(p,QE) # returns E0(Q), I(Q), gamma(Q)
        unsafevoigt_σγ_nc!(σω[i]*ones(rS),rγ,getE(QE)-rE0,vy);
        (modes>1)?conv[i,:]+=sum(wv.*rS.*vy,1):conv[i]+=sum(wv.*rS.*vy)
    end
    conv/=sqrt(pi*pi*pi) # rescale by the overall convolution normalization
end

function convolute{T<:TripleAxis,N}(ff::Delta,p::Any,t::Array{T,N};
        method::Function=get(ff,:method,popoviciRmM),
        RM::Array{ResRM,N}=get(ff,:RM, has(ff,:RM) ? nothing : method.(t)::Array{ResRM,N}),
        kwds...)
    # S must take one or more LatticeVector objects as input and return (vQω,S,γQω), three lists
    # describing the delta functions. Likely the list for each point will be truncated by the
    # function, based on the description of the sample and/or instrument. If any truncation is
    # performed, this function as written requires each point in t to be given-back the same number
    # of delta functions.
    @assert compatible(t,RM) "The passed $(size(t)) TripleAxis is not compatible with the $(size(RM)) resolution information"
    #net=length(t) # number of elements t contains
    verifymultiplicityandreturneltype!(p,t,ff) # ensure that ff.multiplicity and ff.returneltype are set
    modes=multiplicity(ff) # this will throw a NullException if ff.SF.multiplicity isn't set!
    # We will define points in (Q,E) relative to the initial position(s) of the TripleAxis object(s)
    # To do this effectively, we need the vectors defining the crystal lattice orientation(s) of the TripleAxis object(s)
    xv=get4vector(getx(t[:])); yv=get4vector(gety(t[:])); zv=get4vector(getz(t[:])); # Array{Lattice4Vector,1}
    ev=map(x->Lattice4Vector(0,0,0,1.,x.lat),xv); # create vector of energy unit vectors :/
    tQ=get4vector(t[:]) # we also need the current positions as 4-vectors
    #
    (rvQω,rS,rγQω)=ff.SF.S(p,tQ) # rvQω and rγQω are Lattice4Vectors in the same lattice as xv,yv,zv
    @assert compatible(rvQω,rS,rγQω) "Returned elements from S are incompatibly sized"
    #
    # Bragg widths are given by the diagonal elements of the resolution matrix:
    σs=hcat(map(x->1./sqrt.(diag(x)),_getM.(RM[:]))...)' # voigtσγ takes Gaussian widths
    # _getM.(Array{ResMat,N}) returns Array{Array{ResFloat,2},N}
    # map(x->1./sqrt(diag(x)),_getM.(RM)) returns Array{Array{ResFloat,1},N}
    # hcat([]...)' then gives an Nx4 Array{ResFloat,2} with columns for each σᵢ
    σx=σs[:,1];σy=σs[:,2];σz=σs[:,3];σω=σs[:,4]; # slice in case modes>1
    #
    # repeat tQ,xv,y,zv,σx,σy,σz,σω up to size, if necessary
    if modes>1
        isa(modes,Float64)&&(warn("modes is a float?!");println(typeof(modes)))
        repvec=[1;modes]
        tQ=repeat(tQ,outer=repvec)
        xv=repeat(xv,outer=repvec)
        yv=repeat(yv,outer=repvec)
        zv=repeat(zv,outer=repvec)
        ev=repeat(ev,outer=repvec)
        σx=repeat(σx,outer=repvec)
        σy=repeat(σy,outer=repvec)
        σz=repeat(σz,outer=repvec)
        σω=repeat(σω,outer=repvec)
    end
    dQ=rvQω-tQ # the distances to each delta function
    # voigtσγ takes (Gaussian width σ, Lorentzian width γ, position vector) and is normalized properly
    # unsafevoigt! doesn't verify array sizes before attempting to access them via slicing
    conv=rS
    vv=similar(conv) # setup arrays to hold the voigt profiles
    # unsafevoigtσγnc!(σ,γ,x,y) takes arrays of Gaussian widths (σ), Lorentzian widths (γ),
    # distances from 0 (x), and an array to store the resultant calculation (y), and calculates
    # an approximation to the Voigt function for different regions of (~x/σ, ~σ/γ) space using
    # logical slicing without verifying that the arrays have the same size
    #
    # Array{lattice4vector}.*Array{lattice4vector} is necessary because dot() overloading stopped working ?!?
    unsafevoigt_σγ_nc!(σx,xv.*(rγQω.*xv),dQ.*xv,vv);conv.*=vv
    unsafevoigt_σγ_nc!(σy,yv.*(rγQω.*yv),dQ.*yv,vv);conv.*=vv
    unsafevoigt_σγ_nc!(σz,zv.*(rγQω.*zv),dQ.*zv,vv);conv.*=vv
    unsafevoigt_σγ_nc!(σω,ev.*(rγQω.*ev),dQ.*ev,vv);conv.*=vv
    #
    # the version above contains, e.g., xv*rγQω*xv which might make the lorentzian width have the wrong units if rγQω isn't a LatticeTensor
    return conv # which has size(conv)=(length(t),) or (length(t),modes)
end
