# FIXME the following code likely won't work properly unless M,L,K,...=0 (i.e., only point detectors are likely to work)
#_SQωmultiplicity{T<:TripleAxis,N}(p,t::Array{T,N},SQω::Function)=(ret=SQω(p,get4vector(t));nd=ndims(ret); nd>N?size(ret,nd):1)
_SQωmultiplicity{T<:TripleAxis}(p,t::Array{T},SQω::Function)=_SQωmultiplicity(p,t[1],SQω)
_SQωmultiplicity(p,t::TripleAxis,SQω::Function)=length(SQω(p,get4vector(t))) # assume that array output from a single TripleAxis object represents multiple modes.
_SQωreturntype{T<:TripleAxis}(p,t::Array{T},SQω::Function)=_SQωreturntype(p,t[1],SQω)
_SQωreturntype(p,t::TripleAxis,SQω::Function)=eltype(SQω(p,get4vector(t)))

function convolute{T<:TripleAxis,N}(p,t::Array{T,N},SQω::Function;
                                    PQω::Function=((p,x)->ones(ResFloat,size(x)...)),BQω::Function=((p,x)->zeros(ResFloat,size(x)...)),
                                    modes::Integer=_SQωmultiplicity(p,t[1],SQω),returntype::DataType=_SQωreturntype(p,t[1],SQω),
                                    useMonteCarlo::Bool=true,seed::Integer=-1,rng::AbstractRNG=signbit(seed)?MersenneTwister():MersenneTwister(seed),
                                    #n::Integer=5, points::Integer=n^4,
                                    n::Integer=151, points::Integer=4n,
                                    rv::Array{ResFloat,2}=zeros(ResFloat,4,points),
                                    xv::Array{ResFloat,2}=zeros(ResFloat,4,points),
                                    qv::Array{Lattice4Vector{ResFloat},1}=Array{Lattice4Vector{ResFloat}}(points),
                                    resolutionmethod::Function=popoviciRmM,RM::Array{ResRM,N}=resolutionmethod(t)::Array{ResRM,N})
    # SQω takes one or more LatticeVector objects as input and returns S(Q,ω) for the described point(s)
    # PQω takes one or more TripleAxis objects as input and returns scaling prefactors either size(t) or (size(t)..., #modes) large
    # BQω takes one or more TripleAxis objects as input and returns background intensity for t that is size(t) large
    #
    # Neiter PQω nor BQω are convoluted with the resolution.
    # PQω is best used if the same phenomenon has been measured in two experiments that can not be easily added together
    # (e.g., on two different intstruments); then SQω can contain any phenomenon-related scaling parameter(s) and
    # PQω can contain inter-experiment scaling parameters.
    # BQω is best used to describe any background that varies smoothly when compared to the resolution.
    # PQω defaults to the anonymous function x->ones(ResFloat,size(x)...)
    # BQω defaults to the anonymous function x->zeros(ResFloat,size(x)...)
    @assert compatible(t,RM)
    conv=useMonteCarlo
        ?convoluteMonteCarlo(p,t,SQω,modes=modes,returntype=returntype,points=points,rng=rng,rv=rv,xv=xv,qv=qv,RM=RM)
      :convoluteGaussHermite(p,t,SQω,modes=modes,returntype=returntype,points=points,rng=rng,rv=rv,xv=xv,qv=qv,RM=RM)
    conv=applynonconvolutedparts(p,t,conv;PQω=PQω,BQω=BQω,RM=RM)
end
function evaluate{T<:TripleAxis,N}(p,t::Array{T,N},SQω::Function;
                                    PQω::Function=((p,x)->ones(ResFloat,size(x)...)),BQω::Function=((p,x)->zeros(ResFloat,size(x)...)),
                                    resolutionmethod::Function=popoviciRmM,RM::Array{ResRM,N}=resolutionmethod(t)::Array{ResRM,N},kwds...)
# It may be useful to fit a function to data without performing a resolution convolution, perhaps as part of a global parameter search.
# To enable this, evaluate(p,t,S;PQω=P,BQω=B,RM=R,...) evaluates S, P, and B in the same way as convolute(...) but without the convolution.
    conv=SQω(p,t[:])
    conv=applynonconvolutedparts(p,t,conv;PQω=PQω,BQω=BQω,RM=RM)
end
function applynonconvolutedparts{T<:TripleAxis,N}(p,t::Array{T,N},conv;
                            PQω::Function=((p,x)->ones(ResFloat,size(x)...)),BQω::Function=((p,x)->zeros(ResFloat,size(x)...)),
                            resolutionmethod::Function=popoviciRmM,RM::Array{ResRM,N}=resolutionmethod(t)::Array{ResRM,N})
    @assert ndims(conv)<3 "conv should only have size (N,) or (N,M)"
    conv.*=PQω(p,t[:]) # Multiply by the slow prefactor (the .* broadcasting multiplication scales-up PQω as necessary)
    (ndims(conv)==2)&&(conv=sum(conv,2)) # if necessary, sum over all modes; conv should only have size (N,) or (N,M)
    conv=reshape(conv,size(t)) # reshape conv to be the same shape as the input TripleAxis Array
    conv.*=_getR(RM)    # rescale conv by R0 (or Rm)
    conv+=BQω(p,t)      # add on the slow background
end
function convoluteMonteCarlo{T<:TripleAxis,N}(p,t::Array{T,N},SQω::Function;
                                    modes::Integer=_SQωmultiplicity(p,t[1],SQω),returntype::DataType=_SQωreturntype(p,t[1],SQω),
                                    seed::Integer=-1,points::Integer=100,rng::AbstractRNG=signbit(seed)?MersenneTwister():MersenneTwister(seed),
                                    rv::Array{ResFloat,2}=zeros(ResFloat,4,points),
                                    xv::Array{ResFloat,2}=zeros(ResFloat,4,points),
                                    qv::Array{Lattice4Vector{ResFloat},1}=Array{Lattice4Vector{ResFloat}}(points),
                                    resolutionmethod::Function=popoviciRmM,RM::Array{ResRM,N}=resolutionmethod(t)::Array{ResRM,N})
    # the number of dimensions and size of RM should be the same as the size of t
    @assert compatible(t,RM)
    # the size of t might be an K×L×M×N×… Array of TripleAxis object, but is likely just an N-length Vector of TripleAxis objects
    net=length(t) # number of elements t contains
    # The convolution output is stored in conv, which will later be reshaped to the size of t
    conv=(modes>1)?zeros(returntype,net,modes):zeros(returntype,net)
    for i=1:net # loop over the number of elements in t
        # get random points distributed following a normal distribution with σ=1, centered on zero
        randn!(rng,rv);
        # Rescale the random values by eigen-Gaussian-widths and move from the eigen-coordinate system to Popovici's (x,y,z,E) frame, this is x⃗-x⃗ₒ.
        # The eigen values of M are  Mᵢ=1/(σ²ᵢ) so σᵢ=1/sqrt(Mᵢ). Multiplying the random vectors by [σ₁,σ₂,σ₃,σ₄] (elementwise) changes the
        # standard deviation of the distribution appropriately. The eigenbasis of the matrix M is defined by M=VSV⁻¹; where S is
        # a diagonal matrix with elements Mᵢ, and V is a square matrix comprised of column vectors of the eigenvectors of M, vᵢ.
        # A vector, x⃗, can be transformed from the basis of M (the Popovici basis) to the eigenbasis by  v⃗=V⁻¹x⃗; and so a vector
        # in the eigenbasis, v⃗, can be transformed to the basis of M by x⃗=Vv⃗
        rv./=sqrt(eigvals(RM[i])) # uses broadcast to scale-up the singleton dimensions of 1/sqrt(eigvals(RM[i])) to the size of rv
        _multMv!(eigvecs(RM[i]),rv,xv) # _multMv!(M,v) performs v=M*v for the n column vectors in v
        popovici2reciprocal!(rv,t[i],qv) # perform the basis change from Popvici's basis to the Sample reciprocal lattice basis, stores the resultant Lattice4Vectors in qv
        Sx⃗ = SQω(p,get4vector(t[i])+qv)
        # a general MonteCarlo integral is the sum of the function [here S(x⃗)R(x⃗-x⃗ₒ)] evaluated at N randomly chosen points divided by N
        # here we've chosed random points according to R(x⃗-x⃗ₒ), so we only need to take the sum of S(x⃗)
        (modes>1)?conv[i,:]+=sum(Sx⃗,1):conv[i]+=sum(Sx⃗)
    end
    conv./=points       # rescale by the overall convolution normalization (1/#points for MonteCarlo)
end
function convoluteGaussHermite{T<:TripleAxis,N}(p,t::Array{T,N},SQω::Function;
                                    modes::Integer=_SQωmultiplicity(p,t[1],SQω),returntype::DataType=_SQωreturntype(p,t[1],SQω),
                                    n::Integer=151,
                                    r::Array{ResFloat,1}=gausshermiteroots(n),
                                    w::Array{ResFloat,1}=gausshermiteweights(n,r),
                                    xv::Array{ResFloat,2}=zeros(ResFloat,4,4n),
                                    qv::Array{Lattice4Vector{ResFloat},1}=Array{Lattice4Vector{ResFloat}}(4n),
                                    resolutionmethod::Function=popoviciRmM,RM::Array{ResRM,N}=resolutionmethod(t)::Array{ResRM,N})
    @assert compatible(t,RM)
    # the size of t might be an K×L×M×N×… Array of TripleAxis object, but is likely just an N-length Vector of TripleAxis objects
    net=length(t) # number of elements t contains
    # The convolution output is stored in conv, which will later be reshaped to the size of t
    conv=(modes>1)?zeros(returntype,net,modes):zeros(returntype,net)
    # Use the Gauss-Hermite quadrature method to evaluate the convolution. For this we need roots and weights for the 4-directions.
    # Since we use the eigenvectors of M, there's no easy way to single out a direction as being less important -- all are integrated with the same number of points
    # r and w are either passed in or calculated in the function call and describe the roots and weights for each direction;
    # we need to evaluate Sum(w₁w₂w₃w₄S(σ₁r₁,σ₂r₂,σ₃r₃,σ₄r₄)) over the 4D grid for each element of the Array{TripleAxis}
    #rv=hcat(ndgridvecs(r,r,r,r)...)' # 4xlength(r) Array{ResFloat,2}
    #wv=prod(hcat(ndgridvecs(w,w,w,w)...),2) # length(r) Array{ResFloat,1}. there's no need to scale this up to size(SQω()), since .* later will do the equivalent for us
    rv=hcat([r;0r;0r;0r],[0r;r;0r;0r],[0r;0r;r;0r],[0r;0r;0r;r])'
    wv=prod(hcat([w-1;0w;0w;0w],[0w;w-1;0w;0w],[0w;0w;w-1;0w],[0w;0w;0w;w-1])+1,2)
    for i=1:net
        # The function SQω needs to be evaluated at the rescaled Hermitian-root positions along the eigen-vectors of M.
        # The rescaling comes from the convolution being ∫∫∫∫dx⃗ S(x⃗) R(x⃗-x⃗ₒ) where R(x⃗)=Rₒ sqrt(|M|)/(2π)²×exp(-½ x⃗'Mx⃗)
        # ie, ∝∫∫∫∫dx⃗ S(x⃗)×exp(-½ [x⃗-x⃗ₒ]'M[x⃗-x⃗ₒ]), while Gauss-Hermite quadrature is an approximation to an integral of
        # the form ∫∫∫∫dy⃗ S(y⃗)×exp(-½y²). If we replace [x⃗-x⃗ₒ] by the eigen-vectors of M, e⃗,
        # then [x⃗-x⃗ₒ]'M[x⃗-x⃗ₒ] becomes e⃗'Me⃗= e₁M₁₁e₁ + e₂M₂₂e₂ + e₃M₃₃e₃ + e₄M₄₄e₄,
        # and exp(-½ [x⃗-x⃗ₒ]'M[x⃗-x⃗ₒ]) → exp(-½M₁₁e₁²)exp(-½M₂₂e₂²)exp(-½M₃₃e₃²)exp(-½M₄₄e₄²).
        # Each exponential is then independent and we can define yᵢ²=Mᵢᵢeᵢ² to enable use of Gauss-Hermite quadradture.
        # yᵢ=sqrt(Mᵢᵢ)eᵢ → dyᵢ=sqrt(Mᵢᵢ)deᵢ; → eᵢ=yᵢ/sqrt(Mᵢᵢ)
        # Gauss-Hermite quadrature states that we only need to evaluate S(y⃗) at the roots of the Hermite polynomial (of rank n)
        # and take an appropriately weighted sum to estimate ∫∫∫∫dy⃗ S(y⃗)×exp(-½y²).
        # the rᵢ calculated above are the roots along each yᵢ-axis;
        #first rescale them to be along ei and multiply by the eigen-vectors of M to go from eᵢ to [x⃗-x⃗ₒ]
        rv./=sqrt(eigvals(RM[i])) # uses broadcast to scale-up the singleton dimensions of 1/sqrt(eigvals(RM[i])) to the size of rv
        _multMv!(eigvecs(RM[i]),rv,xv)
        # convert from Popovici's coordinate system to the sample reciprocal space, adding the converted x⃗-x⃗ₒ to Q⃗ₒ and evaluating S(Q⃗)
        popovici2reciprocal!(xv,t[i],qv) # perform the basis change from Popvici's basis to the Sample reciprocal lattice basis, stores the resultant Lattice4Vectors in qv
        Sx⃗ = wv.*SQω(p,get4vector(t[i])+qv)
        (modes>1)?conv[i,:]=sum(Sx⃗,1):conv[i]=sum(Sx⃗) # perform the Gauss-Hermite sum
    end
    conv/=pi^2          # rescale by the overall convolution normalization
end
