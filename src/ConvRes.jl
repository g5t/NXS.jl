


# M is only a part of the general 5x5 matrix that defines any conic section in 4D
#      | xx xy xz xE x0 |      |xx xy xz xE|
#      | yx yy yz yE y0 |      |yx yy yz yE|
#  W = | zx zy zz zE z0 |, M = |zx zy zz zE|
#      | Ex Ey Ez EE E0 |      |Ex Ey Ez EE|
#      | x0 y0 y0 E0 -R |
# In general, if det(W) != 0 then the conic section is non-degenerate and the
# det(M) determines the type of conic section that W describes
#   if det(M) < 0, W is a hyperboloid
#   if det(M) = 0, W is a paraboloid
#   if det(M) > 0, W is an ellipsoid (which should always be the case here)
# the remaining elements of W not included in M determine the origin of the conic section,
# and a overall rescaling parameter. Our calculation of M has recentered W such that
# QE4v_0=(0,0,0,0) and pulled out the rescaling parameter as R0 or Rm to give
#      | xx xy xz xE 0 |
#      | yx yy yz yE 0 |
#  W = | zx zy zz zE 0 |
#      | Ex Ey Ez EE 0 |
#      |  0  0  0  0 -1|
# With this formalism, the ellipsoid is described by the equation QE1'*W*QE1=0;
# where QE1 is the momentum-energy-unity 5-vector -- QE1=(Qx,Qy,Qz,E,1) -- or QE'*M*QE=1 where
# QE is the momentum-energy 4-vector -- QE=(Qx,Qy,Qz,E).
# And the elements of M are (at least proportional to) inverse squared-radii of the ellipsoid.
function _calculate_σs{T<:ResRM,N}(RM::Union{Array{T,1},Array{T,N}})
    @assert all(det(RM).>0) "The determiniant of M is not positive, and therefore describes a non-ellipsoid."
    @assert all(trace(RM).>0) "The trace of M must be positive if the resolution ellipsoid is real."
    @assert all(issymmetric(RM)) "Not all resolution matricies in RM are symmetric." # which this approximation requires to be true
    # The following analysis assumes that the z-axis is completely independent in the resolution,
    # which is not strictly true (e.g., if the sample shape is not symmetric about the xy plane)
    # Desribe the matrix by its elements
    # | Mxx Mxy Mxz Mxω |          | Mxx Mxy 0   Mxω |
    # | Mxy Myy Myz Myω |          | Mxy Myy 0   Myω |
    # | Mzx Mzy Mzz Mzω | == M0 ~  | 0   0   Mzz 0   |
    # | Mωx Mωy Mωz Mωω |          | Mωx Mωy 0   Mωω |
    M0xx=_getMij(RM,1,1)
    M0yx=_getMij(RM,2,1); M0yy=_getMij(RM,2,2)
    M0zx=_getMij(RM,3,1); M0zy=_getMij(RM,3,2); M0zz=_getMij(RM,3,3)
    M0ωx=_getMij(RM,4,1); M0ωy=_getMij(RM,4,2); M0ωz=_getMij(RM,4,3); M0ωω=_getMij(RM,4,4);
    # Begin by integrating over z.
    # | M1xx M1xy M1xω |
    # | M1yx M1yy M1yω | == M1
    # | M1ωx M1ωy M1ωω |
    M1xx = M0xx-(M0zx.*M0zx)./M0zz
    M1yx = M0yx-(M0zy.*M0zx)./M0zz; M1yy = M0yy-(M0zy.*M0zy)./M0zz
    M1ωx = M0ωx-(M0ωz.*M0zx)./M0zz; M1ωy = M0ωy-(M0ωz.*M0zy)./M0zz; M1ωω = M0ωω-(M0ωz.*M0ωz)./M0zz
    # Next integrate over energy to get=
    # | M2xx M2xy |
    # | M2yx M2yy | == M2
    M2xx=M1xx-(M1ωx.*M1ωx)./M1ωω
    M2yx=M1yx-(M1ωy.*M1ωx)./M1ωω; M2yy=M1yy-(M1ωy.*M1ωy)./M1ωω
    # Finally, integrate over this matrix to get the single value "M3xx"
    M3xx=M2xx-(M2yx.*M2yx)./M2yy
    # Pull together the Gaussian widths (or FWHM) that describe the resolution ellipsoid
    # 1/sqrt(Mᵢᵢ)=1/sqrt(1/σ²ᵢᵢ)=σᵢᵢ as long as σᵢᵢ>0
    # This is where the indenpendence of z assumption comes into play. If z is not independent, there should be up to three additional widths
    σzz=1./sqrt(M0zz)
    σωω=1./sqrt(M1ωω)
    σωy=-M1ωy./M1ωω./sqrt(M2yy)
    σωx=-(M1ωx./M1ωω - M1ωy./M1ωω.*M2yx./M2yy)./sqrt(M3xx);
    σyy=1./sqrt(M2yy)
    σyx=-M2yx./M2yy./sqrt(M3xx)
    σxx=1./sqrt(M3xx)
    (σxx,σyy,σzz,σωω,σyx,σωx,σωy)
end
# this alternative version is slower (~30%) but, perhaps, more accurate -- and certainly more compact
function _calculate_σs_alternative{T<:ResRM,N}(RM::Union{Array{T,1},Array{T,N}})
    @assert all(det(RM).>0) "The determiniant of M is not positive, and therefore describes a non-ellipsoid."
    @assert all(trace(RM).>0) "The trace of M must be positive if the resolution ellipsoid is real."
    Mxyzω=_getM(RM)                     # (m×n×l×...×)4x4
    Mxyω =_integrateProject(Mxyzω,3)    # (m×n×l×...×)3x3
    Mxy  =_integrateProject(Mxyω, 3)    # (m×n×l×...×)2x2
    Mx   =_integrateProject(Mxy , 2)    # (m×n×l×...×)1x1
    # Pull together the Gaussian widths that describe the resolution ellipsoid
    # the eigen values of M, Mᵢᵢ=1/(σ²ᵢᵢ) so σᵢᵢ=1/sqrt(Mᵢᵢ)
    # _getMij returns the [i,j] elements, the same shape as the outer array M, Array{Matrix{ResFloat},N},
    # e.g. Array{ResFloat,N}. If each M_p is (m×n×l×...×)N×N, then _getMij(M_p,i,j) will be m×n×l×...
    # for future reference: if typeof(A)=Array{ResFloat,2}, create an object of type Array{Array{ResFloat,2},1} via Array{ResFloat,2}[A,A,A,...]
    σxx=1./sqrt(_getMij(Mx   ,1,1))
    σyy=1./sqrt(_getMij(Mxy  ,2,2))
    σzz=1./sqrt(_getMij(Mxyzω,3,3))
    σωω=1./sqrt(_getMij(Mxyω ,3,3))
    σyx=-σxx.*  _getMij(Mxy  ,2,1)./_getMij(Mxy  ,2,2)
    σωy=-σyy.*  _getMij(Mxyω ,3,2)./_getMij(Mxyω ,3,3)
    σωx=-σxx.*( _getMij(Mxyω ,3,1)./_getMij(Mxyω ,3,3) - _getMij(Mxyω ,3,2)./_getMij(Mxyω ,3,3).*_getMij(Mxy  ,2,1)./_getMij(Mxy  ,2,2))
    (σxx,σyy,σzz,σωω,σyx,σωx,σωy)
end
#
convres(S::Function,p::Any,t::TripleAxis,o...;kw...)=convres(S,p,[t];kw...)[1]
function convres{T<:TripleAxis,N,I<:Integer,R<:ResRM,P<:Number}(SQω::Function,p::Array{P},t::Union{Array{T,1},Array{T,N}};
                                             PQω::Function=((p,x)->ones(ResFloat,size(x)...)),
                                             BQω::Function=((p,x)->zeros(ResFloat,size(x)...)),
                                             method::Function=popoviciRmM,
                                             useMonteCarlo::Bool=true,seedMonteCarlo::Integer=-1,
                                             MCpoints::Integer=100,
                                             accuracy::Union{I,Array{I,1}}=10,
                                             RM::Union{Array{R,1},Array{R,N}}=method(t))
    # SQω, PQω, and BQω are function objects (or anonymous functions)
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

    # the size of t might be an K×L×M×N×… Array of TripleAxis object, but is likely just an N-length Vector of TripleAxis objects
    ndt=ndims(t) # number of dimensions of t (guaranteed to be at least one)
    net=prod([size(t)...]) # number of elements t contains

    # calculate the prefactor and resolution matrix, returns an Array of type ResRM,
    # each ResRM has fields :R and :M for R0 (or Rm) and M, respectively
    #RM=method(t) # moved into the function call to allow pre-processing of the resolution matricies
    # the number of dimensions and size of RM should be the same as the size of t
    @assert compatible(t,RM)
    (σxx,σyy,σzz,σωω,σyx,σωx,σωy)=_calculate_σs(RM)
    # check to see how many modes are returned by SQω
    # SQω *should* return an array that is [number of points](x[number of detectors])(x[number of modes])
    # or possibly (NxMxLxKx...)x[number of modes] where N,M,L,K=size(get4vector(t))
    # FIXME the following code likely won't work properly unless M,L,K,...=0 (i.e., only point detectors are likely to work)
    ret=SQω(p,get4vector(t))
    nd=ndims(ret)
    modes=(nd>ndt)?size(ret,nd):1 # OK for N>1-D detectors
    # The convolution output is stored in conv, which will later be reshaped to the size of t
    conv=(modes>1)?zeros(typeof(ret[1]),net,modes):zeros(typeof(ret[1]),net)
    # We will define points in (Q,E) relative to the initial position(s) of the TripleAxis object(s)
    # To do this effectively, we need the vectors defining the crystal lattice orientation(s) of the TripleAxis object(s)
    xv=getx(t); yv=gety(t); zv=getz(t); # Array{Lattice3Vector,N} depending on ndims(t)
    if useMonteCarlo
        mcr = MCpoints; m=1; max_mcr = 3*10^7 # tm of (3e7 x 4) in Float64s is ~915 MB
        (mcr>max_mcr)&&(m=cld(mcr,max_mcr); mcr=cld(mcr,m))
        # calculate the overall convolution normalization
        convn=1/(m*mcr*pi^4) # for any MonteCarlo integration, it's 1/(number of points) here the importance sampling adds a factor of 1/pi^dimensionality to normalize the probability distribution function
        # setup array for MonteCarlo random numbers (four arrays mcr-long should be faster than one mcr×4 array because it avoids slicing in the for loop)
        tx=zeros(ResFloat,mcr)
        ty=zeros(ResFloat,mcr)
        tz=zeros(ResFloat,mcr)
        tω=zeros(ResFloat,mcr)
        # allow for predictability to enable reproducable results (say, for taking numerical derivatives of a function)
        rng=signbit(seedMonteCarlo)?RandomDevice():MersenneTwister(seedMonteCarlo) # create a random-number-generator object just for this function, to avoid reseeding here effecting code elsewhere
        for i=1:net # loop over the number of elements in t
        for j=1:m # in case the MonteCarlo step has been subdivided
            # get random points distributed evenly between (-pi/2,pi/2) and take their tangent
            rand!(rng,tx)
            rand!(rng,ty)
            rand!(rng,tz)
            rand!(rng,tω)
            tx=tan.((tx-0.5)*pi)
            ty=tan.((ty-0.5)*pi)
            tz=tan.((tz-0.5)*pi)
            tω=tan.((tω-0.5)*pi)
            # this gives a distribution in tx,ty,tz,tw that is strongly peaked at 0
            # determine the normalization for each (tx,ty,tw,tz) point
            norm=exp.(-(tx.^2+ty.^2+tz.^2+tω.^2)/2)#.*(1+tx.^2).*(1+ty.^2).*(1+tz.^2).*(1+tω.^2)
            (modes>1)&&(norm=repeat(norm,outer=[1,modes]))
            # Get the QE Lattice4Vector of t[i], calculate the integration grid as an Array of Lattice4Vectors centered on QE
            tv=Lattice4Vector(xv[i]*σxx[i]*tx + yv[i]*(σyy[i]*ty + σyx[i]*tx) + zv[i]*σzz[i]*tz ,σωx[i]*tx + σωy[i]*ty + σωω[i]*tω);
            ret=SQω(p,get4vector(t[i])+tv).*norm # calculate S(Q,ω) and multiply by normalization
            (modes>1)?conv[i,:]+=sum(ret,1):conv[i]+=sum(ret) # perform the integration, possibly in multiple MonteCarlo steps
        end
        end
    else
        (isa(accuracy,ResFloat))&&(accuracy=[accuracy])  # Make accuracy a Vector
        (length(accuracy)==1)&&(accuracy=[accuracy;accuracy;0;accuracy;]) # with four elements
        (length(accuracy)==2)&&(accuracy=[accuracy[1],accuracy[1],accuracy[2],accuracy[1]])
        (length(accuracy)==3)&&(accuracy=[accuracy[1],accuracy[2],0,accuracy[3]])
        # As the code is written now, there are 6+2*modes+1Lattice4Vector+1TripleAxis vectors all of this size at any one time.
        #
        # On a 64bit system a Vector{Lattice4Vector} might be length(Vector{Lattice4Vector})x8bytes long if it is a vector of pointers
        # which Julia might decide to create if it is advantageous (many copies of the same Lattice4Vector in the Vector{Lattice4Vector}?)
        # On my system, Julia chose to make a 9e8-element Vector{Lattice4Vector{Float64}} (created by push!(ing) the same L4V onto an array) out of pointers
        # which then "occupied" ~763 MB.
        # Creating a 9e8-element Vector{Lattice4Vector{Float32}} in the same way promted Julia to instead create a contiguous array totalling 6.7 GB.
        # So calculating overhead isn't so straight forward.

        # Create integration vectors, which are equally spaced points in 4D angle space (-pi/2,pi/2)
        (tx,ty,tw,tz)=ndgrid(_acc2vec(accuracy[1]),_acc2vec(accuracy[2]),_acc2vec(accuracy[3]),_acc2vec(accuracy[4]))
        # and calculate the overall convolution normalization
        convn=prod(_acc2step(accuracy[1:4]))/pi^4 # the sum of the importance sampling distribution function is normalized by multiplying by the step-size (in the angle) and dividing by pi for each dimension
        any(accuracy[1:4].==0)&&(convn*=(0.79788)^sum(accuracy[1:4].==0))
        # take the tangent of each angle to give a strong weighting of the grid near (0,0,0,0), switch from 4D arrays to vectors
        tx=tan.(tx[:]); ty=tan.(ty[:]); tw=tan.(tw[:]); tz=tan.(tz[:])
        # determine the normalization for each (tx,ty,tw,tz) point
        norm=exp.(-0.5*(tx.^2+ty.^2+tw.^2+tz.^2)).*(1+tx.^2).*(1+ty.^2).*(1+tw.^2).*(1+tz.^2)
        (modes>1)&&(norm=repeat(norm,outer=[1,modes]))
        for i=1:net # find a way to enforce parallelization here?
            # Get the QE Lattice4Vector of t[i], calculate the integration grid as an Array of Lattice4Vectors centered on QE
            QE= get4vector(t[i]) + Lattice4Vector(xv[i]*σxx[i]*tx + yv[i]*(σyy[i]*ty + σyx[i]*tx) + zv[i]*σzz[i]*tz, σωx[i]*tx + σωy[i]*ty + σωω[i]*tw)
            ret=SQω(p,QE).*norm # calculate S(Q,ω) and multiply by normalization
            (modes>1)?conv[i,:]=sum(ret,1):conv[i]=sum(ret) # perform the integration
        end
    end
    ndp=ndims(PQω(p,t))
    ((ndp<nd)&(modes==1))&&(warn("PQω returned an Array with $ndp dimensions \n while one with $nd was expected."))
    pref=PQω(p,t[:]) # convert t to a Vector before calculating PQω
    ndp=ndims(pref)
    ((ndp<2)&(modes>1))&&(pref=repeat(pref,outer=[1,cld(modes,size(pref,2))])) # cld(x,y) == ceil(Integer,x/y)
    conv.*=(modes>1)?pref[:,1:modes]:pref # Multiply by the slow prefactor
    (modes>1)&&(conv=sum(conv,2))          # sum over all modes

    conv=reshape(conv,size(t)) # reshape conv to be the same shape as the input TripleAxis Array
    conv.*=_getR(RM)    # rescale conv by R0 (or Rm)
    conv./=sqrt.(det(RM))# divide by the square root of |M|, det is overloaded to return an array the same shape as RM
    conv.*=convn        # rescale by the overall convolution normalization
    conv+=BQω(p,t)      # add on the slow background
end
