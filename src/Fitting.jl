"""
    checkcorrelations(fixed,correlated,scaled,centered)

Checks to ensure that no more than one of each set of correlated values is free,
i.e., that only one of their corresponding entries in `fixed` is `false`. The function
ensures that all correlated values are set to have the same centering and scaling, as well.

The correlation matrix is specified instead as the vector `correlated` with integer
values of `±n` for elements that are correlated, where the `+n` elements act as the
primary values and the `-n` elements are chained to be the same.
`n` represents values `1,2,...,number_correlated`.

As an example, for a system with six parameters representing the area, position, and width
of two Gaussian peaks where the two peak widths should be constrained to be equal and the
positions are enforced to be symmetric about 16, the four input vectors
might have the values

    fixed=[false false false false false false]
    correlated=[0 0 1 0 0 -1]
    scaled=[1 1 1 1 -1 1]
    centered=[0 16 0 0 16 0]

which would result in the returned values

    fixed=[false false false false false true]
    correlated=[0 0 1 0 0 -1]
    scaled=[1 1 1 1 -1 1]
    centered=[0 16 0 0 16 0]

The scaling and centering is subsequently handled by `unscaleandceter` and `rescaleandcenter`
as necessary.
"""
function checkcorrelations(fixed::BitArray,correlated::Vector{I},scaled::Vector{R},centered::Vector{S}) where {I<:Integer,R<:Real,S<:Real}
    @assert compatible(fixed,correlated,scaled,centered)
    for i=1:maximum(correlated) # 0 or n_max
        aci=abs.(correlated).==i
        ci=correlated.==i
        if sum(aci)>1 && sum(ci)>0
            if sum(ci)>1 # if more than one +n exists, fix all but the first.
                ci[findfirst(ci)]=false
                correlated[ci]=-i
            end
            nci=correlated.==-i
            pci=correlated.==i # there should be only one now
            fixed[nci]=true # ensure that correlated positions with negative indicies are fixed
            centered[nci]=centered[pci][1] # set all correlated values to have the same centering
#            scaled[nci]=scaled[pci][1] # set all correlated values to have the same scaling
        else
            warn("A correlation indicator was skipped, used only once, or used only for fixed (-$i) positions; indicator $i has been ignored")
        end
    end
    (fixed,correlated,scaled,centered)
end
checkcorrelations(f::BitArray,correlated::Vector{T},s::Vector,c::Vector) where {T<:Real} = checkcorrelations(f,Int.(correlated),s,c)

"""
    unscaleandcenter(parameters,scaled,centered)

    unscaleandcenter(parameters,fixed,correlated,scaled,centered)

In either form, `unscaleandceter` removes first the centering and then the scaling from the
input parameters and returns the results:

    centered_scaled_parameters = (parameters-centered)./scaled

The function `rescaleandcenter` does the opposite.
"""
unscaleandcenter(z::Vector{T},s::Vector{R},c::Vector{S}) where {T<:Real,R<:Real,S<:Real}=(z-c)./s
unscaleandcenter(z::Vector{T},f::BitArray,r::Vector{I},s::Vector{R},c::Vector{S}) where {T<:Real,I<:Integer,R<:Real,S<:Real}=unscaleandcenter(z,s,c)
"""
    rescaleandcenter(centered_scaled_parameters,scaled,centered)

The opposite of `unscaleandcenter`, applies scaling and centering to return the `parameters`

    parameters= centered_scaled_parameters .* scaled + centered
"""
rescaleandcenter(z::Vector{T},s::Vector{R},c::Vector{S}) where {T<:Real,R<:Real,S<:Real}=z.*s+c
"""
    rescaleandcenter(reduced_centered_scaled_parameters,fixed,correlated,scaled,centered,original_parameters)

In the case that some of the `original_parameters` were fixed and/or correlated to one another,
`rescaleandcenter` takes the `reduced_centered_scaled_parameters`, e.g., from the output of fitting the parameters
of a function to some data, and first expands them to the full size of `original_parameters`
and then utilizes rescaleandcenter! to remove the scaling and centering.
"""
function rescaleandcenter(z::Vector{T},f0::BitArray,r::Vector{I},s::Vector{R},c::Vector{S},z0::Vector{T})  where {T<:Real,I<:Integer,R<:Real,S<:Real}
    z1=similar(z0)
    z1[f0]=z0[f0]
    rescaleandcenter!(z,f0,r,s,c,z1)
    z1
end
"""
    rescaleandcenter!(centered_scaled_parameters,fixed,correlated,scaled,centered,parameters)

Elements of `parameters` which were not kept constant, according to `fixed`, are replaced
by the elements of `centered_scaled_parameters` in place.
Then elements of `parameters` with `-n` entries in `correlated` are replaced by their
corresponding (primary) values, those with `+n` as their entry in `correlated`.
Finally, the scaling is applied to `parameters` and the centering is added.
"""
function rescaleandcenter!(z::Vector{T},f0::BitArray,r::Vector{I},s::Vector{R},c::Vector{S},z1::Vector{T}) where {T<:Real,I<:Integer,R<:Real,S<:Real}
    tmp=copy(z1)
    tmp[.!f0]=z
    tmp-=c # remove the centering here to avoid adding it twice
    for i=1:maximum(r)
        nri=r.==-i
        pri=r.==i
        if sum(pri)==1 && sum(nri)>0
            tmp[nri]=tmp[pri][1]
        end
    end
    tmp.*=s
    tmp+=c
    z1[:]=tmp[:]
end
"""
    splitLMkeys(keyvals)

Check for the presence of `LsqFit.levenberg_marquardt` specific keyword value pairs
in `keyvals`. If any are found, return them as `LMkv` and all remaining as `others`
in the output tuple `(LMkv,others)`
"""
function splitLMkeys(keyvals)
    # check for (and remove) levenberg_marquardt specific keyworks in kwds
    LMkv=[]; others=[]
    LMkeys=[:tolX,:tolG,:maxIter,:lambda,:show_trace]
    for (k,v) in keyvals
        any(k.==LMkeys)?push!(LMkv,(k,v)):push!(others,(k,v))
    end
    (LMkv,others)
end
"""
    parseOptimResult(results,f,g,p0,fcsc)

Parse the result of an Optim.jl fitting to pull out information which we want to keep in a `FitResult`.
`result` should be a `Optim.OptimizationResults`.
"""
function parseOptimResult{F<:Function,G<:Function}(results::OptimBase.OptimizationResults,f::F,g::G,p0::Vector,w,fcsc)
    @assert OptimBase.summary(results)=="Levenberg-Marquardt"
    if !OptimBase.converged(results)
        warn("The iteration limit was reached in $(OptimBase.iterations(results)) steps.")
    end
    p = OptimBase.minimizer(results)
    fitresidual = f(p)      # the weighted difference between the fit function and the data
    residual=fitresidual.*w # the raw difference between fit function and the data
    fitjacobian = g(p)      # the weighted Jacobian for each data point
    jacobian=fitjacobian.*w # the raw Jacobian for each data point
    dof = length(residual) - length(p)
    fit=FitResult(dof, rescaleandcenter(p,fcsc...,p0), ones(Measured,size(p)...), fcsc..., residual, fitresidual, jacobian)
    # dumb hack to include errors in FitResult
    return FitResult(fit.dof,fit.p,errors1sigma(fit),fit.fixed,fit.correlated,fit.scaled,fit.centered,fit.residual,fit.fitresidual,fit.jacobian)
end
function parseOptimResult{F<:Function,G<:Function}(results,f::F,g::G,p0::Vector,w,fcsc,ff::FittingFunction) # hack to include fitting function
    tf=parseOptimResult(results,f,g,p0,w,fcsc)
    FitResult(tf.dof,tf.p,tf.pM,tf.fixed,tf.correlated,tf.scaled,tf.centered,tf.residual,tf.fitresidual,tf.jacobian,ff)
end

# fitting a scan "in-place" should store the fitting results in the scan
nlfit!{F<:Function}(S::F,p0::Vector,s::Scatterd;k...)=(s.fitres=nlfit(S,p0,s;k...))
nlfit!(S::FittingFunction,p0::Vector,s::Scatterd;k...)=(s.fitres=nlfit(S,p0,s;k...))

"""
    nlfit(convfunc,p0,tripleaxis,y,w;convoluted,fixed,correlated,scaled,centered,keywords...)

For a resolution convoluted function `convfunc`, use `LsqFit.jl`'s Levenberg-Marquardt
algorithm to find the minimum of

     χ²(p) = sum(abs2, (convfunc(p,tripleaxis) - y)/w )

by modifying `p`; which may start as a reduced form of `p0` depending on the values of
the optional keyword arguments `fixed`, `correlated`, `scaled`, and `centered`)
See `checkcorrelations`, `rescaleandcenter`, and `unscaleandceter` for details on the behaviour of
these four keywords.
The optional keyword `convoluted` can potentailly speed-up the fitting routine by not
performing the numerical convolution while evaluating ∂χ²/∂pᵢ if `convoluted[i]==false`.
"""
function nlfit{T<:TripleAxis,R<:Number}(F::Convoluted,p0::Vector,t::Array{T,1},y::Array{R,1},w::Array{R,1}=ones(y);
                                    convoluted::BitArray=trues(length(p0)),fixed::BitArray=falses(length(p0)),
                                    correlated::Vector=zeros(Int,length(p0)),scaled::Vector=ones(p0),centered::Vector=zeros(p0),kwds...)
    @assert compatible(t,y)
    @assert compatible(p0,convoluted)
    fcsc=checkcorrelations(fixed,correlated,scaled,centered)
    pv=unscaleandcenter(p0,fcsc...) # remove scaling and centering from initial p vector
    (LMkwds,kwds)=splitLMkeys(kwds)
    f=_makefunction(pv,t,y,w,F,convoluted,fcsc...;kwds...)
    g=_makejacobian(pv,t,y,w,F,convoluted,fcsc...;kwds...)
    results = LsqFit.levenberg_marquardt(f, g, p0[.!fixed]; LMkwds...)
    return parseOptimResult(results,f,g,p0,w,fcsc,F)
end
"""
    nlfit(barefunc,p0,tripleaxis,y,w;convoluted,fixed,correlated,scaled,centered,keywords...)

For a non-resolution-convoluted function `barefunc`, use `Optim.jl`'s Levenberg-Marquardt
algorithm to find the minimum of

     χ²(p) = sum(abs2, (barefunc(p,tripleaxis) - y)/w )

by modifying `p`; which may start as a reduced form of `p0` depending on the values of
the optional keyword arguments `fixed`, `correlated`, `scaled`, and `centered`)
See `checkcorrelations`, `rescaleandcenter`, and `unscaleandceter` for details on the behaviour of
these four keywords.
"""
function nlfit{T<:TripleAxis,R<:Number}(F::FittingFunction,p0::Vector,t::Array{T,1},y::Array{R,1},w::Array{R,1}=ones(y);
        fixed::BitArray=falses(length(p0)),correlated::Vector=zeros(Int,length(p0)),scaled::Vector=ones(p0),centered::Vector=zeros(p0),kwds...)
    @assert compatible(t,y)
    fcsc=checkcorrelations(fixed,correlated,scaled,centered)
    pv=unscaleandcenter(p0,fcsc...) # remove scaling and centering from initial p vector
    (LMkwds,kwds)=splitLMkeys(kwds)
    f=_makefunction(pv,t,y,w,F,fcsc...;kwds...)
    g=_makejacobian(pv,t,y,w,F,fcsc...;kwds...)
    results = LsqFit.levenberg_marquardt(f, g, p0[.!fixed]; LMkwds...)
    return parseOptimResult(results,f,g,p0,w,fcsc,F)
end


function scaled_weights_from_error{T}(δ::Array{T})
    δ=abs.(δ)
    δ0=δ.==0
    all(δ0)&&(return ones(δ)) # all equal weights if no uncertainty
    mδ=minimum(δ[.~δ0])
    w=δ/mδ # relative uncertanties >=1
    w[δ0]=1 # if any have zero uncertainty, treat them as equivalent to the next-most precise
    return w
#    # an alterate approach:
#    δi=δ.>=typemax(T)
#    mδ=maximum(δ[.~δi]) # maximum non-infinite uncertainty
#    w=δ/mδ # relative uncertanties 0<=w<=1
#    w[δi]=1 # infinite uncertainty gets minimum precision, 1
#    w[δ0]=minimum(w[~δ0]) # zero uncertainty gets maximum precision (that is non-zero)
#    return w
end
function scaled_weights_from_weights{T}(win::Array{T})
    win=abs.(win)
    win0=win.==0
    all(win0)&&(return ones(win))
    mwin=minimum(win[.~win0])
    w=win/mwin
    w[win0]=1
    return w
end

"""
    nlfit(func,p0,scatterd;keywords...)
Fit the function `func` to the default y-values in `scatterd` at the instrument conditions stored
in `scatterd.instrument` by varying the parameters `p0`. Observations in `scatterd` with corresponding
`scaterd.mask==false` are ignored. The minimization function is internally defined as
`f(p)= sum(abs2,(func(p,scatterd)-getVal(scatterd.y))/getErr(scatterd.y))` -- that is, raw χ². (where real χ² is this divided by the degrees of freedom)

If `func` is a standard `Function` then `nlfit` accepts an optional keyword boolean argument
`fitinstrument` which is `true` by default.
If `fitinstrument` is `true`, `nlfit` behaves as indicated above. Otherwise, the evaluation of
`func` is modified to `func(p, getVal(scatterd.x) )`.
"""
function nlfit(S::FittingFunction,p0::Vector,s::Scatterd;k...)
    m=.!s.mask # only include non-masked points
    nlfit(S,p0,s.instrument[m],getVal(s,s.y)[m],scaled_weights_from_error(getErr(s,s.y)[m]);k...)
end
function nlfit{F<:Function}(S::F,p0::Vector,s::Scatterd;fitinstrument::Bool=true,k...)
    x=fitinstrument? s.instrument : getVal(s,s.x)
    m=.!s.mask # only include non-masked points
    nlfit(S,p0,x[m],getVal(s,s.y)[m],scaled_weights_from_error(getErr(s,s.y)[m]);k...)
end
function nlfit{T<:Scatterd,F<:Union{FittingFunction,Function}}(S::F,p0::Vector,s::Array{T,1};k...)
    t=s[1].instrument
    y=getVal(s[1],s[1].y)
    w=getErr(s[1],s[1].y)
    m=.!s[1].mask
    for i=2:length(s)
        t=vcat(t,s[i].instrument)
        y=vcat(y,getVal(s[i],s[i].y))
        w=vcat(w,getErr(s[i],s[i].y))
        m=vcat(m,.!s[i].mask)
    end
    nlfit(S,p0,t[m],y[m],scaled_weights_from_error(w[m]);k...)
end
function nlfit{T<:TripleAxis,R<:Number,S<:Function,P<:Function,B<:Function}(SQω::S,p0::Vector,t::Array{T,1},y::Array{R,1},w::Array{R,1}=ones(y);
                                    PQω::P=scale_one,BQω::B=background_zero,fixed::BitArray=falses(length(p0)),
                                    correlated::Vector=zeros(Int,length(p0)),scaled::Vector=ones(p0),centered::Vector=zeros(p0),kwds...)
    @assert compatible(t,y)
    fcsc=checkcorrelations(fixed,correlated,scaled,centered)
    pv=unscaleandcenter(p0,fcsc...) # remove scaling and centering from initial p vector
    (LMkwds,kwds)=splitLMkeys(kwds)
    fitfun=Bare(SQω;P=PQω,B=BQω,kwds...)
    f=_makefunction(pv,t,y,w,fitfun,fcsc...;kwds...)
    g=_makejacobian(pv,t,y,w,fitfun,fcsc...;kwds...)
    results = LsqFit.levenberg_marquardt(f, g, p0[.!fixed]; LMkwds...)
    return parseOptimResult(results,f,g,p0,w,fcsc,fitfun)
end
"""
    nlfit(func,p0,x,y,w;fixed,correlated,scaled,centered,keywords...)

For a regular `Function` `func`, use `Optim.jl`'s Levenberg-Marquardt
algorithm to find the minimum of

     χ²(p) = sum(abs2, (func(p,x) - y)/w )

by modifying `p`; which may start as a reduced form of `p0` depending on the values of
the optional keyword arguments `fixed`, `correlated`, `scaled`, and `centered`)
See `checkcorrelations`, `rescaleandcenter`, and `unscaleandcenter` for details on the behaviour of
these four keywords.

`y` and `w` must be `Vector`s of the same length. `x` will typically be be a `Vector` the same
length as `y` but can be any size as long as the output of `func` is a `Vector` the same length as `y`.
It is possible to create specialized versions of `nlfit` for different dimensionalities of `y` and `x`,
but the possible cases are likely too numerous to make this worthwhile. Instead one can usually
write `func` as a function of multiple dependent variables and then let `x` be a `Matrix`
comprised of column vectors of the dependent variables and to vectorize `y` as well.
"""
function nlfit{T<:Number,R<:Number,F<:Function}(S::F,p0::Vector,x::AbstractArray{T},y::AbstractArray{R,1},w::AbstractArray{R,1}=ones(y);
                                    fixed::BitArray=falses(length(p0)),correlated::Vector=zeros(Int,length(p0)),scaled::Vector=ones(p0),centered::Vector=zeros(p0),kwds...)
    @assert compatible(w,y)
    fcsc=checkcorrelations(fixed,correlated,scaled,centered)
    pv=unscaleandcenter(p0,fcsc...) # remove scaling and centering from initial p vector
    (LMkwds,kwds)=splitLMkeys(kwds)
    f=_makefunction(pv,x,y,w,S,fcsc...;kwds...)
    g=_makejacobian(pv,x,y,w,S,fcsc...;kwds...)
    results = LsqFit.levenberg_marquardt(f, g, p0[.!fixed]; LMkwds...)
    fitfun=Single(S)
    return parseOptimResult(results,f,g,p0,w,fcsc,fitfun)
end


function _checkforRM(ta,kwds)
    # look to see if the resolution Matrix and prefactor have been passed as part of the kwds
    RMfound=false; for kwd in kwds; (kwd[1]==:RM)&&(RMfound=true;resmat=kwd[2]); end
    if !RMfound # look for a specification of the resolution calculation method
        resmeth=popoviciRmM # default to the Popovici approximation normalized to the source flux
        for kwd in kdws; (kdw[1]==:resolutionmethod)&&(resmeth=kwd[2]); end
        resmat=resmeth.(ta)
        push!(kwds,(:RM,resmat)) # do the resolution calculation now if it hasn't been done already
    end
end


# macro and use of centralrule modified from Calculus.jl
macro centralrule(x, e)
    x, e = esc(x), esc(e)
    :($e = cbrt(eps(eltype($x))) * max(one(eltype($x)), abs($x)))
end

function _makeminimizationfunction{T<:TripleAxis,R<:Number,N}(xfull::Vector,ta::Array{T,N},obs::Array{R,N},weights::Array{R,N},
        FF::Convoluted,conved::BitArray,fix::BitArray,corr::Vector,scale::Vector,center::Vector;kwds...)
    F=deepcopy(FF) # avoid issues with adding/replacing keywords
    !has(F,:RM)&&(add(F,:RM,get(F,:method,popoviciRmM).(ta))) # assure that RM is pre-calculated
    MCflag=useMonteCarlo(F) # MonteCarlo is used by by default
    MCflag && has(F,:rng) || add( F,:rng,MersenneTwister( get(F,:seed,rand(UInt32)) ) )
    xfd=similar(xfull); dx=zero(eltype(xfull))
    # now create the function which gets returned. Variables in f(x,grad) other than x and grad are pulled from the wrapper function.
    function f(x::Vector,grad::Union{Void,Vector})
        rescaleandcenter!(x,fix,corr,scale,center,xfull)
        c=convolute(F,xfull,ta) # if nmodes>1, this is size(ta)xModes, otherwise size(ta)
        if !(grad === nothing) && length(grad)>0
            gradptr=pointer(grad)
            MCflag&&(seed=rand(UInt32)) # for the derivatives
            xfd[:]=xfull[:]
            for (j,i) in enumerate(find(.!fix))
                @centralrule xfd[i] dx
                oldx = xfd[i]
                xfd[i]=oldx+dx
                if !conved[i] # this parameter doesn't need to be convoluted
                    f_xplus=applyPB(F,xfd,ta,c)
                    xfd[i]=oldx-dx
                    f_xminus=applyPB(F,xfd,ta,c)
                else
                    # if MonteCarlo is used, reset the random number generator with the same seed for each derivative
                    MCflag&&(replace!(F,:rng,MersenneTwister(seed)))
                    f_xplus=evaluate(S,xfd,ta)
                    xfd[i]=oldx-dx
                    MCflag&&(replace!(F,:rng,MersenneTwister(seed)))
                    f_xminus=evaluate(S,xfd,ta)
                end
                xfd[i]=oldx
                grad[j]= sum((f_xplus-f_xminus)./weights/2/dx)
            end
            @assert gradptr==pointer(grad) "The pointer to grad has been modified.\n"*
            "Only assign elements of grad by grad[i]=new_value or grad[:]=new_values syntax"
        end
        sum(abs2,(applyPB(F,xfull,ta,c)-obs)./weights) # sum(abs2,)=sum(abs2())=sum(abs().^2)
    end
end
_reverseinputs{F<:Function}(f::F)=(reversef(o...)=f(reverse(o)...))
# creates the minimization function f(x⃗,∇x⃗) suitable for use with NLopt.jl
minimizationfunctionNLopt(opt...;kwds...)=_makeminimizationfunction(opt...;kwds...)
# creates the minimization function f(∇x⃗,x⃗) suitable for use with Optim.jl/LsqFit.jl
minimizationfunctionOptim(opt...;kwds...)=_reversinputs(_makeminimizationfunction(opt...;kwds...))

# For use with LsqFit.levenberg_marquardt, following LsqFit.jl (which used a general function from Calculus.jl to create the Jacobian function)

# MetaConvoluted <: FittingFunction is *like* Convoluted but might modify ta as well,
# so it needs to be treated separately :(
function _makejacobian{T<:TripleAxis,R<:Number}(xfull::Vector,ta::Vector{T},obs::Vector{R},weights::Vector{R},
        MFF::MetaConvoluted,conved::BitArray,fix::BitArray,corr::Vector,scale::Vector,center::Vector)
    F=deepcopy(MFF.C)
    # !has(F,:RM)&&(add(F,:RM,get(F,:method,popoviciRmM).(ta))) # assure that RM is pre-calculated # FIXME we can't do this with a MetaConvoluted
    MCflag=useMonteCarlo(F) # MonteCarlo is used by by default
    MCflag && has(F,:rng) || add( F,:rng,MersenneTwister( get(F,:seed,rand(UInt32)) ) )
    xfull=deepcopy(xfull); xfd=similar(xfull); dx=zero(eltype(xfull))
    # now create the function which gets returned.
    function g(x::Vector)
        rescaleandcenter!(x,fix,corr,scale,center,xfull)
        J=zeros(eltype(obs),length(ta),sum(.!fix))
        modta=MFF.M(xfull,ta) # make our modifications to the TripleAxis object(s)
        any(.!conved)&&(c=convolute(F,xfull,modta)) # only evaluate if some parameters aren't convoluted.
        MCflag&&(seed=rand(UInt32))
        xfd[:]=xfull[:]
        for (j,i) in enumerate(find(.!fix))
            @centralrule xfd[i] dx
            oldx = xfd[i]
            xfd[i]=oldx+dx
            if !conved[i] # this parameter doesn't need to be convoluted
                f_xplus=applyPB(F,xfd,modta,c)
                xfd[i]=oldx-dx
                f_xminus=applyPB(F,xfd,modta,c)
            else
                # if MonteCarlo is used, reset the random number generator with the same seed for each derivative
                MCflag&&(replace!(F,:rng,MersenneTwister(seed)))
                f_xplus=evaluate(F,xfd,modta)
                xfd[i]=oldx-dx
                MCflag&&(replace!(F,:rng,MersenneTwister(seed)))
                f_xminus=evaluate(F,xfd,modta)
            end
            xfd[i]=oldx
            J[:,j]= (f_xplus-f_xminus)./weights/2/dx
        end
        return J
    end
end
# No need to define _makefunction just for MetaConvoluted, the default FittingFunction version works fine.

# Convoluted <: FittingFunction
function _makejacobian{T<:TripleAxis,R<:Number}(xfull::Vector,ta::Array{T,1},obs::Array{R,1},weights::Array{R,1},
        FF::Convoluted,conved::BitArray,fix::BitArray,corr::Vector,scale::Vector,center::Vector)
    F=deepcopy(FF) # avoid issues with adding/replacing keywords
    !has(F,:RM)&&(add(F,:RM,get(F,:method,popoviciRmM).(ta))) # assure that RM is pre-calculated
    MCflag=useMonteCarlo(F) # MonteCarlo is used by by default
    MCflag && has(F,:rng) || add( F,:rng,MersenneTwister( get(F,:seed,rand(UInt32)) ) )
    xfull=deepcopy(xfull); xfd=similar(xfull); dx=zero(eltype(xfull))
    # now create the function which gets returned.
    function g(x::Vector)
        rescaleandcenter!(x,fix,corr,scale,center,xfull)
        J=zeros(eltype(obs),length(ta),sum(.!fix))
        any(.!conved)&&(c=convolute(F,xfull,ta)) # only evaluate if some parameters aren't convoluted.
        MCflag&&(seed=rand(UInt32))
        xfd[:]=xfull[:]
        for (j,i) in enumerate(find(.!fix))
            @centralrule xfd[i] dx
            oldx = xfd[i]
            xfd[i]=oldx+dx
            if !conved[i] # this parameter doesn't need to be convoluted
                f_xplus=applyPB(F,xfd,ta,c)
                xfd[i]=oldx-dx
                f_xminus=applyPB(F,xfd,ta,c)
            else
                # if MonteCarlo is used, reset the random number generator with the same seed for each derivative
                MCflag&&(replace!(F,:rng,MersenneTwister(seed)))
                f_xplus=evaluate(F,xfd,ta)
                xfd[i]=oldx-dx
                MCflag&&(replace!(F,:rng,MersenneTwister(seed)))
                f_xminus=evaluate(F,xfd,ta)
            end
            xfd[i]=oldx
            J[:,j]= (f_xplus-f_xminus)./weights/2/dx
        end
        return J
    end
end
function _makefunction{T<:TripleAxis,R<:Number,N}(xfull::Vector,ta::Array{T,N},obs::Array{R,N},weights::Array{R,N},
        FF::Convoluted,conved::BitArray,fix::BitArray,corr::Vector,scale::Vector,center::Vector)
    F=deepcopy(FF)
    !has(F,:RM)&&(add(F,:RM,get(F,:method,popoviciRmM).(ta)))
    xfull=deepcopy(xfull); xfd=similar(xfull); dx=zero(eltype(xfull))
    # now create the function which gets returned.
    function f(x::Vector)
        rescaleandcenter!(x,fix,corr,scale,center,xfull)
        (evaluate(F,xfull,ta)-obs)./weights
    end
end
# Non-convoluted FittingFunctions (and MetaConvoluted, paying the price for using evaluate for all parameters)
function _makejacobian{T<:TripleAxis,R<:Number}(xfull::Vector,ta::Array{T,1},obs::Array{R,1},weights::Array{R,1},
        S::FittingFunction,fix::BitArray,corr::Vector,scale::Vector,center::Vector)
    xfull=deepcopy(xfull); xfd=similar(xfull); dx=zero(eltype(xfull))
    # now create the function which gets returned.
    function g(x::Vector)
        rescaleandcenter!(x,fix,corr,scale,center,xfull)
        J=zeros(eltype(obs),length(ta),sum(.!fix))
        xfd[:]=xfull[:]
        for (j,i) in enumerate(find(.!fix))
            @centralrule xfd[i] dx
            oldx = xfd[i]
            xfd[i]=oldx+dx
            f_xplus=evaluate(S,xfd,ta)
            xfd[i]=oldx-dx
            f_xminus=evaluate(S,xfd,ta)
            xfd[i]=oldx
            J[:,j]= (f_xplus-f_xminus)./weights/2/dx
        end
        return J
    end
end
function _makefunction{T<:TripleAxis,R<:Number,N}(xfull::Vector,ta::Array{T,N},obs::Array{R,N},weights::Array{R,N},
        S::FittingFunction,fix::BitArray,corr::Vector,scale::Vector,center::Vector)
    xfull=deepcopy(xfull); xfd=similar(xfull); dx=zero(eltype(xfull))
    # now create the function which gets returned.
    function f(x::Vector)
        rescaleandcenter!(x,fix,corr,scale,center,xfull)
        (evaluate(S,xfull,ta)-obs)./weights
    end
end
# And for non-SQEFunctions
#function _makejacobian{T<:Number,R<:Number}(xfull::Vector,ta::AbstractArray{T,1},obs::AbstractArray{R,1},weights::AbstractArray{R,1},
#  S::Function,fix::BitArray,corr::Vector,scale::Vector,center::Vector;kwds...)
#    xfull=deepcopy(xfull); xfd=similar(xfull); dx=zero(eltype(xfull))
#    # now create the function which gets returned.
#    function g(x::Vector)
#        rescaleandcenter!(x,fix,corr,scale,center,xfull)
#        J=zeros(eltype(obs),length(ta),sum(.!fix))
#        xfd[:]=xfull[:]
#        for (j,i) in enumerate(find(.!fix))
#            @centralrule xfd[i] dx
#            oldx = xfd[i]
#            xfd[i]=oldx+dx
#            f_xplus=S(xfd,ta;kwds...)
#            xfd[i]=oldx-dx
#            f_xminus=S(xfd,ta;kwds...)
#            xfd[i]=oldx
#            J[:,j]= (f_xplus-f_xminus)./weights/2/dx
#        end
#        return J
#    end
#end
#function _makefunction{T<:Number,R<:Number,N}(xfull::Vector,ta::AbstractArray{T,N},obs::AbstractArray{R,N},weights::AbstractArray{R,N},
#  S::Function,fix::BitArray,corr::Vector,scale::Vector,center::Vector;kwds...)
#    xfull=deepcopy(xfull); xfd=similar(xfull); dx=zero(eltype(xfull))
#    # now create the function which gets returned.
#    function f(x::Vector)
#        rescaleandcenter!(x,fix,corr,scale,center,xfull)
#        (S(xfull,ta;kwds...)-obs)./weights
#    end
#end
function _makejacobian{T<:Number,R<:Number,F<:Function}(xfull::Vector,ta::AbstractArray{T},obs::AbstractArray{R,1},weights::AbstractArray{R,1},
  S::F,fix::BitArray,corr::Vector,scale::Vector,center::Vector;kwds...)
    xfull=deepcopy(xfull); xfd=similar(xfull); dx=zero(eltype(xfull))
    # now create the function which gets returned.
    function g(x::Vector)
        rescaleandcenter!(x,fix,corr,scale,center,xfull)
        J=zeros(eltype(obs),length(obs),sum(.!fix))
        xfd[:]=xfull[:]
        for (j,i) in enumerate(find(.!fix))
            @centralrule xfd[i] dx
            oldx = xfd[i]
            xfd[i]=oldx+dx
            f_xplus=S(xfd,ta;kwds...)
            xfd[i]=oldx-dx
            f_xminus=S(xfd,ta;kwds...)
            xfd[i]=oldx
            J[:,j]= (f_xplus-f_xminus)./weights/2/dx
        end
        return J
    end
end
function _makefunction{T<:Number,R<:Number,F<:Function,N}(xfull::Vector,ta::AbstractArray{T},obs::AbstractArray{R,N},weights::AbstractArray{R,N},
  S::F,fix::BitArray,corr::Vector,scale::Vector,center::Vector;kwds...)
    xfull=deepcopy(xfull); xfd=similar(xfull); dx=zero(eltype(xfull))
    # now create the function which gets returned.
    function f(x::Vector)
        rescaleandcenter!(x,fix,corr,scale,center,xfull)
        (S(xfull,ta;kwds...)-obs)./weights
    end
end
