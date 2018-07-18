export chi2
"""
`FitResult` is an immutable type to contain output from nonlinear least squares fitting.
Fields of each object are:

|  Fieldname   | Description                                                    |
|:------------:|:---------------------------------------------------------------|
|     `dof`    | `Integer` number of degrees of freedom of the fit              |
|     `p`      | `Vector` of output fitting parameters                          |
|     `pM`     | `Vector{Measured}`, output fitting parameters                  |
|              | with uncertainties (1σ by default)                             |
|   `fixed`    | `BitArray`, true for parameters that were fixed during fitting |
| `correlated` | `Vector{Integer}` indicating correlations                      |
|              | between parameters enforced during fitting                     |
|   `scaled`   | `Vector` of relative scaling for correlated parameters         |
|  `centered`  | `Vector` of relative center-point for correlated parameters    |
|  `residual`  | `Vector` of residual at `p`                                    |
|  `jacobian`  | `Matrix` of the Jacobian at `p`                                |
|   `fitfun`   | `FittingFunction` used in the fitting                          |
"""
immutable FitResult{T,R,S,M<:Measured,I<:Integer}
    dof::Integer
    p::Vector{T}
    pM::Vector{M}
    fixed::BitArray
    correlated::Vector{I}
    scaled::Vector{R}
    centered::Vector{S}
    residual::Vector{T}
    fitresidual::Vector{T}
    jacobian::Matrix{T}
    fitfun::FittingFunction
    FitResult{T,R,S,M,I}() where {T,R,S,M<:Measured,I<:Integer}=new() # allow for empty initialization calls to create dummy FitResults in, e.g., TASPd initializaiton
    FitResult{R,L,M,S,U}(a::Integer,b::Vector{R},c::Vector{S},
                     d::BitArray,f::Vector{U},g::Vector{L},
                     h::Vector{M},i::Vector{R},fi::Vector{R},
                     j::Matrix{R},k::FittingFunction) where {R,L,M,S<:Measured,U<:Integer}=new(a,b,c,d,f,g,h,i,fi,j,deepcopy(k))
end
FitResult{T,R,S,M<:Measured,I<:Integer}(a::Integer,b::Vector{T},c::Vector{M},d::BitArray,
          f::Vector{I},g::Vector{R},h::Vector{S},i::Vector{T},fi::Vector{T},j::Matrix{T},
          k::FittingFunction)=FitResult{T,R,S,M,I}(a,b,c,d,f,g,h,i,fi,j,k)
#FitResult{T,M<:Measured,I<:Integer}(a::Integer,b::Vector{T},c::Vector{Measured},d::BitArray,
#          f::Vector{I},g::Vector{T},h::Vector{T},i::Vector{T},j::Matrix{T},
#          k::UInt8,l::Function,m::Function,n::Function)=FitResult{T,M,I}(a,b,c,d,f,g,h,i,j,FittingFunction(k,l,m,n))
#FitResult{T,M<:Measured,I<:Integer}(a::Integer,b::Vector{T},c::Vector{Measured},d::BitArray,
#          f::Vector{I},g::Vector{T},h::Vector{T},i::Vector{T},j::Matrix{T},
#          k::UInt8,l::Function)=FitResult{T,M,I}(a,b,c,d,f,g,h,i,j,FittingFunction(k,l))
FitResult{T,R,S,M<:Measured,I<:Integer}(a::Integer,b::Vector{T},c::Vector{M},d::BitArray,
          f::Vector{I},g::Vector{R},h::Vector{S},i::Vector{T},fi::Vector{T},j::Matrix{T})=FitResult{T,R,S,M,I}(
              a,b,c,d,f,g,h,i,fi,j,Bare((p,x)->ones(size(x)...))) # part of a stupid hack

Base.copy(a::FitResult)=a # since immutable

function showFitResult(io::IO,f::FitResult,compact::Bool=false)
    if hasbeenfit(f)
        compact?(Base.showcompact(io,f.pM)):(Base.show(io,f.pM))
    else
        compact?(Base.print(io,"#undef")):(Base.println(io,"#undef"))
    end
end
Base.show(io::IO,f::FitResult)=showFitResult(io,f,false)
Base.showcompact(io::IO,f::FitResult)=showFitResult(io,f,true)

function (+)(a::FitResult,b::Number)
    @assert hasbeenfit(a);
    FitResult(a.dof,a.p,a.pM,a.fixed,a.correlated,a.scaled,a.centered,
              a.residual,a.fitresidual,a.jacobian,a.fitfun+b)
end
function (-)(a::FitResult,b::Number)
    @assert hasbeenfit(a)
    FitResult(a.dof,a.p,a.pM,a.fixed,a.correlated,a.scaled,a.centered,
              a.residual,a.fitresidual,a.jacobian,a.fitfun-b)
end
function (*)(a::FitResult,b::Number)
    @assert hasbeenfit(a)
    FitResult(a.dof,a.p,a.pM,a.fixed,a.correlated,a.scaled,a.centered,
              a.residual*b,a.fitresidual,a.jacobian,a.fitfun*b)
end
function (/)(a::FitResult,b::Number)
    @assert hasbeenfit(a)
    FitResult(a.dof,a.p,a.pM,a.fixed,a.correlated,a.scaled,a.centered,
              a.residual/b,a.fitresidual,a.jacobian,a.fitfun/b)
end
function (-)(a::FitResult)
    FitResult(a.dof,a.p,a.pM,a.fixed,a.correlated,a.scaled,a.centered,
              -a.residual,a.fitresidual,a.jacobian,-a.fitfun)
end
(+)(a::Number,b::FitResult)=b+a
(-)(a::Number,b::FitResult)=(-b)+a
(*)(a::Number,b::FitResult)=b*a

for z in (:+,:-,:*,:/)
    @eval $z{T<:FitResult}(v::Array{T},b::Number)=broadcast($z,v,b)
    @eval $z{T<:FitResult}(a::Number,v::Array{T})=broadcast($z,a,v)
end

"""
Computes covariance matrix of fit parameters using QR decomposition
"""
function estimate_covar(fit::FitResult)
# compute the covariance matrix from the QR decomposition
Q,R = qr(fit.jacobian)
!isInvertable(R)&&(warn("Upper triangular matrix R was not invertable, expect unphysical error estimates."))
Rinv=isInvertable(R)?inv(R):R
covar = Rinv*Rinv'*chi2(fit)
end

"""
    estimate_errors(fit,alpha=0.68269)
Computes error estimates from `FitResult` `fit` and confidence interval `alpha`.
With alpha left at its default value of 0.68269, 1σ errors are calculated.
"""
function estimate_errors(fit::FitResult, alpha=0.68269)
covar = estimate_covar(fit)
# then the standard errors are given by the sqrt of the diagonal
std_error = sqrt.(diag(covar))
# scale by quantile of the student-t distribution
std_error *= quantile(Distributions.TDist(fit.dof), alpha)
# correlated parameters should have correlated errors:
#Measured(fit.p,sqrt(abs(rescaleandcenter(std_error.^2,fit.fixed,fit.correlated,fit.scaled,0*fit.centered,0*fit.p)))) ### Measured_asym.jl
Measured(fit.p,abs.(rescaleandcenter(std_error.^2,fit.fixed,fit.correlated,fit.scaled,0*fit.centered,0*fit.p)))       ### Measured.jl
end

errors1sigma(fit::FitResult)=estimate_errors(fit,0.68269)
errors2sigma(fit::FitResult)=estimate_errors(fit,0.95449)
errors3sigma(fit::FitResult)=estimate_errors(fit,0.99730)
