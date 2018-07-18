# Plotting.
# Type aliases for plotting, should  be Union{}'s of equivalent plotting types
const Detector0d =  Union{TASPd,C5d{2},DMCd,EIGERd,FLEXXd,IN22d}
const PolarizedDetector0d = C5d{3}
const Detector1d =  MACSd
#const Detector1d = DMCd(?)
#const Detector2d = RITAd

#
# Take advantage of plot(Measured,Measured) to plot 0d detector data
function plot(p::Detector0d,x=p.x,y=p.y;pfit=true,kwargs...)
    xv=getDat(p,x);yv=getDat(p,y);m=p.mask
    any(m) && (errorbar(xv[m],yv[m];masked=true,kwargs...))
    h=Measureds.errorbar(xv[.!m],yv[.!m];kwargs...)
    pfit&&hasbeenfit(p)&&(plotfit(p,x;kwargs...))
    standardAxisLabels(p,x,y)
    return h
end
function plot(p::PolarizedDetector0d,x=p.x,y=p.y;pfit=true,kwargs...)
    nchan=size(p.data,2) # polarized instruments store different channels in the second data dimension
    (cv,fv,mv,lv,fcv)=createcolorlistkeyed(nchan;color=["k"],fillstyle=["left","top","right","bottom"],kwargs...) # if kwargs include color or fillstyle that is used
    xm=getDat(p,x); ym=getDat(p,y); m=p.mask; pfit=pfit&hasbeenfit(p)
    for i=1:nchan
        any(m) && (errorbar(xm[m,i],ym[m,i];masked=true,kwargs...,color=cv[i],fillstyle=fv[i],marker=mv[i],linestyle=lv[i],facecolor=fcv[i]))
        errorbar(xm[.!m,i],ym[.!m,i];kwargs...,color=cv[i],fillstyle=fv[i],marker=mv[i],linestyle=lv[i],facecolor=fcv[i])
#        pfit&&(plotfit(p,x;kwargs...,color=cv[i],fillstyle=fv[i],marker=mv[i],linestyle=lv[i],facecolor=fcv[i],channel=i))
# TODO wite a plotfit that can handle polarized data (this will also require nlfit that can handle polarized data)
    end
    standardAxisLabels(p,x,y)
end
standardAxisLabels(p,x,y)= (xlabel(plainstring(getCol(p,x)));ylabel(plainstring(getCol(p,y))))
function plot(p::Detector1d,y=p.y;kwargs...)
    info("Plotting of Detector1d objects is not yet supported. There may be some logical approach to creating a default Slice.")
end
function plot(p::Scatterd{1};x=p.x[1],y=p.x[2],v=p.y,showcolorbar=true,k...)
    xv=getDat(p,x); yv=getDat(p,y); vv=getDat(p,v)
    if ndims(vv)==2 && all(size(vv).>(1,1)) # we can plot a pcolormesh
        h=pcolormesh(value.(xv),value.(yv),value.(vv);k...)
        xlabel(plainstring(getCol(p,x)));ylabel(plainstring(getCol(p,y)))
        standardAxisLabels(p,x,y)
        if showcolorbar
            c=colorbar(h)
            c["set_label"](plainstring(getCol(p,v)))
            return (h,c)
        end
    elseif ndims(vv)<3
        xv=vec(xv); yv=vec(yv); vv=vec(vv);
        # likely xv or yv will be constant
        plotx=std(xv)>0
        h=Measureds.errorbar( plotx ? xv : yv, vv; k...)
        #pfit&&hasbeenfit(p)&&(plotfit(p,x;kwargs...))
        standardAxisLabels(p,plotx?x:y,v)
    else
        error("A Scatterd{1} object should not return ndim>2 columns!")
    end
    return h
end

function scatter(p::Detector0d,x::AbstractString,y::AbstractString=p.x;c=p.y,showcolorbar=true,kwargs...)
    cv=getDat(p,c); sp=sortperm(cv)
    xv=getDat(p,x)[sp];yv=getDat(p,y)[sp]; m=p.mask[sp]
    cv=cv[sp]
    #any(m) && (errorbar(xv[m],yv[m];masked=true,kwargs...))
    h=Measureds.scatter(xv[.!m],yv[.!m];c=cv[.!m],kwargs...)
    standardAxisLabels(p,x,y)
    if showcolorbar
        cb=colorbar(h)
        cb["set_label"](plainstring(getcolumns(p)[getIdx(p,c)]))
        return (h,cb)
    else
        return (h,nothing)
    end
end

# it would be usefull to have keyword arguments to pass to evaluate and convolute
# TODO: implement kwargs splitting for convolution keywords and plotting keywords?
function plotfit(s::Detector0d,xname=s.x;smooth=0,smoothcriteria=0.25,smoothspread=3,standarddeviations=1,peaks=0,color="k",kwargs...)
    for (i,(k,v)) in enumerate(kwargs)
        (k==:linestyle)&& ( (v=="-")||0==standarddeviations?kwargs[i]=(k,"--"):kwargs[i]=(k,"-") ) # replace linestyle with "-" unless if it is already
    end
    # when plotting a vector of fit objects we need to remove keywords that we can't use
    bkw=[:facecolor,:marker,:fillstyle]; gkw=trues(kwargs)
    for (i,(k,v)) in enumerate(kwargs)
        gkw[i] = !any(k.==bkw)
    end
    kwargs=kwargs[gkw] # prune the troublesome kwargs
    (color,facecolor)=decodecolor(color,"full")


    ~isa(peaks,Vector)&&(peaks=[peaks])
    ~isa(eltype(peaks),Int)&&(peaks=map(Int,peaks))
    if any(fieldnames(s).==:fitres)&&hasbeenfit(s.fitres) # FIXME temporary bad hack
        fun=s.fitres.fitfun
        p=deepcopy(s.fitres.pM) # typeof(pM)=Vector{Measured} #p=s.fitres.param
        1==standarddeviations||(for i=1:length(p); p[i]=Measured(value(p[i]),standarddeviations^2*variance(p[i])); end)
        peaks=peaks[0.<peaks.<=length(p)] # restrict to valid peaks
        npks=length(peaks)
        if npks>1
            for i=1:npks
                thisp=deepcopy(p)
                thisp[peaks[vcat(1:i-1,i+1:npks)]]=0 # zero the other peaks

                (x,y)=callsmooth(fun,thisp,s,xname;repeat=smooth,criteria=smoothcriteria,spread=smoothspread) # x is the plotting x defined by xname
                plot(x,y;kwargs...)
            end
        end
        # callsmooth evaluates fun and intelligently intercalates points to produce a smooth(er) output, if requested
        (x,y)=callsmooth(fun,p,s,xname;repeat=smooth,criteria=smoothcriteria,spread=smoothspread) # x is the plotting x defined by xname
        plot(x,y;color=color,kwargs...)
    end
    if any(fieldnames(s).==:fits)&&length(s.fits)>0 # FIXME temporary bad hack
        for i=1:length(s.fits)
            if hasbeenfit(s.fits[i])
                fun=s.fits[i].fitfun; p=deepcopy(s.fits[i].pM)
                1==standarddeviations||(for i=1:length(p); p[i]=Measured(value(p[i]),standarddeviations^2*variance(p[i])); end)

                peaks=peaks[0.<peaks.<=length(p)] # restrict to valid peaks
                npks=length(peaks)
                if npks>1
                    for i=1:npks
                        thisp=deepcopy(p)
                        thisp[peaks[vcat(1:i-1,i+1:npks)]]=0 # zero the other peaks

                        (x,y)=callsmooth(fun,thisp,s,xname;repeat=smooth,criteria=smoothcriteria,spread=smoothspread) # x is the plotting x defined by xname
                        plot(x,y;kwargs...)
                    end
                end

                (x,y)=callsmooth(fun,p,s,xname;repeat=smooth,criteria=smoothcriteria,spread=smoothspread)
                plot(x,y;kwargs...)
            end
        end
    end
end

# Simple plotting of vectors of Scatterd objects FIXME this only makes sense for 0d detectors or cuts through Nd detectors
for pfn in (:plot,:plotfit)
    @eval begin
        function $pfn{T<:Scatterd}(v::Array{T},o...;kwargs...)
            (cv,fv,mv,lv,fcv)=createcolorlistkeyed(v;kwargs...)
            h=map((x,y,z,a,b,c)->$pfn(x,o...;kwargs...,color=y,fillstyle=z,marker=a,linestyle=b,facecolor=c),v,cv,fv,mv,lv,fcv)
        end
    end
end

function callsmooth(f::Function,p::Array,xin::AbstractArray;repeat::Integer=1,criteria::Real=0.01,spread::Integer=3)
    # ensure that the input x are sorted (and a normal vector)
    x=sort(collect(xin))
    # first evaluate f at the input x
    y=f(p,x)
    ## here we can loop through multiple times
    for no=1:repeat
        idxm=calcidxmat(x,y,criteria,spread)
        intx=vec(mean(x[idxm],2))
        inty=f(p,intx) # evalute f at the new x positions
        pushvecs!(x,intx,y,inty) # push the new points onto the existing vectors
        (x,y)=sortvecs(x,y)
    end
    return (x,y)
end
function callsmooth(f::Single,p::Array,s::Detector0d,xname;repeat::Integer=1,criteria::Real=0.01,spread::Integer=3)
    x=getVal(s,s.x)
    # first evaluate f at the input x
    y=f(p,x)
    # grab also the plotting x
    px=getVal(s,xname)
    # and sort all three
    (px,x,y)=sortvecs(px,x,y) # sorts the three vectors based on px
    ## here we can loop through multiple times
    for no=1:repeat
        idxm=calcidxmat(px,y,criteria,spread)
        intx=vec(mean(x[idxm],2))
        intp=vec(mean(px[idxm],2))
        inty=f(p,intx) # evalute f at the new x positions
        pushvecs!(x,intx,y,inty,px,intp) # push the new points onto the existing vectors
        (px,x,y)=sortvecs(px,x,y) # sorts the three vectors based on px
    end
    return (px,y)
end
function callsmooth(f::FittingFunction,p::Array,s::Detector0d,xname;repeat::Integer=1,criteria::Real=0.01,spread::Integer=3)
    x=copy(s.instrument) # the "x" values are the instrument configurations stored in s -- copy to avoid modifying
    y=f(p,x) # evaluate f at the input x
    # try to be clever about updating plotx:
    isa(xname,Number) && (xname=s.columns[xname]) # look-up the name in the columns list
    if     xname=="qh" || xname == "h"
        getx=geth
    elseif xname=="qk" || xname == "k"
        getx=getk
    elseif xname=="ql" || xname == "l"
        getx=getl
    elseif xname=="en" || xname == "e"
        getx=getE
    elseif xname[1]=='a' # (character comparison) x should be one of the angles
        col=parse(Int,xname[2])
        getx=x->getangle(x,col)/pi*180 # angles are stored in radian but always(?) used in degrees
    else
        info("Unknown plotting x-column, $xname. You must extend `callsmooth` to enable smoothed plotting for $xname with a `FittingFunction`.")
        getx=x->getDat(s,xname)
        repeat=0
    end
    idxm=hcat(1:length(x)-1,2:length(x))        # double up
    intx=vec( x[idxm[:,1]]/2 + x[idxm[:,2]]/2 ) # on all
    pushvecs!(x,intx,y,f(p,intx))               # points first
    px=getx(x)
    (px,x,y)=sortvecs(px,x,y) # sorts the three vectors based on px
    for no=1:repeat
        idxm=calcidxmat(px,y,criteria,spread)
        intx=vec( x[idxm[:,1]]/2 + x[idxm[:,2]]/2 )
        length(intx)==0&&(status(:hint,"callsmooth does nothing for zero-curvature functions");break)
        pushvecs!(x,intx,y,f(p,intx)) # push the new x and y points onto the existing vectors
        px=getx(x)
        (px,x,y)=sortvecs(px,x,y) # sorts the three vectors based on px
    end
    return (px,y)
end
function pushvecs!{T1,T2}(x::Vector{T1},a::Vector{T2})
    @assert T2<:T1
    for i=1:length(a); push!(x,a[i]); end
end
function pushvecs!{T1,T2,R1,R2}(x::Vector{T1},a::Vector{T2},y::Vector{R1},b::Vector{R2})
    @assert (T2<:T1)&(R2<:R1)&(length(x)==length(y))&(length(a)==length(b))
    for i=1:length(a); push!(x,a[i]); push!(y,b[i]); end
end
function pushvecs!{T1,T2,R1,R2,S1,S2}(x::Vector{T1},a::Vector{T2},y::Vector{R1},b::Vector{R2},z::Vector{S1},c::Vector{S2})
    @assert (T2<:T1)&(R2<:R1)&(S2<:S1)&(length(x)==length(y)==length(z))&(length(a)==length(b)==length(c))
    for i=1:length(a); push!(x,a[i]); push!(y,b[i]); push!(z,c[i]); end
end
"""
    calcidxmat(x,y,criteria)

Calculate the index matrix based on the local curvature of y=f(x) and the passed criteria.
"""
function calcidxmat{T<:Measured}(x::Vector,y::Vector{T},criteria::Real,spread::Integer=3)
    yr = (value(y)-minimum(value(y)))/(maximum(value(y))-minimum(value(y)))
    dx = centerdiff(x)
    dx = ((dx-minimum(dx))/(maximum(dx)-minimum(dx))).^2
    κ0 = dx.*yr.*abs.(localcurvature(x,value(y))) # absolute local curvature
    κm = dx.*yr.*abs.(localcurvature(x,value(y)-negativeuncertainty(y)))
    κp = dx.*yr.*abs.(localcurvature(x,value(y)+positiveuncertainty(y)))
    κ0/=maximum(κ0); κm/=maximum(κm); κp/=maximum(κp)
    i=find((κm.>criteria).|(κ0.>criteria).|(κp.>criteria)) # find points of high curvature in any of the three lines
    m=indexvectortomatrix(i,spread) # build up 2*spread*Nx2 index matrix, without duplicate rows
    #return m[ all(length(x).>=m.>0,2), :]
    m[squeeze(all(length(x).>=m.>0,2),2),:]
end
function calcidxmat{T<:Real}(x::Vector,y::Vector{T},criteria::Real,spread::Integer=3)
    yr = (y-minimum(y))/(maximum(y)-minimum(y))
    dx = centerdiff(x)
    dx = ((dx-minimum(dx))/(maximum(dx)-minimum(dx))).^2
    c = dx.*yr.*abs.(localcurvature(x,y)) # favor moderate curvature with few points over high curvature with many points
    c./=maximum(c)
    i=find(c.>criteria) # find the N indicies of y=f(x) with large local curvature
    m=indexvectortomatrix(i,spread) # build up 2*spread*Nx2 index matrix, without duplicate rows
    #return m[ all(length(x).>=m.>0,2), :] # ensure m is within bounds for x
    m[squeeze(all(length(x).>m.>0,2),2),:]
end
function indexvectortomatrix{T<:Integer}(i::Vector{T},pm::Integer=1)
    n=length(i)
    out=Array{T}(2*pm*n,2)
    for j=1:pm
        out[(1:n)+(2*(j-1)*n  ),:]=hcat(i-j  ,i-j+1)
        out[(1:n)+(2*(j-1)*n+n),:]=hcat(i+j-1,i+j  )
    end
    return unique(out,1)
end

"""
    sortvecs(x,y,[z])

Sort two or three vectors based on the sorting permutation of `x`.
"""
function sortvecs(x::Vector,y::Vector,z::Vector)
    @assert length(x)==length(y)==length(z)
    s=sortperm(x)
    return (x[s],y[s],z[s])
end
function sortvecs(x::Vector,y::Vector)
    @assert length(x)==length(y)
    s=sortperm(x)
    return (x[s],y[s])
end
"""
    centerdiff(y,[x=1:length(y)])

Compute the central finite difference of the vector `y`. If `x` is included in the call,
the central finite difference approximation to the derivative is caluclated, defined as:

    ∂y/∂x[i] = (y[i+1]-y[i-1])/(x[i+1]-x[i-1])

for 1<i<length(y); for the two end points the forward- and backwards- finite difference
derivatives are caluclated

    ∂y/∂x[1] = (y[2]-y[1]/(x[2]-x[1])
    ∂y/∂x[end] = (y[end]-y[end-1])/(x[end]-x[end-1])
"""
function centerdiff(y::Vector,x=collect(1:length(y));l=length(y))
    @assert length(x)==length(y)==l>1
    f=(y[2]-y[1])/(x[2]-x[1])
    b=(y[l]-y[l-1])/(x[l]-x[l-1])
    2==l && (return [f;b])
    [f;(y[3:l]-y[1:l-2])./(x[3:l]-x[1:l-2]);b]
end
function forwarddiff(y::Vector,x=collect(1:length(y));l=length(y))
    @assert length(x)==length(y)==l>1
    [(y[2:l]-y[1:l-1])./(x[2:l]-x[1:l-1]); (y[l]-y[l-1])/(x[l]-x[l-1])]
end
function backwarddiff(y::Vector,x=collect(1:length(y));l=length(y))
    @assert length(x)==length(y)==l>1
    [(y[2]-y[1])/(x[2]-x[1]);(y[2:l]-y[1:l-1])./(x[2:l]-x[1:l-1])]
end

"""
    localcurvature(x,y)

Calculate the local curvature of y=f(x) for the passed x and y vectors using the
central finite differences approximation to the derivative and second derivative.
"""
function localcurvature(x::Vector,y::Vector)
    @assert length(x)==length(y)
    1==length(x) && (return [1]) # not defined for a single point
    ∂ =centerdiff(y,x) #  first derivative
    ∂²=centerdiff(∂,x) # second derivative
    κ = ∂²./(1+∂.^2).^(3/2) # keep sign information in curvature
end
