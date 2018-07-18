sliceplotdimension(s::Slice)=sum(nobins(s.spec).>1) # how many bin specifications have more than one bin

function setupplot(s::Slice)
    nbins=nobins(s.spec) # the number of sliced bins in all slicing directions
    plotable=find(nbins.>1) # we need more than one bin to use a slice direction as a plot axis
    singlets=(find(nbins.==1)...) # this *should* be !plotable, but since we're going straight to indicies, it needs to be calculated separately
    nplotable=length(plotable)
    ndd=ndims(s.data)
    colnames=name.(getcolumns(s))
    return (plotable,singlets,ndd,colnames)
end
function _plot1d(plotfun::Function,s::Slice,plotable::Vector{I},singlets::Tuple,ndd::Integer,colnames::Vector{T};
              v::AbstractString=s.source.y,x::AbstractString="",vmin::Real=-Inf,vmax::Real=Inf,kwargs...) where {I<:Integer,T<:AbstractString}
    @assert 1==length(plotable)
    xbs=s.spec[plotable[1]]
    isempty(x) && (x=xbs.col)
    xidx=findfirst(colnames.==x)
    yidx=findfirst(colnames.==v)
    xval=squeeze(slicedim(s.data,ndd,xidx),singlets) # slicedim now removes integer indexed dimensions (ndd only)
    yval=squeeze(slicedim(s.data,ndd,yidx),singlets) # we still need to squeeze-out the singleton bins for plotting
    p=sortperm(xval)
    xval=xval[p]
    yval=yval[p]
    m=isfinite.(yval)
    #any(!,m) && (errorbar(xval[.!m],yval[.!m];masked=true,kwargs...))
    h=plotfun(xval[m],yval[m];kwargs...)
    yextrm=value.(extrema(yval[m])) # value to ensure we only pass floats to python/matplotlib
    vmin==-Inf && (vmin=yextrm[1])
    vmax==+Inf && (vmax=yextrm[2])
    setp(gca(),ylim=[vmin,vmax])
    xlabel(plainstring(s.columns[xidx]));ylabel(plainstring(s.columns[yidx]))
    return(h)
end
function _plot2d(plotfun::Function,centersoredges::Symbol,s::Slice,plotable::Vector{I},singlets::Tuple,ndd::Integer,colnames::Vector{T};
              v::AbstractString=s.source.y,x::AbstractString="",y::AbstractString="",vmin::Real=-Inf,vmax::Real=Inf,linewidth=0.1,
              cmap=cmap_inferno,norm::Nullable{Function}=Nullable{Function}(),showcolorbar=true,denancycles::Int=1,kwargs...) where {I<:Integer,T<:AbstractString}
    @assert 2==length(plotable)
    if isa(s.source.x,Array) && length(s.source.x)>1
        isempty(x) && (x=s.source.x[1])
        isempty(y) && (y=s.source.x[2])
    end
    xbs=s.spec[plotable[1]]
    ybs=s.spec[plotable[2]]
    vecx=false
    vecy=false
    isempty(x) && (x=xbs.col)
    isempty(y) && (y=ybs.col)
    xidx=findfirst(colnames.==x) # also used below for labeling
    yidx=findfirst(colnames.==y)
    if x==xbs.col # we can use the binned edges directly
      xval= :edges==centersoredges ? getedges(xbs) : getcenters(xbs)
      vecx=true
    else # we need to do some extra work
      xcen=value.(squeeze(slicedim(s.data,ndd,xidx),singlets)) # (N,M) centers
      denangrid!(xcen,denancycles)
      xval= :edges==centersoredges ? grid_centers_to_corners( xcen ) : xcen # (N+1,M+1) corners
    end
    if y==ybs.col
      yval= :edges==centersoredges ? getedges(ybs) : getcenters(ybs)
      vecy=true
    else
      ycen=value.(squeeze(slicedim(s.data,ndd,yidx),singlets))
      denangrid!(ycen,denancycles)
      yval= :edges==centersoredges ? grid_centers_to_corners( ycen ) : ycen
    end
    xmin,xmax=extrema(xval[isfinite.(xval)])
    xrng=[xmin,xmax]+[-1,1]*(xmax-xmin)/100;
    ymin,ymax=extrema(yval[isfinite.(yval)])
    yrng=[ymin,ymax]+[-1,1]*(ymax-ymin)/100
    if vecx&&vecy
        (xval,yval)=ndgrid(xval,yval)
    elseif vecx
        xval=repeatouter(xval,[1,size(yval,2)])
    elseif vecy
        yval=repeatouter(yval,[1,size(xval,1)])'
    end
    zidx=findfirst(colnames.==v)
    zval=value.(squeeze(slicedim(s.data,ndd,zidx),singlets))
    @assert ndims(xval)==ndims(yval)==ndims(zval)==2 "Something went wrong with squeezing/slicing. x: $(size(xval)), y: $(size(yval)), z: $(size(zval))"
    isnull(norm) || (zval=get(norm).(zval))
    plotint=zval[isfinite.(zval)] # use isfinite instead of .!isnan to remove +/-Inf too
    if !isempty(plotint)
        # vmin,vmax=extrema(intrange)
        vmin==-Inf && (vmin=minimum(plotint))
        vmax==+Inf && (vmax=maximum(plotint))
        # now xval and yval are (N+1,M+1) while zval is (N,M); so we can use pcolormesh
        h=plotfun(xval,yval,zval;vmin=vmin,vmax=vmax,cmap=cmap,linewidth=linewidth,kwargs...)
        xlabel(plainstring(s.columns[xidx]));ylabel(plainstring(s.columns[yidx]))
        setp(gca(),"xlim",xrng,"ylim",yrng)
        return showcolorbar ? colorbar(h,s,v) : h
    end
end
function _plot3d(plotfun::Function,centersoredges::Symbol,s::Slice,plotable::Vector{I},singlets::Tuple,ndd::Integer,colnames::Vector{T};
              v::AbstractString=s.source.y,x::AbstractString="",y::AbstractString="",z::AbstractString="",vmin::Real=-Inf,vmax::Real=Inf,independentscales=false,
	      linewidth=0.1,
              cmap=cmap_inferno,norm::Nullable{Function}=Nullable{Function}(),showcolorbar=true,denancycles::Int=1,kwargs...) where {I<:Integer,T<:AbstractString}
    @assert 3==length(plotable)
    if isa(s.source.x,Array) && length(s.source.x)>1
        isempty(x) && (x=s.source.x[1])
        isempty(y) && (y=s.source.x[2])
    end
    xbs=s.spec[plotable[1]]
    ybs=s.spec[plotable[2]]
    zbs=s.spec[plotable[3]]
    vecx=false
    vecy=false
    vecz=false
    isempty(x) && (x=xbs.col)
    isempty(y) && (y=ybs.col)
    isempty(z) && (z=zbs.col)
    xidx=findfirst(colnames.==x) # also used below for labeling
    yidx=findfirst(colnames.==y)
    zidx=findfirst(colnames.==z)
    if x==xbs.col # we can use the binned edges directly
      xval= :edges==centersoredges ? getedges(xbs) : getcenters(xbs)
      vecx=true
      xmin,xmax=extrema(xval)
    else # we need to do some extra work
      xcen=value.(squeeze(slicedim(s.data,ndd,xidx),singlets)) # (N,M,L) centers
      xmin,xmax=extrema(xcen[isfinite.(xcen)])
    end
    if y==ybs.col
      yval= :edges==centersoredges ? getedges(ybs) : getcenters(ybs)
      vecy=true
      ymin,ymax=extrema(yval)
    else
      ycen=value.(squeeze(slicedim(s.data,ndd,yidx),singlets))
      ymin,ymax=extrema(ycen[isfinite.(ycen)])
    end
    if z==zbs.col
        zcen= getcenters(zbs) # always centers for z-value
        vecz=true
    else
        zcen=value.(squeeze(slicedim(s.data,ndd,zidx),singlets))
    end
    xrng=[xmin,xmax]+[-1,1]*(xmax-xmin)/100
    yrng=[ymin,ymax]+[-1,1]*(ymax-ymin)/100
    vecx&&vecy && ( (xmat,ymat)=ndgrid(xval,yval) )

    vidx=findfirst(colnames.==v)
    vval=value.(squeeze(slicedim(s.data,ndd,vidx),singlets))

    if !independentscales
        plotint=vval[isfinite.(vval)]
        vmin==-Inf && (vmin=minimum(plotint))
        vmax==+Inf && (vmax=maximum(plotint))
    else
        vminin=vmin
        vmaxin=vmax
    end

    number_of_z_slices=size(vval,3) # <- this is bad (earlier we should ensure s.data has the right permutation for this
    #nvert=floor(Int,sqrt(number_of_z_slices))
    nvert=1 # number_of_z_slices
    nhorz=cld(number_of_z_slices,nvert)
    (fig,ax)=subplots(nvert,nhorz,sharex=true,sharey=true)

    hc=[]
    for i=1:number_of_z_slices
        axes(ax[i])
        if !vecx
          thisxcen=xcen[:,:,i]
          denangrid!(thisxcen,denancycles)
          xmat= :edges==centersoredges ? grid_centers_to_corners( thisxcen ) : thisxcen # (N+1,M+1) corners
          vecy && (ymat=repeat(outer,yval,[1,size(xmat,1)])')
        end
        if !vecy
          thisycen=ycen[:,:,i]
          denangrid!(thisycen,denancycles)
          ymat= :edges==centersoredges ? grid_centers_to_corners( thisycen ) : thisycen
          vecx && ( xmat=repeat(outer,xval,[1,size(ymat,2)]) )
        end
        if vecz
            zval=zcen[i]
        else
            thiscen=zcen[:,:,i]
            denangrid!(thiscen,denancycles)
            zval=mean(thiscen)
        end
        vmat=vval[:,:,i]

        @assert ndims(xmat)==ndims(ymat)==ndims(vmat)==2 "Something went wrong with squeezing/slicing. x: $(size(xmat)), y: $(size(ymat)), z: $(size(vmat))"
        isnull(norm) || (vmat=get(norm).(vmat))
        plotint=vmat[isfinite.(vmat)] # use isfinite instead of .!isnan to remove +/-Inf too
        if !isempty(plotint)
            independentscales && (vmin = vminin==-Inf ? minimum(plotint) : vminin)
            independentscales && (vmax = vmaxin==+Inf ? maximum(plotint) : vmaxin)
            # now xval and yval are (N+1,M+1) while zval is (N,M); so we can use pcolormesh
            h=plotfun(xmat,ymat,vmat;vmin=vmin,vmax=vmax,cmap=cmap,linewidth=linewidth,kwargs...)
            mod(i,nhorz)==1 && ylabel(plainstring(s.columns[yidx]));
            nhorz*(nvert-1)<i && xlabel(plainstring(s.columns[xidx]))
            title(name(s.columns[zidx])*" = $zval "*unit(s.columns[zidx]))
            setp(gca(),"xlim",xrng,"ylim",yrng)
            push!(hc, showcolorbar && (independentscales || mod(i,nhorz)==0) ? colorbar(h,s,v) : (h,nothing) )
        end
    end
    return hc
end
function plot(s::Slice;kwargs...)
    (plotable,singlets,ndd,colnames)=setupplot(s)
    nplotable=length(plotable)
    if 1==nplotable
        _plot1d(errorbar,s,plotable,singlets,ndd,colnames;kwargs...)
    elseif 2==nplotable
        _plot2d(pcolormesh,:edges,s,plotable,singlets,ndd,colnames;kwargs...)
    elseif 3==nplotable
        _plot3d(pcolormesh,:edges,s,plotable,singlets,ndd,colnames;kwargs...)
    elseif 0==nplotable
        error("At least one plotable slice specification is necessary.")
    else
        error("I don't know how to handle >2 dimensional plots. Maybe write something like `sliceomatic`?")
    end
end
plot1d(s::Slice;kwargs...)=_plot1d(errorbar,s,setupplot(s)...;kwargs...)
plot2d(s::Slice;kwargs...)=_plot2d(pcolormesh,:edges,s,setupplot(s)...;kwargs...)
contour(s::Slice;showcolorbar::Bool=false,kwargs...)=_plot2d(contour,:centers,s,setupplot(s)...;showcolorbar=showcolorbar,kwargs...)

function waterfallplot(s::Slice; v::AbstractString=s.source.y,x1::AbstractString="",x2::AbstractString="",vmin::Real=-Inf,vmax::Real=Inf)
    @assert 2==sliceplotdimension(s) "A waterfallplot requires exactly two independent dimensions, not $(sliceplotdimension(s))"
    (plotable,singlets,ndd,colnames)=setupplot(s)

    x1bs=s.spec[plotable[1]]
    x2bs=s.spec[plotable[2]]
    vecx1=false
    vecx2=false
    isempty(x1) && (x1=x1bs.col)
    isempty(x2) && (x2=x2bs.col)
    x1idx=findfirst(colnames.==x1) # also used below for labeling
    x2idx=findfirst(colnames.==x2)

    if x1==x1bs.col # we can use the binned centers directly
      x1val=getcenters(x1bs)
      vecx1=true
    else # we need to do some extra work
      x1val=value.(squeeze(slicedim(s.data,ndd,x1idx),singlets)) # (N,M) centers
    end
    if x2==x2bs.col
      x2val=getcenters(x2bs)
      vecx2=true
    else
      x2val=value.(squeeze(slicedim(s.data,ndd,x2idx),singlets))
    end
    x1min,x1max=extrema(x1val[isfinite.(x1val)])
    x1rng=[x1min,x1max]#+[-1,1]*(x1max-x1min)/100;
    x2min,x2max=extrema(x2val[isfinite.(x2val)])
    x2rng=[x2min,x2max]#+[-1,1]*(x2max-x2min)/100
    if vecx1&&vecx2
        (x1val,x2val)=ndgrid(x1val,x2val)
    elseif vecx1
        x1val=repeatouter(x1val,[1,size(x2val,2)])
    elseif vecx2
        x2val=repeatouter(x2val,[1,size(x1val,1)])'
    end
    zidx=findfirst(colnames.==v)
    zval=value.(squeeze(slicedim(s.data,ndd,zidx),singlets))
    @assert ndims(x1val)==ndims(x2val)==ndims(zval)==2 "Something went wrong with squeezing/slicing. x: $(size(x1val)), y: $(size(x2val)), z: $(size(zval))"
    #isnull(norm) || (zval=get(norm).(zval))
    plotint=zval[isfinite.(zval)] # use isfinite instead of .!isnan to remove +/-Inf too
    if !isempty(plotint)
        # vmin,vmax=extrema(intrange)
        vmin==-Inf && (vmin=minimum(plotint))
        vmax==+Inf && (vmax=maximum(plotint))
        # now x1val, y2val and zval is (N,M); so we can plot in 3D
        h=plot3D(vec(x1val),vec(x2val),vec(zval),kwargs...)
        xlabel(plainstring(s.columns[xidx]));ylabel(plainstring(s.columns[yidx]));zlabel(plainstring(s.columns[zidx]))
        setp(gca(),"xlim",xrng,"ylim",yrng,"zlim",[vmin,vmax])
        return h
    end
end

import PyPlot2TikZ: colorbar
function colorbar(h::PyPlot2TikZ.ColoredSurfaces,s::Slice,v::AbstractString=s.source.y)
    c=colorbar(h)
    c["set_label"](plainstring(s.columns[findfirst(name.(s.columns).==v)]))
    return (h,c)
end
