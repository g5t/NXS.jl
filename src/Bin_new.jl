function bin!(p::Scatterd{0},bs::Vector{T};
              averagecounters::Bool=true,collapsedetectors::Bool=false,
              fixbinnedvariance::Bool=false,removeemptybins::Bool=true,
              refillinstrument::Bool=true,timer::Vector{Float64}=[time()],
              outputlevel=0,useccode::Bool=false
              ) where {T<:BinSpec}
    outputlevel>0 && calltimer(timer,"start of bin! function")
    bs=replaceBinSteps(p,bs) # BinStep objects don't contain enough information, so replace them in the vector
    nbs=length(bs)
    edges=map(getedges,bs)
    cents=map(getcenters,bs)
    lcens=map(length,edges)-1
    colidxs=map(x->getIdx(p,x.col),bs)
    colvals=map(x->getVal(p,x.col),bs) # this *should* be OK for Scatterd{D>0}
    oldsize=size(p.data) # this *should not* be OK for Scatterd{D>0}
    outputlevel>0 && calltimersub(timer,"allocated inital arrays")
    pntidxs=determinebin.(edges,colvals) # first determine subscript indicies # XXX A julia parallel version of this was slower
    outputlevel>0 && calltimersub(timer,"subscript indicies determined")
    # calculate the linear index from the subscript indicies e.g., for size(p.data)=(P,X,Y,N)=(100,10,5)
    # we're indexing into just the (P,X,Y) part
    #   (3,2,4) -> (4-1)*100*10 + (2-1)*100 + 3 = 3103
    binspan=[1;cumprod(lcens)]
    pntidx=pntidxs[1]
    for i=2:length(edges); pntidx+=(pntidxs[i]-1)*binspan[i]; end
    outputlevel>0 && calltimersub(timer,"converted subscripts to linear indicies")
    # Check if any of the subscript indicies are out-of-bounds (and should be excluded)
    oob=falses(oldsize[1:end-1]...)
    for i=1:length(edges)
        setindex!(oob,true,pntidxs[i].<1) # can combine as (pntidxs[i].<1)|(pntidxs[i].>lcens[i])
        setindex!(oob,true,pntidxs[i].>lcens[i]) # but make sure the parentheses are preserved!
    end
    pntidx[oob]=zero(eltype(pntidx)) # out-of-bound pixels should be set to the 0th bin
    outputlevel>0 && calltimersub(timer,"found and removed out-of-bound linear indicies")

    newsize=collapsedetectors?(prod(lcens),oldsize[end]):(prod(lcens),oldsize[2:end]...) # if colloapsing, go down to 2D.
    binno=zeros(typeof(1),newsize...)
    nd=length(oldsize) # the number of dimensions of the datablock (the dimension which will contain columns)
    outputlevel>0 && calltimersub(timer,"allocated bin data and total-contributing-pixel arrays")
    if collapsedetectors
        if useccode
          bindata=reduce_pixels_c!(getarray(p.data),pntidx,binno,newsize)
        else
          bindata=reduce_pixels!(getarray(p.data),pntidx,binno,newsize)
        end
        outputlevel>0 && calltimersub(timer,"reduced pixels into bins using "*(useccode?"C":"julia")*" code")
        if fixbinnedvariance # then colvals is used again later
            detdims=(2:nd-1...) # the detector dimensions are 2:nd-1
            for j=1:length(colvals); colvals[j]=mean(colvals[j],detdims); end
            outputlevel>0 && calltimersub(timer,"found mean of column values in detector dimensions")
        end
        # for consistency we need to increase the dimensionality of bindata/binno from 2 to nd (from MxN to Mx1x1...x1xN)
        newsize=([newsize[1];ones(eltype(newsize),nd-2);newsize[2]]...)
        bindata=reshape(bindata,newsize)
        binno=reshape(binno,newsize)
        outputlevel>0 && calltimersub(timer,"reshaped the reduced binned data up to input dimensionality")
    else
        if useccode
          bindata=distribute_pixels_c!(getarray(p.data),pntidx,binno,newsize)
        else
          bindata=distribute_pixels!(getarray(p.data),pntidx,binno,newsize) # fills binno and bindata
        end
        outputlevel>0 && calltimersub(timer,"distributed pixels into bins using "*(useccode?"C":"julia")*" code")
    end
    filled=sum(binno,nd).>0 # (Nbins1*Nbins2*...,[NdetX,[NdetY,[...]]],1)
    outputlevel>0 && calltimersub(timer,"found which bins have contributing pixels")
    if !averagecounters
        ctcolon=(find(matchBA(name.(getcolumns(p)),[p.counters;p.timers]))-1)*prod(newsize[1:end-1]) #
        # find(filled) gives the linear indicies of (Prod(Nbins),[NdetX,[NdetY,[...]]]) which have intensity
        # ctcolon adds the last dimension on for the positions of the counters and timers
        binno[find(filled).+ctcolon']=one(eltype(binno))
        outputlevel>0 && calltimersub(timer,"set contributing pixels counts to 1 for counters and timers")
    end
    # take the average of each binned value -- will produce Inf for (bin,detector[s]) with no contributing detector (pixels)
    bindata./=binno #   parallel_norm!(bindata,binno) # ~30% slower than simple approach
    outputlevel>0 && calltimersub(timer,"divided binned data by contributing pixel numbers")
    if fixbinnedvariance # FIXME this is *𝐯𝐞𝐫𝐲* slow?!
        for j=1:nbs
            for i=1:newsize[1]
                # calculate standard deviations for combined data in each column
                thisbin=pntidx.==i # Bool, (Npts,[NdetX,[NdetY,...]])
                thisval=value.(slicedim(slicedim(bindata,nd,colidxs[j]),1,i)) # size (1,[NdetX,[NdetY,...]])
                thiscen=colvals[j][thisbin] # already normal numbers, since getVal() returns such. size (<=Npts,[NdetX,[Ndet,...]])
                thisstd=sum(abs2,2thiscen.-2thisval,1)./sum(abs2,thisbin,1) # (1,[NdetX,[NdetY,...]])
                #println("thisval: ",size(thisval)," thisstd: ",size(thisstd))
                thisdat=Measured(thisval,thisstd)
                thisidx=[i;repeat([Colon()],outer=[nd-2]);colidxs[j]]
                #println("thisdat: ",size(thisdat)," bindata: ",size(bindata)," thisidx: ",thisidx)
                setindex!(bindata, thisdat, thisidx...)
            end
        end
        outputlevel>0 && calltimersub(timer,"replaced variance of binned data by deviation from their mean")
    end
    # for the binned column(s), replace their uncertainty with the bin-width(s)
    colons=repeat([Colon()],outer=[nd-1])
    repdets=[1;[newsize[2:end-1]...]]
    for j=1:nbs
        outrep=deepcopy(lcens)
        outrep[j]=1
        #
        thisval=value.(slicedim(bindata,nd,colidxs[j])) # (Nbins1*Nbins2*...,[NdetX,[NdetY,...]])
        thisvar=abs2.(diff(edges[j])/2) # (Nbinsj,)
        thisvar=repeatouter(thisvar,ones(Int64,nbs)) # (Nbinsj,1,1,...,1)
        thisvar=permutedims(thisvar,(circshift(1:nbs,j-1)...)) # (1,...,1,Nbinsj,1,...,1) # size(thisvar,j)==Nbinsj
        thisvar=repeatouter(thisvar,outrep) # (Nbins1,Nbins2,...)
        thisvar=vec(thisvar) # (Nbins1*Nbins2*... ,)
        thisvar=repeatouter(thisvar,repdets) # (Nbins1*Nbins2*...,[NdetX,[NdetY,[...]]])
        # the next one-liner can replace the previous six lines at the expense of readability
        #thisvar=repeatouter( vec(repeatouter( permutedims( repeatouter(thisvar,ones(Int64,nbs)) ,(circshift(1:nbs,j-1)...)) ,outrep)) ,repdets)
        setindex!(bindata, Measured(thisval,thisvar),[colons;colidxs[j]]...)
    end
    outputlevel>0 && calltimersub(timer,"replaced variance of bin-determining data by their bin widths")
    if removeemptybins
        notempty=falses(newsize[1])
        for i=1:newsize[1]; notempty[i]=any(filled[i,:]); end # as a reminder, filled is based on binno==0 not bindata==0 (which is desired)
        necolon=sub2ind(newsize,ones(Int64,prod(newsize[2:nd])),ind2sub(newsize[2:nd],1:prod(newsize[2:nd]))...)-1 # give all 2nd to ndth elements from 1st dimension element
        bindata=reshape(bindata[find(notempty).+necolon'],(sum(notempty),newsize[2:end]...)) # keep only bins with at least one contributing point
        outputlevel>0 && calltimersub(timer,"removed empty bins")

    end
    colcolon=(0:newsize[nd]-1)*prod(newsize[1:nd-1]) # for a linear index into the first nd-1 dimensions, colcolon provides the linear indexing for the "columns" dimension
    bindata[find(any(isinf.(bindata),nd)).+colcolon']=NaN
    outputlevel>0 && calltimersub(timer,"Replaced ±Inf with NaN")
    p.mask=squeeze(any(isnan.(bindata),nd),nd) # detectors which don't contribute to a particular bin have NaN intensity
    outputlevel>0 && calltimersub(timer,"Set mask from NaN values")
    p.data=NArray(bindata,name.s(p.data)) # save the binned data as a NArray, column names remain unchanged
    outputlevel>0 && calltimersub(timer,"Created output NArray")
    refillinstrument && (fillinstrument!(p); outputlevel>0 && calltimersub(timer,"Instrument (re)filled") )# take the nuclear approach and rebuild the instrument array :/ FIXME
    outputlevel>0 && calltimer(timer,"finished bin!")
    return p
end


function bin!(p::Scatterd{1},bs::Vector{T};
              averagecounters::Bool=true,collapsedetectors::Bool=false,
              fixbinnedvariance::Bool=false,removeemptybins::Bool=true,
              refillinstrument::Bool=true,timer::Vector{Float64}=[time()],
              outputlevel=0,useccode::Bool=false
              ) where {T<:BinSpec}
    outputlevel>0 && calltimer(timer,"start of bin! function")
    bs=replaceBinSteps(p,bs) # BinStep objects don't contain enough information, so replace them in the vector
    nbs=length(bs)
    edges=map(getedges,bs)
    cents=map(getcenters,bs)
    lcens=map(length,edges)-1
    colidxs=getIdx.(p,bs) # Vector{Int}(nbs)
    colvals=getVal.(p,bs) # Vector{Array{Measured}(size(p.detector)[1:end-1])}(nbs)
    datsize=size(p.data) # (N_points, N_0D_columns)
    detsize=size(p.detector) # (N_Points, Dim_1:M..., N_MD_columns)
    fullsize=(detsize[1:end-1]...,datsize[end]+detsize[end]) # the full block would be (N_points, Dim_1:M..., N_0D_columns+N_MD_columns)

    outputlevel>0 && calltimersub(timer,"allocated inital arrays")
    pntidxs=determinebin.(edges,colvals) # first determine subscript indicies # XXX A julia parallel version of this was slower
    outputlevel>0 && calltimersub(timer,"subscript indicies determined")
    # calculate the linear index from the subscript indicies e.g., for size(p.data)=(P,X,Y,N)=(100,10,5)
    # we're indexing into just the (P,X,Y) part
    #   (3,2,4) -> (4-1)*100*10 + (2-1)*100 + 3 = 3103
    binspan=[1;cumprod(lcens)]
    pntidx=pntidxs[1]
    for i=2:length(edges); pntidx+=(pntidxs[i]-1)*binspan[i]; end
    outputlevel>0 && calltimersub(timer,"converted subscripts to linear indicies")
    # Check if any of the subscript indicies are out-of-bounds (and should be excluded)
    oob=falses(fullsize[1:end-1]...) # XXX this is detsize[1:end-1]
    for i=1:length(edges)
        setindex!(oob,true,pntidxs[i].<1) # can combine as (pntidxs[i].<1)|(pntidxs[i].>lcens[i])
        setindex!(oob,true,pntidxs[i].>lcens[i]) # but make sure the parentheses are preserved!
    end
    pntidx[oob]=zero(eltype(pntidx)) # out-of-bound pixels should be set to the 0th bin
    outputlevel>0 && calltimersub(timer,"found and removed out-of-bound linear indicies")

    # collapsing the detector bank of a Scatterd{>0} might not be a great idea :/ FIXME
    newsize=collapsedetectors?(prod(lcens),fullsize[end]):(prod(lcens),fullsize[2:end]...) # if collapsing, go down to 2D.
    binno=zeros(typeof(1),newsize...)
    nd=length(fullsize) # the number of dimensions of the datablock (the dimension which will contain columns)
    outputlevel>0 && calltimersub(timer,"allocated bin data and total-contributing-pixel arrays")
    if collapsedetectors
        if useccode
          bindata=reduce_pixels_c!(p.data,pntidx,binno,newsize)
        else
          bindata=reduce_pixels!(p.data,pntidx,binno,newsize)
        end
        outputlevel>0 && calltimersub(timer,"reduced pixels into bins using "*(useccode?"C":"julia")*" code")
        if fixbinnedvariance # then colvals is used again later
            detdims=(2:nd-1...) # the detector dimensions are 2:nd-1
            for j=1:length(colvals); colvals[j]=mean(colvals[j],detdims); end
            outputlevel>0 && calltimersub(timer,"found mean of column values in detector dimensions")
        end
        # for consistency we need to increase the dimensionality of bindata/binno from 2 to nd (from MxN to Mx1x1...x1xN)
        newsize=([newsize[1];ones(eltype(newsize),nd-2);newsize[2]]...)
        bindata=reshape(bindata,newsize)
        binno=reshape(binno,newsize)
        outputlevel>0 && calltimersub(timer,"reshaped the reduced binned data up to input dimensionality")
    else
        if useccode
          bindata=distribute_pixels_c!(p.data,pntidx,binno,newsize)
        else
          bindata=distribute_pixels!(p.data,pntidx,binno,newsize) # fills binno and bindata
        end
        outputlevel>0 && calltimersub(timer,"distributed pixels into bins using "*(useccode?"C":"julia")*" code")
    end
    filled=sum(binno,nd).>0 # (Nbins1*Nbins2*...,[NdetX,[NdetY,[...]]],1)
    outputlevel>0 && calltimersub(timer,"found which bins have contributing pixels")
    if !averagecounters
        ctcolon=(find(matchBA(name.(p.columns),[p.counters;p.timers]))-1)*prod(newsize[1:end-1]) #
        # find(filled) gives the linear indicies of (Prod(Nbins),[NdetX,[NdetY,[...]]]) which have intensity
        # ctcolon adds the last dimension on for the positions of the counters and timers
        binno[find(filled).+ctcolon']=one(eltype(binno))
        outputlevel>0 && calltimersub(timer,"set contributing pixels counts to 1 for counters and timers")
    end
    # take the average of each binned value -- will produce Inf for (bin,detector[s]) with no contributing detector (pixels)
    # if allparallel
    #   parallel_norm!(bindata,binno) # ~30% slower than simple approach
    # else
      bindata./=binno
    # end
    outputlevel>0 && calltimersub(timer,"divided binned data by contributing pixel numbers")
    if fixbinnedvariance # FIXME this is *𝐯𝐞𝐫𝐲* slow?!
        for j=1:nbs
            for i=1:newsize[1]
                # calculate standard deviations for combined data in each column
                thisbin=pntidx.==i # Bool, (Npts,[NdetX,[NdetY,...]])
                thisval=value.(slicedim(slicedim(bindata,nd,colidxs[j]),1,i)) # size (1,[NdetX,[NdetY,...]])
                thiscen=colvals[j][thisbin] # already normal numbers, since getVal() returns such. size (<=Npts,[NdetX,[Ndet,...]])
                thisstd=sum(abs2,2thiscen.-2thisval,1)./sum(abs2,thisbin,1) # (1,[NdetX,[NdetY,...]])
                #println("thisval: ",size(thisval)," thisstd: ",size(thisstd))
                thisdat=Measured(thisval,thisstd)
                thisidx=[i;repeat([Colon()],outer=[nd-2]);colidxs[j]]
                #println("thisdat: ",size(thisdat)," bindata: ",size(bindata)," thisidx: ",thisidx)
                setindex!(bindata, thisdat, thisidx...)
            end
        end
        outputlevel>0 && calltimersub(timer,"replaced variance of binned data by deviation from their mean")
    end
    # for the binned column(s), replace their uncertainty with the bin-width(s)
    colons=repeat([Colon()],outer=[nd-1])
    repdets=[1;[newsize[2:end-1]...]]
    for j=1:nbs
        outrep=deepcopy(lcens)
        outrep[j]=1
        #
        thisval=value.(slicedim(bindata,nd,colidxs[j])) # (Nbins1*Nbins2*...,[NdetX,[NdetY,...]])
        thisvar=abs2.(diff(edges[j])/2) # (Nbinsj,)
        thisvar=repeatouter(thisvar,ones(Int64,nbs)) # (Nbinsj,1,1,...,1)
        thisvar=permutedims(thisvar,(circshift(1:nbs,j-1)...)) # (1,...,1,Nbinsj,1,...,1) # size(thisvar,j)==Nbinsj
        thisvar=repeatouter(thisvar,outrep) # (Nbins1,Nbins2,...)
        thisvar=vec(thisvar) # (Nbins1*Nbins2*... ,)
        thisvar=repeatouter(thisvar,repdets) # (Nbins1*Nbins2*...,[NdetX,[NdetY,[...]]])
        # the next one-liner can replace the previous six lines at the expense of readability
        #thisvar=repeatouter( vec(repeatouter( permutedims( repeatouter(thisvar,ones(Int64,nbs)) ,(circshift(1:nbs,j-1)...)) ,outrep)) ,repdets)
        setindex!(bindata, Measured(thisval,thisvar),[colons;colidxs[j]]...)
    end
    outputlevel>0 && calltimersub(timer,"replaced variance of bin-determining data by their bin widths")
    if removeemptybins
        notempty=falses(newsize[1])
        for i=1:newsize[1]; notempty[i]=any(filled[i,:]); end # as a reminder, filled is based on binno==0 not bindata==0 (which is desired)
        necolon=sub2ind(newsize,ones(Int64,prod(newsize[2:nd])),ind2sub(newsize[2:nd],1:prod(newsize[2:nd]))...)-1 # give all 2nd to ndth elements from 1st dimension element
        p.data=reshape(bindata[find(notempty).+necolon'],(sum(notempty),newsize[2:end]...)) # keep only bins with at least one contributing point
        #fullnotempty=repeatouter(notempty,[1;newsize[2:end]...])
        #p.data=reshape(bindata[fullnotempty],(sum(notempty),newsize[2:end]...)) # keep only bins with at least one contributing point
        outputlevel>0 && calltimersub(timer,"removed empty bins while saving binned data to output")
    else
        p.data=bindata
        outputlevel>0 && calltimersub(timer,"saved binned data for output")
    end
    colcolon=(0:newsize[nd]-1)*prod(newsize[1:nd-1]) # for a linear index into the first nd-1 dimensions, colcolon provides the linear indexing for the "columns" dimension
    p.data[find(any(isinf.(p.data),nd)).+colcolon']=NaN
    outputlevel>0 && calltimersub(timer,"Replaced ±Inf with NaN")
    p.mask=squeeze(any(isnan.(p.data),nd),nd) # detectors which don't contribute to a particular bin have NaN intensity
    outputlevel>0 && calltimersub(timer,"Set mask from NaN values")
    refillinstrument && (fillinstrument!(p); outputlevel>0 && calltimersub(timer,"Instrument (re)filled") )# take the nuclear approach and rebuild the instrument array :/ FIXME
    outputlevel>0 && calltimer(timer,"finished bin!")
    return p
end
