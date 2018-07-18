export BinSpec,BinStep,BinInt,BinRange,BinEdges,BinCenters
"""
A `BinSpec` details how to combine an N-dimensional block of data.
The properties of any `BinSpec` must in some way define which piece of the
data block should be used as a binning variable and also must define the
bins into which the data is collected.

The simplest, `BinStep`, contains the name of the data to use as the binning criteria
and the size of each bin. `BinRanges` adds the minimum and maximum bin center values.
The `BinCenters` type contains a `Vector` of the bin centers, allowing non-uniform bin sizes.
The `BinEdges` type also contains a `Vector` allowing for non-uniform bin sizes,
but it is more precise than a list of bin centers as all center-information can be
calculated from edge-information but the reverse is not true.
"""
abstract type BinSpec{T} end
# Base.start(a::BinSpec)=1

# TODO add cyclic binning type that has a periodic boundary condidition (to bin angular ranges?)
immutable BinStep{T} <: BinSpec{T}
    col::AbstractString
    step::T
    function BinStep{R}(a::AbstractString,b::Vector{R}) where R
        @assert 1==length(b)
        return new(a,b[1])
    end
end
immutable BinInt{T} <: BinSpec{T}
    col::AbstractString
    start::T
    stop::T
    BinInt{R}(a::AbstractString,b1::R,b2::R) where R=new(a,minimum([b1,b2]),maximum([b1,b2]))
    function BinInt{R}(a::AbstractString,b::Vector{R}) where R
        @assert 2==length(b)
        return new(a,minimum(b),maximum(b))
    end
end
# Base.length(a::BinInt)=2
# Base.done(a::BinInt,i)=i>2
# Base.next(a::BinInt,i)=( i==1?a.start:i==2?a.stop:NaN, i+1 )
Base.first(a::BinInt)=a.start
Base.last(a::BinInt)=a.stop
Base.extrema(a::BinInt)=(a.start,a.stop)

BinInt(p::Scatterd,col::AbstractString)=BinInt(col,extrema(getVal(p,col))...) # automatic complete integration over a column
immutable BinRange{T} <: BinSpec{T}
    col::AbstractString
    start::T
    step::T
    stop::T
    function BinRange{R}(a::AbstractString,b1::R,b2::R,b3::R) where R
        ss=[b1,b3]
        ssdiff=abs(b1-b3)
        st=abs(b2)
        ssdiff<st && (st=ssdiff)
        sv=minimum(ss):st:maximum(ss) # to ensure an integer number of bins
        return new(a,first(sv),st,last(sv))
    end
    function BinRange{R}(a::AbstractString,b::Vector{R}) where R
        @assert 3==length(b)
        ss=b[[1,3]]
        ssdiff=abs(ss[1]-ss[2])
        st=abs(b[2])
        ssdiff<st && (st=ssdiff)
        sv=minimum(ss):st:maximum(ss) # to ensure an integer number of bins
        return new(a,first(sv),st,last(sv))
    end
end
#Base.length(a::BinRange) = ceil(Int, 1+(a.stop-a.start)/a.step )
#Base.done(a::BinRange,i)= i > length(a)
#Base.next(a::BinRange,i)=( a.start + (i-1)*a.step , i+1 )
Base.first(a::BinRange)=a.start
Base.last(a::BinRange)=a.stop
Base.extrema(a::BinRange)=(a.start,a.stop)


function BinRange(p::Scatterd,bs::BinStep)
    col=getVal(p,bs.col)
    BinSpec(bs.col,[minimum(col),bs.step,maximum(col)])
end
immutable BinEdges{T} <: BinSpec{T}
    col::AbstractString
    edges::Vector{T}
end
Base.first(a::BinEdges)=Base.first(a.edges)
Base.last(a::BinEdges)=Base.last(a.edges)
Base.extrema(a::BinEdges)=Base.extrema(a.edges)

immutable BinCenters{T} <: BinSpec{T}
    col::AbstractString
    centers::Vector{T}
end
Base.first(a::BinCenters)=Base.first(a.centers)
Base.last(a::BinCenters)=Base.last(a.centers)
Base.extrema(a::BinCenters)=Base.extrema(a.centers)

function BinSpec{T<:Real}(s::AbstractString,a::Union{Range{T},AbstractArray{T,1},T})
    spec=collect(a) # now guaranteed to be a dense Array
    len=length(spec)
    R=eltype(spec)
    1==len && (return BinStep{R}(s,spec))
    2==len && (return BinInt{R}(s,spec))
    3==len && (return BinRange{R}(s,spec))
    return BinEdges{R}(s,spec)
end
BinSpec(s::AbstractString,a::Real)=BinSpec(s,[a]) # allow for easy BinStep creation

# simple(r) creation of BinSpec objects.
# We already have ("h",[-1,0.1,1]) being interpreted at BinSpec("h",[-1,0.1,1])
# in the calls to, e.g., bin!; but this doesn't allow for ("h",[-1,0.1,1])+h0
# if bin"h -1 0.1 1" is turned into BinSpec("h",[-1,0.1,1]) then our range-modifying
# overloading will work.
export @bin_str
macro bin_str(binstring)
    s=split(replace(binstring,['[',',',']'],' '))
    BinSpec(shift!(s),[parse.(s)...])
end


function BinSpec{T<:Real}(p::Scatterd,col::AbstractString,a::Union{Range{T},AbstractArray{T,1},T})
  spec=collect(a)
  len=length(spec)
  R=eltype(spec)
  1==len && (return BinRange(p,BinStep(col,spec))) # auto convert to BinRange
  2==len && (return BinInt{R}(col,spec))
  3==len && (return BinRange{R}(col,spec))
  return BinEdges{R}(col,spec)
end
BinSpec(p::Scatterd,col::AbstractString)=BinInt(p,col)
BinSpec(p::Scatterd,col::AbstractString,step)=BinRange(p,BinStep(col,step))
function BinSpec(p::Scatterd,col::AbstractString,start,stop)
  (a,b)=promote(start,stop)
  BinInt{typeof(a)}(col,[a,b])
end
function BinSpec(p::Scatterd,col::AbstractString,one,two,three)
  (a,b,c)=promote(one,two,three)
  BinRange{typeof(a)}(col,[a,b,c])
end


# overloading addition and subtraction for shifting ranges seems like a good idea
for y in (:+,:-)
    @eval $y{R,T}(a::BinInt{T},n::R)=BinInt{promote_type(T,R)}(a.col,$y(a.start,n),$y(a.stop,n))
    @eval $y{R,T}(a::BinRange{T},n::R)=BinRange{promote_type(T,R)}(a.col,$y(a.start,n),a.step,$y(a.stop,n))
    @eval $y{R,T}(a::BinEdges{T},n::R)=BinEdges{promote_type(T,R)}(a.col,$y(a.edges,n))
    @eval $y{R,T}(a::BinCenters{T},n::R)=BinCenters{promote_type(T,R)}(a.col,$y(a.centers,n))
end

showbinspec(io::IO,a::BinStep   ,compact::Bool=false)=Base.print(io," δ"*a.col*"=$(a.step)")
showbinspec(io::IO,a::BinInt    ,compact::Bool=false)=Base.print(io," "*a.col*"=[$(a.start),$(a.stop)]")
showbinspec(io::IO,a::BinRange  ,compact::Bool=false)=Base.print(io," "*a.col*"=($(a.start):$(a.step):$(a.stop))")
function showbinspec(io::IO,a::BinEdges  ,compact::Bool=false)
    Base.print(io," "*a.col)
    if compact
        edgelist=sort(a.edges)
        if 4 >= length(edgelist)
            Base.print(io,"=[",join(edgelist,","),"]")
        else
            Base.print(io,"=[",join(edgelist[1:2],","),",…,",join(edgelist[end-1:end],","),"]")
        end
    else
        Base.print(io," edges=$(a.edges)")
    end
end
function showbinspec(io::IO,a::BinCenters,compact::Bool=false)
    Base.print(io," "*a.col)
    if compact
        centerslist=sort(a.centers)
        if 4 >= length(centerslist)
            Base.print(io,"=(",join(centerslist,","),")")
        else
            Base.print(io,"=(",join(centerslist[1:2],","),",…,",join(centerslist[end-1:end],","),")")
        end
    else
        Base.print(io," centers=$(a.centers)")
    end
end

Base.show(io::IO,a::BinSpec)=showbinspec(io,a,false)
Base.showcompact(io::IO,a::BinSpec)=showbinspec(io,a,true)

export isbinning
isbinning(a::BinInt)=false
isbinning(a::BinSpec)=true
isbinning{T<:BinSpec}(a::Array{T})=map(isbinning,a)

export nobins
nobins(a::BinInt)=1
nobins(a::BinRange)=length(a.start:a.step:a.stop)
nobins(a::BinEdges)=length(a.edges)-1
nobins(a::BinCenters)=length(a.centers)
nobins(a::BinSpec)=throw(TypeError(:nobins,"BinND.jl needs more information",Union{BinInt,BinRange,BinEdges,BinCenters},typeof(a)))
nobins{T<:BinSpec}(a::Array{T})=map(nobins,a)


getedges(bi::BinInt)=[bi.start,bi.stop]
getedges(br::BinRange)=(cen=getcenters(br); vcat(cen-br.step/2,cen[end]+br.step/2) )
getedges(be::BinEdges)=be.edges
getedges(bc::BinCenters)=(d=diff(bc.centers); vcat(bc.centers[1]-d[1]/2, bc.centers[1:end-1]+d/2, bc.centers[end]+d[end]/2))
getedges(bs::BinSpec)=throw(TypeError(:edges,"Not enough information to calculate range",Union{BinRange,BinEdges,BinCenters},typeof(bs)))
export getedges
getcenters(bi::BinInt)=[(bi.start+bi.stop)/2]
getcenters(br::BinRange)=collect(br.start:br.step:br.stop)
getcenters(be::BinEdges)=be.edges[1:end-1]+diff(be.edges)/2
getcenters(bc::BinCenters)=bc.centers
getcenters(bs::BinSpec)=throw(TypeError(:centers,"Not enough information to calculate range",Union{BinRange,BinEdges,BinCenters},typeof(bs)))
export getcenters

getcolumn(bs::BinSpec)=bs.col

# Detector0d, Detector1d, Detector2d, ...
# original Detector0d binning still exists in v0.4 copy of ScatteringInstruments
bin!(p::Scatterd,bst::Tuple...;kwds...)=bin!(p,[BinSpec(p,t...) for t in bst];kwds...)
bin!(p::Scatterd,bs::BinSpec...;kwds...)=bin!(p,[bs...];kwds...)
function fmttime(time)
    @sprintf("%8.4f",time)
end
function calltimer(timer,msg="timer")
    push!(timer,time())
    status(:timer,fmttime(timer[end]-timer[1])*": "*msg)
end
function calltimersub(timer,msg="timer")
    push!(timer,time())
    status(:timer,fmttime(timer[end]-timer[1])*": ("*fmttime(timer[end]-timer[end-1])*") "*msg)
end
function replaceBinSteps{T<:BinSpec}(p::Scatterd,bs::Vector{T})
    nbs=length(bs)
    if any(x->isa(x,BinStep),bs) # we need to completely replace bs :(
        bintypes=map(typeof,bs)
        toreplace=find(x->x<:BinStep,bintypes)
        for i in toreplace
            bintypes[i]=typeof(BinRange(p,bs[i]))
        end
        newbs=Array{Union{bintypes...}}(nbs)
        for i=1:nbs
            if !any(i.==toreplace)
                newbs[i]=bs[i]
            end
        end
        for i in toreplace
            newbs[i]=BinRange(p,bs[i])
        end
    else
        newbs=bs
    end
    return newbs
end
Base.bin(p::Scatterd,o...;k...)=(a=copy(p);bin!(a,o...;k...);a)
function bin!(p::Scatterd{0},bs::Vector{T};
              averagecounters::Bool=true,collapsedetectors::Bool=false,
              fixbinnedvariance::Bool=false,removeemptybins::Bool=true,
              excludemaskedpoints::Bool=false,
              refillinstrument::Bool=true,timer::Vector{Float64}=[time()],
              outputlevel=0,useccode::Bool=false
              ) where {T<:BinSpec}
    outputlevel>0 && calltimer(timer,"start of bin! function")
    bs=replaceBinSteps(p,bs) # BinStep objects don't contain enough information, so replace them in the vector
    nbs=length(bs)
    edges=map(getedges,bs)
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

    excludemaskedpoints && ( pntidx[p.mask]=zero(eltype(pntidx)) ) # points with mask==true should be set to 0th bin

    newsize=collapsedetectors?(prod(lcens),oldsize[end]):(prod(lcens),oldsize[2:end]...) # if colloapsing, go down to 2D.
    binno=zeros(typeof(1),newsize...)
    nd=length(oldsize) # the number of dimensions of the datablock (the dimension which will contain columns)
    outputlevel>0 && calltimersub(timer,"allocated bin data and total-contributing-pixel arrays")
    if collapsedetectors
        if useccode
        #   bindata=reduce_pixels_c!(getarray(p.data),pntidx,binno,newsize)
          bindata=reduce_pixels_c!(p.data,pntidx,binno,newsize)
        else
        #   bindata=reduce_pixels!(getarray(p.data),pntidx,binno,newsize)
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
        #   bindata=distribute_pixels_c!(getarray(p.data),pntidx,binno,newsize)
          bindata=distribute_pixels_c!(p.data,pntidx,binno,newsize)
        else
        #   bindata=distribute_pixels!(getarray(p.data),pntidx,binno,newsize) # fills binno and bindata
          bindata=distribute_pixels!(p.data,pntidx,binno,newsize)
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
    p.data=bindata # save the binned data as a NArray, column names remain unchanged
    outputlevel>0 && calltimersub(timer,"Stored output Array in data field")
    refillinstrument && (fillinstrument!(p); outputlevel>0 && calltimersub(timer,"Instrument (re)filled") )# take the nuclear approach and rebuild the instrument array :/ FIXME
    outputlevel>0 && calltimer(timer,"finished bin!")
    return p
end
bin!(p::Scatterd{D},bs::Vector{T};k...) where {D,T} = error("bin! of a Scatterd{>0} object does not exist yet.")

function repeatouter(a::AbstractArray,outer::Vector)
    lot=length(outer)
    nda=ndims(a)
    @assert all(outer.>0)
    out=copy(a)
    tmp=copy(a)
    for i=1:lot
        for j=2:outer[i] # only repeat for >1 times
            out=cat(i,out,tmp)
        end
        tmp=copy(out)
    end
    lot>nda&&2>outer[lot]&&(out=cat(lot,out)) # promote the shape to higher dimension if necessary
    return out
end

function determinebin{T}(edg::Union{Range{T},AbstractArray{T,1}},val,msk=falses(val))
    nd=ndims(val)
    edg=permutedims(cat(nd+1,collect(edg)),circshift(1:nd+1,-1))# move vector into nd+1 dimension (1,1,...,1,length(edg))
    idx=squeeze(sum(val.>=edg,nd+1),nd+1)
    idx[msk]=0
    return idx
end
# XXX a parallel version of determinebin was slower than the implementation above.

function mynormbins(q::SharedArray)
    idx=indexpids(q)
    0==idx && (return 1:0)
    nchunks=length(procs(q))
    splits=[round(Int,s) for s in linspace(0,length(q),nchunks+1)]
    splits[idx]+1:splits[idx+1]
end
function parallel_norm_chunk!(a::SharedArray,b::SharedArray)
    for i in mynormbins(a)
        a[i]=a[i]/b[i]
    end
end
function parallel_norm!(a::SharedArray,b::SharedArray)
    @assert size(a)==size(b)
    @sync begin
        for p in procs(a)
            @async remotecall_wait(parallel_norm_chunk!,p,a,b)
        end
    end
end
function parallel_norm!{T<:Measureds.Measured,R<:Measureds.Measured}(a::Array{T},b::Array{R})
  (ia,ib)=Measureds.all_sym_or_asym(a,b) # creates new matricies of a single bits type
  _bitstype_parallel_norm!(ia,ib)
  a[:]=ia[:]
end
function parallel_norm!{T<:Measureds.Measured}(a::Array{T},b::Array)
    ia=Measureds.all_sym_or_asym(a)
    _bitstype_parallel_norm!(ia,b)
    a[:]=ia[:]
end
parallel_norm!{T,R<:Measureds.Measured}(a::Array{T},b::Array{R})=_bitstype_parallel_norm!(a,Measureds.all_sym_or_asym(b))
parallel_norm!{T}(a::Array{T},b::Array{T})=_bitstype_parallel_norm!(a,b)
function _bitstype_parallel_norm!(a::Array,b::Array)
    shareda=convert(SharedArray,a)
    parallel_norm!(shareda,convert(SharedArray,b))
    a[:]=shareda[:]
end

function distribute_pixels_c!{I<:Integer,T<:Measureds.Measured}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
  # a special version for an array of the abstract Measured type
  if allSymmetric(source)
    (sval,svar)=libnxs_dpp2!(value.(source),variance.(source),idx,totl,sinksize)
    return MeasuredSymmetric(sval,svar)
  else
    (sval,svap,svam)=libnxs_dpp3!(value.(source),positivevariance.(source),negativevariance.(source),idx,totl,sinksize)
    return MeasuredAsymmetric(sval,svap,svam)
  end
end
function distribute_pixels_c!{I<:Integer,T<:Measureds.MeasuredSymmetric}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
  (sval,svar)=libnxs_dpp2!(value.(source),variance.(source),idx,totl,sinksize)
  return MeasuredSymmetric(sval,svar)
end
function distribute_pixels_c!{I<:Integer,T<:Measureds.MeasuredAsymmetric}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
  (sval,svap,svam)=libnxs_dpp3!(value.(source),positivevariance.(source),negativevariance.(source),idx,totl,sinksize)
  return MeasuredAsymmetric(sval,svap,svam)
end
function distribute_pixels_c!{I<:Integer,T<:Float64}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
  sink=libnxs_dpp1(source,idx,totl,sinksize)
end


function reduce_pixels_c!{I<:Integer,T<:Measureds.Measured}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
  # a special version for an array of the abstract Measured type
  if allSymmetric(source)
    (sval,svar)=libnxs_rpp2!(value.(source),variance.(source),idx,totl,sinksize)
    return MeasuredSymmetric(sval,svar)
  else
    (sval,svap,svam)=libnxs_rpp3!(value.(source),positivevariance.(source),negativevariance.(source),idx,totl,sinksize)
    return MeasuredAsymmetric(sval,svap,svam)
  end
end
function reduce_pixels_c!{I<:Integer,T<:Measureds.MeasuredSymmetric}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
  (sval,svar)=libnxs_rpp2!(value.(source),variance.(source),idx,totl,sinksize)
  return MeasuredSymmetric(sval,svar)
end
function reduce_pixels_c!{I<:Integer,T<:Measureds.MeasuredAsymmetric}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
  (sval,svap,svam)=libnxs_rpp3!(value.(source),positivevariance.(source),negativevariance.(source),idx,totl,sinksize)
  return MeasuredAsymmetric(sval,svap,svam)
end
function reduce_pixels_c!{I<:Integer,T<:Float64}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
  sink=libnxs_rpp1!(source,idx,totl,sinksize)
end

function reduce_pixels_g!{I<:Integer,T<:Measureds.Measured}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
  # a special version for an array of the abstract Measured type
  if allSymmetric(source)
    (sval,svar)=libnxs_rpg2!(value(source),variance(source),idx,totl,sinksize)
    return MeasuredSymmetric(sval,svar)
  else
    (sval,svap,svam)=libnxs_rpg3!(value(source),positivevariance(source),negativevariance(source),idx,totl,sinksize)
    return MeasuredAsymmetric(sval,svap,svam)
  end
end
function reduce_pixels_g!{I<:Integer,T<:Measureds.MeasuredSymmetric}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
  (sval,svar)=libnxs_rpg2!(value(source),variance(source),idx,totl,sinksize)
  return MeasuredSymmetric(sval,svar)
end
function reduce_pixels_g!{I<:Integer,T<:Measureds.MeasuredAsymmetric}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
  (sval,svap,svam)=libnxs_rpg3!(value(source),positivevariance(source),negativevariance(source),idx,totl,sinksize)
  return MeasuredAsymmetric(sval,svap,svam)
end
function reduce_pixels_g!{I<:Integer,T<:Float64}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
  sink=libnxs_rpg1!(source,idx,totl,sinksize)
end

function distribute_pixels_common{I<:Integer,T}(source::AbstractArray{T},idx::AbstractArray{I},totl::AbstractArray{I},sink::AbstractArray{T})
  dsource=size(source); didx=size(idx); dtotl=size(totl); dsink=size(sink)
  @assert dsource[2:end]==dsink[2:end]
  @assert dtotl==dsink
  @assert didx==dsource[1:end-1]
  nd=length(dsource)
  nsrcpix=prod(dsource[1:nd-1])
  nc=dsource[end]
  sourcecolon=(0:nc-1)*prod(dsource[1:nd-1])
  sinkcolon=(0:nc-1)*prod(dsink[1:nd-1])
  return (sourcecolon,sinkcolon,nsrcpix,nc)
end

function mysinkbins(q::SharedArray)
    idx=indexpids(q)
    0==idx && (return 1:0)
    nchunks=length(procs(q))
    splits=[round(Int,s) for s in linspace(0,size(q,1),nchunks+1)]
    splits[idx]+1:splits[idx+1]
end
function distribute_pixels!{I<:Integer,T<:Measureds.Measured}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
  # The Measureds.Measured type is no longer a bitstype (but MeasuredSymmetric and MeasuredAsymmetric are)
  # convert(SharedArray,x) only works for x<:bitstype so we need to make sure our Array is a singular bitstype Array
  _bitstype_distribute_pixels!(Measureds.all_sym_or_asym(source),idx,totl,sinksize) # creates new matrix of a single bits type :(
end
distribute_pixels!{I<:Integer,T}(s::Array{T},i::Array{I},t::Array{I},k)=_bitstype_distribute_pixels!(s,i,t,k)
# Before this splitting of distribute_pixels! into a type-depedenent wrapper and an internal
# SharedArray creator, the d_p!{I,T<:Measured} version would stack overflow by repeatedly calling itself.
function _bitstype_distribute_pixels!{I<:Integer,T}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
    sharedsink=SharedArray{T}(sinksize, init = S -> S[Base.localindexes(S)] = 0)
    sharedtotl=convert(SharedArray,totl)
    sharedidx=convert(SharedArray,idx)
    sharedsource=convert(SharedArray,source)
    distribute_pixels_parallel!(sharedsource,sharedidx,sharedtotl,sharedsink)
    totl[:]=sharedtotl[:];
    sink=sdata(sharedsink)
    finalize(sharedsource); finalize(sharedidx); finalize(sharedtotl); finalize(sharedsink);
    return sink
end
function distribute_pixels_parallel!{I<:Integer,T}(source::SharedArray{T},idx::SharedArray{I},totl::SharedArray{I},sink::SharedArray{T})
    (sourcecolon,sinkcolon,nosrcpix,nocols)=distribute_pixels_common(source,idx,totl,sink)
    @sync begin
        for p in procs(sink)
            @async remotecall_wait(distribute_pixels_chunk!,p,source,sourcecolon,idx,totl,sink,sinkcolon)
        end
    end
end
function distribute_pixels_chunk!(src,srccln,idx,ttl,snk,snkcln)
    didx=size(idx)
    dsnk=size(snk)
    for i in mysinkbins(snk)
        thisidx=find(idx.==i) # linear indexing of matching pixels (i,j,k)->i+a(j-1)+ab(k-1)
        # subidx=ind2sub(didx,thisidx) # their subinidicies ( [i1,i2,...], [j1,j2,...], ... )
        # snkidx=i*ones(subidx[1]) # [i,i,...] the size of [i1,i2,...]
        # for j=2:length(subidx) # for the remaining subinidicies
        #   snkidx+=(subidx[j]-1)*prod(dsnk[1:j-1]) # find the linear indexing into snk
        # end
        # # the following is too simplistic (?) since we're only looping over the first dimension with i
        srcidx=srccln .+ thisidx' # (ncols, length(thisidx))
        snkidx=snkcln .+ i
        for j=1:length(thisidx)
            snk[snkidx]+=src[srcidx[:,j]]
            ttl[snkidx]+=1
        end
    end
end


function reduce_pixels_common{I<:Integer,T}(source::AbstractArray{T},idx::AbstractArray{I},totl::AbstractArray{I},sink::AbstractArray{T})
  dsource=size(source); didx=size(idx); dtotl=size(totl); dsink=size(sink); nd=length(dsource); nsrcpix=prod(dsource[1:nd-1]); nc=dsource[nd]
  @assert length(dsink)==2 "Reducing requires that sink be (Nbins,Ncolumns) not $dsink"
  @assert dsink[2]==nc "Source and sink must have equal Ncolumns, $(nc)≠$(dsink[2])"
  @assert dtotl==dsink "The total and sink arrays must have equivalent sizes"
  @assert didx==dsource[1:nd-1] "The index array must be (Npoints,[NdetectorsX,[NdetectorsY,...]])"
  sourcecolon=(0:nc-1)*prod(dsource[1:nd-1])
  sinkcolon=(0:nc-1)*dsink[1] # sink is (Nbins,Ncolumns), so we want [0,1,2,3...]*Nbins for the last-dimension "colon" operator
  return (sourcecolon,sinkcolon,nsrcpix,nc)
end

function reduce_pixels!{I<:Integer,T<:Measureds.Measured}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
  # The Measureds.Measured type is no longer a bitstype (but MeasuredSymmetric and MeasuredAsymmetric are)
  # convert(SharedArray,x) only works for x<:bitstype so we need to make sure our Array is a singular bitstype Array
  _bitstype_reduce_pixels!(Measureds.all_sym_or_asym(source),idx,totl,sinksize) # creates a new matrix of a single bits type :(
end
reduce_pixels!{I<:Integer,T}(s::Array{T},i::Array{I},t::Array{I},k)=_bitstype_reduce_pixels!(s,i,t,k)
function _bitstype_reduce_pixels!{I<:Integer,T}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
    sharedsink=SharedArray{T}(sinksize, init = S -> S[Base.localindexes(S)] = 0)
    sharedtotl=convert(SharedArray,totl)
    sharedidx=convert(SharedArray,idx)
    sharedsource=convert(SharedArray,source)
    reduce_pixels_parallel!(sharedsource,sharedidx,sharedtotl,sharedsink)
    sink=sdata(sharedsink) # pull-back the SharedArray to just an Array
    totl[:]=sharedtotl[:]
    finalize(sharedsource); finalize(sharedidx); finalize(sharedtotl); finalize(sharedsink)
    return sink
end
function reduce_pixels_parallel!{I<:Integer,T}(source::SharedArray{T},idx::SharedArray{I},totl::SharedArray{I},sink::SharedArray{T})
    (sourcecolon,sinkcolon,nosrcpix,nocolumns)=reduce_pixels_common(source,idx,totl,sink)
    @sync begin
        for p in procs(sink)
            @async remotecall_wait(reduce_pixels_chunk!,p,source,sourcecolon,idx,size(idx),totl,sink,size(sink),sinkcolon)
        end
    end
end
function reduce_pixels_chunk!(src,srccln,idx,didx,ttl,snk,dsnk,snkcln)
    for i in mysinkbins(snk)
        thisidx=find(idx.==i) # linear indexing of matching pixels (i,j,k)->i+a(j-1)+ab(k-1)
        srcidx=srccln .+ thisidx' # (ncols, length(thisidx))
        snkidx=snkcln .+ i # (ncols,) << all go into a single sink bin for certain
        for j=1:length(thisidx)
            snk[snkidx]+=src[srcidx[:,j]]
        end
        ttl[snkidx]=length(thisidx) # since all go into a single bin, the total contribution is easy
    end
end
