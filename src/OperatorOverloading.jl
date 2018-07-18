function guessnormalization(a::Scatterd{0})
    ds=size(a.data)
    nd=length(ds)-1 # getVal returns an ndims(a.data)-1 dimensioned array
    torc=[a.counters;a.timers]
    z=zeros(length(torc))
    # Using the standard deviation of the difference caused problems
    # when two scans of different monitor counts were combined.
    # Here use the mean derivative, so single-point jumps in normalization (hopefully)
    # don't foul up the selection.
    for i=1:length(z)
        tcv=getVal(a,torc[i])
        for j=nd:-1:2; tcv=squeeze(sum(tcv,j),j); end # sum over detector dimensions
        z[i]=mean(abs.(diff(tcv)./tcv[1:end-1])) # z is the mean per-point derivative of the counts over the whole detector(s)
    end
    m=indmin(z) # indmin (and minimum) thankfully ignore NaN values
    ( torc[m] , minimum(getVal(a,torc[m])) )
end
function guessnormalization(a::Scatterd)
    torc=[a.counters;a.timers]
    z=zeros(length(torc))
    for i=1:length(z)
        tcv=value.( get_data_or_detector_column(a,torc[i]) )
        d=ntuple(i->i+1,ndims(tcv)-1) # (2,3,...)
        isempty(d) || ( tcv=squeeze(sum(tcv,d),d) )  # sum over detector dimensions
        z[i]=mean(abs.(diff(tcv)./tcv[1:end-1])) # z is the mean per-point derivative of the counts over the whole detector(s)
    end
    m=indmin(z) # indmin (and minimum) thankfully ignore NaN values
    ( torc[m] , minimum(value.(get_data_or_detector_column(a,torc[m]))) )
end

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
# XXX vvvv This code needs to be rewritten to allow for >0 dimensional detectors vvvv XXX
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
function combineRepeatedPoints(a::Detector0d,k...) # FIXME a.x doesn't necessarily exist for >0 dimensional detectors
    #TODO instead of picking a bin based on within-error ax, pick a bin by comparing
    #     all non counter/monitors/timers? This would likely also require treating temperature
    #     and field devices separately.
    norminfo=guessnormalization(a)
    ax=getDat(a,a.x)
    ii=fuzzyindexin(ax,ax) # equivalent (within error) points have the same ii value
    # fuzzyindexin returns binno for Measured values
    uii=unique(ii)
    binno=Array{Int}(size(ax,1))
    for i=1:length(ax)
        binno[i] = findfirst(uii.==ii[i]) # there should be only 1 match
    end
    @assert maximum(binno)==length(uii)
    nbins=maximum(binno)
    outdata=Array{eltype(a.data)}(nbins,size(a.data,(2:ndims(a.data))...))
    numdata=Array{Int}(size(outdata)...)
    for i=1:nbins
        outdata[i,:]=sum(a.data[binno.==i,:],1)
        numdata[i,:]=sum(binno.==i)
    end
    # average everything *except* the y-column(s) and any monitors/timers/counters
    numdata[:,matchBA(name.(getcolumns(a)),vcat(a.y,a.counters,a.timers))]=1
    acopy=deepcopy(a)
    acopy.data=outdata./numdata;
    acopy.instrument=combineRepeatedInstruments(a.instrument,binno,nbins;k...)
    acopy.mask=falses(nbins); # the mask must be the same length as size(data,1)
    normalize!(acopy,norminfo...)
    return acopy
end
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX
# XXX ^^^^ This code needs to be rewritten to allow for >0 dimensional detectors ^^^^ XXX
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX

# Basic operator overloading
function *{T<:Scatterd}(a::T,b::T)
    nda=ndims(a.data)
    ndb=ndims(b.data)
#    @assert typeof(a)==typeof(b) "Both Scatterd objects must be of the same type"
    @assert nda==ndb "Data blocks must have the same dimensionality"
    @assert size(a.data,2:nda...)==size(b.data,2:ndb...) "Only first dimension size can differ between data blocks"
    @assert all(name.(getcolumns(a)).==name.(getcolumns(b))) "Column names must match for concatination"
    p=T(a) # FIXME this only keeps the field information of a, which is a problem.
    p.filename=a.filename*"\n"*b.filename
    p.command=a.command*"\n"*b.command
    p.header=[a.header;b.header] # this might be problematic
    p.data=cat(1,a.data,b.data)
    any(fieldnames(p).==:detector) && (p.detector=cat(1,a.detector,b.detector))
    p.mask=cat(1,a.mask,b.mask)
    p.instrument=cat(1,a.instrument,b.instrument)
    (a.definst==b.definst)||warn("combining two scans with different default instruments")
    return p
end
*{T<:Scatterd,R<:Scatterd}(::T,::R)=throw(TypeError(:*,"OperatorOverloading.jl",T,R))
(+){T<:Scatterd}(a::T,b::T)=(p=a*b;combineRepeatedPoints(p))
(^)(a::Scatterd,n::Integer)=(for i=2:n; a=a*a; end, a) # this is of dubious use

#(-){T<:Scatterd}(a::T,b::T)=a+(-b) # FIXME Subtraction can not be this simple.
# XXX using such a simplistic subtraction adds monitor counts, which is usually not what one
# XXX intends when determining the difference between two scans.
# XXX Instead, it is important that both Scatterd objects are first normalized in the same way
# XXX and then only the y values should be subtracted.
function (-){T<:Scatterd}(a::T,b::T)
    (ann,anv)=guessnormalization(a) # name and value of normalization for a
    (bnn,bnv)=guessnormalization(b) # name and value of normalization for b
    (ann==bnn)||warn("Different normalization columns in input, $ann and $bnn")
    bcopy=deepcopy(((ann==bnn)&&(anv==bnv))?b:normalize(b,ann,anv))
    # flip the intensity column
    setDat(bcopy,bcopy.y,-getDat(bcopy,bcopy.y))
    # and any fitting function intensity, if present
#    hasbeenfit(bcopy)&&(bcopy.fitres=-bcopy.fitres)

    # Both a and bcopy have anv counts in column ann. When they are combined this will
    # be a problem. Replace each with half.
    acopy=deepcopy(a)

    # the columns present in both objects must match for subtraction to work.
    # remove any columns present in only one object
    acol=name.(acopy.columns)
    bcol=name.(bcopy.columns)
    icol=intersect(acol,bcol)
    
    exta=[!any(x.==icol) for x in acol]
    extb=[!any(x.==icol) for x in bcol]
    any(exta) && delCol!(acopy,acol[exta])
    any(extb) && delCol!(bcopy,bcol[extb])

    # both objects now have the same columns, but their order might be different
    permuteCols(acopy, sortperm(name.(acopy.columns)) )
    permuteCols(bcopy, sortperm(name.(bcopy.columns)) )

    ax=getDat(acopy,acopy.x)
    bx=getDat(bcopy,acopy.x) # in case they don't have the same default x
    ii=fuzzyindexin(ax,bx) # same size/shape as ax with 0 anywhere ax is not in bx and a largest positive integer if ax is in bx
    any(ii.==0)&&warn("Extra points in positive object will be ignored. Maybe implement interpolation?")
    axbin=zeros(Int,size(ax))
    bxbin=zeros(Int,size(bx))
    j=1
    for i=1:length(ax) # the number of positive points
        if ii[i]>0 # this point has a matching negative point
            allmatching=map(x->isapprox(bx[ii[i]],x),bx) # it's possible multiple bx points are within error of this ax point
            axbin[i]=j
            bxbin[allmatching]=j
            j+=1
        end
    end
    nbins=maximum(axbin)
    outdata=Array{eltype(acopy.data)}(nbins,size(acopy.data,(2:ndims(acopy.data))...))
    for i=1:nbins
        tbxbin=bxbin.==i
        outdata[i,:]=sum(vcat(acopy.data[axbin.==i,:],sum(bcopy.data[tbxbin,:],1)/sum(tbxbin)),1)
    end
    # average everything *except* the y-column(s)
    outdata[:,.!matchBA(name.(getcolumns(acopy)),acopy.y)]/=2
    acopy.data=outdata;
    acopy.mask=a.mask[axbin.>0]
    return acopy
end

# Scatterd {+,-,*,/} Number
for y in [:+,:-,:*,:/]
    @eval $y(a::Scatterd,b::Number)=(c=deepcopy(a);setDat(c,c.y,$y(getDat(c,c.y),b)); hasbeenfit(c)&&(c.fitres=$y(c.fitres,b));c)
    #@eval $y{T<:Scatterd}(a::Array{T},b::Number)=(out=Array{T}(size(a)); for i=1:length(a); out[i]=$y(a[i],b); end; out)
end
# commutative operators
for y in [:+,:*]
    @eval $y(b::Number,a::Scatterd)=$y(a,b)
end
# I think this is OK. Only the y-values should be negated (not all counters) since negative monitor counts would be problematic
(-)(c::Scatterd)=(a=deepcopy(c);setDat(a,a.y,-getDat(a,a.y));hasbeenfit(a) && (a.fitres=-a.fitres);a)
(-)(b::Number,c::Scatterd)=b+(-c) # of dubious use

prod(a::Scatterd)=deepcopy(a);
function prod{T<:Scatterd{0}}(a::AbstractArray{T})
  ndm=map(x->ndims(x.data),a)
  nd=ndm[1]
  @assert all(x->nd==x,ndm) "Data blocks must all have the same dimensionality"
  szs=map(x->size(x.data,2:nd...),a)
  sz=szs[1]
  @assert all(x->sz==x,szs) "Only first dimension size can differ between data blocks"
  cns=map(x->name.(getcolumns(x)),a)
  cn=cns[1]
  samecolumns = [all(x.==cns[1]) for x in cns]
  if any(!,samecolumns)
      cantbefixed=false
      while !cantbefixed && (badidx=findfirst(.!samecolumns))>0
          # length(cns[1])==length(cns[badidx]) || (cantbefixed=true; continue)
          all(sort(cns[1]).==sort(cns[badidx])) || (cantbefixed=true; continue)
          # p1=sortperm(cns[1]) # cns[1][p1] puts cns[1] in order
          pb=sortperm(cns[badidx]) # cns[badidx][pb] puts cns[badidx] in the *same* order
          p1rev=sortperm(sortperm(cns[1])) # cns[1][p1][p1rev] === cns[1]
          # so cns[badidx][pb][p1rev] === cns[1] as well! And [pb[p1rev]] === [pb][p1rev]
          permuteCols(a[badidx],pb[p1rev])
          samecolumns[badidx]=true
      end
  end
  @assert all(samecolumns) "Column names must match for concatination"
  n=length(a) # number of input Scatterd objects
  m=prod(map(x->size(x.data,1),a)) # total number of points
  p=T(a[1])# FIXME this only keeps the field information of a, which is a problem.
  p.filename= join(map(x->x.filename,a),"\n")
  p.command= join(map(x->x.command,a),"\n")
  p.header=vcat(map(x->x.header,a)...)
  p.mask=cat(1,map(x->x.mask,a)...)
  p.data=cat(1,map(x->x.data,a)...)
  any(fieldnames(p).==:instrument) && ( p.instrument=cat(1,map(x->x.instrument,a)...) )
  return p
end
sum{T<:Scatterd}(a::AbstractArray{T};k...)=(p=prod(a);combineRepeatedPoints(p;k...))

# Need to define the next two lines (corretly) for the third and fourth to work:
+{T<:Scatterd,R<:Number}(::Range{T}, ::Range{R})=throw(InexactError())
+{T<:Scatterd,R<:Number}(::Range{R}, ::Range{T})=throw(InexactError())
+{T<:Scatterd,R<:Number}(::Range{T}, ::AbstractArray{R})=throw(InexactError())
+{T<:Scatterd,R<:Number}(::AbstractArray{R}, ::Range{T})=throw(InexactError())
+{T<:Scatterd,R<:Number}(a::AbstractArray{T},b::Range{R})=a+collect(b)
+{T<:Scatterd,R<:Number}(a::Range{R},b::AbstractArray{T})=collect(a)+b
+{T<:Scatterd,R<:Number}(a::AbstractArray{T},b::AbstractArray{R})=(size(a)==size(b) && (c=deepcopy(a); for i=1:length(a); c[i]=a[i]+b[i]; end; return c))
+{T<:Scatterd,R<:Number}(b::AbstractArray{R},a::AbstractArray{T})=(a.+b)
