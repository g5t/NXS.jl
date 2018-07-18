copy_for_slice(a::Scatterd{0},cols)=Base.copy(a;columns=cols)

copy_for_slice(a::Scatterd,cols::Vector{T}) where {T<:Column} = copy_for_slice(a,name.(cols))
function copy_for_slice(a::Scatterd,cols::Vector{T}) where {T<:AbstractString}
    out=zerotype(a)() # All Scatterd{>0} objects *should* have an equivalent Scatterd{0} type, which has an empty constructor
    fn=fieldnames(a) # all fieldnames
    ft=typeof.(getfield.(a,fn))
    as=[t<:AbstractString for t in ft]
    for i in find(as); setfield!(out,fn[i],identity(getfield(a,fn[i]))); end
    fn=fn[.!as]
    hnms=[:motors,:counters,:timers,:data,:detector] # to be dealt with separately
    mnms=[:fitres,:fits,:definst]
    hflg=matchBA(fn,hnms)
    mflg=matchBA(fn,mnms)
    easy=find(.!(hflg.|mflg))
    medi=find(mflg)
    for i in easy; setfield!(out,fn[i],copy(getfield(a,fn[i]))); end
    for i in medi; setfield!(out,fn[i],deepcopy(getfield(a,fn[i]))); end
    # now the hard part(s)
    datcols=name.(getnames(a.data))
    detcols=name.(getnames(a.detector))
    # ic was the columns of the input object
    oc=cols[isCol.(a,cols)]   # those columns which have been requested and that exist
    length(oc)>0 || info("None of the specified columns exist.")
    ocdat=a.data[matchBA(datcols,oc)]     # (Npoints,                Ndatcolsout)
    ocdet=a.detector[matchBA(detcols,oc)] # (Npoints,NdetX,NdetY,...,Ndetcolsout)
    dd=ndims(ocdet)
     # scale-up the 0D output data to match the ND detector data
    ocdat=repeat(ocdat,outer=(1,size(ocdet,2:dd-1...)...)) # repeat is overloaded to append the named (last) dimesion as 1 in outer
    newdata=cat(dd,ocdat,ocdet)
    newcols=getnames(newdata); newdata=getarray(newdata); # XXX bad :(
    setfield!(out,:columns,newcols)
    newtimers=copy(a.timers[matchBA(a.timers,oc)])
    newmotors=copy(a.motors[matchBA(a.motors,oc)])
    newcounters=copy(a.counters[matchBA(a.counters,oc)])
    n=Any[newmotors,newcounters,newtimers,newdata]
    for i=1:length(n); setfield!(out,hnms[i],n[i]); end
    return out
end


function Base.copy(a::Scatterd{0};columns=name.(getcolumns(a)))
    # FIXME: this doesn't necessarily appropriately copy vector properties
    # One could use deepcopy to ensure all outputs are divorced from all inputs
    # but that approach is *very* slow, so stick with copy for now ¯\_(ツ)_/¯
    out=typeof(a)() # All Scatterd objects have empty constructors
    fn=fieldnames(a) # get all fieldnames.
    ft=[typeof(getfield(a,f)) for f in fn] # since we can't copy AbstractString types anymore :(
    as=[t<:AbstractString for t in ft]
    for i in find(as); setfield!(out,fn[i],identity(getfield(a,fn[i]))); end
    fn=fn[.!as]
    hnms=[:data,:columns,:motors,:counters,:timers]
    mnms=[:fitres,:definst] # medium difficulty fields
    #enms=[:instrument,...] # FIXME the instrument field isn't actually easy to copy
    hflg=matchBA(fn,hnms) # a BitArray to indicate which of fn are hard to copy
    mflg=matchBA(fn,mnms)
    easy=find(.!(hflg.|mflg)) # the indicies of fn which are easy to copy
    medi=find(mflg) # indicies for fields which are moderately easy to copy
    for i in easy; setfield!(out,fn[i],copy(getfield(a,fn[i]))); end
    for i in medi; setfield!(out,fn[i],deepcopy(getfield(a,fn[i]))); end
    # now the hard part(s)
    ic=getcolumns(a) # those columns from the object
    oc=columns[isCol.(a,columns)]   # those columns which have been requested and that exist
    length(oc)>0 || info("None of the specified columns exist.")
    oci=getIdx.(a,oc) # the column indicies of the remaining column(s)
    ocn=name.(ic[oci]) # the name(s) of the remaining column(s)
    oki=falses(size(ic)...)
    oki[oci]=true # oki is true for columns which we want to keep

    newdata=slicedim(a.data,ndims(a.data),oki) # this makes a copy
    newcolumns=copy(a.columns[oki])
    newtimers=copy(a.timers[matchBA(a.timers,ocn)])
    newmotors=copy(a.motors[matchBA(a.motors,ocn)])
    newcounters=copy(a.counters[matchBA(a.counters,ocn)])
    n=Any[newdata,newcolumns,newmotors,newcounters,newtimers]
    for i=1:length(hnms)
        #od[hnms[i]]=Array{eltype(n[i]),ndims(n[i])}[n[i]]
        setfield!(out,hnms[i],n[i])
    end
    return out
end

const IntStrCol =  Union{Integer,AbstractString,Column}

# Column (Motor or Counter) verification:
isCol(p::Scatterd{0},i::Integer;k...)= i>0 && (i<=size(p.data,ndims(p.data)))
isCol(p::Scatterd,   i::Integer;k...)= i>0 && i<=size(p.data,ndims(p.data))+size(p.detector,ndims(p.detector))
isCol(p::Scatterd{0},n::AbstractString;add::Bool=false) = any(name.(getcolumns(p)) .== n) || ( add && checkandcalculate(p,n) )
isCol(p::Scatterd{D},n::AbstractString;add::Bool=false) where D = ( hasname(p.data,n) || (D>0 && hasname(p.detector,n)) ) || ( add && checkandcalculate(p,n) )
function checkandcalculate(p::Scatterd{D},n::AbstractString) where D
    # there is a p.instrument field that contains
    # information about the instrument from which Ei, Ef, E, a1, a2, a3, a4, a5, a6, h, k, l, ki, kf
    # can be calculated. for these values if there isn't already a same-named column
    # calculate the values, add them to p.data and p.columns and return true
    canbecalculated=["a1","a2","a3","a4","a5","a6","ei","ef","ki","kf","phi","h","k","l","e","incoi_h","incoi_k","incoi_l","incof_h","incof_k","incof_l","q"]
    withunitsof=["°","°","°","°","°","°","meV","meV","Å⁻¹","Å⁻¹","radian","rlu","rlu","rlu","meV","rlu","rlu","rlu","rlu","rlu","rlu","Å⁻¹"]
    idx=findfirst(canbecalculated.==n)
    idx > 0 || (return false)
    if      7 >idx; newcol=getangle(p.instrument,idx)
    elseif  7==idx; newcol=getEi(p.instrument)
    elseif  8==idx; newcol=getEf(p.instrument)
    elseif  9==idx; newcol=getki(p.instrument)
    elseif 10==idx; newcol=getkf(p.instrument)
    elseif 11==idx; newcol=getphi(p.instrument)
    elseif 12==idx; newcol=geth(p.instrument)
    elseif 13==idx; newcol=getk(p.instrument)
    elseif 14==idx; newcol=getl(p.instrument)
    elseif 15==idx; newcol=getE(p.instrument)
    elseif 16==idx; newcol=get_incoi_h(p.instrument)
    elseif 17==idx; newcol=get_incoi_k(p.instrument)
    elseif 18==idx; newcol=get_incoi_l(p.instrument)
    elseif 19==idx; newcol=get_incof_h(p.instrument)
    elseif 20==idx; newcol=get_incof_k(p.instrument)
    elseif 21==idx; newcol=get_incof_l(p.instrument)
    elseif 22==idx; newcol=norm.(getQ(p.instrument))
    else; return false; # this should probably be an error instead
    end
    7>idx && (newcol*=180/pi) # angles are stored in TripleAxis objects as radians; we typically want degrees for, e.g., plotting
    colname=Column(n;unit=withunitsof[idx])
    if D==0
        p.data=cat(ndims(p.data),p.data,newcol) # allow for cat(3,...) etc.
        p.columns=vcat(p.columns,colname)
        #addname(p.data,newcol,colname)
    else
        isdat = size(newcol)==size(p.data,1:ndims(p.data)...)
        isdet = size(newcol)==size(p.detector,1:ndims(p.detector)...)
        xor(isdat,isdet) || error("ambiguity in calculated column, it matches neither or both the data and detector field shapes")
        isdat?addname!(p.data,newcol,colname):addname!(p.detector,newcol,colname)
    end
    return true
end

getCol(p::Scatterd,i)= getcolumns(p)[getIdx(p,i)]

# Convert between column name and column index (string to int)
getIdx(p::Scatterd,i::Integer)= isCol(p,i) ? i : error("There is no column $i")
getIdx(p::Scatterd{0},n::AbstractString) = isCol(p,n) ? findfirst(name.(getcolumns(p)) .== n) : error("There is no column $n")
function getIdx(p::Scatterd,n::AbstractString)
    isCol(p,n) || error("There is no column $n")
    dati=SimpleNamedArrays.findindex(p.data,n)
    dati > 0 && (return dati)
    Base.size(p.data,ndims(p.data))+SimpleNamedArrays.findindex(p.detector,n)
end
# Functions to return Value, Variance, and Errors of any column
# {get|set}{Val|Var|Err} need to be updated if multi-dimensional detectors are to be supported
# allow for the possibility of getting multiple columns at once:
# getDat for Scatterd{N>0} requires logic for whether the requested column is part of the data block or detector block
function get_data_or_detector_column(p::Scatterd,i::Integer)
    @assert isCol(p,i) "$p \n has no data or detector column with index $i"
    get_data_or_detector_column(p,getcolumns(p)[i])
end
function get_data_or_detector_column(p::Scatterd,i)
    hasname(p.data,i) ? p.data[i] : hasname(p.detector,i) ? p.detector[i] : error("$p has no data or detector column with name $i")
end
function getDat(p::Scatterd,i::Integer)
    @assert isCol(p,i) "$p \n has no data or detector column with index $i"
    getDat(p,getcolumns(p)[i])
end
function getDat(p::Scatterd,i) # need to check p.data and p.detector. scale-up p.data to size of p.detector columns.
    isdata = hasname(p.data,i)
    dat=get_data_or_detector_column(p,i)
    if isdata
        detcoldim=1:ndims(p.detector)-1
        dat=repeat( dat,outer=fld.( Base.size(p.detector,detcoldim...),Base.size(dat,detcoldim...) ) )
    end
    return dat
end
function getDat(p::Scatterd,i::Vector) # need to check p.data and p.detector. scale-up p.data to size of p.detector columns.
    error("getDat Scatterd{D} is not yet defined for a vector of column indicies")
end
# and for Scatterd{N==0}
getDat(p::Scatterd{0},ni::Vector)= slicedim(p.data,ndims(p.data),getIdx.(p,ni))
getDat(p::Scatterd{0},i::Integer) = slicedim(p.data,ndims(p.data),i)
getDat(p::Scatterd{0},i) = slicedim(p.data,ndims(p.data),getIdx(p,i))
#
# legacy (useless?) functions:
getVal(p::Scatterd,i) = value.(getDat(p,i))
getVar(p::Scatterd,i) = variance.(getDat(p,i))
getErr(p::Scatterd,i) = uncertainty.(getDat(p,i))

# Functions to set
setDat(p::Scatterd{0},i::Integer,v) = isCol(p,i) ? setindex!(p.data,v,repeat([Colon()],outer=[ndims(p.data)-1])...,i) : throw(BoundsError())
function setDat(p::Scatterd,i::Integer,v,newcolname=isCol(p,i)?getcolumns(p)[i]:nothing) # need to check p.data and p.detector. potentially scale-up v to size of p.data or p.detector columns.
    @assert isCol(p,i) "Can not set the data for $i as it is not present in\n $p"
    oldcolname=getcolumns(p)[i]
    if hasname(p.data,oldcolname)
        replacename!(p.data,oldcolname,newcolname,v)
    elseif hasname(p.detector,oldcolname)
        replacename!(p.detector,oldcolname,newcolname,v)
    else # This shouldn't be possible thanks to assertion above
        error("Can not set non-existent $i for $p.")
    end
    return nothing
end
setVal(p::Scatterd,i::Integer,v) = setDat(p,i,v) # for consistency
setDat(p::Scatterd,n::AbstractString,v) = setDat(p,getIdx(p,n),v)
setVal(p::Scatterd,n::AbstractString,v) = setVal(p,getIdx(p,n),v)

# Overloaded functions to access and set elements of an object, obj, via obj[index] and obj[index]=newvalue, respectively
# These utilize getVal and setVal.
getindex(p::Scatterd,n) = getDat(p,n)
getindex(p::Scatterd,n,m::Integer...) = getindex(getDat(p,n),m...)
# setindex!(p::Scatterd,v,n) = setDat(p,n,v)
# setindex!(p::Scatterd,v,n,m::Integer...) = (vold=getDat(p,n);vold[m...]=v;setDat(p,n,vold))

# individual column add|del|move|exch(ange)
function addCol{T<:Measured}(p::Scatterd{0},v::Array{T},name::Column=Column("addedCol"),no::Integer=length(p.columns)+1)
    nd=ndims(p.data)
    ndv=ndims(v)
    repouter=ones(Int64,nd)
    matcheddim=false
    if 1==length(v) # scale up a single value
        repouter[1:nd-1]=[size(p.data,(1:nd-1)...)...]
        p.data=cat(nd,p.data, repeat([v[1]],outer=repouter))
        matcheddim=true
            # first case works unless v is a vector with size (length,) vs size(p.data,(1:1)...) returning an integer
    elseif (ndv==nd-1 && size(p.data,(1:nd-1)...)==size(v)) || (1==ndv&&2==nd&&size(p.data,1)==size(v,1))
        # the shape of v is appropriate to be pushed on
        p.data=cat(nd,p.data,v)
        matcheddim=true
    elseif ndv==nd-2 # v has one less dimension than (npts,ndet...) so it might be possible to scale it up
        for i=1:nd-1
            chkdims=vcat(1:i-1,i+1:nd-1) # [2,3,4,...], [1,3,4,...], [1,2,3,...], ...
            if size(p.data,chkdims...)==size(v)
                repouter[i]=size(p.data,i)
                p.data=cat(nd, p.data, repeat(v,outer=repouter) )
                matcheddim=true
            end
            matcheddim && break
        end
    elseif ndv==nd-3 # v has two less dimensions than (npts,ndet...) so it might be possible to scale it up
        for i=1:nd-2,j=1:nd-1
            chkdims=vcat(1:i-1,i+1:j-1,j+1:nd-1) # [3,4,...],[2,4,...],[2,3,...],...,[1,4,...],[1,3,...]
            if size(p.data,chkdims...)==size(v)
                repouter[[i,j]]=[size(p.data,i,j)...]
                p.data=cat(nd, p.data, repeat(v,outer=repouter) )
                matcheddim=true
            end
            matcheddim && break
        end
    elseif ndv==nd-4 # v has three less dimensions than (npts,ndet...) so it might be possible to scale it up
        for i=1:nd-3,j=1:nd-2,k=1:nd-1
            chkdims=vcat(1:i-1,i+1:j-1,j+1:k-1,k+1:nd-1) # [3,4,...],[2,4,...],[2,3,...],...,[1,4,...],[1,3,...]
            if size(p.data,chkdims...)==size(v)
                repouter[[i,j,k]]=[size(p.data,i,j,k)...]
                p.data=cat(nd, p.data, repeat(v,outer=repouter) )
                matcheddim=true
            end
            matcheddim && break
        end
    else
        error("Ambiguous size mismatch between p [size $(size(p.data))] and v [size $(size(v))]")
    end
    matcheddim || error("No up-scaling of v [size $(size(v))] is compatible with p [size $(size(p.data))]")
    push!(p.columns,name)
    moveCol(p,length(p.columns),no)
end
addCol{T<:Measured}(p::Scatterd,v::Array{T},name::AbstractString,o...)=addCol(p,v,Column(name),o...) # let just the new column name be specified
addCol(p::Scatterd,v::AbstractArray,o...)=addCol(p,collect(v),o...)
addCol(p::Scatterd,v::Array,o...)=addCol(p,Measured(v),o...)

function delCol!(p::Scatterd{0},v::AbstractArray) # don't specify Array type, let isCol/getIdx type dispatching do its job
    v=v[map(x->isCol(p,x),v)] # we can only remove columns that exist
    length(v)>0 || info("None of the specified columns exist.")
    ci=map(x->getIdx(p,x),v) # the column indicies of the remaining column(s)
    cn=p.columns[ci] # the name(s) of the remaining column(s)
    ki=trues(size(p.columns)...)
    ki[ci]=false # ki is true for columns which we want to keep
    tmpdata=slicedim(p.data,ndims(p.data),ki) # slice out just those columns that we'll keep
    sum(ki)==1 ? p.data=cat(ndims(p.data),tmpdata) : p.data=tmpdata # if only keeping one column, make sure the data dimensionality stays the same
    p.columns=getcolumns(p)[ki]
    cnn=name.(cn) # p.columns is now contains type Column instead of String. timers, motors, counters, all contain String instead
    p.timers=p.timers[.!matchBA(p.timers,cnn)]
    p.motors=p.motors[.!matchBA(p.motors,cnn)]
    p.counters=p.counters[.!matchBA(p.counters,cnn)]
    return p
end
delCol!(p::Scatterd,n)=delCol(p,[n]) # if n isn't an Array, make it one
delCol{T<:Scatterd}(p::T,o...) = (d=T(p); delCol!(d,o...); d)
function moveCol(p::Scatterd{0},orig,dest=length(p.columns))
    n=getIdx(p,orig); m=getIdx(p,dest);
    incr=(n<m)?-1:1
    # repeatedly exchange columns to reach desired state
    # using the origin column as a pseudo buffer
    # (slow but effective?)
    while n!=m
        exchCol(p,n,m)
        m+=incr
    end
    p
end
function exchCol(p::Scatterd{0},orig,dest=length(p.columns))
    n=getIdx(p,orig);
    m=getIdx(p,dest);
    d=ndims(p.data)
    siv=repeat( [Colon()], outer=[d-1] )
    datacopy=copy( slicedim(p.data,d,m) )
    colcopy=deepcopy(p.columns[m]);
    setindex!(p.data,copy( slicedim(p.data,d,n) ), vcat(siv,m)...)  # move origin to destination
    p.columns[m]=deepcopy(p.columns[n]);
    setindex!(p.data,datacopy,vcat(siv,n)...) # copy back destination to origin
    p.columns[n]=colcopy;
    p
end

function permuteCols(p::Scatterd{0},perm::Vector{T}) where T<:Integer
    N=length(p.columns)
    otN=1:N
    @assert all(sort(perm).==otN) "A permutation specification must consist of 1:length(columns), {likely permuted}"
    s=sortperm(perm) # the permutations that put perm in order
    loopcounter=0
    while any(s.!=otN) && loopcounter < N+1
        for orig=otN
            if orig!=s[orig] # only do something if we need to
                dest=s[orig] # the destination is where the permutation says it should go
                p=exchCol(p,orig,dest) # exchange this column with the one that should be here
                # and exchange the sorting permutations' entries
                s[orig]=s[dest] # move the destination's sorting permutation index to the origin's location
                s[dest]=dest # add the destination's index to its location (as it's made it where it wanted to go)
            end
        end
        loopcounter+=1
    end
    loopcounter<N+1 || warn("permuteCols did not finish permuting after $loopcounter iterations. Improve the algorithm.")
    return p
end

#getcolumns(a::Scatterd{0})=copy(a.columns)
getcolumns(a::Scatterd{0})=a.columns # getnames(a.data);
getcolumns(a::Scatterd)=vcat(getnames(a.data),getnames(a.detector))
