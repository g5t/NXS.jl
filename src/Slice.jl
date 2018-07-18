export Slice,slice
"""
The `Slice` type is a container for `Scatterd` object data which has been
treated by the `slice` function.

| Fieldname | Description |
|:---------:|:------------|
| spec | The vector of binning specifications that produced this `Slice` |
| data | The binned and reshaped data from a `Scatterd` object |
| columns | A subset of the columns field of the source `Scatterd` object |
| source | The `Scatterd` object from which the `Slice` was created |

"""
type Slice{N,T}
    spec::Vector{BinSpec}              # the slicing specifications that made this Slice object
    data::Array{T,N}                   # the binned and re-shaped data; a 1D slice has 2D data to hold column information
    columns::Vector{Column}
    source::Scatterd                   # the full source object from which the slice came (should hopefully be a pointer)
end
function Slice{M,R,B<:BinSpec,A<:Column}(s::Vector{B},d::Array{R,M},c::Vector{A},src::Scatterd)
    #@assert sum(isbinning(s))==M-1 # in order to preserve the column information for each bin, a 2D slice has 3D data
    @assert (nobins(s)...)==(M>2?size(d,1:M-1...):(size(d,1),)) "Number of bins in spec vector, $((nobins(s)...)), do not match sizes of data array, $(size(d,1:M-1...))"
    @assert length(c)==size(d,M) "Passed $(length(cl)) column names, expected $(size(d,M))."
#    new(s,d,c,src)
    Slice{M,R}(s,d,c,src)
end
getcolumns(a::Slice)=copy(a.columns)

function Base.copy{N,T}(s::Slice{N,T})
    spec=deepcopy(s.spec)
    data=deepcopy(s.data)
    cols=deepcopy(s.columns)
    source=s.source # avoid copying the source object
    Slice{N,T}(spec,data,cols,source)
end

specstring(a::Slice)=strip(join(["$x" for x in a.spec]))

function Base.show(io::IO,p::Slice)
    #for x in p.spec; Base.showcompact(io,x);end; #Base.println(io)
    Base.showcompact(io,specstring(p))
    sd=join(["$x" for x in size(p.data)],"×")
    Base.print(io," "*"data (",sd,"): "*join(p.columns,' '))
end
Base.show(io::IO,::MIME"text/plain",p::Slice)=Base.print(io,"$(typeof(p)):\n",p)


"""
    slice(s::Scatterd, spec::Vector{BinSpec}; removeintdims=false)
The function `slice` takes a `Scatterd` object, a vector of `BinSpec` specifications, and
optionally a `Bool` to indicate if integrated dimensions should be removed from the
returned `Slice` object.

See help for `BinSpec` for details on the available binning specifications.
See help for `Slice` for details on the returned object type.
"""
function slice{T<:BinSpec}(s::Scatterd,spec::Vector{T};
                           removeintdims::Bool=false,
                           allcolumns::Bool=false,
                           extracolumns::Vector=isa(s.x,Vector)?s.x:[],
                           timer::Vector{Float64}=[time()],
                           outputlevel=0,verbose=false,
                           useccode::Bool=false)
    verbose && (outputlevel+=1)
    outputlevel>0 && calltimer(timer,"start slice")
    columns=allcolumns?name.(getcolumns(s)):typeof(s.y)[s.y;extracolumns] # typeof to protect against columns being Array{Any}
    nss=length(spec)
    for i=1:nss; isa(spec[i],BinStep) && (spec[i]=BinRange(s,spec[i])); end
    # make sure we're keeping the BinSpec cols
    bincols=[x.col for x in spec]
    for i=1:nss; any(columns.==bincols[i]) || push!(columns,bincols[i]); end
    outputlevel>0 && calltimersub(timer,"bin specification vector checked and column names added to columns vector")
    (normcol,normval)=guessnormalization(s) # guess the normalization values before binning
    any(columns.==normcol)||push!(columns,normcol) # keep the normalization column in addition to those specified in `columns`
    outputlevel>0 && calltimersub(timer,"intensity normalization column name determined"*(allcolumns?"":" and added to columns vector"))
    # since bin! overwrites the data block, we need to copy s
    p=copy_for_slice(s,columns) # to cut down on memory usage only specified columns
    outputlevel>0 && calltimersub(timer,"input Scatterd object copied with reduced columns")
    # Use the binning code to collect slice. Flags needed to change default behavior.
    bin!(p,spec;averagecounters=false,collapsedetectors=true,removeemptybins=false,excludemaskedpoints=true,refillinstrument=false,timer=timer,outputlevel=outputlevel,useccode=useccode)
    normalize!(p,normcol,normval) # normalize the counters and timers based on pre-bin guess
    outputlevel>0 && calltimersub(timer,"binned data normalized")
    # The detectors were combined during binning. p.data *should be* (prod(nobins(spec)),1,...,1,nocols)
    # With the first column being linearized, we can simply reshape the array wihtout squeezing first
    @assert size(p.data,1)==prod(nobins(spec)) "Wrong number of bins?! $(size(p.data,1))!=$(prod(nobins(spec)))"
    nobinscols=[nobins(spec);size(p.data,ndims(p.data))]
    slicedata=reshape(p.data,(nobinscols...)) # reshape from (N*M*L...,C) to (N,M,L,...,C) #TODO think about not doing this?
    outputlevel>0 && calltimersub(timer,"binned data reshaped to slice dimensions")
    if removeintdims
        singletons=nobinscols.==1
        spec=copy(spec[!singletons])
        cutdims=1:ndims(slicedata)
        cutdims=cutdims[singletons]
        slicedata=squeeze(slicedata,(cutdims...))
        outputlevel>0 && calltimersub(timer,"integrated slice dimensions removed from output slice data")
    end
    return Slice(spec,slicedata,p.columns,s) # return a Slice object with the binning specifications, sliced data, and Scatterd object
end
slice(s::Scatterd,t::Tuple...;k...)=slice(s,[BinSpec(s,x...) for x in t];k...)
slice(s::Scatterd,a::BinSpec...;k...)=slice(s,[a...];k...)

function slice(s::Slice,newspec::Vector{T};k...) where T<:BinSpec
    oldspec = s.spec
    newcols = getcolumn.(newspec)
    old2keep = [all(x.!=newcols) for x in getcolumn.(oldspec)] # those old columns which are not in newcols
    #
    spec = vcat(oldspec[old2keep],newspec)
    cols = getcolumn.(spec) # all names for binned columns in new slice
    extr = [all(x.!=cols) for x in name.(s.columns)] # possible extra columns to keep in ouput
    return any(extr) ? slice(s.source, spec; extracolumns=name.(s.columns)[extr]) : slice(s.source, spec)
end
slice(s::Slice,t::Tuple...;k...)=slice(s,[BinSpec(s.source,x...) for x in t];k...)
slice(s::Slice,a::BinSpec...;k...)=slice(s,[a...];k...)

# http://matplotlib.org/examples/color/colormaps_reference.html
# http://matplotlib.org/api/pyplot_summary.html?highlight=colormaps#matplotlib.pyplot.colormaps
import PyPlot: cm_get_cmap
cmap_viridis=cm_get_cmap("viridis")
cmap_inferno=cm_get_cmap("inferno")
cmap_plasma =cm_get_cmap("plasma")
cmap_magma  =cm_get_cmap("magma")
cmap_spectral=cm_get_cmap("nipy_spectral") # cmap spectral has been depricated :/
for i in (:cmap_viridis,:cmap_inferno,:cmap_plasma,:cmap_spectral,:cmap_magma)
    @eval $i[:set_under]("none")
    @eval $i[:set_over]("gray")
end

function grid_centers_to_corners{T}(x::AbstractVector{T})
  dx=diff(x)/2
  x=vcat(x[1,:]-dx[1,:],x+vcat(dx,dx[end,:]))
  return x
end
function grid_centers_to_corners{T}(x::AbstractArray{T,2})
  dx=diff(x,1)/2
  x=vcat( x[1:1,:]-dx[1:1,:], x+vcat(dx,dx[end:end,:]) )
  dx=diff(x,2)/2
  x=hcat( x[:,1:1]-dx[:,1:1], x+hcat(dx,dx[:,end:end]) )
  return x
end
function grid_centers_to_corners{T,N}(x::AbstractArray{T,N})
  f=[1:1,repeat([Colon()],outer=N-1)...]
  l=[N:N,repeat([Colon()],outer=N-1)...]
  for i=1:N
    d=diff(x,i)/2
    fi=circshift(f,i-1)
    fl=circshift(l,i-1)
    x=cat(i, view(x,fi...)-view(d,fi...), x + cat(i,d,view(d,fl...)))
  end
  return x
end
function gridneighbors{T}(x::AbstractArray{T,2},i::Integer)
  (N,M)=size(x)
  @assert 0<i<=N*M "The passed integer $i is not a valid linear index for a ($N,$M) matrix"
  if 1==i
    n=[1,N,N+1]
  elseif N==i
    n=[-1,N,N-1]
  elseif N*M == i
    n=[-N-1,-N,-1]
  elseif N*M - N + 1 == i
    n=[-N,-N+1,1]
  else
    i1=(i-1)%N+1
    i2=fld(i-1,N)+1
    if 1<i1<N && 1<i2<M
      n=[N+1,N,N-1,1,-1,-N+1,-N,-N-1]
    elseif 1<i1<N
      n=vcat([1,-1], 1==i2? [N+1,N,N-1] : [-N+1,-N,-N-1])
    elseif 1<i2<M
      n=vcat([N,-N], 1==i1? [-N+1,1,N+1] : [-N-1,-1,N-1])
    else
      error("Can not determine neighbors for ($N,$M) matrix at index $i -> [$i1,$i2]")
    end
  end
  return i.+n
end
function denangrid!{T}(x::Array{T,2},cycles::Int=1)
  for cycle=1:cycles
    !any(isnan,x) && return # do nothing if we have no NaNs
    # first, try to be clever assuming one direction will have <diff(x,i)^2>=0
    n1=.!isnan.(vec(sum(x,2))) # vec( (N,1) ) -> (N,)
    n2=.!isnan.(vec(sum(x,1))) # vec( (1,M) ) -> (M,)
    nnx=x[n1,n2] # no rows or columns that have *any* NaN values
    if !isempty(nnx)
      idx=find(isnan,x)
      (r,c)=ind2sub(size(x),idx)
      if sum(abs,diff(nnx,1))==0 # columns are all one value
        colvals=map( i->all(isnan,x[:,i])?NaN:mean(x[:,i][.!isnan.(x[:,i])]) ,1:size(x,2) )
        x[idx]=colvals[c]
      elseif sum(abs,diff(nnx,2))==0 # rows are all one value
        rowvals=map( i->all(isnan,x[i,:])?NaN:mean(x[i,:][.!isnan.(x[i,:])]) ,1:size(x,1) )
        x[idx]=rowvals[r]
      end
    end
    # if we're unlucky and <diff(x,i)^2>!=0 for any i, or a row/column was all NaN, keep going
    stillnan=find(isnan,x)
    !isempty(stillnan) && (oldx=deepcopy(x)) # to avoid screwing up complete missing rows/columns in the <diff(x,i)^2>=0 case
    for i in stillnan
      n=oldx[gridneighbors(x,i)]
      n=n[.!isnan.(n)] # avoid other NaN values
      !isempty(n) && (x[i]=mean(n))
    end
  end
end

# move plotting of slices to its own
include("SlicePlotting.jl")

import Base: log, log10, log1p
for x in (:log,:log10,:log1p)
    # @eval $x(s::Slice)=_function_on_slice_intensity($x,s)
    @eval $x(s::Slice,o...)=_function_on_slice($x,s,o...)
end
_function_on_slice(fn::Function,sl::Slice,col::AbstractString)=_function_on_slice(fn,sl,[col])
function _function_on_slice{N,S<:AbstractString}(fn::Function,sl::Slice{N},
         cols::AbstractVector{S}=name.(getcolumns(sl))[matchBA(name.(getcolumns(sl)),vcat(sl.source.counters,sl.source.timers))])
  out=copy(sl)
  npix=prod(size(out.data,1:N-1...))
  ofst=-npix
  indx=Array{Int}(length(cols)*npix)
  for col in cols
    c=findfirst(name.(getcolumns(out)).==col)
    indx[(1:npix)+(ofst+=npix)]=(1:npix)+(c-1)*npix
    out.columns[c]=addfunction(out.columns[c],fn) # modify the way that the column name is printed
  end
  out.data[indx]=fn.(out.data[indx])
  return out
end
