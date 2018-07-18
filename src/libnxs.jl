const libnxs = joinpath(Pkg.dir("NXS"),"deps","usr","lib","libnxs.so")
function maxfrequency{T<:Integer}(i::Array{T})
  d=Dict{T,T}(); o=one(T);
  for x in i; x>0 && (d[x]=haskey(d,x)?d[x]+o:o);end
  #return maximum(values(d))
  # if d is big, [fun](values(d)) is much faster than [fun](collect(values(d)))
  return (sort(collect(keys(d))),length(keys(d)),maximum(values(d)))
end
function targeted_sinks{T<:Integer}(idx::Array{T})
  d=Dict{T,T}(); o=one(T);
  for i in idx; i>0 && (d[i]=haskey(d,i)?d[i]+o:o); end # only keep track of valid destination indicies (>0)
  #return (sort(collect(keys(d))),length(keys(d))) # (targets,number_of_targets)
  return ( collect(keys(d)) ,length(keys(d)) ) # (targets,number_of_targets)
end

""" Wrapper for OpenMP C code to distribute pixels from one source into one sink. See `distribute_pixels`."""
function libnxs_dpp1!{I<:Integer,T<:Float64}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
  sink=zeros(T,sinksize...)
  (sourcecolon,sinkcolon,nosrcpix,nocolumns)=distribute_pixels_common(source,idx,totl,sink)
  (sinkwithcontrib,nosinkwithcontrib,contribpersink,maxcontrib)=maxfrequency(idx)
  ccall( (:dpp1,libnxs), Void,
         (Ptr{Cdouble}, Clonglong, Ptr{Clonglong},
          Ptr{Clonglong}, Clonglong, Ptr{Clonglong},
          Ptr{Cdouble}, Clonglong, Ptr{Clonglong},
          Clonglong),
        source, size(source,1), collect(sourcecolon), idx, maxcontrib, totl, sink, size(sink,1), collect(sinkcolon), nocolumns)
  return sink # totl and sink are modified curing ccall
end
""" Wrapper for OpenMP C code to distribute pixels from two sources into two sinks using the same asignment indicies. See `distribute_pixels`. """
function libnxs_dpp2!{I<:Integer,T<:Float64}(src1::Array{T},src2::Array{T},idx::Array{I},totl::Array{I},sinksize)
  @assert size(src1)==size(src2)
  snk1=zeros(T,sinksize...); snk2=zeros(T,sinksize...)
  (sourcecolon,sinkcolon,nosrcpix,nocolumns)=distribute_pixels_common(src1,idx,totl,snk1)
  (sinkwithcontrib,nosinkwithcontrib,contribpersink,maxcontrib)=maxfrequency(idx)
  ccall( (:dpp2,libnxs), Void,
         (Ptr{Cdouble}, Ptr{Cdouble}, Clonglong, Ptr{Clonglong},
          Ptr{Clonglong}, Clonglong, Ptr{Clonglong},
          Ptr{Cdouble}, Ptr{Cdouble}, Clonglong, Ptr{Clonglong},
          Clonglong),
         src1, src2, size(src1,1), collect(sourcecolon), idx, maxcontrib, totl, snk1, snk2, size(snk1,1), collect(sinkcolon), nocolumns)
  return (snk1,snk2) # totl, snk1 and snk2 are modified curing ccall
end
""" Wrapper for OpenMP C code to distribute pixels from three sources into three sinks using the same asignment indicies. See `distribute_pixels`."""
function libnxs_dpp3!{I<:Integer,T<:Float64}(src1::Array{T},src2::Array{T},src3::Array{T},idx::Array{I},totl::Array{I},sinksize)
  @assert size(src1)==size(src2)==size(src3)
  snk1=zeros(T,sinksize...); snk2=zeros(T,sinksize...); snk3=zeros(T,sinksize...);
  (sourcecolon,sinkcolon,nosrcpix,nocolumns)=distribute_pixels_common(srcval,idx,totl,snkval)
  (sinkwithcontrib,nosinkwithcontrib,contribpersink,maxcontrib)=maxfrequency(idx)
  ccall( (:dpp3,libnxs), Void,
        (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Clonglong, Ptr{Clonglong},
         Ptr{Clonglong}, Clonglong, Ptr{Clonglong},
         Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Clonglong, Ptr{Clonglong},
         Clonglong),
        src1, src2, src3, size(src1,1), collect(sourcecolon), idx, maxcontrib, totl, snk1, snk2, snk3, size(snk1,1), collect(sinkcolon), nocolumns)
  return (snk1,snk2,snk3) # totl, snk1, snk2 and snk3 are modified curing ccall
end

""" Wrapper for OpenMP C code to reduce pixels from one source into one sink. See `reduce_pixels`."""
function libnxs_rpp1!{I<:Integer,T<:Float64}(source::Array{T},idx::Array{I},totl::Array{I},sinksize)
  sink=zeros(T,sinksize...)
  (sourcecolon,sinkcolon,nosrcpix,nocolumns)=reduce_pixels_common(source,idx,totl,sink)
  (targetsinks,notargets)=targeted_sinks(idx)
  ccall( (:rpp1,libnxs), Void,
         (Ptr{Cdouble}, Ptr{Clonglong}, Clonglong, Clonglong,
          Ptr{Clonglong}, Ptr{Clonglong},
          Ptr{Cdouble}, Ptr{Clonglong},
          Ptr{Clonglong}, Clonglong),
         source, collect(sourcecolon), nosrcpix, nocolumns, idx-1, totl, sink, collect(sinkcolon), targetsinks-1, notargets)
  return sink # totl and sink are modified curing ccall
end
""" Wrapper for OpenMP C code to reduce pixels from two sources into two sinks using the same asignment indicies. See `reduce_pixels`."""
function libnxs_rpp2!{I<:Integer,T<:Float64}(src1::Array{T},src2::Array{T},idx::Array{I},totl::Array{I},sinksize)
  @assert size(src1)==size(src2)
  snk1=zeros(T,sinksize...); snk2=zeros(T,sinksize...)
  (sourcecolon,sinkcolon,nosrcpix,nocolumns)=reduce_pixels_common(src1,idx,totl,snk1)
  (targetsinks,notargets)=targeted_sinks(idx)
  ccall( (:rpp2,libnxs), Void,
         (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Clonglong}, Clonglong, Clonglong,
          Ptr{Clonglong}, Ptr{Clonglong},
          Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Clonglong},
          Ptr{Clonglong}, Clonglong),
         src1, src2, collect(sourcecolon), nosrcpix, nocolumns, idx-1, totl, snk1, snk2, collect(sinkcolon), targetsinks-1, notargets)
  return (snk1,snk2) # totl, snk1 and snk2 are modified curing ccall
end
""" Wrapper for OpenMP C code to reduce pixels from three sources into three sinks using the same asignment indicies. See `reduce_pixels`."""
function libnxs_rpp3!{I<:Integer,T<:Float64}(src1::Array{T},src2::Array{T},src3::Array{T},idx::Array{I},totl::Array{I},sinksize)
  @assert size(src1)==size(src2)==size(src3)
  snk1=zeros(T,sinksize...); snk2=zeros(T,sinksize...); snk3=zeros(T,sinksize...);
  (sourcecolon,sinkcolon,nosrcpix,nocolumns)=reduce_pixels_common(src1,idx,totl,snk1)
  (targetsinks,notargets)=targeted_sinks(idx)
  ccall( (:rpp3,libnxs), Void,
         (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Clonglong}, Clonglong, Clonglong,
          Ptr{Clonglong}, Ptr{Clonglong},
          Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Clonglong},
          Ptr{Clonglong}, Clonglong),
         src1, src2, src3, collect(sourcecolon), nosrcpix, nocolumns, idx-1, totl, snk1, snk2, snk3, collect(sinkcolon), targetsinks-1, notargets)
  return (snk1,snk2,snk3) # totl, snk1, snk2 and snk3 are modified curing ccall
end

# the CUDA versions of the reduce_pixels functions need the sizes of the source
# and sink arrays in order to allocate and copy data to the graphics card.
# since they already need nocolumns and nosrcpix, passing nosinkpix is enough
# extra information to calculate all array sizes
""" Wrapper for CUDA C code to reduce pixels from one source into one sink. See `reduce_pixels`."""
function libnxs_rpg1!{I<:Integer,T<:Float64}(src1::Array{T},idx::Array{I},totl::Array{I},sinksize)
  snk1=zeros(T,sinksize...); nosnkpix=prod(sinksize[1:end-1]);
  (sourcecolon,sinkcolon,nosrcpix,nocolumns)=reduce_pixels_common(src1,idx,totl,snk1)
  (targetsinks,notargets)=targeted_sinks(idx)
  ccall( (:rpg1,libnxs), Void,
         (Ptr{Cdouble}, Ptr{Clonglong}, Clonglong, Clonglong,
          Ptr{Clonglong}, Ptr{Clonglong},
          Ptr{Cdouble}, Ptr{Clonglong}, Clonglong,
          Ptr{Clonglong}, Clonglong),
         src1, collect(sourcecolon), nosrcpix, nocolumns, idx-1, totl, snk1, collect(sinkcolon), nosnkpix, targetsinks-1, notargets)
  return snk1 # totl and snk1 are modified curing ccall
end
""" Wrapper for CUDA C code to reduce pixels from two sources into two sinks using the same asignment indicies. See `reduce_pixels`."""
function libnxs_rpg2!{I<:Integer,T<:Float64}(src1::Array{T},src2::Array{T},idx::Array{I},totl::Array{I},sinksize)
  @assert size(src1)==size(src2)
  snk1=zeros(T,sinksize...); snk2=zeros(T,sinksize...); nosnkpix=prod(sinksize[1:end-1]);
  (sourcecolon,sinkcolon,nosrcpix,nocolumns)=reduce_pixels_common(src1,idx,totl,snk1)
  (targetsinks,notargets)=targeted_sinks(idx)
  ccall( (:rpg2,libnxs), Void,
         (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Clonglong}, Clonglong, Clonglong,
          Ptr{Clonglong}, Ptr{Clonglong},
          Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Clonglong}, Clonglong,
          Ptr{Clonglong}, Clonglong),
         src1, src2, collect(sourcecolon), nosrcpix, nocolumns, idx-1, totl, snk1, snk2, collect(sinkcolon), nosnkpix, targetsinks-1, notargets)
  return (snk1,snk2) # totl, snk1 and snk2 are modified curing ccall
end
""" Wrapper for CUDA C code to reduce pixels from three sources into three sinks using the same asignment indicies. See `reduce_pixels`."""
function libnxs_rpg3!{I<:Integer,T<:Float64}(src1::Array{T},src2::Array{T},src3::Array{T},idx::Array{I},totl::Array{I},sinksize)
  @assert size(src1)==size(src2)==size(src3)
  snk1=zeros(T,sinksize...); snk2=zeros(T,sinksize...); snk3=zeros(T,sinksize...); nosnkpix=prod(sinksize[1:end-1]);
  (sourcecolon,sinkcolon,nosrcpix,nocolumns)=reduce_pixels_common(src1,idx,totl,snk1)
  (targetsinks,notargets)=targeted_sinks(idx)
  ccall( (:rpg3,libnxs), Void,
         (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Clonglong}, Clonglong, Clonglong,
          Ptr{Clonglong}, Ptr{Clonglong},
          Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Clonglong}, Clonglong,
          Ptr{Clonglong}, Clonglong),
         src1, src2, src3, collect(sourcecolon), nosrcpix, nocolumns, idx-1, totl, snk1, snk2, snk3, collect(sinkcolon), nosnkpix, targetsinks-1, notargets)
  return (snk1,snk2,snk3) # totl, snk1, snk2 and snk3 are modified curing ccall
end
