import Base: log, log10, log1p
for x in (:log,:log10,:log1p)
    @eval $x(s::Detector0d,o...)=_function_on_data($x,s,o...)
end
_function_on_data(fn::Function,sl::Detector0d,col::AbstractString)=_function_on_slice(fn,sl,[col])
function _function_on_data(fn::Function,d::Detector0d,cols::AbstractVector{S}=name.(getcolumns(d))[matchBA(name.(getcolumns(d)),vcat(d.counters,d.timers))]) where S<:AbstractString
  out=copy(d)
  for col in cols
    c=findfirst(name.(getcolumns(out)).==col)
    setDat(out,col,fn.(getDat(out,col)))
    out.columns[c]=addfunction(out.columns[c],fn) # modify the way that the column name is printed
  end
  return out
end
