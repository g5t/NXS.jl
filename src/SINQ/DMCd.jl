function generalPDLoad(filename::AbstractString,normfile::AbstractString="")
  isfile(filename) || error("Expected the path to a file, got '$filename' instead")
  sb=split(basename(filename),".")
  ext=length(sb)>1?sb[end]:""
  "hdf"==ext && error("Loading HDF files has not yet been implemented. Sorry")
  # now that we know the file doesn't have the HDF extension, try treating it
  # as a regular text file.
  lines=readlines(filename)
  # try to determine the filetype from the header
  if "DMC"==lines[1][1:3] # this is a normalized file from DMC
    sl=split(lines[1],",")
    title=strip(length(sl)>1?join(sl[2:end],","):"")
    header=deepcopy(lines[1:3])
    columns=[Column("a4",unit="°"),Column("cnts",unit="counts"),Column("m1",unit="counts")]
    par=Dict{AbstractString,Any}()
    lp=map(strip,split(replace(lines[2],r"=\s*","="),",")) # preserve the space in the Date field

    # lines[3] always contains a4_min, a4_step, a4_max, mon_counts
    # and optionally contails further key=value parameter pairs after a comma
    vals_pars=split(lines[3],",")
    (a4min,a4step,a4max,mon)=(map(parse,split(vals_pars[1]))...) # this should *always* work
    length(vals_pars)>1&&(lp=vcat(lp,map(strip,vals_pars[2:end]))) # add the extra parameters
    # the counts start from line[4] and are usually in 10 columns, but we don't
    # need to restrict ourselves to that case
    a4vals=collect(a4min:a4step:a4max)
    counts=Array{Float64}(length(a4vals)) # could be non-integer due to pre-normalization
    errors=similar(counts)
    monits=mon*ones(counts)
    j=1; nopts=length(a4vals); nolines=length(lines)
    i=4
    while j<nopts && i<=nolines && isa(map(x->parse(x,raise=false),split(lines[i])),Array)
      v=map(parse,split(lines[i]))
      n=length(v)
      counts[j:j+n-1]=v
      j+=n
      i+=1
    end
    j=1
    while j<nopts && i<=nolines && isa(map(x->parse(x,raise=false),split(lines[i])),Array)
      v=map(parse,split(lines[i]))
      n=length(v)
      errors[j:j+n-1]=v
      j+=n
      i+=1
    end
    while i<nolines
      lp=vcat(lp, map(strip,split(lines[i+=1],";")) ) # parameters at the end are separated by ; instead of ,
      header=vcat(header,lines[i])
    end
    for p in lp
      sp=split(p,"=")
      par[sp[1]]=isa(parse(sp[2],raise=false),Number)?parse(sp[2]):sp[2]
    end
    haskey(par,"dT")&&haskey(par,"T")&&(par["T"]=Measured(par["T"],par["dT"]^2);delete!(par,"dT"))
    haskey(par,"Date")&&(par["Date"]=DateTime(par["Date"],"'Y-m-d H:M:S'"))

    data=hcat( Measured(a4vals,1e-4), Measured(counts,errors.^2), Measured(monits,0) )

    return (title,header,columns,par,data)
  else
    error("Unknown instrument with first line\n\t$(lines[1])")
  end
end

#abstract Scatterd
#abstract Neutrond <: Scatterd
#abstract TripleAxisd <: Neutrond
#abstract SINQd <: TripleAxisd
type DMCd <: Diffractometerd{0}
# DMC is a banana-detector double-axis neutron diffractometer at SINQ
    filename::AbstractString
    command::AbstractString
    header::Array{AbstractString,1}
    columns::Array{Column,1}
    motors::Array{AbstractString,1}
    counters::Array{AbstractString,1}
    timers::Array{AbstractString,1}
    parameters::Dict{AbstractString,Any}
    data::Array{Measured,2}
    mask::BitArray
    instrument::Array{TwoAxis,1}
    definst::Function
    fitres::FitResult
    fits::Vector{FitResult}
    x::AbstractString
    y::AbstractString
    #t::Int

    DMCd()=new()
    DMCd(d::DMCd)=deepcopy(d)

    function DMCd(filename::AbstractString,INSTdef::Function=defineDMC)
        this=new()
        this.filename=filename
        (this.command, this.header, this.columns, this.parameters, this.data)=generalPDLoad(filename)
        columnnames=name.(this.columns) # this.columns is of type Array{Column,1}
        # alist=matchBA(columnnames,[r"^a[1-4]$",r"^[s,m,a]g[u,l]$"]) # angles (a1-a4,tilts)
        clist=matchBA(columnnames,"cnts").|matchBA(columnnames,r"^m[1-4]$") # counter and monitors
        # dlist=matchBA(columnnames,r"^[s,m,a]t[u,l]$") # translations
        # elist=matchBA(columnnames,r"^e[i,f,n]?$") # matches ei,ef,en,e
        # flist=matchBA(columnnames,"mf") # magnetic field is always(?) "mf"
        # klist=matchBA(columnnames,[r"^t[t,s]$",r"^temp?$"]) # non-exhaustive list of possible temperature "motors"
        # rlist=matchBA(columnnames,r"^q?[h,k,l]$") # matches qh,qk,ql,h,k,l
        tlist=matchBA(columnnames,"time")
        # wlist=matchBA(columnnames,r"^k[i,f]$") # ki, kf
        this.counters=columnnames[clist]
        this.timers=columnnames[tlist]
        this.motors=columnnames[.!(clist.|tlist)]
        # To start, all points should be unmasked (for fitting/plotting)
        this.mask=falses(size(this.data,1))
        # Determining default x and y columns is somewhat straightforward.
        this.y="cnts" # Default y should be "cnts" from generalPDLoad.
        this.x="a4" # and will (almost) always be a4
        this.definst=INSTdef
        fillinstrument!(this)
        this.fitres=FitResult{Float16,Float16,Float16,Measured,Int8}() # as a small(ish) placeholder
        this.fits=FitResult[] # initialize an empy array for the fits
        return this
    end
end
DMCd{T<:AbstractString}(fn::AbstractArray{T},INSTdef::Function=defineDMC)=map(x->DMCd(x,INSTdef),fn)
DMCload(o...)=DMCd(SINQdataPath("dmc",o...))

include("DMCd2TripleAxis.jl")
