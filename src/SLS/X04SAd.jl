function generalXYEload(filename::AbstractString)
    isfile(filename) || error("Expected the path to a file, got '$filename' instead")
    lines=filter(r"^\s*[0-9\.]+\s*[0-9\.]+\s*[0-9\.]+\s*$",readlines(filename)) # keep only lines which are NNN.NNN NNN.NNN NNN.NNN
    no=length(lines)
    x=Array{Float64}(no); y=similar(x); e=similar(y);
    for i=1:no
        xye=split(lines[i])
        x[i]=parse(xye[1])
        y[i]=parse(xye[2])
        e[i]=parse(xye[3])
    end
    return (x,y,e)
end

function X04SAparameters(filename::AbstractString)
    isfile(filename) || error("Expected a valid file path")
    isfile(filename*"_PAR")&&(filename*="_PAR")
    lines=readlines(filename)
    validlines=ismatch.(r"^#",lines) # invalid lines should be concatinated with their preceeding lines
    invalidlinenos=find(!,validlines)
    if !isempty(invalidlinenos) && all(invalidlinenos.>1)
        lines[invalidlinenos-1].*=lines[invalidlinenos]
        lines=lines[validlines]
    end
    lines=strip.(lines,Set(['#',' ']))
    outdict=Dict{AbstractString,Any}()
    for pair in lines
      sp=split(pair,":")
      length(sp)>2 && (sp=[sp[1];join(sp[2:end],":")]) # we only want to split off the key

      ssp=strip(sp[2],['"',' '])
      p=parse(ssp,raise=false)
      key=replace(strip(sp[1]),' ','_')
      if isa(p,Number)
        outdict[key]=p
      elseif isa(p,Expr) && :tuple == p.head && isa(eval(p),NTuple)
        outdict[key]=[eval(p)...]
      elseif isa(p,Expr) && :vect == p.head
        outdict[key]=eval(p)
      else # punt, just store the string
        outdict[key]=ssp
      end
    end

    if haskey(outdict,"Temp._device_X04SA-ES2-CRST")
        tempvals=split(outdict["Temp._device_X04SA-ES2-CRST"])
        outdict["Temp._device_X04SA-ES2-CRST"]=tempvals[1]
        outdict["Temp"]=Measured(parse(tempvals[2]),parse(tempvals[4])^2)
        outdict["No_Temps"]=parse(tempvals[5])
    end
    return outdict
end



#abstract Scatterd
#abstract XRayd <: Scatterd
#abstract XRayPowderd <: Xrayd
#abstract X04SAd <: XrayPowderd
type X04SAd # <: ???
# X04SA is a x-ray powder diffractometer at the SLS, but might be other things too?
    filename::AbstractString
    columns::Array{Column,1}
    motors::Array{AbstractString,1}
    counters::Array{AbstractString,1}
    parameters::Dict{AbstractString,Any}
    data::Array{Measured,2}
    mask::BitArray
    #instrument::Array{TwoAxis,1}
    #definst::Function
    fitres::FitResult
    fits::Vector{FitResult}
    x::AbstractString
    y::AbstractString

    X04SAd()=new()
    X04SAd(d::X04SAd)=deepcopy(d)

    function X04SAd(filename::AbstractString)
        this=new()
        this.filename=filename
        (tth,int,err)=generalXYEload(filename)
        nopts=length(tth)
        thedata=hcat(Measured(tth,1e-8),Measured(int,err.^2))
        this.data=Array{Measured}(nopts,2)
        this.data[:]=thedata[:]

        this.parameters=X04SAparameters(filename*"_PAR")

        this.columns=Column(["2θ","Intensity"])
        columnnames=name.(this.columns) # this.columns is of type Array{Column,1}
        clist=matchBA(columnnames,"Intensity")
        this.counters=columnnames[clist]
        this.motors=columnnames[.!clist]
        # To start, all points should be unmasked (for fitting/plotting)
        this.mask=falses(size(this.data,1))
        # Determining default x and y columns is somewhat straightforward.
        this.y="Intensity" # Default y should be "cnts" from generalPDLoad.
        this.x="2θ"
        this.fitres=FitResult{Float16,Float16,Float16,Measured,Int8}() # as a small(ish) placeholder
        this.fits=FitResult[] # initialize an empy array for the fits
        return this
    end
end
X04SAd{T<:AbstractString}(fn::AbstractArray{T})=map(x->X04SAd(x),fn)
# X04SAload(o...)=X04SAd(SLSdataPath("x04sa",o...))

function plot(p::X04SAd,x=p.x,y=p.y;pfit=true,kwargs...)
    xv=p.data[:,1];yv=p.data[:,2];m=p.mask
    any(m) && (Measureds.errorbar(xv[m],yv[m];masked=true,kwargs...))
    h=Measureds.errorbar(xv[.!m],yv[.!m];kwargs...)
    #pfit&&length(p.fits)>0&&(plotfit(p,x;kwargs...))
        xlabel(plainstring(p.columns[1]))
        ylabel(plainstring(p.columns[2]))
    return h
end
