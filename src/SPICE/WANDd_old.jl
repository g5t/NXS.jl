#abstract Scatterd
#abstract Neutrond <: Scatterd
type WANDd <: Diffractometerd
# WAND is a banana-detector double-axis neutron diffractometer at HFIR
    filename::AbstractString
    command::AbstractString
    header::Array{AbstractString,1}
    columns::Array{Column,1}
    motors::Array{AbstractString,1}
    counters::Array{AbstractString,1}
    timers::Array{AbstractString,1}
    parameters::Dict{AbstractString,Any}
    data::Array{Measured,3}
    mask::BitArray
    #instrument::Array{TwoAxis,2}
    definst::Function
    fitres::FitResult
    fits::Vector{FitResult}
    x::AbstractString
    y::AbstractString
    #t::Int

    WANDd()=new()
    WANDd(d::WANDd)=deepcopy(d)

    function WANDd(filename::AbstractString,INSTdef::Function=defineWAND)
        this=new()
        this.filename=filename
        (this.command, this.header, thecolumns, spicedata, this.parameters) = SPICEload(filename)
        # reshape the WAND data from SPICE format to 1-D detector format
        # do this now while the data are Float64 as it is *much* quicker than if the data are Measureds.Measured
        (this.columns,spicedata)=SPICE2WAND(thecolumns,spicedata)
        # spicedata is an Array{Float64} but we want Array{Measured}, so do some conversion fun
        this.data=Array{Measured}(size(spicedata)); this.data[:]=spicedata[:] # promote to Array{Measured}


        columnnames=name.(this.columns) # this.columns is of type Array{Column,1}
        alist=matchBA(columnnames,[r"^[mas][12]$",r"^sg[ul]$",r"^marc"])
        clist=matchBA(columnnames,["detector","monitor"]) # the column "mcu" *is* a unit and should not be treated as a counter
        dlist=matchBA(columnnames,[r"^st[ul]$",r"^slit_[bflprt]+$"]) # translations
        # elist=matchBA(columnnames,r"^e[i,f,n]?$") # matches ei,ef,en,e
        # flist=matchBA(columnnames,"mf") # magnetic field is always(?) "mf"
        klist=matchBA(columnnames,r"^(vti|sample|temp?)[_ab]*$") # non-exhaustive list of temperature columns
        # rlist=matchBA(columnnames,r"^q?[h,k,l]$") # matches qh,qk,ql,h,k,l
        tlist=matchBA(columnnames,r"^c?time$")
        # wlist=matchBA(columnnames,r"^k[i,f]$") # ki, kf
        this.counters=columnnames[clist]
        this.timers=columnnames[tlist]
        this.motors=columnnames[.!(clist.|tlist)]
        # Add units and uncertainties to known columns
        cc=(Colon(),Colon())
        any(clist) && (countingerror!(view(this.data,cc...,clist));     addunit!(view(this.columns,clist),"counts"))
        any(tlist) && ( replaceerror!(view(this.data,cc...,tlist),1e-2);addunit!(view(this.columns,tlist),"seconds"))
        any(alist) && ( replaceerror!(view(this.data,cc...,alist),1e-2);addunit!(view(this.columns,alist),"°"))
        any(dlist) && ( replaceerror!(view(this.data,cc...,dlist),1e-2);addunit!(view(this.columns,dlist),"m";power=-3))
        any(klist) && ( replaceerror!(view(this.data,cc...,klist),1e-3);addunit!(view(this.columns,klist),"K"))

        # To start, all points should be unmasked (for fitting/plotting)
        this.mask=falses(size(this.data,1))
        # Determining default x and y columns is somewhat straightforward.
        this.y="detector" # Default y should be "cnts" from generalPDLoad.
        this.x="s2" # and will (almost) always be s2
        this.definst=INSTdef
        #fillinstrument!(this)
        this.fitres=FitResult{Float16,Measured,Int8}() # as a small(ish) placeholder
        this.fits=FitResult[] # initialize an empty array for the fits
        return this
    end
end
WANDd{T<:AbstractString}(fn::AbstractArray{T},INSTdef::Function=defineWAND)=map(x->WANDd(x,INSTdef),fn)
WANDload(o...)=WANDd(SPICEdataPath("wand",o...))

const WANDa4offset=0:0.2:127.8;
"""
    SPICE2WAND!(columns,data)
The SPICE data format was designed for 0D detector instruments while WAND has a 1D detector.
To shoehorn WAND data into a SPICE data file some (constant) information is discarded and all detector pixels are stored in a single "point".
This works for data storage but is undesirable for data visualization and manipulation.
This function transforms the 2D SPICE data block into a 3D WAND data block restoring information which was stripped in SPICE.
"""
function SPICE2WAND(cols::Vector{T},data::Array{R,2}) where {T<:Column,R}
    (npt,ncol)=size(data)

    cn=name.(cols)
    chans=matchBA(cn,r"^chan[0-9]+$") # true for columns which are the channels
    keepasis=.!(chans.|matchBA(cn,["twotheta","detector"]))
    ndet=sum(chans)
    chanidx=Array{Int64}(ndet)
    for i=1:ndet; chanidx[i]=findfirst(matchBA(cn,"chan$i")); end
    @assert all(!iszero,chanidx) "Could not find channel index for chan$(findfirst(chanidx.==0))"
    ttcol=findfirst(matchBA(cn,"twotheta"))
    newdata=Array{eltype(data)}(npt,ndet,sum(keepasis))
    extdata=Array{eltype(data)}(npt,ndet,2) # detector, twotheta
    for j=1:npt,i=1:ndet
        newdata[j,i,:]=data[j,keepasis]
        extdata[j,i,1]=data[j,chanidx[i]]
        extdata[j,i,2]=data[j,ttcol]+WANDa4offset[i]
    end

    retdata=cat(3,newdata,extdata)
    retcols=vcat(cols[keepasis],Column(["detector","s2"]))
    return (retcols,retdata)
end

include("WANDd2TwoAxis.jl")
