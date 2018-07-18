#abstract Scatterd
#abstract Neutrond <: Scatterd
type WANDd <: Diffractometerd{1}
# WAND is a banana-detector double-axis neutron diffractometer at HFIR
    filename::AbstractString
    command::AbstractString
    header::Array{AbstractString,1}
    motors::Array{AbstractString,1}
    counters::Array{AbstractString,1}
    timers::Array{AbstractString,1}
    parameters::Dict{AbstractString,Any}
    data::NArray{Measured,2,Column}
    detector::NArray{Measured,3,Column}
    mask::BitArray
    instrument::Array{TwoAxisBanana,1}
    definst::Function
    fitres::FitResult
    fits::Vector{FitResult}
    x::Vector{AbstractString}
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
        (pointcol,pointdat,detectorcol,detectordat)=SPICE2WAND(thecolumns,spicedata)
        mpd=Array{Measured}(size(pointdat)); mpd[:]=pointdat[:] # promote to Array{Measured}
        mdd=Array{Measured}(size(detectordat)); mdd[:]=detectordat[:] # promote to Array{Measured}
        cntrrgx=["detector","monitor"] # the column "mcu" *is* a unit and should not be treated as a counter
        cntrunt="counts"
        tmrrgx=r"^c?time$"
                 #  angles                                translations                      temperatures                  timers
        knwnrgx=([r"^[mas][12]$",r"^sg[ul]$",r"^marc"], [r"^st[ul]$",r"^slit_[bflprt]+$"], r"^(vti|sample|temp?)[_ab]*$",tmrrgx)
        knwnerr=(1e-2,1e-2,1e-3,1e-2)
        knwnunt=("°","m","K","seconds")
        knwnpwr=(0,-3,0,0)
        for (c,d,cc) in zip( (pointcol,detectorcol), (mpd,mdd), ( (Colon(),), (Colon(),Colon()) ) )
            n=name.(c)
            m=matchBA(n,cntrrgx)
            any(m) && (countingerror!(view(d,cc...,m)); addunit!(view(c,m),cntrunt))
            for (r,e,u,p) in zip(knwnrgx,knwnerr,knwnunt,knwnpwr)
                m=matchBA(n,r)
                any(m) && (replaceerror!(view(d,cc...,m),e); addunit!(view(c,m),u;power=p) )
            end
        end
        allnames=name.(vcat(pointcol,detectorcol))
        cntrs=matchBA(allnames,cntrrgx)
        timrs=matchBA(allnames,tmrrgx)
        this.counters=allnames[cntrs]
        this.timers=allnames[timrs]
        this.motors=allnames[.!(cntrs.|timrs)]

        this.data=NArray(mpd,pointcol)
        this.detector=NArray(mdd,detectorcol)

        # To start, all points should be unmasked (for fitting/plotting)
        this.mask=falses(size(this.data,1)) # FIXME this should *probably* be based on the detector Array
        # Determining default x and y columns is somewhat straightforward.
        this.y="detector" # Default y is detector.
        this.x=["s1","s2"]
        this.definst=INSTdef
        fillinstrument!(this)
        this.fitres=FitResult{Float16,Float16,Float16,Measured,Int8}() # as a small(ish) placeholder
        this.fits=FitResult[] # initialize an empty array for the fits
        return this
    end
end
WANDd{T<:AbstractString}(fn::AbstractArray{T},INSTdef::Function=defineWAND)=map(x->WANDd(x,INSTdef),fn)
WANDload(o...)=WANDd(SPICEdataPath("wand",o...))

"""
    SPICE2WAND!(columns,data)
The SPICE data format was designed for 0D detector instruments while WAND has a 1D detector.
To shoehorn WAND data into a SPICE data file some (constant) information is discarded and all detector pixels are stored in a single "point".
This works for data storage but is undesirable for data visualization and manipulation.
This function transforms the 2D SPICE data block into a 3D WAND data block restoring information which was stripped in SPICE.
"""
function SPICE2WAND(cols::Vector{T},data::Array{R,2}) where {T<:Column,R}
    WANDa4offset=0:0.2:127.8;
    (npt,ncol)=size(data)
    cn=name.(cols)
    chans=matchBA(cn,r"^chan[0-9]+$") # true for columns which are the channels
    keepasis=.!(chans.|matchBA(cn,["twotheta","detector"]))
    ndet=sum(chans)
    chanidx=Array{Int64}(ndet)
    for i=1:ndet; chanidx[i]=findfirst(matchBA(cn,"chan$i")); end
    @assert all(!iszero,chanidx) "Could not find channel index for chan$(findfirst(chanidx.==0))"
    ttcol=findfirst(matchBA(cn,"twotheta"))
    detectordata=Array{R}(npt,ndet,2) # intensity and 2θ
    for j=1:npt,i=1:ndet
        detectordata[j,i,1]=data[j,chanidx[i]]
        detectordata[j,i,2]=data[j,ttcol]+WANDa4offset[i]
    end
    detectorcolumns=[Column("detector",unit="counts"),Column("s2",unit="°")]
    return (cols[keepasis],data[:,keepasis],detectorcolumns,detectordata)
end

# WAND is a banana-detector double-axis neutron diffractometer at HFIR
# WANDd0 holds data from a WANDd object in a single data block, this
# requires replicating 0-D information for each point and is very memory
# hungry as a result. This type should probably only be used by Slice, which
# removes as much unneccessary 0-D data as possible
type WANDd0 <: Diffractometerd{0}
    filename::AbstractString
    command::AbstractString
    header::Array{AbstractString,1}
    columns::Array{Column,1}
    motors::Array{AbstractString,1}
    counters::Array{AbstractString,1}
    timers::Array{AbstractString,1}
    parameters::Dict{AbstractString,Any}
    data::Array{Measured,3}
    #data::NArray{Measured,3,Column} # -- hopefully this works
    mask::BitArray
    #instrument::Array{TwoAxis,2}
    definst::Function
    fitres::FitResult
    fits::Vector{FitResult}
    x::AbstractString
    y::AbstractString
    #t::Int

    WANDd0()=new()
    WANDd0(d::WANDd)=deepcopy(d)
end
zerotype(::WANDd)=WANDd0

include("WANDd2TwoAxis.jl")
