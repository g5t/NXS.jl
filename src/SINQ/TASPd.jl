TASPload(o...)=TASPd(SINQdataPath("tasp",o...))

#abstract Scatterd
#abstract Neutrond <: Scatterd
#abstract TripleAxisd <: Neutrond
#abstract SINQd <: TripleAxisd
type TASPd <: SINQd{0}
# TASP is a single-detector triple-axis neutron spectrometer at SINQ
    filename::AbstractString
    command::AbstractString
    header::Array{AbstractString,1}
    columns::Array{Column,1}
    motors::Array{AbstractString,1}
    counters::Array{AbstractString,1}
    timers::Array{AbstractString,1}
    parameters::Dict{AbstractString,Any}
    variables::Dict{AbstractString,Any}
    zeros::Dict{AbstractString,Any}
    data::Array{Measured,2}
    mask::BitArray
    instrument::Array{TripleAxis,1}
    definst::Function
    fitres::FitResult
    fits::Vector{FitResult}
    x::AbstractString
    y::AbstractString
    #t::Int

    TASPd()=new()
    TASPd(d::TASPd)=deepcopy(d)

    function TASPd(filename::AbstractString,TASdef::Function=defineTASP)
        this=new()
        this.filename=filename
        isfile(filename) || (warn("The file $filename does not exist; returning uninitialized structure");return)
        (this.command, this.header, this.columns, dataval, flatcone, this.variables, this.parameters, this.zeros)=ILLTASLoad(filename)
        isempty(flatcone) || warn("TASP flatcone information ignored. Please extend my code to match new hardware capabilities.")
        columnnames=name.(this.columns) # this.columns is of type Array{NamedArrays.Column,1}
        #(Note, map(Measured,dataval) would return Array{MeasuredSymmetric} and we need the abstract Array{Measured})
        thedata=Array{Measured}(size(dataval)...); thedata[:]=dataval[:] # promote dataval to Array{Measured}
        # setup Regex, uncertainties, unit-name and unit-power for known columns
        tmrrgx=r"^c?time$"
        reup=[ ([r"^a[1-6]$",r"^[sma]g[ul]$"],1e-2,  "°", 0), # angles
               ( r"^[sma]t[ul]$"             ,1e-2,  "m",-3), # translations
               ([r"^t[ts]$",r"^temp?$"]      ,1e-3,  "K", 0), # temperatures
               (r"^e[ifn]?$"                 ,5e-4, "eV",-3), # energies
               (r"^mf$"                      ,1e-3,  "T", 0), # magnetic field(s)
               (r"^q?[hkl]$"                 ,1e-4,"rlu", 0), # Qhkl, hkl
               (r"^k[if]$"                   ,1e-4,"Å⁻¹", 0), # k
               (tmrrgx                       ,1e-2,  "s", 0)  # timers
               ]
        # handle counter columns separately since the variance is the counts
        cntrs=matchBA(columnnames,[r"^cnts$",r"^m[1-4]$"])        # detector and monitors
        any(cntrs) && (countingerror!(view(thedata,Colon(),cntrs)); addunit!(view(this.columns,cntrs),"counts"))
        # loop through the sets of Regex, uncertainty, Unit, Power for each known column name
        for (r,e,u,p) in reup
            m=matchBA(columnnames,r)
            any(m) && (replaceerror!(view(thedata,Colon(),m),e); addunit!(view(this.columns,m),u;power=p) )
        end
        # Now deal with the special column types
        timrs=matchBA(columnnames,tmrrgx)
        this.counters=columnnames[cntrs]
        this.timers=columnnames[timrs]
        this.motors=columnnames[.!(cntrs.|timrs)]
        # And stash the Array{Measured} data in the object
        this.data=thedata
        # To start, all points should be unmasked (for fitting/plotting)
        this.mask=falses(size(this.data,1))
        # Determining default x and y columns is somewhat straightforward.
        this.y="cnts" # Default y should be "cnts" on TASP.
        this.x=guessscanvariableILLTAS(this.header,this.command,this.parameters,columnnames,this.data) # Default x can be determined from the entered command:
        this.definst=TASdef
        fillinstrument!(this)
        this.fitres=FitResult{Float16,Float16,Float16,Measured,Int8}() # as a small(ish) placeholder
        this.fits=FitResult[] # initialize an empty array for the fits
        return this
    end
end
TASPd{T<:AbstractString}(fn::AbstractArray{T},TASdef::Function=defineTASP)=map(x->TASPd(x,TASdef),fn)


"""
    TASP_stick_points(scan_numbers::Vector{Integer},stick_positions::Vector,[year::Integer])

Loads the input scans from their `scan_numbers` and `year` (defaulting to the current year),
then uses `one_from_many` to combine the scans with the `stick_postions` added as a new
column, which is set to the default independent variable before all scans are combined.
"""
function TASP_stick_points{I<:Integer}(scannumbers::AbstractVector{I},stickpos::AbstractVector,year::Integer=Dates.year(now()))
    @assert length(scannumbers)==length(stickpos)
    scans=TASPload(scannumbers,year)
    stickname=addunit(Column("stick"),"°")
    return one_from_many(scans,stickpos,stickname)
end

"""
    one_from_many(scans::Vector{Scatterd},values::Vector,named)

If many scans have been performed while a parameter that is not recorded in the scan
was changed, this function can be used to add the parameter value at each scan as an
additional `named` column (::String, or ::Column), change the default independent
variable to this parameter, and combine all scans into a single `Scatterd` output.
"""
function one_from_many{S<:Scatterd}(scans::AbstractVector{S},values::AbstractVector,newcol::Column=Column("newcol"))
    @assert length(scans)==length(values)
    colname=name.(newcol)
    for i=1:length(scans)
        scans[i]=addCol(scans[i],[values[i]],newcol)
        scans[i].x=colname
    end
    return prod(scans)
end
one_from_many(s::AbstractVector,v::AbstractVector,newcol::String)=one_from_many(s,v,Column(newcol))

export TASP_stick_points,one_from_many

include("TASPd2TripleAxis.jl")
