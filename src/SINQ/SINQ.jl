""" A non-exhaustive list of file extensions used by various instruments at SINQ. Please extend as necessary."""
function SINQextension(instrument::AbstractString)
    ext=".dat"
    "eiger"==instrument && (ext=".scn")
    "rita2"==instrument && (ext=".hdf")
    return ext
end

"""
    SINQinstPath(instrument,number,year=[this year])
`SINQinstPath` returns the path to the local `instrument` copy of a specified
scan `number` for a given `year` (or this year, if omitted).
Since the local copy of a data file for a running scan is more complete than
its AFS copy, this allows plotting scans in progress if running on the instrument
computer.
"""
function SINQinstPath(instrument::AbstractString,no::Integer,year::Integer=Dates.year(now()))
    instrument=lowercase(instrument)
    file=joinpath("/home",instrument,"data","$year",pad(fld(no,1000),3),instrument*"$(year)n"*pad(no,6)*SINQextension(instrument))
end

"""
    SINQafsPath(instument,number,year=[this year])
By default `SINQafsPath` returns the path to the AFS copy of a specified scan
`number` for a given `year` (this year, if omitted).

By modifying `NXS.config["/sinq/afsroot"]` one can instead point this function
to a non-AFS copy of their data. For example, if you have scans 4023 through
4875 collected by TASP in 2016 stored in the folder
`C:\\data\\SINQ\\2016\\tasp\\004\\` you would need to put
`NXS.config["/sinq/datapath/"]="C:\\data\\SINQ"` in a script before loading those
scans with, e.g., `TASPload(4023:4875,2016)`.
"""
function SINQafsPath(instrument::AbstractString,no::Integer,year::Integer=Dates.year(now()))
    # This function determines the (now) standard locations for SICS instrument data files at SINQ
    # There is a high probability that it does not work for all instrument/year combinations.
    # For example, TASP first has datafiles in the AFS drive from 2000, but does not use the
    # current layout scheme until 2002 with 2003 being the first full-year of the current scheme.
    # This function can be extended to comply with earlier schemes, or a new function can be defined.
    # At this point in time, without a defined need for this functionality, the effort to extend
    # this function to all instrument/year cases seems wasteful and will not be done.
    instrument=lowercase(instrument)
    file=joinpath(config["/sinq/afsroot"],"$year",instrument,pad(fld(no,1000),3),instrument*"$(year)n"*pad(no,6)*SINQextension(instrument))
end
# For non-AFS (also non-PSI?) systems, allow the root location of
# all SINQ AFS data to be changed at runtime.
config["/sinq/afsroot"]=joinpath("/afs","psi.ch","project","sinqdata")

"""
    SINQdataPath(intrument,number,year)
Returns the full path to the specified scan `number` from the given `instrument` in
the given `year` (or the current year if omitted). If `julia` is running on the
instrument server this path points to the local filesystem copy, otherwise the
AFS path is returned (with the possibility to modify the AFS root, see the help for `SINQafsPath`)
"""
function SINQdataPath(inst::AbstractString,no::Integer,year::Integer=Dates.year(now()))
    hdt=map(lowercase,split(Base.gethostname(),"."))
    if length(hdt)==3 && "ch" == hdt[3] && "psi"==hdt[2] && lowercase(inst)==hdt[1]
        fdp=SINQinstPath(inst,no,year)
    else
        fdp=SINQafsPath(inst,no,year)
    end
    return fdp
end
# Dot notation "should" eliminate the need for this definition
# SINQdataPath(inst::AbstractString,no::AbstractVector{T},year::Integer=Dates.year(now())) where T === SINQdataPath.(inst,no,year) ?
SINQdataPath{T<:Integer}(inst::AbstractString,no::AbstractVector{T},year::Integer=Dates.year(now()))=map(x->SINQdataPath(inst,x,year),no)


function SINQload(instrument::AbstractString,o...)
    instrument=lowercase(instrument)
    #(ismatch(r"^rita",instrument))&&(RITAload(o...))
    (ismatch(r"^eiger",instrument))&&(EIGERload(o...))
    (ismatch(r"^tasp",instrument))&&(TASPload(o...))
end

# There are multiple problems with using Measured angle/QE values here. Drop them for the time being.
"""
    fillTASinstrument!(::SINQd)
can be used to fill the instrument field for any triple-axis instrument that refers
to the six requisite angles as A1,A2,A3,A4,A5,A6 and the components of (Q,E) as
qh,qk,ql,en.
"""
fillTASinstrument!(a::SINQd)=_fillTASinstrument!(a,["a1","a2","a3","a4","a5","a6"],["qh","qk","ql","en"])
# If we ever want to load in, e.g, SPiCE triple-axis data we can define
#  fillTASinstrument!(a::SPiCEd)=_fillTASinstrument!(a,["m1","m2","s1","s2","a1","a2"],["h","k","l","e"])
# We can define a similar function for any type of Triple axis instrument
function _fillTASinstrument!(a::TripleAxisd,anglecols::Vector{T},qecols::Vector{R}) where {T<:AbstractString,R<:AbstractString}
    @assert length(anglecols)==6 "A triple-axis needs six named angle columns to be defined"
    @assert length(qecols)==4 "A triple-axis needs four named columns to define a (Q,E) point"
    oneTAS=a.definst(a)
    nopt=size(a.data,1) # number of points in the scan object
    if mean(isCol.(a,anglecols)) >= mean(isCol.(a,qecols)) # use angles or (Q,E), preferring angles.
        headAs=getMultiple(a.variables,anglecols,0.) # the values in the header
        As=zeros(nopt,length(anglecols)) # without a Type, defaults to system-Float
        for i=1:length(anglecols)
            As[:,i]= isCol(a,anglecols[i])?getVal(a,anglecols[i]):headAs[i]
        end
        inst=setangles(oneTAS,permutedims(As/180*pi,(2,1))) # SINQ data files store degrees, TripleAxis objects want radians
    else
        headqes=getMultiple(a.variables,qecols,0.)
        QE=zeros(nopt,length(qecols))
        for i=1:length(qecols)
            QE[:,i]=isCol(a,qecols[i])?getVal(a,qecols[i]):headqes[i]
        end
        QE4v=LatticeVector(permutedims(QE,(2,1)),oneTAS.sample.crystal.rlat)
        inst=gotoQE(oneTAS,QE4v)
    end
    a.instrument=inst
end

include("TASPd.jl")
include("DMCd.jl")
include("EIGERd.jl")
include("RITA2.jl")
