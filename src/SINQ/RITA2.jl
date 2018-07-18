include("RITA2d.jl")
include("CAMEAd.jl")
function RITA2(filename,command,header,columns,dataval,variables,parameters,zeros,counters,timers,motors,reup,defINST)
    columnnames=name.(columns)
    #(Note, map(Measured,dataval) would return Array{MeasuredSymmetric} and we need the abstract Array{Measured})
    thedata=Array{Measured}(size(dataval)...); thedata[:]=dataval[:] # promote dataval to Array{Measured}
    cntrs=matchBA(columnnames,counters)
    any(cntrs) && (countingerror!(view(thedata,Colon(),cntrs)); addunit!(view(columns,cntrs),"counts"))
    # loop through the sets of Regex, uncertainty, Unit, Power for each known column name
    for (r,e,u,p) in reup
        m=matchBA(columnnames,r)
        any(m) && (replaceerror!(view(thedata,Colon(),m),e); addunit!(view(columns,m),u;power=p) )
    end
    # To start, all points should be unmasked (for fitting/plotting)
    mask=falses(size(thedata,1))
    # Determining default x and y columns is somewhat straightforward.
    y="cnts" # Default y should be "cnts" on RITA2.
    # Default x can be determined from the entered command:
    x=guessscanvariableILLTAS(header,command,parameters,columnnames,thedata)
    RITA2d(filename,command,header,columns,motors,counters,timers,parameters,variables,zeros,thedata,mask,x,y,defINST)
end
function CAMEA(filename,command,header,columns,dataval,flatcone,variables,parameters,zeros,counters,timers,motors,reup,defINST)
    # here we combine the data and flatcone information
    #(pointcol,pointdat,detectorcol,detectordat)=flat2CAMEA(columns,vcat(counters,timers),dataval,flatcone)
    (pointcol,pointdat,detectorcol,detectordat)=flat2CAMEA(columns,dataval,flatcone)
    push!(counters,"intensity")
    #
    mpd=Array{Measured}(size(pointdat)); mpd[:]=pointdat[:] # promote to Array{Measured}
    mdd=Array{Measured}(size(detectordat)); mdd[:]=detectordat[:] # promote to Array{Measured}
    cntrrgx=["cnts","m1","m2","intensity"]; cntrunt="counts"
    for (c,d,cc) in zip( (pointcol,detectorcol), (mpd,mdd), ( (Colon(),), (Colon(),Colon(),Colon()) ) )
        n=name.(c)
        m=matchBA(n,cntrrgx)
        any(m) && (countingerror!(view(d,cc...,m)); addunit!(view(c,m),cntrunt))
        for (r,e,u,p) in reup
            m=matchBA(n,r)
            any(m) && (replaceerror!(view(d,cc...,m),e); addunit!(view(c,m),u;power=p) )
        end
    end
    data=NArray(mpd,pointcol)
    detector=NArray(mdd,detectorcol)
    mask=falses(size(detector,1:3...)) # include all points, all a4 channels, all energy channels

    y="intensity"
    x1=guessscanvariableILLTAS(header,command,parameters,data)
    x1=="a3" || status(:hint,"Expected CAMEA scan variable a3, instead found $x1")
    x=vcat(x1,"a4")

    CAMEAd(filename,command,header,motors,counters,timers,parameters,variables,zeros,data,detector,mask,x,y,defINST)
end
function RITA2d(filename::AbstractString,defINST::Function=defineRITA2)
    isfile(filename) || (warn("The file $filename does not exist; returning nothing");return)
    (command, header, columns, dataval, flatcone, variables, parameters, zeros)=ILLTASLoad(filename)
    columnnames=name.(columns) # columns is of type Array{NamedArrays.Column,1}
    # handle counter columns separately since the variance is the counts
    cntrs=matchBA(columnnames,[r"^cnts$",r"^m[1-4]$"])        # detector and monitors
    # Now deal with the special column types
    tmrrgx=r"^c?time$"
    timrs=matchBA(columnnames,tmrrgx)
    counters=columnnames[cntrs]
    timers=columnnames[timrs]
    motors=columnnames[.!(cntrs.|timrs)]
    # setup Regex, uncertainties, unit-name and unit-power for known columns
    reup=[ ([r"^a[1-6]$",r"^[sma]g[ul]$",r"^xx$"],1e-2,  "°", 0), # angles
           ( r"^[sma]t[ul]$"             ,1e-2,  "m",-3), # translations
           ([r"^t[ts]$",r"^temp?$"]      ,1e-3,  "K", 0), # temperatures
           (r"^e[ifn]?$"                 ,7e-3, "eV",-3), # energies
           (r"^mag$"                     ,1e-3,  "T", 0), # magnetic field(s)
           (r"^q?[hkl]$"                 ,2e-3,"rlu", 0), # Qhkl, hkl
           (r"^k[if]$"                   ,2e-3,"Å⁻¹", 0), # k
           (tmrrgx                       ,1e-2,  "s", 0)  # timers
           ]
    if isempty(flatcone)
        return RITA2(filename,command,header,columns,dataval,variables,parameters,zeros,counters,timers,motors,reup,defINST)
    else
        return CAMEA(filename,command,header,columns,dataval,flatcone,variables,parameters,zeros,counters,timers,motors,reup,defINST)
    end
end
RITA2d{T<:AbstractString}(fn::AbstractArray{T},TASdef::Function=defineRITA2)=map(x->RITA2d(x,TASdef),fn)

RITA2path(no::Integer,path="/hzb/RITA2/data")=path*"/"*pad(no,6)
RITA2path(no::AbstractArray,o...;k...)=RITA2path.(no,o...;k...)

include("RITA2d2TripleAxis.jl")
