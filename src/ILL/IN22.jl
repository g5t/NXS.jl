include("IN22d.jl")
#include("multiIN22d.jl")
function singleIN22(filename,command,header,columns,dataval,variables,parameters,zeros,counters,timers,motors,reup,defINST)
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
    y="cnts" # Default y should be "cnts" on IN22.
    # Default x can be determined from the entered command:
    x=guessscanvariableILLTAS(header,command,parameters,columnnames,thedata)
    IN22d(filename,command,header,columns,motors,counters,timers,parameters,variables,zeros,thedata,mask,x,y,defINST)
end
#function multiIN22(filename,command,header,columns,dataval,flatcone,variables,parameters,zeros,counters,timers,motors,reup,defINST)
#    # here we combine the data and flatcone information
#    #(pointcol,pointdat,detectorcol,detectordat)=flat2multiIN22(columns,vcat(counters,timers),dataval,flatcone)
#    (pointcol,pointdat,detectorcol,detectordat)=flat2multiIN22(columns,dataval,flatcone)
#    push!(counters,"intensity")
#    #
#    mpd=Array{Measured}(size(pointdat)); mpd[:]=pointdat[:] # promote to Array{Measured}
#    mdd=Array{Measured}(size(detectordat)); mdd[:]=detectordat[:] # promote to Array{Measured}
#    cntrrgx=["cnts","m1","m2","intensity"]; cntrunt="counts"
#    for (c,d,cc) in zip( (pointcol,detectorcol), (mpd,mdd), ( (Colon(),), (Colon(),Colon(),Colon()) ) )
#        n=name.(c)
#        m=matchBA(n,cntrrgx)
#        any(m) && (countingerror!(view(d,cc...,m)); addunit!(view(c,m),cntrunt))
#        for (r,e,u,p) in reup
#            m=matchBA(n,r)
#            any(m) && (replaceerror!(view(d,cc...,m),e); addunit!(view(c,m),u;power=p) )
#        end
#    end
#    data=NArray(mpd,pointcol)
#    detector=NArray(mdd,detectorcol)
#    mask=falses(size(detector,1:3...)) # include all points, all a4 channels, all energy channels
#
#    y="intensity"
#    x1=guessscanvariableILLTAS(header,command,parameters,data)
#    x1=="a3" || status(:hint,"Expected multiIN22 scan variable a3, instead found $x1")
#    x=vcat(x1,"a4")
#
#    multiIN22d(filename,command,header,motors,counters,timers,parameters,variables,zeros,data,detector,mask,x,y,defINST)
#end
function IN22d(filename::AbstractString,defINST::Function=defineIN22)
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
    reup=[ ([r"^a[1-6]$",r"^g[ulam]$"],1e-2,  "°", 0), # angles
           ( r"^[sma]t[ul]$"             ,1e-2,  "m",-3), # translations
           ([r"^t[ts]$",r"^temp?$"]      ,1e-3,  "K", 0), # temperatures
           (r"^e[ifn]?$"                 ,1.5e-2, "eV",-3), # energies
           (r"^mag$"                     ,1e-3,  "T", 0), # magnetic field(s)
           (r"^q?[hkl]$"                 ,2e-3,"rlu", 0), # Qhkl, hkl
           (r"^k[if]$"                   ,2e-3,"Å⁻¹", 0), # k
           (tmrrgx                       ,1e-2,  "s", 0)  # timers
           ]
    if isempty(flatcone)
        return singleIN22(filename,command,header,columns,dataval,variables,parameters,zeros,counters,timers,motors,reup,defINST)
    else
        error("This program was not aware that IN22 could use a flatcone analyzer!")
        #return multiIN22(filename,command,header,columns,dataval,flatcone,variables,parameters,zeros,counters,timers,motors,reup,defINST)
    end
end
IN22d{T<:AbstractString}(fn::AbstractArray{T},TASdef::Function=defineIN22)=map(x->IN22d(x,TASdef),fn)

IN22path(no::Integer,path=".")=joinpath(path,pad(no,6))
IN22path(no::AbstractArray,o...;k...)=IN22path.(no,o...;k...)

include("IN22d2TripleAxis.jl")
