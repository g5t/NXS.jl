

function ICELoad(filename::AbstractString)
    # Loads ICE formatted file "filename"
    # Returns the tuple (command, header, columns, values, variables, parameters, zeros)
    # where command is the executed string
    #       header contains all commented lines
    #       columns contains the column names
    #       values is a matrix of the data block
    #       variables is a dictionary of motor values (only non-scanned values are correct)
    #       parameters is a dictionary of sample/instrument parameters
    #       zeros is a dictionary of motor soft zeros
    lines=readlines(filename)
    # pre-data lines contain some useful information
    # however, the #Columns, #DetectorEfficiencies, and #ScanDescr lines should be dealt with independently
    plines=filter(s->(s[1]=='#'&&s[1:8]!="#Columns"&&s[1:10]!="#DetectorE"&&s[1:6]!="#ScanD"),lines)
    # turn the remaining pre-data lines into dictionary keys and values
    parameters=Dict{AbstractString,AbstractString}()
    paramregex=r"#([a-zA-Z]+).?\s*(.*)"
    all(ismatch.(paramregex,plines)) || error("The parameter Regex needs to be updated")
    for i=1:length(plines)
        rmc=match(paramregex,plines[i]).captures
        push!(parameters,lowercase(rmc[1])=>rmc[2])
    end
    # switch all names to lowercase for easier(?) later matching via lowercase()
    # split the string into an Array{SubString{AbstractString},1} at whitespace, via split()
    columns=split(filter(s->s[1:8]=="#Columns",lines)[1])[2:end]
    # any line without a leading # is a data block line
    dlines=filter(s->s[1]!='#', lines)
    # keep all lines starting with # as the header information
    hlines=filter(s->s[1]=='#', lines)
    # define an anonymous function to create parse the datablock, one line at a time
    myreplace(x,y,z)=(for (a,b) in zip(y,z); x=replace(x,a,b); end; x)
    gl2v(g,l) = map(x->parse(Float64,x),split(strip(myreplace(l,["BEO","BE","PG","IN","OUT","N/A"],["-3.","-2.","-1.","1.","-0.","NaN"])))[g])'
    bad=falses(columns)
    #for bn in ("AColMon","BColMon","BeFilMon","MgFilMon","PgFilMon","MCFX","cfxbe","cfxhopg","cfxmgf","FLIP","HKL")
    for bn in ("FLIP","HKL") # the replacement "fixes" the other bad columns, though perhaps not usefully
        bad.|=columns.==bn
    end
    good=.!bad
    columns=columns[good]
    values=Array{Float64}(length(dlines),length(columns))
    for i=1:length(dlines)
        values[i,:]=gl2v(good,dlines[i])
    end
    # the detector efficiencies are stored in the #DetectorEfficiencies line:
    deff=split(filter(s->s[1:10]=="#DetectorE",lines)[1])[2:end] # an Array{AbstractString,1} of [DetectorName]=[Float]
    counters=Dict{AbstractString,Float64}()
    for i=1:length(deff)
        kv=split(deff[i],"=")
        push!(counters,lowercase(kv[1])=>parse(Float64,kv[2]))
    end
    command=join(split(strip(filter(s->s[1:6]=="#ScanD",lines)[1]))[2:end]," ") # this is sort of the scan command
    # filter repeated columns
    b=uniqueBA(columns) # b is true for columns which should be kept, false otherwise
    columns=columns[b] # use BitArray slicing to remove extra column names
    values=values[:,b] # and their associated columns in values matrix
    return (command, Column(lowercase.(columns)), values, parameters, counters, hlines)
end # ICELoad

const macsa4offset=-76:8:76 # these are the a4 offsets for the 20 analyzers relative to "KIDNEY"
"""
MACS is a multi-analyzer cold neutron spectrometer at NCNR, NIST. It has 20 double-analyzers
at constant scattering-angle offset from each other (Δ2θ=8°). Each double-analyzer channel has
two detectors, one labeled "diff" placed in transmission behind the first analyzer and
the second labeled "spec" in the Bragg condition of the second analyzer.

MACS can be run in a triple-axis mode where one of the 20 detectors is positioned at a specified
(Q,E) point for each of a series of such points. The chosen detector can change during a scan to
overcome motor limits with an integer "ptai" specifying which analyzer channel was utilized for
each point in the scan.
This mode is used, e.g., while aligning a sample in the instrument. It is likely not used
during the remainder of an experiment, as doing so would be a gross misuse of the instrument.
As such, `MACSd` ignores the triple-axis mode columns and assumes that the mapping-mode is always used.
If one were so inclined, it would be straightforward to implement a triple-axis mode version of MACSd.
"""
type MACSd <: Neutrond{0} # FIXME this should be {1}
    filename::AbstractString
    command::AbstractString
    header::Array{AbstractString,1}
    columns::Array{Column,1}
    motors::Array{AbstractString,1}
    counters::Array{AbstractString,1}
    timers::Array{AbstractString,1}
    parameters::Dict{AbstractString,AbstractString}
    #data::Array{Union{MeasuredSymmetric,MeasuredAsymmetric},2} ### Measured_asym.jl
    data::Array{Measured,3}                                     ### Measured.jl
    mask::BitArray
    instrument::Array{TripleAxis,2}
    definst::Function
    fitres::FitResult
    x::Vector{AbstractString}
    y::AbstractString

    MACSd()=new()
    function MACSd(filename::AbstractString,MACSdef::Function=defineMACS)
        this=new()
        this.filename=filename
        (this.command, this.columns, macsdata, this.parameters, counterdict,this.header)=ICELoad(filename)
        # XXX                                              XXX
        # XXX Now deal with some of the MACS peculiarities XXX
        # XXX                                              XXX
        # XXX 1) remove PTAI, SPEC, DIFF, H, K, L, E, Qx, Qy, Qz columns since they all correspond to the triple-axis mode
        keep=.!matchBA(name.(this.columns),["ptai","spec","diff","h","k","l","e","ef","qx","qy","qz","a4","a5","a6"])
        macsdata=macsdata[:,keep] # keep only the non-tasmode columns
        this.columns=this.columns[keep]
        #
        columnnames=name.(this.columns)
        cntrnames=collect(keys(counterdict))
        isthetimer=cntrnames.=="time"
        this.timers =cntrnames[  isthetimer]
        this.counters=cntrnames[.!isthetimer]
        arecounters=matchBA(columnnames,cntrnames)
        this.motors=columnnames[.!arecounters]
        # filter for negative time/monitor/counts # sometimes MACS records -1 for various counters and then repeats the same point with real counts
        macsdata=macsdata[ squeeze(minimum(macsdata[:,arecounters],2).>=0, 2) ,:]

        # now that we've filtered the raw data as much as possible, convert to a Measured Array
        #(Note, map(Measured,dataval) would return Array{MeasuredSymmetric} and we need the abstract Array{Measured})
        thedata=Array{Measured}(size(macsdata)); thedata[:]=macsdata[:] # promote macsdata to Array{Measured}
        # we need to insert uncertainty for counters at this point *BEFORE* adjusting for detector efficiencies
        # so we might as well do everything else too
        clist=matchBA(columnnames,this.counters)
        tlist=matchBA(columnnames,this.timers)
        alist=matchBA(columnnames,[r"^a[1-6]$",r"^monblade",r"smpl[ul]tilt"]).|matchBA(columnnames,["basesampletheta","kidney"]) # angles
        dlist=matchBA(columnnames,r"^smpl[xyz]$") # translations
        elist=matchBA(columnnames,r"^e[if]?$") # matches ei,ef,e
        #flist=matchBA(columnnames,"mf") # magnetic field is always(?) "mf"
        klist=matchBA(columnnames,[r"^t[ts]$",r"^temp?$"]) # non-exhaustive list of possible temperature "motors"
        rlist=matchBA(columnnames,r"^[hkl]$") # h,k,l
        wlist=matchBA(columnnames,[r"^q[xyz]$",r"^k[if]$"]) #matches qx,qy,qz, ki, kf
        # countingerror! inserts variance equal to the counts; replaceerror! inserts symmetric variance equal to the passed error squared
        any(clist) && (countingerror!(view(thedata,Colon(),clist));      addunit!(view(this.columns,clist),"counts") )
        any(tlist) && ( replaceerror!(view(thedata,Colon(),tlist),1e-2); addunit!(view(this.columns,tlist),"seconds") )
        any(rlist) && ( replaceerror!(view(thedata,Colon(),rlist),1e-4); addunit!(view(this.columns,rlist),"rlu") )
        any(elist) && ( replaceerror!(view(thedata,Colon(),elist),5e-4); addunit!(view(this.columns,elist),"eV";power=-3) )
        any(alist) && ( replaceerror!(view(thedata,Colon(),alist),1e-2); addunit!(view(this.columns,alist),"°") )
        any(dlist) && ( replaceerror!(view(thedata,Colon(),dlist),1e-2); addunit!(view(this.columns,dlist),"m";power=-3) )
        any(klist) && ( replaceerror!(view(thedata,Colon(),klist),1e-3); addunit!(view(this.columns,klist),"K") )
        # any(flist) && ( replaceerror!(view(thedata,Colon(),flist),1e-4); addunit!(view(this.columns,flist),"T") )
        any(wlist) && ( replaceerror!(view(thedata,Colon(),wlist),1e-4); addunit!(view(this.columns,wlist),"Å⁻¹") )

        # now that true counting uncertainties are calculated,
        # normalize the counter intensities according to counterdict's efficiency values
        for k in this.counters # this.counters is a subset of the counterdict keys
            thedata[ : , find(columnnames.==k) ] ./= counterdict[k]
        end
        (npt,ncol)=size(thedata)
        ndet=length(macsa4offset)
        # Up to here macsdata/thedata are 2D. the following reshapes into a 3D matrix
        newdata=Array{eltype(thedata)}(npt,ndet,ncol-3*ndet) # 3*ndet -> spec##, diff##, analyzertheta##; 00<##<=20
        assd=Array{eltype(thedata)}(npt,ndet,5) #diff, spec, a4, a5, a6
        notdet=.!matchBA(columnnames,[r"^spec",r"^diff",r"^analyzertheta"])
        @assert sum(!,notdet)==3ndet
        speccol=Array{Int64}(ndet)
        diffcol=similar(speccol)
        analcol=similar(speccol)
        kdnycol=findfirst( matchBA(columnnames,"kidney") )
        for i=1:ndet
            speccol[i]=findfirst( matchBA(columnnames,["spec"*pad(i,2)]) )
            diffcol[i]=findfirst( matchBA(columnnames,["diff"*pad(i,2)]) )
            analcol[i]=findfirst( matchBA(columnnames,["analyzertheta"*pad(i,2)]) )
        end
        for j=1:npt,i=1:ndet
            newdata[j,i,:]=thedata[j,notdet]
            assd[j,i,1]=thedata[j,diffcol[i]]              # diffraction detector
            assd[j,i,2]=thedata[j,speccol[i]]              # spectroscopy detector
            assd[j,i,3]=thedata[j,kdnycol]+macsa4offset[i] # a4
            assd[j,i,4]=thedata[j,analcol[i]]              # a5
        end
        assd[:,:,5]=2assd[:,:,4] # a6 is forced to be 2*a5 for each analyzer
        # finish reshapeing the data and in doing so removed a large number of columns
        this.data=cat(3,newdata,assd)
        # now we can reduce the columns field to just those that remain:
        this.columns=vcat(this.columns[notdet],addunit.(Column.(["diff","spec","a4","a5","a6"]),["counts","counts","°","°","°"]))
        this.x=["a3","a4"]
        this.y="spec" #plot the Spectroscopic detectors by default

        # To start, all points should be unmasked (for fitting/plotting)
        this.mask=falses(size(this.data,1,2))

        rmcnt=matchBA(this.counters,[r"^spec[0-9]+",r"^diff[0-9]+"])
        rmmot=matchBA(this.motors,r"^analyzertheta[0-9]+")
        this.counters=vcat(this.counters[.!rmcnt],["diff","spec"])
        this.motors=vcat(this.motors[.!rmmot],["a4","a5","a6"])

        this.definst=MACSdef
        fillinstrument!(this) # fills in this.instrument with TripleAxis objects
        # use the TripleAxis utilities to calculate h,k,l, Ef, E?

        this.fitres=FitResult{Float16,Float16,Float16,Measured,Int8}() # as a (smallish) placeholder
        return this
    end
end
MACSd{T<:AbstractString}(fn::AbstractArray{T},o...;k...)=map(x->MACSd(x,o...;k...),fn)
