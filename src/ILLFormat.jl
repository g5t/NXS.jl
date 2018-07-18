
function ILLTASkeyLines2Dict{T<:AbstractString}(lines::AbstractArray{T},key::AbstractString)
    lines=map(strip,lines)
    outdict=Dict{AbstractString,Any}()
    keyregex=Regex("^$key(.*)") # captures the rest of the line
    keptvals=[match(keyregex,l).captures[1] for l in lines[find(ismatch.(keyregex,lines))]]
    keptvals=split(join(replace.(keptvals," ",""),","),",",keep=false) # Array{String}, e.g., ["DM=3.35","SM=-1",...]
    splitvals=split.(lowercase.(keptvals),"=") # [["dm","3.35"],["sm","-1"],...]
    for (dictkey,val) in [tuple(x...) for x in splitvals]
        try
            pval=parse(val)
            ptype=typeof(pval)
            if ptype<:AbstractFloat
                ival=round(Integer,pval)
                outdict[dictkey]= abs(pval-ival)<=eps(pval) ? ival : pval
            elseif ptype<:Number
                outdict[dictkey]= pval
            else
                outdict[dictkey]= val
            end
        catch problem
            isa(problem,ParseError) ? (outdict[dictkey]=NaN) : rethrow(problem)
        end
    end
    return outdict
end

function stringlines_to_matrix(l::AbstractArray{T},d::DataType=Float64,m::Integer=0) where T<:AbstractString
    n=length(l)
    ndims(l)>0 && n>0 || (return Array{d}(n,m))
    m<1 && (m=length(split(strip(l[1]))))
    mok = length.(split.(strip.(l))) .== m
    if all(mok)
        out = hcat([parse.(d,x) for x in split.(strip.(l))]...)'
    else
        out = stringlines_to_matrix(l[mok],d,m)
    end
    size(out)==(n,m) || status(:warn,"Expected to make Matrix of size ($n,$m) instead made $(size(out))?!")
    return out::Array{d,2}
end
stringlines_to_matrix(l::AbstractArray{Union{T,Void}},d::DataType=Float64,m::Integer=0) where T = Array{d}(0,0)

function ILLTASLoad(filename::AbstractString)
    # Loads ILL TAS formatted file "filename"
    # Returns the tuple (command, header, columns, values, variables, parameters, zeros)
    # where command is the executed string
    #       header contains all commented lines
    #       columns contains the column names
    #       values is a matrix of the data block
    #       variables is a dictionary of motor values (only non-scanned values are correct)
    #       parameters is a dictionary of sample/instrument parameters
    #       zeros is a dictionary of motor soft zeros
    isfile(filename) || return
    lines=readlines(filename)

    isComment(s) = s[6]==':'
    comments=filter(isComment, lines)

    # Find the command string which was executed
    commandline=findlast(ismatch.(r"^COMND:",comments)) # there should be only one, but take the last incase
    command=lowercase(strip(join(split(comments[commandline],":")[2:end],":"))) # join with : incase ":" is in command somehow


    # determine position of "DATA_:" key
    p=findlast(ismatch.(r"^DATA_:",lines))
    0==p && error("ILLTASLoad only works with ILL-TAS formatted files which have a line `DATA_:`")

    # "column" names comprise the next line.
    # strip off leading and trailing whitespace plus newline chacter with strip()
    # switch all names to lowercase for easier(?) later matching via lowercase()
    # split the string into an Array{SubString{AbstractString},1} at whitespace, via split()
    columns=split(lowercase(strip(lines[p+1])))

    # In ILL TAS format files the data block immediately follows the "columns" line
    # and may continue to the end of the file.
    alldatalines=lines[p+2:end]
    # If this is a flatcone data file there will be lines that being with "flat:"
    # followed by lines of " endflat" which will mess up my number parsing.
    freg=r"^flat:\s*(.*)"
    flatlines=[match(freg,x).captures[1] for x in alldatalines[ismatch.(freg,alldatalines)]]
    # Regular point lines will always(?) start with white space followed by an integer
    preg=r"^\s*([0-9]+.*)"
    pointlines=[match(preg,x).captures[1] for x in alldatalines[ismatch.(preg,alldatalines)]]
    #
    values=stringlines_to_matrix(pointlines,Float64,length(columns))
    size(values,1)>0 || status(:info,"No point lines with $(length(columns)) columns in file $filename")

    flat=stringlines_to_matrix(flatlines)

    # filter repeated columns
    b=uniqueBA(columns) # b is true for columns which should be kept, false otherwise
    b.&= columns .!= "pnt" # also only keep non "pnt" columns
    columns=columns[b] # use BitArray slicing to remove extra column names
    values=values[:,b] # and their associated columns in values matrix

    # the returned columns variable should be of type Array{Column,1} which is a wrapper
    # class for the column name and its units.
    columns=Column(columns) # let the individual types set their units

    #
    #keys=["INSTR","EXPNO","USER","LOCAL","FILE","DATE","TITLE","COMND","POSQE","STEPS","PARAM","VARIA","ZEROS"]
    variables=merge(ILLTASkeyLines2Dict(comments,"VARIA:"),ILLTASkeyLines2Dict(comments,"POSQE:"))
    parameters=merge(ILLTASkeyLines2Dict(comments,"PARAM:"),ILLTASkeyLines2Dict(comments,"STEPS:"))
    zeros=ILLTASkeyLines2Dict(comments,"ZEROS:")

    return (command, comments, columns, values, flat, variables, parameters, zeros)
end # ILLTASLoad

function guessscanvariableILLTAS(header::Array{R},command::AbstractString,parameters::Dict,data::NArray) where R<:AbstractString
    guessscanvariableILLTAS(header,command,parameters,name.(names(data)),array(data))
end
function guessscanvariableILLTAS(header::Array{R},command::AbstractString,parameters::Dict,columnnames::Array{S},data::Array) where {R<:AbstractString,S<:AbstractString}
    # FIXME this ignores the possibility that a word in the command could have a "d" and not be a step-size (e.g., a motor named "td")
    steps=ILLTASkeyLines2Dict(header,"STEPS:") # some instruments have d[motorname]=xxx in the STEPS: keys, others just have [motorname]=xxx :/
    #xhead=[keys(steps)...][ indmax( abs.([values(steps)...]) ) ] # if two motors have the same step size the one this picks will be the first, but Dicts aren't ordered so we don't know which one.
    xhead=[keys(steps)...][  abs.([values(steps)...]) .> 0 ] # instead keep an array of all motors with steps
    # some instruments (TASP) include d[qh,qk,ql] while others (FLEXX) only keep [qh,qk,ql]
    # remove the extra 'd' to aid matching, if it's present
    all(x->x[1]=='d',xhead) && (xhead=map(x->x[2:end],xhead))

    x=""
    firstd=findfirst(split(command,"").=="d")
    if firstd>0
        post=split(command[firstd:end])# the first d[motorname] in the command is post[1]
        while isempty(x) && length(post)>=4 # we need at least `[motorname] [value] d[motorname] [value]` to parse
            if post[1]=="dqh" # this could be 'dqh N N N N', 'dqh N N N', or 'dqh N dqk N dql N'
                hkle=["qh","qk","ql","en"]
                dhkle=getMultiple(parameters,"d".*hkle,0.) # pull step sizes from the file header
                if isapprox(sum(abs2,dhkle),0)
                    m=matchBA(columnnames,hkle)
                    sum(m)==4 || warn("some of $hkle are not present in columns")
                    Δhkle=[sum(diff(value.(data[:,matchBA(columnnames,x)]))) for x in hkle]
                    x=hkle[indmax(abs.(Δhkle))]
                else
                    xv=findfirst(abs.(dhkle).>0)
                    xv>0 && (x=hkle[xv])
                end
            else
                tc=post[1][2:end]
                abs(mean(diff(data[:,findfirst(columnnames.==tc)])))>0 && (x=tc)
            end
            if isempty(x) # remove motor/step specifications until the next d[motorname] is found
                # join the Array{AbstractString} with spaces, split the resultant AbstractString into individual
                # characters, find the first "d", and keep only from there on, then split on spaces
                jpost=join(post[2:end]," ")
                dpos=findfirst(split(jpost,"").=="d")
                post= (dpos>0)? (split(jpost[dpos:end])) : []
            end
        end
    end
    if isempty(x)
        x=xhead[1]
    else
        any(x.==xhead) || status(:info,"The determined step motor $x differs from header recorded step motor(s) $xhead")
    end
    isempty(x) && (warn("The scanned variable in  $command could not be determined."))
    return x
end
