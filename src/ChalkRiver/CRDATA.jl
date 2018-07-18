#module CRDATA
#using FixedSizeArrays

rRunNo=r"^run\s*([0-9]+)"
type AxF{N,T<:AbstractString,R<:Integer}
    lines::Vector{T}
    runNos::NTuple{N,R}
    runPos::NTuple{N,R}
end
function AxF{T<:AbstractString}(lines::Vector{T})
    pos=find(map(x->ismatch(rRunNo,x),lines))
    nos=[parse(match(rRunNo,lines[x]).captures[1]) for x in pos]
    AxF(lines,(nos...),(pos...))
end
AxF(f::AbstractString)=isfile(f)?AxF(filter(!isempty,lowercase.(readlines(f)))):error("$f is not a file")# make everything lowercase ¯\_(ツ)_/¯


function getlines{N,T,R}(a::AxF{N,T,R},start::R,stop::R)
    @assert 0<start<=length(a.lines) && start<=stop<=length(a.lines)
    a.lines[start:stop]
end
getlines{N,T,R}(a::AxF{N,T,R},n::NTuple{2,R})=getlines(a,n[1],n[2]-1)
function getlines{N,T,R}(a::AxF{N,T,R},runno::Integer)
    idx=findfirst(x->x==runno,a.runNos) # this is likely always overkill and could be idx=runno :-/
    0==idx && (return Array{typeof(a.lines)}(0)) # run not found
    getlines(a, idx<N?a.runPos[idx+(0:1)]:(a.runPos[idx],length(a.lines)+1) )
end


abstract type CR end
type CRone{N} <: CR
    aof::AxF{N}
end
type CRtwo{N,M} <: CR
    aof::AxF{N}
    acf::AxF{M}
end
function showCR{N}(io::IO,a::CRone{N},compact::Bool=false)
    if compact
        Base.print(io,"CR:",N)
    else
        Base.print(io,"Chalk River data with $N scans")
    end
end
function showCR{N,M}(io::IO,a::CRtwo{N,M},compact::Bool=false)
    if compact
        Base.print(io,"CR:$N,$M")
    else
        Base.print(io,"Chalk River data with $N main and $M supplement scans")
    end
end
Base.show{T<:CR}(io::IO,a::T)=showCR(io,a,false)
Base.showcompact{T<:CR}(io::IO,a::T)=showCR(io,a,true)

function loadCR(aoffile::AbstractString)
    CRone(AxF(aoffile))
end
function loadCR(aoffile::AbstractString,acffile::AbstractString)
    CRtwo(AxF(aoffile),AxF(acffile))
end

getlines(a::CRone,o...)=getlines(a.aof,o...)
function getlines(a::CRtwo,o...)
    ol=getlines(a.aof,o...)
    cl=getlines(a.acf,o...)
    # As far as I can tell, the AOF file scan should have one more line than the ACF file scan
    # and the first 12 lines of both scans (the header) should be identical.
    if length(cl)==length(ol)-1 && length(cl)>12 && all(ol[1:12].==cl[1:12])
        head=ol[1:12] # hold onto the header from AOF
        ol=ol[13:end]
        cl=vcat(ol[1],cl[13:end]) # append the missing column names for the ACF file
        (o_ncl,o_ncn,o_lpp,o_nrl,o_ndl,o_np)=determineCRstructure(ol)
        (c_ncl,c_ncn,c_lpp,c_nrl,c_ndl,c_np)=determineCRstructure(cl)
        # in the only case that I'm aware of, o_ncl=c_ncl+1, but it may be possible for this to be more complicated
        @assert o_np==c_np "the main and supplement scans do not have the same number of points"
        @assert o_nrl==c_nrl "the main and supplement scans do not have equal number of repeated lines"
        columnlines=vcat(ol[1:o_ncl],cl[2:c_ncl]) # 2:c_ncl to avoid the appended missing column names
        ol=ol[o_ncl+1:end]
        cl=cl[c_ncl+1:end]
        # Ensure that we calculate the number of data lines in case a last partially-repeated
        # point is to be dropped. Recalling that the ACF file repeats the first line of the AOF file for each point
        datalines=Array{eltype(ol)}(o_lpp*o_np+(c_lpp-1)*c_np)
        ld=0;lo=0;lc=0;
        for i=1:o_np
            datalines[ld+=1]=ol[lo+=1];
            @assert cl[lc+=1]==datalines[ld] "mismatching lines! \n$(cl[lc])\n\t!=\n$(datalines[ld])"
            for j=1:o_nrl+1 # nrl is just the number of *repeats*, there is always 1 line per point
                for k=2:o_ncl; datalines[ld+=1]=ol[lo+=1];end
                for k=2:c_ncl; datalines[ld+=1]=cl[lc+=1];end
            end
        end
        out=vcat(head,columnlines,datalines)
    else
        out=ol
    end
    return out
end

function determineCRstructure{T<:AbstractString}(a::Vector{T})
    columnlineflag=map(x->ismatch(r"^point",x),a) # logical vector, true where line starts with Point
    no_column_lines=sum(columnlineflag) # a simple scan has 1 "Point" line
    no_column_names=map(x->length(split(replace(a[x],r"=\s+","="))),find(columnlineflag))-1 # number of column names in each column line, excluding "Point"
    point1flag=map(x->ismatch(r"^1\s",x),a) # true for every line that is "1 ..."
    # if one or more line(s) has the point and second column run together then point1flag may be, e.g., {FF...F}TTFT{FF...F}
    # we can check for this if diff(point1flag) has more than one +1 or -1 value
    dp1f=diff(vcat(false,point1flag)) # pad against TTFTFF
    if sum(dp1f.==1)>1||sum(dp1f.==-1)>1 # we may need to check both in case of FTTFT, FTTFTF
        point1flag[findfirst(dp1f.==-1):findlast(dp1f.==1)]=true # fill in the intermediate falses
    end
    lines_per_point=sum(point1flag) # should only be more than 1 if no_column_lines>1
    firstpoint=1; idx=findlast(x->ismatch(r"^1\s",x),a)
    if lines_per_point<=0
        idx=findlast(columnlineflag)
        firstpoint=parse(split(a[idx+=1])[1])
        lines_per_point=1
        samepoint=true
        while idx<length(a)&&samepoint
            samepoint=firstpoint==parse(split(a[idx+=1])[1])
            samepoint && (lines_per_point+=1)
        end
        !samepoint && (idx-=1) # if the last line is not the same point, step back for concatination check
#        info("Scan starts at point $firstpoint, with $lines_per_point lines per point")
    end
    # check to ensure we haven't missed any lines of the first point. (this happens if the temperature is in the second column and above 999.9999 K
    if idx<length(a) # idx points to the first line *after* the firstpoint lines
        no_digits=floor(Int,log10(firstpoint))+1 # 1-9>1, 10-99>2, etc.
        checkagainst=firstpoint%10 # we only need to check against the ones place
        issamepoint=true
        while issamepoint && idx<length(a)
            issamepoint=checkagainst==parse(Int,a[idx+=1][no_digits]) # the character to check. if 1-9 check that; if 10-99 check the second digit only
            issamepoint&&(lines_per_point+=1)
        end
    end
    times_repeated=no_column_lines>1?fld(lines_per_point-1,no_column_lines-1)-1:lines_per_point-1 # might be zero
    no_data_lines=length(a)-findlast(columnlineflag)
    no_points=fld(no_data_lines,(no_column_lines-1)*times_repeated+no_column_lines)
    float_no_points=no_data_lines/((no_column_lines-1)*times_repeated+no_column_lines)
    if float_no_points>no_points
        warn("Unfinished repeated point dropped")
    elseif float_no_points<no_points
        error("Something has gone drastically wrong.")
    end
    return (no_column_lines,no_column_names,lines_per_point,times_repeated,no_data_lines,no_points)
end

function parseCRMon{T<:AbstractString}(mons::Vector{T})
   mondict=Dict{AbstractString,Any}()
   idx=find(map(x->ismatch(r"^mon",x),mons))
   isempty(idx)&&(warn("""no monitors found in '$(join(mons," "))'""");return mondict)
   slc=hcat(idx[1:end],vcat(idx[2:end]-1,length(mons)))
   base=1
   mrgx=r"(mon[0-9]+)\[(.*)\]=([0-9]+)?\*([0-9]+)"
   for i=1:length(idx)
       mtc=match(mrgx,join(mons[slc[i,1]:slc[i,2]]))
       isa(mtc.captures[3],Void)||(base=parse(mtc.captures[3]))
       mondict[mtc.captures[1]*"/name"]=mtc.captures[2]
       mondict[mtc.captures[1]*"/value"]=base*parse(mtc.captures[4])
   end
   return mondict
end
function parseCRTemp{T<:AbstractString}(tmps::Vector{T})
    tmpdict=Dict{AbstractString,Int}()
    trgx=r"([a-z])sensr\(([0-9]+)\)"
    for i=1:length(tmps)
        mtc=match(trgx,tmps[i])
        isa(mtc,Void)?info("No temperature sensor in $(tmps[i])"):( tmpdict[ "temp/"*mtc.captures[1]*"/sensor" ] = parse(mtc.captures[2]) )
    end
    return tmpdict
end
function parseCRbragg(line::AbstractString)
    brgdict=Dict{AbstractString,Any}()
    brgx=r"([a-zA-Z]+)x\s*([a-zA-Z0-9]+)\s*\[\s*([0-9\.]+)\]\#\s*[0-9]*"
    mtchs=matchall(brgx,line)
    for i=1:length(mtchs)
        tmtch=match(brgx,mtchs[i])
        brgdict[tmtch.captures[1]*"/name"]=tmtch.captures[2]
        brgdict[tmtch.captures[1]*"/d"]=parse(tmtch.captures[3])
    end
    return brgdict
end

#parseCRdataline(line)=map(x->parse(Float64,x),split(line)[2:end]) # cut off the Point column

function parseCRdataline(line,no)
    sl=split(line)
#    print_with_color(:green,STDOUT,line,"\n")
    while length(sl)!=no+1
#        print_with_color(:green,STDOUT,"$(no+1)!=$(length(sl)): ")
#        print_with_color(:yellow,STDOUT,join(sl," "),"\n")
        if length(sl)<no+1 # look for two numbers that have merged due to fixed-width fields
            npr=map(x->tryparse(Float64,x),sl) # Array{Nullable{Float64},1} with #NULL for un-parseable values
            isn=map(isnull,npr)
            if any(isn)
                thenull=find(isn)
                for i=1:length(thenull)
                    j=thenull[i]
                    # now we start guessing
                    # some fields have 4 decimal places, these seem most likely to run into each other
                    ssl=matchall(r"[+-]?[0-9]+\.?[0-9]{0,4}",sl[j])
                    if all(x->~isnull(tryparse(Float64,x)),ssl)
                        sl=vcat(sl[1:j-1],ssl,sl[j+1:end])
                        thenull+=(length(ssl)-1)
                    end
                end
            else # no unparsable values, but still too few.
                # the most likely culprit (that I've seen) is a temperature above 999.9999 K in the
                # second column running into the point index column. Check for this first.
                rbigtemp=r"([0-9]+)([+-]?[0-9]{4}\.?[0-9]{0,4})"
                if ismatch(rbigtemp,sl[1])
                    btm=match(rbigtemp,sl[1])
                    sl=vcat(btm.captures[1],btm.captures[2],sl[2:end])
                else # maybe the longest string needs to be cut? This is a bad idea, generally
#                    info("$(sl[1]) did not match $rbigtemp")
                    lens=map(length,sl)
                    j=findfirst(lens.==maximum(lens))
                    ssl=matchall(r"[+-]?[0-9]+\.?[0-9]{0,4}",sl[j]) # e.g., "-7.000012394949" -> "-7.0000" "12394949"
                    if 1>=length(ssl) # two integers run together?
                        info("Treating $(sl[j]) as two integers run together")
                        h=fld(lens[j],2)
                        ssl=[sl[j][1:h],sl[j][h+1:end]]
                    end
                    if all(x->~isnull(tryparse(Float64,x)),ssl)
                        sl=vcat(sl[1:j-1],ssl,sl[j+1:end])
                    end
                end
            end
        else # too many columns?!
            error("too many values in Chalk River data line")
        end
    end
    map(x->parse(Float64,x),sl[2:end]) # always cut off the "point" column
end

function CRTASLoad(crdata::CR,run::Integer)
    # pulls CR TAS formatted scan number [run] from CR data structure
    # parses AOF (and ACF files) to return
    # the tuple ("file",command, header, columns, values, variables, parameters, zeros)
    # where command is the executed string
    #       header contains all non-data-block lines
    #       columns contains the column names
    #       values is a matrix of the data block
    #       variables is a dictionary of motor values (only non-scanned values are correct)
    #       parameters is a dictionary of sample/instrument parameters
    #       zeros is a dictionary of motor soft zeros
    #       "file" is a pseudo filename, e.g., C5A:BFMA13#150

    runlines=getlines(crdata,run) # combines AOF and ACF files, if present in crdata
    isempty(runlines) && (warn("Run $run not found in $crdata"); return)
    # getlines(CR) returns map(lowercase,getlines(CR.AOF&CR.ACF))

    # The first line in a CR run is *always*
    #   Run [0-9]+ Seq [0-9]+ Rec [0-9]+ File [0-9A-Z:]+ Date `d-u-y` `H:M:S.s`
    #
    # which are all (possibly) useful, stash them in the parameters dictionary
    parameters=Dict{AbstractString,Any}()
    sl=split(runlines[1])
    parameters[sl[1]]=parse(sl[2])
    parameters[sl[3]]=parse(sl[4])
    parameters[sl[5]]=parse(sl[6])
    parameters[sl[7]]=sl[8]
    parameters[sl[9]]=DateTime(join(sl[10:11],"T"),"d-u-yTH:M:S.s")
    # the second line is *always*(?)
    #  Mode [A-Z] Arm [0-9]+ Npts [0-9]+ Mon1\[\s*[0-9A-Z]*\]= [0-9]+ * [0-9]+ Mon2\[\s*[0-9A-Z]+\]=* [0-9]+
    sl=split(runlines[2])
    parameters[sl[1]]=parse(sl[2])
    parameters[sl[3]]=parse(sl[4])
    parameters[sl[5]]=parse(sl[6])
    # the third lines contains two more monitors (always?)
    parameters=merge(parameters,parseCRMon(vcat(sl[7:end],split(runlines[3]))))
    # the fourth line has some (floating point) parameters with equal signs
    sl=split(replace(runlines[4],"="," "))
    for i=1:2:length(sl); parameters[sl[i]]=parse(sl[i+1]); end
    # same with line 5
    sl=split(replace(runlines[5],"="," "))
    for i=1:2:length(sl); parameters[sl[i]]=parse(sl[i+1]); end
    # line 6 contains temperature information
    #   Temp: FX    Csensr(1)         Rsensr(3)         Temp=    1.0000; MAG =    152.0"
    sl=split(replace(replace(replace(runlines[6],"="," "),":"," "),";"," "))
    parameters["tempmode"]=sl[2]
    parameters=merge(parameters,parseCRTemp(sl[3:4])) # this could give issues if more or less sensors
    for i=5:2:length(sl); parameters[sl[i]]=parse(sl[i+1]); end
    # line 7: "Monox HE111   [ 3.43470]# 0     Anax HE111   [ 3.43470]# 0"
    parameters=merge(parameters,parseCRbragg(runlines[7]))
    # line 8 (probably) just has "Drv : VFHE=   3.702"
    sl=split(replace(join(split(runlines[8],":")[2:2:end]),"="," "))
    for i=1:2:length(sl); parameters[sl[i]]=parse(sl[i+1]); end
    # line 9 (might be) "Helm: HHCONST=  0.00 HHIHF=  4.000 MFIELD= 152.000 IVFT=  4.00 IVFB=  4.00"
    sl=split(replace(join(split(runlines[9],":")[2:2:end]),"="," "))
    for i=1:2:length(sl); parameters[sl[i]]=parse(sl[i+1]); end
    # 10: "Reserved for collimator, filter, specimen and detector data."
    # 11: "Use NAME command to store text on this line."
    # 12: "Use COMMENT command to store text on this line."

    # now, finally, we start to get variable data
    # first as one or more lines starting with "Point" that give us column names
    # to make life extra fun, the columns are wrapped *and* each point can be measured more than
    # once (with only parameters on the wrapped lines changing!
    line=13 # this should *always* be the first Column names line

    # number of: column lines, column names per column line, lines per point, repetitions per point, total data lines, points
    (nCL,nCN,nLPP,nRep,nDL,nPts)=determineCRstructure(runlines)
    # nLPP should only be more than one if nCL is greater than one as well
    1==nCL&&nLPP>1&&warn("CRTASLoad may not work for more than one line per point with only one column line.")

    nEPL=nCL-1 # number of extra Point lines
    columns=split(replace(replace(join(runlines[line+(0:nEPL)]," "),"point",""),r"=\s*","="))
    # now to start dealing with actual data

    zeroline=line+nEPL # advance line to just before the start of the data

    nCol=length(columns) # number of columns (excluding Point column)
    claimedPoints=get(parameters,"npts",nPts)
    if claimedPoints<nPts
        warn("Run header indicates $claimedPoints points but we seem to have $nPts points")
    elseif claimedPoints>nPts
        sl=split(runlines[end])
        lastnumber=parse(sl[1])
        stoppedearly=lastnumber==nPts
        nosteps=1==nPts&&( fld(claimedPoints-1,2)+1==lastnumber )
        stoppedearly||nosteps|| warn("Run header indicates $claimedPoints points but we only have $nPts points")
    end
    if nRep>0 # the more complicated case:
        values=Array{Float64}(nPts,1+nRep,nCol)
        l=0
#        if 1==nEPL
#            for i=1:nPts
#                firstvals=parseCRdataline(runlines[zeroline+(l+=1)])
#                for j=1:nRep+1
#                    remvals=parseCRdataline(runlines[zeroline+(l+=1)])
#                    @assert length(firstvals)+length(remvals)==nCol
#                    values[i,j,:]=vcat(firstvals,remvals)
#                end
#            end
#        else # more than one line per point :/
            for i=1:nPts
                firstvals=parseCRdataline(runlines[zeroline+(l+=1)],nCN[1])
                for j=1:nRep+1
                    remvals=parseCRdataline(runlines[zeroline+(l+=1)],nCN[2])
                    for k=2:nEPL;remvals=vcat(remvals,parseCRdataline(runlines[zeroline+(l+=1)],nCN[k+1]));end
                    @assert length(firstvals)+length(remvals)==nCol
                    values[i,j,:]=vcat(firstvals,remvals)
                end
            end
#        end
    else
        values=Array{Float64}(nPts,nCol)
        if 0==nEPL
            for i=1:nPts
                values[i,:]=parseCRdataline(runlines[zeroline+i],nCN[1]) #TODO make sure nCN is always a vector
            end
        else
            l=0
            for i=1:nPts
                firstvals=parseCRdataline(runlines[zeroline+(l+=1)],nCN[1])
                remvals=parseCRdataline(runlines[zeroline+(l+=1)],nCN[2])
                for k=2:nEPL;remvals=vcat(remvals,parseCRdataline(runlines[zeroline+(l+=1)],nCN[k+1]));end
                @assert length(firstvals)+length(remvals)==nCol
                values[i,:]=vcat(firstvals,remvals)
            end
        end
    end

#    @assert all(uniqueBA(columns)) "Some columns are repeated?!\n $columns"
    # filter repeated columns
    uc=uniqueBA(columns) # b is true for columns which should be kept, false otherwise
    columns=columns[uc] # use BitArray slicing to remove extra column names
    values=nRep>0?values[:,:,uc]:values[:,uc] # and their associated columns in values matrix

    # the returned columns variable should be of type Array{Column,1} which is a wrapper
    # class for the column name and its units.
    columns=Column(columns) # let the individual types set their units

    file=get(parameters,"file","crt:UNKNWN")*"#$run"
    command="""rec $(get(parameters,"rec",0)) seq $(get(parameters,"seq",0))"""
    header=runlines[1:zeroline] # everything before the data block
    return (file, command, header, columns, values, parameters)
end # CRTASLoad

export CR,loadCR,getlines,determineCRstructure,CRTASLoad
#end # module
