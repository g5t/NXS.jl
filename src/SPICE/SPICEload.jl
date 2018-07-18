# SPICE data files contain a number of "comment" lines as a header, a number of
# uncommented data lines, and a number of more commented lines as a footer.
# The distinction between header and footer is unimportant, it only exists within
# the data files due to the unavailability of footer information before the scan
# is finished.
# Within Scatterd objects this distinction does not exist.

function SPICEcomments2dict(lines::AbstractArray{T}) where T<:AbstractString
  lines=strip.(lines,Set(['#',' '])) # remove the comment # and its following space
  ewe = ismatch.(r"^[^=]*=$",lines) # true for lines that End With Equal signs (and don't have preceeding =)
  fewe=find(ewe)
  lines[fewe].*=lines[fewe+1] # concatinate the next line(s) for line(s) that end in =
  lines=lines[circshift(.!ewe,1)] # keep only lines with an = within them, circshift to get rid of the line *following* the one which ended with an =
  outdict=Dict{AbstractString,Any}()
  outdict["completed"]=false
  for pair in lines
    sp=split(pair,"=")
    length(sp)>2 && (sp=[sp[1];join(sp[2:end],"=")]) # we only want to split off the key
    if length(sp)==1
      sp=split(pair)
      if "completed."==sp[end]
        outdict["completed"]=true
        outdict["completed_at"]=join(sp[1:3],' ')
      else
        error("Unknown equal-less line: $pair")
      end
    else
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
  end
  return outdict
end

function addSPICEpeaks!(ubfile::AbstractString,dict::Dict{T}) where T<:AbstractString
  lines=filter(r"^.+$",readlines(ubfile)) #read the file, filtering out empty lines
  rsec=r"^\[([a-zA-Z]+)\]$"
  sectionheaders=find(ismatch.(rsec,lines)) # matches, e.g., [UBMode] and [Data]
  sections=Array{AbstractString}(length(sectionheaders))
  for i=1:length(sectionheaders)
      sections[i]=match(rsec,lines[sectionheaders[i]]).captures[1] # guaranteed to work
  end
  # We need to know which UBMode was in use (1 or 2. does 3 exist?)
  ubmodeheader=sectionheaders[findfirst("UBMode".==sections)]
  modeline=split(lines[ubmodeheader+1],"=")
  ubmode=parse("Mode"==modeline[1]?modeline[2]:"-1",raise=false)
  (isa(ubmode,Number)&&(3>ubmode>0))||(return)

  dataheader=findfirst("Data".==sections)
  datalines=lines[sectionheaders[dataheader]+1:sectionheaders[dataheader+1]-1]
  datadict=SPICEcomments2dict(datalines)
  if 1==ubmode
    # there should be "ScatteringPlaneVectors" and "Peak1" in the [Data] section
    haskey(datadict,"ScatteringPlaneVectors")&&haskey(datadict,"Peak1")||(return)
    uv=datadict["ScatteringPlaneVectors"]
    p1=datadict["Peak1"][1:3]
    p2=sum(abs,cross(p1,uv[1:3]))>0?uv[1:3]:uv[4:6]
  elseif 2==ubmode
    haskey(datadict,"Peak1")&&haskey(datadict,"Peak2")||(return)
    p1=datadict["Peak1"][1:3]
    p2=datadict["Peak2"][1:3]
  else # this shouldn't be possible
    return
  end
  dict["orientx"]=p1
  dict["orienty"]=p2
  return
end


function SPICEload(filename::AbstractString)
    isfile(filename) || return

    lines=readlines(filename)

    comments=filter(x->x[1]=='#', lines)
    databloc=filter(x->x[1]!='#', lines)

    parameters=SPICEcomments2dict(comments)
    # a bad hack to pull in UBconf information:
    #    this all assumes that we're dealing with a copy of the experiment folder created by SPICE
    #    this *will* fail if that assumption is wrong
    if haskey(parameters,"ubconf")
        expdir=splitdir(splitdir(filename)[1])[1] # should be /.../expXXX with first INST_expXXXX_scanYYYY.dat then Datafiles cut away
        ubconf=joinpath(expdir,"UBConf",parameters["ubconf"])
        isfile(ubconf) && addSPICEpeaks!(ubconf,parameters)
    end


    command=get(parameters,"command","unknown")
    columns=split(lowercase( parameters["col_headers"] ) ) # error if col_headers doesn't exist

    values=Array{Float64}(length(databloc),length(columns))
    for i=1:length(databloc); values[i,:]=parse.(Float64,split(databloc[i])); end

    # filter out repeated columns
    b=uniqueBA(columns)
    b.&=columns.!="pt." # the point number *is* (excluding writing errors) the index, so cut it out
    columns=columns[b]
    values=values[:,b]

    columns=Column(columns) # wrap the column names into Column objects

    return (command,comments,columns,values,parameters)
end
