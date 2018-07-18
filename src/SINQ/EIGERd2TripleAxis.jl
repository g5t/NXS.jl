# Sample information, collimations, and orientation can be read from the datafile.
# A more accurate TripleAxis object can be created when loading data files by first
# setting appropriate keys in the EIGERsetup dictionary, e.g., from your script files:
#
#         using NXS
#         NXS.config["eiger/analyzer/curvature/horizontal"]=NXS.flatBragg
#         scanobject=EIGERd(filename)
#
# would define TripleAxis objects for the points in the scan with flat analyzers (since
# the EIGER analyzer is only horizontally focusing to begin with)
#
# If one forgets to update alf[1-4] in SICS but changes the collimators, the correct
# value(s) can be forced into the TripleAxis objects by, e.g.,
#
#        using NXS
#        NXS.config["eiger/collimators/monochromator/sample"]=40; # arc-minutes
#        NXS.config["eiger/collimators/sample/analyzer"]=80; # arc-minutes
#        scanobject=EIGERd(filename)
#
# if the collimators in place were -/40'/80'/- but the SICSserver had [any other values]
#
# A last alternative is to create your own equivalent defineEIGER() function and pass it
# to the EIGERd constructor, e.g.,
#
# using NXS
# function mydefineEIGER()
# ...
# end
# defaultEIGERscanobject = EIGERd(filename)
# modifiedEIGERscanobject= EIGERd(filename,mydefineEIGER)

# config and get/check setup functions defined in Config.jl

function defaultconfigEIGER()
    src=Struct(); mono=Struct(); moni=Struct(); ana=Struct();
    det=Struct(); dist=Struct(); coll=Struct(); eiger=Struct()
    src[:width]= 80.  #mm
    src[:height]=150. #mm
    #
    mono[:mosaic]=pi/360. # 30' in radian
    mono[:width]=15*20.+ 14*1. # 15×20mm crystals wide + 14×1mm gaps
    mono[:height]=9*20.+  8*1. #  9×20mm crystals tall +  8×1mm gaps
    mono[:depth]=2. # 2mm thick crystals
    mono["/curvature/horizontal"]=horizontalyCurvedBragg
    mono["/curvature/vertical"]=  verticallyCurvedBragg
    #
    moni[:width]=  50. #mm
    moni[:height]=160. #mm
    #
    ana[:mosaic]=pi/270. # 40' in radian
    ana[:width]=7*25.+6*1. # 7×25mm crystals plus 6×1mm gaps wide
    ana[:height]=3*50. # 3×50mm crystals tall
    ana[:depth]=2. # 2mm thick
    ana["/curvature/horizontal"]=horizontalyCurvedBragg
    ana["/curvature/vertical"]=flatBragg
    #
    det[:width] = 20. #mm, most common width, though variable with often-used values 10, 20, 30, 50 mm
    det[:height]=160. #mm, height of standard detector mask
    #
    dist[:source_monochromator] =3250. # virtual source
    dist[:source_monochromator] =9400. # guide entrance window
    dist[:monochromator_sample] =2530.
    dist[:sample_analyzer]      =1090.
    dist[:analyzer_detector]    =1100.
    dist[:monochromator_monitor]=1500.
    #
    coll[:source_monochromator]=coll[:monochromator_sample]=coll[:sample_analyzer]=coll[:analyzer_detector]=800.
    #
    eiger[:source]=src; eiger[:monochromator]=mono; eiger[:monitor]=moni;
    eiger[:analyzer]=ana;  eiger[:detector]=det;       eiger[:distances]=dist;
    eiger["/collimators/horizontal"]=coll; eiger["/collimators/vertical"]=coll
    return eiger
end
reset_config_EIGER()=config[:eiger]=defaultconfigEIGER()
reset_config_EIGER()

# A general function to create a TripleAxis object representative of a typical EIGER setup.
function defineEIGER(a::EIGERd)
    deft=defaultconfigEIGER()
    eiger=config[:eiger]
    #
    src=get(eiger,:source,Struct())
    if !isa(src,Source)
        def=deft[:source]
        src=Source(get(src,:width,def[:width]),get(src,:height,def[:height]))
    end
    mono=get(eiger,:monochromator,Struct())
    if !isa(mono,Bragg)
        def=deft[:monochromator]
        tau   =get(mono,:tau   ,2pi/get(a.parameters,"dm",-1))
        named =get(mono,:name  ,tauname(tau))
        sense =get(mono,:sense ,get(a.parameters,"sm",-1))
        mosaic=get(mono,:mosaic,def[:mosaic])
        if isa(mosaic,Struct)
            horz=get(mosaic,"/horizontal",def[:mosaic])
            vert=get(mosaic,"/vertical"  ,def[:mosaic])
            mosaic=[horz,vert]
        end
        isa(mosaic,Array) || (mosaic=[mosaic,mosaic])
        width =get(mono,:width,def[:width])
        height=get(mono,:height,def[:height])
        depth =get(mono,:depth,def[:depth])
        extent=[width,height,depth] # overall w×h×t mm³; will be converted to second moments automatically
        cfh=get(mono,"/curvature/horizontal",def["/curvature/horizontal"])
        cfv=get(mono,"/curvature/vertical"  ,def["/curvature/vertical"]  )
        mono=Bragg(tau,named,sense,mosaic,extent,[cfh,cfv])
    end
    moni=get(eiger,:monitor,Struct())
    if !isa(moni,Detector)
        def=deft[:monitor]
        moni=Detector(get(moni,:width,def[:width]),get(moni,:height,def[:height]))
    end
    sample=get(eiger,:sample,Struct())
    if !isa(sample,Sample)
        o1=get(sample,"/orient/x",getMultiple(a.parameters,["ax","ay","az"],1))
        o2=get(sample,"/orient/y",getMultiple(a.parameters,["bx","by","bz"],1))
        if haskey(sample,:crystal) && isa(sample[:crystal],Crystal)
            crystal=sample[:crystal]
        else
            ld=get(sample,"/lattice/lengths",getMultiple(a.parameters,["as","bs","cs"],2pi))
            la=get(sample,"/lattice/angles", getMultiple(a.parameters,["aa","bb","cc"],90.))
            crystal=Crystal(ld...,la...)
        end
        mosaic=get(sample,:mosaic,Struct())
        isa(mosaic,Struct) && (mosaic=[get(mosaic,:horizontal,0.),get(mosaic,:vertical,0.)])
        isa(mosaic,Number) && (mosaic=[mosaic,mosaic])
        sense=get(sample,:sense,get(a.parameters,"ss",1))
        if haskey(sample,:geometry)&&isa(sample[:geometry],Lattices.SampleMesh)
            invAbsLen=get(sample,:inverse_absorption_length,0.)
            sample=Sample(crystal,o1,o2,sample[:geometry],mosaic,sense,invAbsLen)
        else
            shape=get(sample,:shape,Struct())
            if isa(shape,Struct)
                xx=get(shape,:x,100/12.) # <x*x> for a 10 mm cube
                yy=get(shape,:y,100/12.) # <y*y> for a 10 mm cube
                zz=get(shape,:z,100/12.) # <z*z> for a 10 mm cube
                shape=[xx,yy,zz]
            end
            isa(shape,Vector) && (shape=diagm(shape))
            shape=convert.(SecondMoment,shape) # a non-op if eltype(shape)==SecondMoment already
            sample=Sample(crystal,o1,o2,shape,mosaic,sense)
        end
    end
    analyzer=get(eiger,:analyzer,Struct())
    if !isa(analyzer,Bragg)
        def=deft[:analyzer]
        tau   =get(analyzer,:tau   ,2pi/get(a.parameters,"da",-1))
        named =get(analyzer,:name  ,tauname(tau))
        sense =get(analyzer,:sense ,get(a.parameters,"sa",-1))
        mosaic=get(analyzer,:mosaic,def[:mosaic])
        if isa(mosaic,Struct)
            horz=get(mosaic,"/horizontal",def[:mosaic])
            vert=get(mosaic,"/vertical"  ,def[:mosaic])
            mosaic=[horz,vert]
        end
        isa(mosaic,Array) || (mosaic=[mosaic,mosaic])
        width =get(analyzer,:width,def[:width])
        height=get(analyzer,:height,def[:height])
        depth =get(analyzer,:depth,def[:depth])
        extent=[width,height,depth] # overall w×h×t mm³; will be converted to second moments automatically
        cfh=get(analyzer,"/curvature/horizontal",def["/curvature/horizontal"])
        cfv=get(analyzer,"/curvature/vertical"  ,def["/curvature/vertical"]  )
        analyzer=Bragg(tau,named,sense,mosaic,extent,[cfh,cfv])
    end
    det=get(eiger,:detector,Struct())
    if !isa(det,Detector)
        def=deft[:monitor]
        det=Detector(get(det,:width,def[:width]),get(det,:height,def[:height]))
    end
    dist=get(eiger,:distances,Struct())
    if isa(dist,Struct)
        def=deft[:distances]
        s2m=get(dist,:source_monochromator,def[:source_monochromator])
        m2s=get(dist,:monochromator_sample,def[:monochromator_sample])
        s2a=get(dist,:sample_analyzer,def[:sample_analyzer])
        a2d=get(dist,:analyzer_detector,def[:analyzer_detector])
        m2m=get(dist,:monochromator_monitor,def[:monochromator_monitor])
        dist=[s2m,m2s,s2a,a2d,m2m]
    end
    isa(dist,Array) || error("config[\"/eiger/distances\"] should be either a Struct or an array of the five distances not $(typeof(dist))")
    isa(dist,Vector)||(dist=vec(dist))
    hcoll=get(eiger,"/collimators/hoizontal",Struct())
    if isa(hcoll,Struct)
        def=deft["/collimators/horizontal"]
        α1=get(hcoll,:source_monochromator,def[:source_monochromator])
        α2=get(hcoll,:monochromator_sample,def[:monochromator_sample])
        α3=get(hcoll,:sample_analyzer,def[:sample_analyzer])
        α4=get(hcoll,:analyzer_detector,def[:analyzer_detector])
        hcoll=[α1,α2,α3,α4]
    end
    isa(hcoll,Array)||error("config[\"/eiger/collimators/horizontal\"] should be either a struct or an array of the four collimations not $(typeof(hcoll))")
    isa(hcoll,Vector)||(hcoll=vec(hcoll))
    vcoll=get(eiger,"/collimators/vertical",Struct())
    if isa(vcoll,Struct)
        def=deft["/collimators/vertical"]
        β1=get(vcoll,:source_monochromator,def[:source_monochromator])
        β2=get(vcoll,:monochromator_sample,def[:monochromator_sample])
        β3=get(vcoll,:sample_analyzer,def[:sample_analyzer])
        β4=get(vcoll,:analyzer_detector,def[:analyzer_detector])
        vcoll=[β1,β2,β3,β4]
    end
    isa(vcoll,Array)||error("config[\"/eiger/collimators/vertical\"] should be either a struct or an array of the four collimations not $(typeof(vcoll))")
    isa(vcoll,Vector)||(vcoll=vec(vcoll))
    # determine the magnitude of the fixed incident or final energy
    # look for a value in fixed/k, fixed/kf, or fixed/ki defaulting to value in scan file
    f=get(eiger,:fixed,Struct())
    fixedk=get(f,:k,get(f,:kf,get(f,:ki,get(a.parameters,"kfix",2.662))))
    # allow for overriding fixed-k by a defined fixed-E, otherwise use fixed-k to calculate fixed-E
    fixede=get(f,:ei,get(f,:ef,get(f,:energy,ħ²2mₙ*fixedk^2)))

    # determine if ki or kf is fixed:
    fx=1!=get(a.parameters,"fx",1) # a.parameters[fx]==1 says kf-fixed; otherwise ki-fixed (backwards from TASP)
    # so, if this fx is true, then the TASMAD file says ki is fixed
    fi=haskey(f,:ki)&&isa(f[:ki],Real) # ki-fixed if true <--- possible override flags
    ff=haskey(f,:kf)&&isa(f[:kf],Real) # kf-fixed if true <-/
    # only override if (true,false) or (false,true) not (true,true) or (false,false)
    is_ki_fixed= xor(fi,ff) ? fi : fx # if xor(fi,ff); fi; else; fx; end

    tas=TripleAxis(src,mono,moni,sample,analyzer,det,dist,hcoll,vcoll,fixede,is_ki_fixed)
    return tas
end

fillinstrument!(a::EIGERd)=fillTASinstrument!(a)
