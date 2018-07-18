# A general function to create a TripleAxis object representative of a typical IN22 setup.
function defaultconfigIN22()
    src=Struct(); mono=Struct(); moni=Struct(); ana=Struct();
    det=Struct(); dist=Struct(); coll=Struct(); in22=Struct()
    src[:width]=60.
    src[:height]=125.
    #
    mono[:mosaic]=pi/270. # 40' in radian
    mono[:width]=140.
    mono[:height]=120.
    mono[:depth]=2. # 2mm thick crystals
    mono["/curvature/horizontal"]=flatBragg
    mono["/curvature/vertical"]=verticallyCurvedBragg
    #
    moni[:width]=70.
    moni[:height]=160.
    #
    ana[:mosaic]=pi/360. # 30' in radian
    ana[:width]=15*13.+14*1. #
    ana[:height]=125.
    ana[:depth]=2. # 2mm thick
    ana["/curvature/horizontal"]=horizontalyCurvedBragg
    ana["/curvature/vertical"]=flatBragg
    #
    det[:width] = 49. #mm, most common width, though variable with often-used values 10, 20, 30, 50 mm
    det[:height]=150. #mm, height of standard detector mask
    #
    dist[:source_monochromator] =1750.
    dist[:monochromator_sample] =1700.
    dist[:sample_analyzer]      =1150.
    dist[:analyzer_detector]    =1100.
    dist[:monochromator_monitor]=1000.
    #
    coll[:source_monochromator]= -2.
    coll[:monochromator_sample]=coll[:sample_analyzer]=coll[:analyzer_detector]=80.
    #
    in22[:source]=src; in22[:monochromator]=mono; in22[:monitor]=moni;
    in22[:analyzer]=ana;  in22[:detector]=det;       in22[:distances]=dist;
    in22["/collimators/horizontal"]=coll; in22["/collimators/vertical"]=coll
    return in22
end
reset_config_IN22()=config[:in22]=defaultconfigIN22()
reset_config_IN22()

# A general function to create a TripleAxis object representative of a typical IN22 setup.
function defineIN22(a::IN22d)
    deft=defaultconfigIN22()
    in22=config[:in22]
    #
    src=get(in22,:source,Struct())
    if !isa(src,Source)
        def=deft[:source]
        src=Source(get(src,:width,def[:width]),get(src,:height,def[:height]))
    end
    mono=get(in22,:monochromator,Struct())
    if !isa(mono,Bragg)
        def=deft[:monochromator]
        tau   =get(mono,:tau   ,2pi/get(a.parameters,"dm",3.355))
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
    moni=get(in22,:monitor,Struct())
    if !isa(moni,Detector)
        def=deft[:monitor]
        moni=Detector(get(moni,:width,def[:width]),get(moni,:height,def[:height]))
    end
    sample=get(in22,:sample,Struct())
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
    analyzer=get(in22,:analyzer,Struct())
    if !isa(analyzer,Bragg)
        def=deft[:analyzer]
        tau   =get(analyzer,:tau   ,2pi/get(a.parameters,"da",3.355))
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
    det=get(in22,:detector,Struct())
    if !isa(det,Detector)
        def=deft[:monitor]
        det=Detector(get(det,:width,def[:width]),get(det,:height,def[:height]))
    end
    dist=get(in22,:distances,Struct())
    if isa(dist,Struct)
        def=deft[:distances]
        s2m=get(dist,:source_monochromator,def[:source_monochromator])
        m2s=get(dist,:monochromator_sample,def[:monochromator_sample])
        s2a=get(dist,:sample_analyzer,def[:sample_analyzer])
        a2d=get(dist,:analyzer_detector,def[:analyzer_detector])
        m2m=get(dist,:monochromator_monitor,def[:monochromator_monitor])
        dist=[s2m,m2s,s2a,a2d,m2m]
    end
    isa(dist,Array) || error("config[\"/in22/distances\"] should be either a Struct or an array of the five distances not $(typeof(dist))")
    isa(dist,Vector)||(dist=vec(dist))
    hcoll=get(in22,"/collimators/hoizontal",Struct())
    if isa(hcoll,Struct)
        def=deft["/collimators/horizontal"]
        α1=get(hcoll,:source_monochromator,def[:source_monochromator])
        α2=get(hcoll,:monochromator_sample,def[:monochromator_sample])
        α3=get(hcoll,:sample_analyzer,def[:sample_analyzer])
        α4=get(hcoll,:analyzer_detector,def[:analyzer_detector])
        hcoll=[α1,α2,α3,α4]
    end
    isa(hcoll,Array)||error("config[\"/in22/collimators/horizontal\"] should be either a struct or an array of the four collimations not $(typeof(hcoll))")
    isa(hcoll,Vector)||(hcoll=vec(hcoll))
    vcoll=get(in22,"/collimators/vertical",Struct())
    if isa(vcoll,Struct)
        def=deft["/collimators/vertical"]
        β1=get(vcoll,:source_monochromator,def[:source_monochromator])
        β2=get(vcoll,:monochromator_sample,def[:monochromator_sample])
        β3=get(vcoll,:sample_analyzer,def[:sample_analyzer])
        β4=get(vcoll,:analyzer_detector,def[:analyzer_detector])
        vcoll=[β1,β2,β3,β4]
    end
    isa(vcoll,Array)||error("config[\"/in22/collimators/vertical\"] should be either a struct or an array of the four collimations not $(typeof(vcoll))")
    isa(vcoll,Vector)||(vcoll=vec(vcoll))
    # determine the magnitude of the fixed incident or final energy
    # look for a value in fixed/k, fixed/kf, or fixed/ki defaulting to value in scan file
    f=get(in22,:fixed,Struct())
    fixedk=get(f,:k,get(f,:kf,get(f,:ki,get(a.parameters,"kfix",2.55))))
    # allow for overriding fixed-k by a defined fixed-E, otherwise use fixed-k to calculate fixed-E
    fixede=get(f,:ei,get(f,:ef,get(f,:energy,ħ²2mₙ*fixedk^2)))

    # determine if ki or kf is fixed:
    fx=2!=get(a.parameters,"fx",2) # a.parameters[fx]==2 says kf-fixed; otherwise ki-fixed
    # so, if this fx is true, then the TASMAD file says ki is fixed
    fi=haskey(f,:ki)&&isa(f[:ki],Real) # ki-fixed if true <--- possible override flags
    ff=haskey(f,:kf)&&isa(f[:kf],Real) # kf-fixed if true <-/
    # only override if (true,false) or (false,true) not (true,true) or (false,false)
    is_ki_fixed= xor(fi,ff) ? fi : fx # if xor(fi,ff); fi; else; fx; end

    tas=TripleAxis(src,mono,moni,sample,analyzer,det,dist,hcoll,vcoll,fixede,is_ki_fixed)
    return tas
end

fillinstrument!(a::IN22d)=_fillTASinstrument!(a,["a1","a2","a3","a4","a5","a6"],["qh","qk","ql","en"])
