function defaultconfigFLEXX()
    src=Struct(); mono=Struct(); moni=Struct(); ana=Struct();
    det=Struct(); dist=Struct(); coll=Struct(); flexx=Struct()
    src[:width]=60.
    src[:height]=125.
    #
    mono[:mosaic]=pi/270. # 40' in radian
    mono[:width]=300.
    mono[:height]=140.
    mono[:depth]=2. # 2mm thick crystals
    mono["/curvature/horizontal"]=horizontalyCurvedBragg
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
    dist[:sample_analyzer]      =1400.
    dist[:analyzer_detector]    = 460. # and variable
    dist[:monochromator_monitor]=1000.
    #
    coll[:source_monochromator]= -2.
    coll[:monochromator_sample]=coll[:sample_analyzer]=coll[:analyzer_detector]=80.
    #
    flexx[:source]=src; flexx[:monochromator]=mono; flexx[:monitor]=moni;
    flexx[:analyzer]=ana;  flexx[:detector]=det;       flexx[:distances]=dist;
    flexx["/collimators/horizontal"]=coll; flexx["/collimators/vertical"]=coll
    return flexx
end
reset_config_FLEXX()=config[:flexx]=defaultconfigFLEXX()
reset_config_FLEXX()

function defaultconfigmultiFLEXX()
    src=Struct(); mono=Struct(); moni=Struct(); ana=Struct();
    det=Struct(); dist=Struct(); coll=Struct(); flexx=Struct()
    src[:width]=60.
    src[:height]=125.
    #
    mono[:mosaic]=pi/270. # 40' in radian
    mono[:width]=300.
    mono[:height]=140.
    mono[:depth]=2. # 2mm thick crystals
    mono["/curvature/horizontal"]=horizontalyCurvedBragg
    mono["/curvature/vertical"]=verticallyCurvedBragg
    #
    moni[:width]=70.
    moni[:height]=160.
    #
    ana[:mosaic]=pi/450. # 0.4° in radian
    ana[:width]=20. #
    ana[:height]=60.
    ana[:depth]=2. # 2mm thick
    ana["/curvature/horizontal"]=horizontalyCurvedBragg
    ana["/curvature/vertical"]=flatBragg
    #
    det[:width] = 25.
    det[:height]= 50.
    #
    dist[:source_monochromator] =1750.
    dist[:monochromator_sample] =1700.
    dist[:sample_analyzer]      = NaN  # dist[3] = [1050 mm, 1220 mm, 1387 mm, 1552 mm and 1732 mm] depending on the energy channel
    dist[:analyzer_detector]    = 400. # and fixed
    dist[:monochromator_monitor]=1000.
    #
    coll[:source_monochromator]= -2.
    coll[:monochromator_sample]=coll[:sample_analyzer]=coll[:analyzer_detector]=80.
    #
    flexx[:source]=src; flexx[:monochromator]=mono; flexx[:monitor]=moni;
    flexx[:analyzer]=ana;  flexx[:detector]=det;       flexx[:distances]=dist;
    flexx["/collimators/horizontal"]=coll; flexx["/collimators/vertical"]=coll
    return flexx
end
reset_config_multiFLEXX()=config[:multiflexx]=defaultconfigmultiFLEXX()
reset_config_multiFLEXX()

# A general function to create a TripleAxis object representative of a typical FLEXX setup.

defineFLEXX(a::FLEXXd)=     defineFLEXXcommon(a,config[:flexx],     defaultconfigFLEXX()     )
defineFLEXX(a::multiFLEXXd)=defineFLEXXcommon(a,config[:multiflexx],defaultconfigmultiFLEXX())
function defineFLEXXcommon(a,srct::Struct,deft::Struct)
    #
    src=get(srct,:source,Struct())
    if !isa(src,Source)
        def=deft[:source]
        src=Source(get(src,:width,def[:width]),get(src,:height,def[:height]))
    end
    mono=get(srct,:monochromator,Struct())
    if !isa(mono,Bragg)
        def=deft[:monochromator]
        tau   =get(mono,:tau   ,2pi/get(a.parameters,"dm",3.3541627))
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
    moni=get(srct,:monitor,Struct())
    if !isa(moni,Detector)
        def=deft[:monitor]
        moni=Detector(get(moni,:width,def[:width]),get(moni,:height,def[:height]))
    end
    sample=get(srct,:sample,Struct())
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
    analyzer=get(srct,:analyzer,Struct())
    if !isa(analyzer,Bragg)
        def=deft[:analyzer]
        tau   =get(analyzer,:tau   ,2pi/get(a.parameters,"da",3.3541627))
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
    det=get(srct,:detector,Struct())
    if !isa(det,Detector)
        def=deft[:monitor]
        det=Detector(get(det,:width,def[:width]),get(det,:height,def[:height]))
    end
    dist=get(srct,:distances,Struct())
    if isa(dist,Struct)
        def=deft[:distances]
        s2m=get(dist,:source_monochromator,def[:source_monochromator])
        m2s=get(dist,:monochromator_sample,def[:monochromator_sample])
        s2a=get(dist,:sample_analyzer,def[:sample_analyzer])
        a2d=get(a.variables,"td",def[:analyzer_detector]) # this is variable and recorded
        m2m=get(dist,:monochromator_monitor,def[:monochromator_monitor])
        dist=[s2m,m2s,s2a,a2d,m2m]
    end
    isa(dist,Array) || error("distances should be either a Struct or an array of the five distances not $(typeof(dist))")
    isa(dist,Vector)||(dist=vec(dist))
    hcoll=get(srct,"/collimators/hoizontal",Struct())
    if isa(hcoll,Struct)
        def=deft["/collimators/horizontal"]
        α1=get(hcoll,:source_monochromator,def[:source_monochromator])
        α2=get(hcoll,:monochromator_sample,def[:monochromator_sample])
        α3=get(hcoll,:sample_analyzer,def[:sample_analyzer])
        α4=get(hcoll,:analyzer_detector,def[:analyzer_detector])
        hcoll=[α1,α2,α3,α4]
    end
    isa(hcoll,Array)||error("collimators/horizontal should be either a struct or an array of the four collimations not $(typeof(hcoll))")
    isa(hcoll,Vector)||(hcoll=vec(hcoll))
    vcoll=get(srct,"/collimators/vertical",Struct())
    if isa(vcoll,Struct)
        def=deft["/collimators/vertical"]
        β1=get(vcoll,:source_monochromator,def[:source_monochromator])
        β2=get(vcoll,:monochromator_sample,def[:monochromator_sample])
        β3=get(vcoll,:sample_analyzer,def[:sample_analyzer])
        β4=get(vcoll,:analyzer_detector,def[:analyzer_detector])
        vcoll=[β1,β2,β3,β4]
    end
    isa(vcoll,Array)||error("collimators/vertical should be either a struct or an array of the four collimations not $(typeof(vcoll))")
    isa(vcoll,Vector)||(vcoll=vec(vcoll))
    # determine the magnitude of the fixed incident or final energy
    # look for a value in fixed/k, fixed/kf, or fixed/ki defaulting to value in scan file
    f=get(srct,:fixed,Struct())
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

fillinstrument!(a::FLEXXd)=_fillTASinstrument!(a,["a1","a2","a3","a4","a5","a6"],["qh","qk","ql","en"])
function fillinstrument!(a::multiFLEXXd)
    oneFLEXX=a.definst(a)

    # try to grab default values for a1-a6; qh,qk,ql,en; dqh,dqk,dql,den
    nopt,noa4,noen=size(a.detector,1,2,3) # should be N,31,5 for multiFLEXX

    angles=["a1","a2","a3","a4","a5","a6"] # angles necessary to define a triple axis position
    qenames=["qh","qk","ql","en"]

    if mean(isCol.(a,angles)) >= mean(isCol.(a,qenames))
        headAs=getMultiple(a.variables,angles,0.) # the values in the header
        As=zeros(nopt,noa4,noen,length(angles))
        for i=1:length(angles)
            As[:,:,:,i]= isCol(a,angles[i])?getVal(a,angles[i]):headAs[i]
        end
        inst=setangles(oneFLEXX,permutedims(As/180*pi,(4,1,2,3)))
    else
        headqes=getMultiple(a.variables,qenames,0.)
        QE=zeros(nopt,noa4,noen,length(qenames))
        for i=1:length(qenames)
            QE[:,:,:,i]=isCol(a,qenames[i])?getVal(a,qenames[i]):headqes[i]
        end
        QE4v=LatticeVector(permutedims(QE,(4,1,2,3)),oneFLEXX.sample.crystal.rlat)
        inst=gotoQE(oneFLEXX,QE4v)
    end
    # now we can modify all TripleAxis objects to have the correct multiFLEXX
    # sample-analyzer distances (and analyzer radius of curvatures?)
    sals=[1050, 1220, 1387, 1552, 1732.] # mm, from the multiFLEXX paper
    for n=1:noen, i in inst[:,:,n]; i.arms[3]=sals[n]; end

    a.instrument=inst
end
