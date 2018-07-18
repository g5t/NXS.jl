function defaultconfigC5()
    src=Struct(); mono=Struct(); moni=Struct(); ana=Struct();
    det=Struct(); dist=Struct(); coll=Struct(); c5=Struct()
    src[:width]=60.32 # mm, == 2.375 in
    src[:height]=155.575 # mm, == 6.125 in
    #
    # FIXME C5 has exchangable monos with all parameters changing
    # these values are for the Heusler:
    mono[:mosaic]=2pi/1000. # 0.36° in radian
    mono[:width]=135.636
    mono[:height]=100.076
    mono[:depth]=2. # 2mm thick crystals
    mono["/curvature/horizontal"]=flatBragg
    mono["/curvature/vertical"]=verticallyCurvedBragg
    #
    moni[:width]=50.
    moni[:height]=160.
    #
    # FIXME same like for the mono
    ana[:mosaic]=2pi/1000. # 0.36° in radian
    ana[:width]=135.636
    ana[:height]=100.076
    ana[:depth]=2. # 2mm thick
    ana["/curvature/horizontal"]=horizontalyCurvedBragg
    ana["/curvature/vertical"]=flatBragg
    #
    det[:width] = 38.1 # mm, ==1.5 in
    det[:height]=127.  # mm, ==5   in
    #
    dist[:source_monochromator] =6654.8
    dist[:monochromator_sample] =1784.35
    dist[:sample_analyzer]      =1016.0
    dist[:analyzer_detector]    = 260.35
    dist[:monochromator_monitor]=1000.
    #
    coll[:source_monochromator]= -2.
    coll[:monochromator_sample]=coll[:sample_analyzer]=coll[:analyzer_detector]=600.
    #
    c5[:source]=src; c5[:monochromator]=mono; c5[:monitor]=moni;
    c5[:analyzer]=ana;  c5[:detector]=det;       c5[:distances]=dist;
    c5["/collimators/horizontal"]=coll; c5["/collimators/vertical"]=coll
    return c5
end
reset_config_C5()=config[:c5]=defaultconfigC5()
reset_config_C5()

# A general function to create a TripleAxis object representative of a typical C5 setup.
function defineC5(a::C5d)
    deft=defaultconfigC5()
    c5=config[:c5]
    #
    src=get(c5,:source,Struct())
    if !isa(src,Source)
        def=deft[:source]
        # C5 has an eliptical source
        src=ElliSource(get(src,:width,def[:width]),get(src,:height,def[:height]))
    end
    mono=get(c5,:monochromator,Struct())
    if !isa(mono,Bragg)
        def=deft[:monochromator]
        tau   =get(mono,:tau   ,2pi/get(a.parameters,"mono/d",-1))
        named =get(mono,:name  ,get(a.parameters,"mono/name",tauname(tau)))
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
    moni=get(c5,:monitor,Struct())
    if !isa(moni,Detector)
        def=deft[:monitor]
        moni=Detector(get(moni,:width,def[:width]),get(moni,:height,def[:height]))
    end
    sample=get(c5,:sample,Struct())
    if !isa(sample,Sample)
        o1=get(sample,"/orient/x",[0,0,1.])
        o2=get(sample,"/orient/y",[1,0,0.])
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
    analyzer=get(c5,:analyzer,Struct())
    if !isa(analyzer,Bragg)
        def=deft[:analyzer]
        tau   =get(analyzer,:tau   ,2pi/get(a.parameters,"ana/d",-1))
        named =get(analyzer,:name  ,get(a.parameters,"ana/name",tauname(tau)))
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
    det=get(c5,:detector,Struct())
    if !isa(det,Detector)
        def=deft[:monitor]
        det=Detector(get(det,:width,def[:width]),get(det,:height,def[:height]))
    end
    dist=get(c5,:distances,Struct())
    if isa(dist,Struct)
        def=deft[:distances]
        s2m=get(dist,:source_monochromator,def[:source_monochromator])
        m2s=get(dist,:monochromator_sample,def[:monochromator_sample])
        s2a=get(dist,:sample_analyzer,def[:sample_analyzer])
        a2d=get(dist,:analyzer_detector,def[:analyzer_detector])
        m2m=get(dist,:monochromator_monitor,def[:monochromator_monitor])
        dist=[s2m,m2s,s2a,a2d,m2m]
    end
    isa(dist,Array) || error("config[\"/c5/distances\"] should be either a Struct or an array of the five distances not $(typeof(dist))")
    isa(dist,Vector)||(dist=vec(dist))
    hcoll=get(c5,"/collimators/hoizontal",Struct())
    if isa(hcoll,Struct)
        def=deft["/collimators/horizontal"]
        α1=get(hcoll,:source_monochromator,def[:source_monochromator])
        α2=get(hcoll,:monochromator_sample,def[:monochromator_sample])
        α3=get(hcoll,:sample_analyzer,def[:sample_analyzer])
        α4=get(hcoll,:analyzer_detector,def[:analyzer_detector])
        hcoll=[α1,α2,α3,α4]
    end
    isa(hcoll,Array)||error("config[\"/c5/collimators/horizontal\"] should be either a struct or an array of the four collimations not $(typeof(hcoll))")
    isa(hcoll,Vector)||(hcoll=vec(hcoll))
    vcoll=get(c5,"/collimators/vertical",Struct())
    if isa(vcoll,Struct)
        def=deft["/collimators/vertical"]
        β1=get(vcoll,:source_monochromator,def[:source_monochromator])
        β2=get(vcoll,:monochromator_sample,def[:monochromator_sample])
        β3=get(vcoll,:sample_analyzer,def[:sample_analyzer])
        β4=get(vcoll,:analyzer_detector,def[:analyzer_detector])
        vcoll=[β1,β2,β3,β4]
    end
    isa(vcoll,Array)||error("config[\"/c5/collimators/vertical\"] should be either a struct or an array of the four collimations not $(typeof(vcoll))")
    isa(vcoll,Vector)||(vcoll=vec(vcoll))
    # determine the magnitude of the fixed incident or final energy
    # look for a value in fixed/k, fixed/kf, or fixed/ki defaulting to value in scan file
    f=get(c5,:fixed,Struct())
    fixedk=get(f,:k,get(f,:kf,get(f,:ki,get(a.parameters,"kfix",-1))))
    # allow for overriding fixed-k by a defined fixed-E, otherwise use fixed-k to calculate fixed-E
    fixede=fixedk>0?ħ²2mₙ*fixedk^2:get(f,:ei,get(f,:ef,get(f,:energy,get(a.parameters,"eprim",3.52)*1e12*value(_planck))))

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

# There are multiple problems with using Measured values here. Drop them for the time being.
function fillinstrument!{N}(a::C5d{N})
    oneC5=a.definst(a)
    # C5 (hopefully) always records angles (2tm,psi,phi,2ta) and Q,E-related variables (zeta,eta,nu)
    # we can determine (Q,E) from (zeta,eta,nu) using the orienting vectors and (a,r) in the parameters dict
    # but this seems dumb, since we should have all angles
    a2=getVal(a,"2tm")
    a3=getVal(a,"psi")
    a4=getVal(a,"phi")
    a6=getVal(a,"2ta")
    a1=a2/2
    a5=a6/2
    # each of a1-a6 is an N-1-dimensional array
    As=cat(N,a1,a2,a3,a4,a5,a6) # concatinate them in the N+1 dimension
    As=permutedims(As,circshift(1:N,1)) # setangles *requires* that a1-a6 are in the first dimension
    c5=setangles(oneC5,As/180*pi) # TripleAxis angles must be in radian
    a.instrument=c5
end
