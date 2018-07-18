function defaultconfigRITA2()
    src=Struct(); mono=Struct(); moni=Struct(); ana=Struct();
    det=Struct(); dist=Struct(); coll=Struct(); rita2=Struct()
    src[:width]=30.
    src[:height]=120.
    #
    mono[:mosaic]=pi/270. # 40' in radian
    mono[:width]=3*50. # 3×50mm crystals wide
    mono[:height]=5*25.+4*1. # 5×25mm crystals + 4×1mm gaps tall
    mono[:depth]=2. # 2mm thick crystals
    mono["/curvature/horizontal"]=flatBragg
    mono["/curvature/vertical"]=verticallyCurvedBragg
    #
    moni[:width]=50.
    moni[:height]=160.
    #
    ana[:mosaic]=pi/270. # 40' in radian
    ana[:width]=7*25.+6*1. # 7×25mm crystals plus 6×1mm gaps wide
    ana[:height]=3*50. # 3×50mm crystals tall
    ana[:depth]=2. # 2mm thick
    ana["/curvature/horizontal"]=horizontalyCurvedBragg
    ana["/curvature/vertical"]=flatBragg
    #
    det[:width] = 30. #mm, most common width, though variable with often-used values 10, 20, 30, 50 mm
    det[:height]=160. #mm, height of standard detector mask
    #
    dist[:source_monochromator] =1750.
    dist[:monochromator_sample] =1700.
    dist[:sample_analyzer]      =1400.
    dist[:analyzer_detector]    =1000.
    dist[:monochromator_monitor]=1000.
    #
    coll[:source_monochromator]=coll[:monochromator_sample]=
         coll[:sample_analyzer]=coll[:analyzer_detector]=800.
    #
    rita2[:source]=src; rita2[:monochromator]=mono; rita2[:monitor]=moni;
    rita2[:analyzer]=ana;  rita2[:detector]=det;       rita2[:distances]=dist;
    rita2["/collimators/horizontal"]=coll; rita2["/collimators/vertical"]=coll
    return rita2
end
reset_config_RITA2()=config[:rita2]=defaultconfigRITA2()
reset_config_RITA2()

# FIXME very little of this is right
function defaultconfigCAMEA()
    camea=Struct(); ana=Struct(); det=Struct()
    #
    ana[:mosaic]=pi/180. # 60' in radian. value from CAMEA paper
    ana[:width]=20.
    ana[:height]=60.
    ana[:depth]=2. # 2mm thick
    ana["/curvature/horizontal"]=flatBragg
    ana["/curvature/vertical"]=horizontalyCurvedBragg # radius is like horizontal or like vertical for horizontal scattering?
    #
    det[:width] =25.
    det[:height]=50.
    #
    camea[:analyzer]=ana; camea[:detector]=det;
    return camea
end
reset_config_CAMEA()=config[:camea]=defaultconfigCAMEA
reset_config_CAMEA()


# A general function to create a TripleAxis object representative of a typical RITA2 setup.
function defineRITA2common(a::Union{RITA2d,CAMEAd})
    deft=defaultconfigRITA2()
    rita2=config[:rita2]
    #
    src=get(rita2,:source,Struct())
    if !isa(src,Source)
        def=deft[:source]
        src=Source(get(src,:width,def[:width]),get(src,:height,def[:height]))
    end
    mono=get(rita2,:monochromator,Struct())
    if !isa(mono,Bragg)
        def=deft[:monochromator]
        tau   =get(mono,:tau   ,2pi/get(a.parameters,"dm",3.3541627))
        named =get(mono,:name  ,tauname(tau))
        sense =get(mono,:sense ,get(a.parameters,"sm",+1))
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
    moni=get(rita2,:monitor,Struct())
    if !isa(moni,Detector)
        def=deft[:monitor]
        moni=Detector(get(moni,:width,def[:width]),get(moni,:height,def[:height]))
    end
    sample=get(rita2,:sample,Struct())
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

    dist=get(rita2,:distances,Struct())
    if isa(dist,Struct)
        def=deft[:distances]
        s2m=get(dist,:source_monochromator,def[:source_monochromator])
        m2s=get(dist,:monochromator_sample,def[:monochromator_sample])
        s2a=get(dist,:sample_analyzer,def[:sample_analyzer])
        a2d=get(dist,:analyzer_detector,def[:analyzer_detector])
        m2m=get(dist,:monochromator_monitor,def[:monochromator_monitor])
        dist=[s2m,m2s,s2a,a2d,m2m]
    end
    isa(dist,Array) || error("config[\"/rita2/distances\"] should be either a Struct or an array of the five distances not $(typeof(dist))")
    isa(dist,Vector)||(dist=vec(dist))
    hcoll=get(rita2,"/collimators/hoizontal",Struct())
    if isa(hcoll,Struct)
        def=deft["/collimators/horizontal"]
        α1=get(hcoll,:source_monochromator,def[:source_monochromator])
        α2=get(hcoll,:monochromator_sample,def[:monochromator_sample])
        α3=get(hcoll,:sample_analyzer,def[:sample_analyzer])
        α4=get(hcoll,:analyzer_detector,def[:analyzer_detector])
        hcoll=[α1,α2,α3,α4]
    end
    isa(hcoll,Array)||error("config[\"/rita2/collimators/horizontal\"] should be either a struct or an array of the four collimations not $(typeof(hcoll))")
    isa(hcoll,Vector)||(hcoll=vec(hcoll))
    vcoll=get(rita2,"/collimators/vertical",Struct())
    if isa(vcoll,Struct)
        def=deft["/collimators/vertical"]
        β1=get(vcoll,:source_monochromator,def[:source_monochromator])
        β2=get(vcoll,:monochromator_sample,def[:monochromator_sample])
        β3=get(vcoll,:sample_analyzer,def[:sample_analyzer])
        β4=get(vcoll,:analyzer_detector,def[:analyzer_detector])
        vcoll=[β1,β2,β3,β4]
    end
    isa(vcoll,Array)||error("config[\"/rita2/collimators/vertical\"] should be either a struct or an array of the four collimations not $(typeof(vcoll))")
    isa(vcoll,Vector)||(vcoll=vec(vcoll))
    # determine the magnitude of the fixed incident or final energy
    # look for a value in fixed/k, fixed/kf, or fixed/ki defaulting to value in scan file
    f=get(rita2,:fixed,Struct())
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

    (src,mono,moni,sample,dist,hcoll,vcoll,fixede,is_ki_fixed)
end

function defineRITA2(a::RITA2d)
    deft=defaultconfigRITA2()
    rita2=config[:rita2]
    (src,mono,moni,sample,dist,hcoll,vcoll,fixede,is_ki_fixed)=defineRITA2common(a)
    analyzer=get(rita2,:analyzer,Struct())
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
    det=get(rita2,:detector,Struct())
    if !isa(det,Detector)
        def=deft[:monitor]
        det=Detector(get(det,:width,def[:width]),get(det,:height,def[:height]))
    end
    tas=TripleAxis(src,mono,moni,sample,analyzer,det,dist,hcoll,vcoll,fixede,is_ki_fixed)
    return tas
end

function defineRITA2(a::CAMEAd)
    deft=defaultconfigCAMEA()
    camea=config[:camea]
    (src,mono,moni,sample,dist,hcoll,vcoll,fixede,is_ki_fixed)=defineRITA2common(a)
    dist[4]=50. # FIXME look this up again
    analyzer=get(camea,:analyzer,Struct())
    if !isa(analyzer,Bragg)
        def=deft[:analyzer]
        tau   =get(analyzer,:tau   ,2pi/get(a.parameters,"da",3.3541627))
        named =get(analyzer,:name  ,tauname(tau))
        sense =get(analyzer,:sense ,get(a.parameters,"sa",+1)) # FIXME ## modify resolution calculations to include out-of-plane scattering? sa=+1 is up?
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
    det=get(rita2,:detector,Struct())
    if !isa(det,Detector)
        def=deft[:monitor]
        det=Detector(get(det,:width,def[:width]),get(det,:height,def[:height]))
    end
    tas=TripleAxis(src,mono,moni,sample,analyzer,det,dist,hcoll,vcoll,fixede,is_ki_fixed)
    return tas
end

fillinstrument!(a::RITA2d)=fillTASinstrument!(a) # simple for "normal" RITA-2

function fillinstrument!(a::CAMEAd)
    oneInst=a.definst(a)

    # the number of a4 bins is probably 104 (8 segments, 13 tubes per segment)
    # but if they're arrayed radially the outer tubes won't be illuminated
    # since segmentΔa4 = [3.82,4.14,4.42,4.73,4.95,5.19,5.41,5.60] #°,
    # from the CAMEA paper. So what a4 value gets assigned to each tube?
    #
    # the number of energy bins could be as small as 8 (one bin per analyzer)
    # or as high as 452 (904 mm active length tubes, 2 mm bins in hardware)
    # XXX the actual details of *how* the tube is split into M-bins could be
    #     extremely important, but lets not worry about that now.
    nopt,noa4,noen=size(a.detector,1,2,3) # should be N,13*8,8≤M≤452 for CAMEA

    minangles = ["a3","a4","a6"] # with only these angles we *can* define a CAMEA triple axis position
    angles=["a1","a2","a3","a4","a5","a6"] # angles necessary to define a triple axis position
    qenames=["qh","qk","ql","en"]

    if mean(isCol.(a,minangles)) >= mean(isCol.(a,qenames))
        headAs=getMultiple(a.variables,angles,0.) # the values in the header
        As=zeros(nopt,noa4,noen,length(angles))
        for i=[1:4;6] # treat a5 separately
            As[:,:,:,i]= isCol(a,angles[i]) ? getVal(a,angles[i]): headAs[i]
        end
        As[:,:,:,5]= isCol(a,angles[5]) ? getVal(a, angles[5]) : As[:,:,:,6]/2;
        inst=setangles(oneInst,permutedims(As/180*pi,(4,1,2,3)))
    else
        headqes=getMultiple(a.variables,qenames,0.)
        QE=zeros(nopt,noa4,noen,length(qenames))
        for i=1:length(qenames)
            QE[:,:,:,i]=isCol(a,qenames[i])?getVal(a,qenames[i]):headqes[i]
        end
        QE4v=LatticeVector(permutedims(QE,(4,1,2,3)),oneInst.sample.crystal.rlat)
        inst=gotoQE(oneInst,QE4v)
    end
    # now we can modify all TripleAxis objects to have the correct CAMEA
    # sample-analyzer distances, analyzer-detector distances,
    # horizontal widths, curvatures(?)
    sals=[930,994,1056,1120,1183,1247,1312,1379.] # mm, from the CAMEA paper
    # The reported analyzer-detector distances are for the center analyzer only.
    # adls=[707,702, 700, 701, 703, 709, 717, 727.] # mm, from the CAMEA paper
    # Instead we should calculate the distance from a6 and the difference in
    # in height between the analyzers and detectors (I think 700 mm).
    analyzer_detector_height = 700.

    # we can create more-accurate analyzer objects for each of the 8 types
    # of analyzers.
    named="pg(002)"
    tau = gettau(named)
    #
    analyzer_efs=[3.21,3.38,3.58,3.80,4.05,4.33,4.64,5.01] # from the CAMEA paper
    analyzer_kfs=sqrt.(analyzer_efs/ħ²2mₙ)
    #analyzer_rot=asin.(tau./2analyzer_kfs) # rotation from horizontal plane
    #
    sense = +1 # assuming we eventually figure out the modifications to Popovici's method, and that +1 in this case is vertical scattering
    mosaic = pi/180*[1,1] # 60' in radian; from the CAMEA paper
    widths = [72.0,81.8,91.6,102.4,112.2,119.1,128.1,139.2] # mm , from the CAMEA paper
    heights = 2*tand(2)*sals # CAMEA paper says ±2°; this calculation isn't right since the analyzers are rotated
    #height ./ = sin.(analyzer_rot) # correct for the rotation angle
    heights .*= 2analyzer_kfs/tau # correct for the rotation angle; avoid sin.(asin.(...)) construction
    depth = 1. # mm, CAMEA paper
    cfh= flatBragg;
    cfv = rowlandCurvedBragg; #
    analyzers=[ Bragg(tau,named,sense,mosaic,[w,h,depth],[cfh,cfv]) for (w,h) in zip(widths,heights) ]

    mod(noen,8)==0 || status(:warn,"$noen energy bins is not divisible by 8!?")
    noen_per_analyzer = noen >> 3 # hopefully mod(noen,8) == 0
    @assert 0<noen_per_analyzer<57 "How did you manage to split $noen energy bins between CAMEA's 8 analyzers?"
    analyzer_index = cld.(1:noen, noen_per_analyzer) # if mod(noen,8)>0 some analyzer_index values will be >8, which will be a problem.
    # worry about that later ¯\_(ツ)_/¯
    for n=1:noen, i in inst[:,:,n]
        i.ana = analyzers[analyzer_index[n]]
        i.arms[3] = sals[analyzer_index[n]]
        i.arms[4] = analyzer_detector_height/sin( getangle(i,6) ) # sin() since TripleAxis angles are stored in radian
    end

    a.instrument=inst
end

# # the a4 values input so far are the segment centers?
# segmentΔa4 = [3.82,4.14,4.42,4.73,4.95,5.19,5.41,5.60] #°, from the CAMEA paper
