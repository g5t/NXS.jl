
# function absorption_correction(t::TripleAxis)
#     isa(t.sample,Lattices.GeometrySample) || (warn("A GeometrySample is needed to calculate an absorption correction"), return 1)
#     Q=getQ(t)
#     Qxa=atan2(dot(Q,gety(t)),dot(Q,getx(t))) # the angle between Q and x
#     Qka=-getphi(t) # the angle between Q and ki
#     rot=basisrotation3xy(Qxa+Qka) # we need to change description axes of the sample shape to have x along ki
#     ingeom = t.sample.geometry
#     t.sample.geometry = Lattices.transform(rot,ingeom)
#     invA = Lattices.absorption_correction(t.a[4],t.sample) # TODO a.a is in radians?
#     a.sample.geometry = ingeom # rotate the geometry back
#     return invA
# end

# function absorption_correction_angles(t::TripleAxis)
#     Q=getQ(t)
#     Qxa=atan2(dot(Q,gety(t)),dot(Q,getx(t))) # the angle between Q and x
#     Qka=-getphi(t) # the angle between Q and ki
#     (Qxa,Qka,t.a[4])
# end
# absorption_correction_angles(t::Neutrond)=absorption_correction_angles.(t.instrument)

function absorption_correction(t::TripleAxis;rQx=1,rQk=-1,kwds...)
    isa(t.sample,Lattices.GeometrySample) || (warn("A GeometrySample is needed to calculate an absorption correction"), return 1)
    Q=getQ(t)
    Qxa=rQx*atan2(dot(Q,gety(t)),dot(Q,getx(t))) # the angle between Q and x
    Qka=rQk*getphi(t) # the angle between Q and ki
    rot=basisrotation3xy(Qxa+Qka) # we need to change description axes of the sample shape to have x along ki
    Lattices.absorption_correction(t.a[4],t.sample.μ,Lattices.transform(rot,t.sample.geometry);kwds...) # TODO a.a is in radians?
end

function correct_absorption!(a::Neutrond;kwds...)
    invA = absorption_correction.(a.instrument;kwds...)
    setDat(a,a.y,invA.*a[a.y])
    # hasbeenfit(a) && (a.fitres= invA*a.fitres) # this can't work since multiplying the fit result changes only a global scale
    a.fits = typeof(a.fits)(0) # remove any fit information if present
    return a
end
correct_absorption(a::Neutrond;kwds...)=correct_absorption!(deepcopy(a);kwds...)

function perform_absorption_correction(a::Scatterd;kwds...)
    error("Calculating the absorption correction is only supported thus far for \n
           neutron scattering, but could be implemented for other scattering processes.")
end
