["""
A conic section in four dimensions is described by a 5×5 matrix `W`,

    \ \ ⎛xx xy xz xE x0⎞\n\t\ \ ⎜yx yy yz yE y0⎟\n\tW=⎜zx zy zz zE z0⎟\n\t\ \ ⎜Ex Ey Ez EE E0⎟\n\t\ \ ⎝x0 y0 z0 E0 -R⎠

For the conic section described by `W` to be non-degenerage, the determinant of `W` must be finite.
The elements `x0`, `y0`, `z0`, and `E0` of `W` determine the origin of the conic section,
while `R` defines an overall scaling parameter.
If a transformation is performed such that the conic section is centered at the origin and the
scaling parameter is `1`, then `W` becomes

    \ \ ⎛xx xy xz xE  0⎞\n\t\ \ ⎜yx yy yz yE  0⎟\n\tW=⎜zx zy zz zE  0⎟\n\t\ \ ⎜Ex Ey Ez EE  0⎟\n\t\ \ ⎝ 0  0  0  0 -1⎠

And the conic section is described by the equation `QE1'*W*QE1=0`
where `QE1` is the momentum-energy-unity 5-vector -- `QE1=[Qx Qy Qz E 1]`.

In the restricted case of a conic section centered at the origin with unity scaling, we can define
a useful submatrix `M` as

    \ \ ⎛xx xy xz xE⎞\n\tM=⎜yx yy yz yE⎟\n\t\ \ ⎜zx zy zz zE⎟\n\t\ \ ⎝Ex Ey Ez EE⎠

and note that `QE'*M*QE=1` where `QE` is the momentum-energy 4-vector -- `QE=(Qx,Qy,Qz,E)` --
is sufficient to describe the conic section.
The type of conic section described by `M` is determined by its determinate with

| `det(M)` | conic section |
|:--------:|:-------------:|
| - | hyperboloid |
| 0 | paraboloid |
| + | ellipsoid |

The eigenvalues of `M` are the inverse squared-radii of the ellipsoid along principal axes given
by the eigenvectors of `M`.
"""]

# A TripleAxis object contains all relevant information about:
#   the monochromator (scattering sense, tau, d spacing, mosaic widths, reflectivity, shape information, and radii of curvature)
#   the sample (scattering sense, crystal lattice, orientation relative to the sample table, mosaic widths, shape information,
#   the analyzer (scattering sense, tau, d spacing, mosaic widths, reflectivity, shape information, and radii of curvature)
#   the complete instrument (all angles, horizontal and vertical collimations, fixed ki or kf, fixed energy value, and distances between elements)
# some quantities might need to be calculated; functions already exist for ki,kf,Q,E, angle(ki,Q)=calcphi, angle(kf,Q)=pi-a4-calcphi

"""
`ResFloat` is an alias for floating point machine precision.
It can be swapped to `Float32` to save memory on a 64-bit system,
or to `Float64` to increase precision on a 32-bit system.
"""
const ResFloat = typeof(1.)
"""`ResMatrix` is a type alias for `Array{ResFloat,2}`"""
const ResMatrix = Array{ResFloat,2}
"""
The type `ResRM` is a composite type to describe the resolution of a triple-axis spectrometer
at a single point.

| Fieldname | Description |
|:---------:|:------------|
| R | The resolution prefactor |
| M | The resolution matrix |
| eigvals | eigenvalues of M |
| eigvecs | eigenvectors of M |

The `eigvals` and `eigvecs` should not be computed prior to constructing a `ResRM` object as
an internal constructor will use the fact that the resolution matrix is (likely) Hermitian to
calculate the eigenvalues and eigenvectors.
"""
type ResRM
    R::ResFloat
    M::ResMatrix
    eigvals::Array{ResFloat,1}
    eigvecs::Array{ResFloat,2}
    function ResRM(a::ResFloat,b::ResMatrix)
        @assert all([size(b)...].==[4,4]) "A resolution matrix must be 4x4"
        @assert isapprox(b,b') "The resolution matrix isn't symmetric (that seems unlikely)"
        M=Hermitian((b+b')/2) # average b and it's transpose because Hermitian only takes the upper-diagonal elements when creating a Hermitian Matrix
        eg=eigfact(M)
        new(a,b,eg.values,eg.vectors)
    end
end
#ResRM(a,b,discarded...)=ResRM(a,b) # Throw away extra arguments. This is probably a bad idea in general, but maybe OK here (used by *R0M() and *RmM()).
ResRM{T<:Number}(a::ResFloat,b::ResMatrix,c::Array{T},d::Array{T})=ResRM(a,b)
ResRM{T<:Number}(a::ResFloat,b::ResMatrix,c::Array{T},d::Array{T},f::Array{T})=ResRM(a,b)

Base.eigfact(a::ResRM)=(a.eigvals,a.eigvecs)
Base.eigvals(a::ResRM)=a.eigvals
Base.eigvecs(a::ResRM)=a.eigvecs

function which_isE(M::ResRM)
    isE=broadcast(isapprox,abs.(M.eigvecs[3:3,:]),1)
    sum(isE)!=1 && warn("Eigen vector along energy can not be determined for $(M.eigvecs)")
    return isE
end
function hkleigvecs(M::ResRM)
    reshape(M.eigvecs[ ([0,0,1,0].+which_isE(M)).==0 ],(3,3))
end
function hkleigvals(M::ResRM)
    M.eigvals[.!vec(which_isE(M))]
end
#hkleigvecs(M::ResRM)=reshape(M.eigvecs[BitArray([1 1 0 1; 1 1 0 1; 0 0 0 0; 1 1 0 1])],(3,3))
#hkleigvals(M::ResRM)=M.eigvals[[1,2,4]]


function resolutionBraggWidths(RM::ResRM)
    FWHMxyEz=_σ2FWHM(1./sqrt(RM.M[(1:4)+4*(0:3)])) # M is in inverse squared Gaussian widths
    return FWHMxyEz # FWHM of a Bragg peak in [Qx,Qy,E,Qz] -- sample orientation depenent
end
function resolutionPhononWidth(RM::ResRM,n::Vector) # n should be the 4-vector normal to a plane "phonon"
    @assert 4==length(n)
    t=diagm(ones(ResFloat,4))
    t[3,:]=n # the rows/columns of RM.M are [Qx,Qy,En,Qz]
    s=inv(t)
    Mp=s'*RM.M*s
    # reduce to 1x1 (  reduce to 2x2 (  reduce to 3x3 (    )  )  ) pull single value out of Array
    # integrate over the three Qi directions, their order is not important but it *is* important
    # to not integrate out the energy axis -- 4,2,1 <-> Qz,Qy,Qx
    _σ2FWHM(1/_integrateProject(_integrateProject(_integrateProject(Mp,4),2),1)[1])
end
function resolutionVanadiumWidth(RM::ResRM) # equivalent to a "phonon" sheet perpendicular to energy
    # integrate over the three Qi directions, their order is not important but it *is* important
    # to not integrate out the energy axis -- 4,2,1 <-> Qz,Qy,Qx
    _σ2FWHM(1/sqrt(_integrateProject(_integrateProject(_integrateProject(RM.M,4),2),1)[1]))
end
for x in (:resolutionBraggWidths,:resolutionVanadiumWidth)
    @eval $x{T<:ResRM}(v::Array{T})=map($x,v)
end
resolutionPhononWidth{T<:ResRM}(v::Array{T},n::Vector)=map(x->resolutionPhononWidth(x,n),v)
function resolutionPhononWidth{T<:ResRM,R<:Vector}(v::Array{T},nv::Array{R})
    #@assert compatible(v,nv)
    broadcast(resolutionPhononWidth,v,nv)
end
export resolutionBraggWidths,resolutionPhononWidth,resolutionVanadiumWidth

function Base.show(io::IO,a::ResRM)
    println(io,"R : $(a.R)")
    Base.show(io,a.M)
#    println(io,"\neigen values:")
#    Base.show(io,a.eigvals)
#    println(io,"\neigen vectors:")
#    Base.show(io,a.eigvecs)
end

_getR(a::ResRM)=a.R;
_getM(a::ResRM)=a.M;
#_getM{T<:ResRM}(a::Array{T,1})=(v=zeros(ResFloat,4,4,length(a)); for i=1:length(a); v[:,:,i]=a[i].M; end; v)
#_getM{T<:ResRM}(a::Array{T,2})=((x,y)=size(a);v=zeros(ResFloat,4,4,x,y); for j=1:y; for i=1:x; v[:,:,i,j]=a[i,j].M; end; end; v)
_getMij(a::ResRM,i::Integer,j::Integer)=(@assert (i<5)&(i>0)&(j<5)&(j>0); a.M[i,j])
_getMij{T<:ResRM}(a::Array{T,1},i::Integer,j::Integer)=(@assert (i<5)&(i>0)&(j<5)&(j>0); la=length(a); v=zeros(ResFloat,la); for k=1:la; v[k]=a[k].M[i,j]; end; v)
_getMij{T<:ResRM,N}(a::Array{T,N},i::Integer,j::Integer)=(@assert (i<5)&(i>0)&(j<5)&(j>0); sa=size(a); v=zeros(ResFloat,sa...); for k=1:prod(sa...); v[k]=a[k].M[i,j]; end; v)

#Base.issymmetric(a::ResRM)=Base.issymmetric(a.M)
Base.issymmetric(a::ResRM)=(sum(a.M-transpose(a.M))<1e-10) # use a fuzzy comparison in case of rounding errors
Base.issymmetric{T<:ResRM}(a::Array{T,1})=map(x->Base.issym(x),a)
Base.issymmetric{T<:ResRM,N}(a::Array{T,N})=map(x->Base.issym(x),a)
Base.det(a::ResRM)=Base.det(a.M)
Base.det{T<:ResRM}(a::Array{T,1})=map(x->Base.det(x),a)
Base.det{T<:ResRM,N}(a::Array{T,N})=map(x->Base.det(x),a)
Base.trace(a::ResRM)=Base.trace(a.M)
Base.trace{T<:ResRM}(a::Array{T,1})=map(x->Base.trace(x),a)
Base.trace{T<:ResRM}(a::Array{T,2})=map(x->Base.trace(x),a) # necessary because trace is typically only defined for Array{Number,2}
Base.trace{T<:ResRM,N}(a::Array{T,N})=map(x->Base.trace(x),a)

["""
What follows is a translitteration of ResMat.m by A. Zheludev initially as part of ResLib.
The translated version here includes bug-fixes released as part of ResLibCal by E. Farhi

ResLib used second-moments of an object's size -- thus allowing for,
e.g., a slab shaped sample rotated about an arbitrary axis; ResLibCal on
the other hand has done away with this convention, instead using FWHM everywhere
and therefore removing an extra factor
`0.424660900144 = FWHM2RMS`
from the `CONVERT1` factor used by ResLib, `CONVERT1=0.4246609*pi/60/180;`.
Instead, ResLibCal uses `CONVERT1=pi/60/180;` which converts arc-minutes to radians
and `CONVERT2=2.072;` which is a less precise representation of my
`const ħ²2mₙ=2.072142 # ħ²/(2mₙ) = 2.072142 meV Å²`

`coopernathansR0M`/`coopernathansRmM`/`popoviciR0M`/`popoviciRmM` all calculate and return
`R0`, the resolution prefactor (converted to incident flux or monitor counts)
Plus the resolution matrix given by

    \ \ ⎛xx xy xz xE⎞\n\tM=⎜yx yy yz yE⎟\n\t\ \ ⎜zx zy zz zE⎟\n\t\ \ ⎝Ex Ey Ez EE⎠
"""]

"""
A single function to calculate `R₀` and `M` independent of which resolution approximation is made.
"""
function _calcR0M(t::TripleAxis,R0::ResFloat,Cov::ResMatrix)
    # passed-in R0 is just sqrt(det(F)/det(H+...)) part
    R0*=pi^2*getki(t)^3*getkf(t)^3*cot(t.a[2]/2)cot(t.a[6]/2)/(4sin(t.a[2]/2)sin(t.a[6]/2)) # Popovici equations (9), and (5); but missing mono- and analyzer angular reflectivities
    if !all(getηs(t.sample).==0) # Account for sample mosaic
        M=inv(Cov) # The resolution matrix is the inverse of the covariance matrix (Popovici)
        # S. A. Werner & R. Pynn, J. Appl. Phys. 42, 4736, (1971), eq 19
        nQηs²=norm(getQ(t))*getηs(t.sample).^2 # getηs returns [ηH,ηH] as Gaussian widths, so there's no need to convert from FWHM
        #R0/=sqrt((1+(nQηs²[2]*M[2,2]))*(1+(nQηs²[1]*M[3,3]))) # this identically matches Zheludev equation 12, which appears to contain a typo
        R0/=sqrt(1+(nQηs²[1]*M[2,2]))
        R0/=sqrt(1+(nQηs²[2]*M[3,3]))
        Cov[2,2]+=nQηs²[1]
        Cov[3,3]+=nQηs²[2]
    end
    M=inv(Cov) # The resolution matrix is the inverse of the covariance matrix (Popovici)
    # Transform prefactor to Chesser-Axe normalization
    R0*=sqrt(det(M))/(2pi)^2
    # Include kf/ki part of cross section
    R0*=getkf(t)/getki(t)
#    if (t.ana.kinref!=0)&&(t.ana.size[3]!=0)
#        # not fully tested correction for analyzer reflectivity (though it seems to work)
#        # ResLib used separate analyzer thickness field, why not just use depth?
#        αs=getαs(t); # horizontal divergences expressed as Gaussian widths
#        #dth=(collect(1:201)./200).*sqrt(2log(2))*((αs[4]<αs[3])?αs[4]:αs[3]) # linspace(0,HWHM,200)
#        # linspace(0,1,N) is slower than (1:N+1)/N by a factor of 4 if no collection is done.
#        # Collecting reverses this, with collect(linspace(...)) being ~30% faster than collect(1:N+1)/N
#        dth=linspace(0,sqrt(2log(2))*(αs[4]<αs[3]?αs[4]:αs[3]),200) # linspace(0,HWHM_collimation,200)
#        wdth=exp(-dth.^2/2/getηh(t.ana)^2) # Gaussian due to analyzer mosaic
#        rdth=1./(1+1./(t.ana.size[3]*t.ana.kinref*(gettau(t.ana)/2)/sqrt(getkf(t)^2-(gettau(t.ana)/2)^2)*wdth/getηh(t.ana)/sqrt(2*pi)))
#        R0*=sum(rdth)/sum(wdth)
#    end
    return (R0,M)
end
"""
Implement the Cooper-Nathans method to approximate the resolution function of a `TripleAxis` object.
"""
function _innercoopernathansR0M(t::TripleAxis)
    A,B=_makeAB(ResFloat);C=_makeC(ResFloat);F,G=_makeFG(ResFloat)
        _fillAB!(t,A,B);    _fillC!(t,C);        _fillFG!(t,F,G)
    # A,B,C,F,G *should* be real matricies, so the distinction between
    # conjugate transposition (') and plain transposition (.') probably doesn't matter
    H=G+C'*F*C       # Popovici Eq 8
    #Ninv=A*inv(H)*A' # Cooper-Nathans (Popovici Eq 10)
    #Cov=B*Ninv*B'   # Popovici Eq 3
    Cov=B*A*inv(H)*A'*B'
    #Cov=B*A/H*A'*B'# A*inv(H) is (perhaps) less accurate than A/H
    #R0=sqrt(det(F)/det(H))    # Cooper-Nathans (Popovici Eq 5 and 9)
    R0=sqrt(det(F)/det(H)) # F is 4x4, H is 8x8, |aN|=a^n|N| if N is nxn. ResLibCal is missing (8log(2))^-2 because it uses FWHM
    R0::ResFloat,M::ResMatrix=_calcR0M(t,R0,Cov)
    # Popovici's formalism has x2 (at the sample) bisecting the sample scattering angle 2θₛ.
    # this likely isn't the case for our sample coordinate system, so rotate it
    Q=getQ(t)
    ra=atan2(dot(Q,gety(t)),dot(Q,getx(t)))
    rot4=basisrotation4xy(ra) # rotate the coordinate system by ra in the (Qx,Qy) plane
    M=rot4'*M*rot4 # rotate like ResLib's ResMatS
    return (R0,M,C,F,G) # Matricies C,F,G are used to convert to monitor counts
end
"""
    coopernathansR0M(t::TripleAxis)
returns a `ResRM` object containing the source-normalized resolution prefactor, `R₀`,
and the Cooper-Nathans approximation for the resolution matrix, `M`, describing a 4-D ellipsoid.
"""
coopernathansR0M(t::TripleAxis)=ResRM(_innercoopernathansR0M(t)...)
"""
    coopernathansRmM(t::TripleAxis)
returns a `ResRM` object containing the monitor-normalized resolution prefactor, `Rₘ`,
and the Cooper-Nathans approximation for the resolution matrix, `M`, describing a 4-D ellipsoid.
"""
coopernathansRmM(t::TripleAxis)=ResRM(convertcoopernathansR02Rm(t,_innercoopernathansR0M(t)...)...)
#@vectorize_1arg TripleAxis coopernathansR0M
#@vectorize_1arg TripleAxis coopernathansRmM
"""
Implement the Popovici method to approximate the resolution function of a `TripleAxis` object.
"""
function _innerpopoviciR0M(t::TripleAxis)
    # The Popovici method includes sample shape effects by considering the mean
    # squared displacements of sample geometry. If the sample is not isotropic
    # then for each point the *shape* must be rotated to account for the
    # physical rotation of the sample to reach a particular (Q,E) point.
    #
    # We then need to rotate the description of the shape into a coordinate
    # system where x is along Q∥, y is along Q⟂, and z is (still) vertical.
    # Further rotations of the description axes will be performed while
    # calculating the S matrix, following Popovici's formalism with x2 bisecting
    # the sample scattering angle 2θs.
    #
    # Instead of rotating the sample shape by a3 and then rotating the axes
    # of that description by a second angle, we can accomplish both by only
    # rotating the axes by the angle between Q and orient-1 (sample x). This is
    # appropriate because the sample orienting vectors and Q are *always*
    # described in the same reciprocal lattice.
    Q=getQ(t)
    ra=atan2(dot(Q,gety(t)),dot(Q,getx(t))) # this the angle between sample-x and Q
    rot3=basisrotation3xy(ra) # rotate from a coordinate system with x along orient1 by ra to put x along Q
    inputshape = t.sample.shape
    t.sample.shape=rot3*inputshape*rot3.' # this correctly changes the axes of inputshape to have x along Q.
    A,B=_makeAB(ResFloat);F,G=_makeFG(ResFloat);D,S,T=_makeDST(ResFloat)
        _fillAB!(t,A,B);      _fillFG!(t,F,G);        _fillDST!(t,D,S,T)
    t.sample.shape=inputshape # now that S is calculated put the sample shape back to where it started
    # A,B,F,G,D,S,T *should* be real matricies, so the distinction between
    # conjugate transposition (') and plain transposition (.') probably doesn't matter
    H=inv(D*inv(S+T.'*F*T)*D.')
    #Ninv=A*inv(H+G)*A' # Popovici equation 20
    #Cov=B*Ninv*B'     # Popovici equation 3
    Cov=B*A*inv(H+G)*A.'*B'
    #R0=sqrt(det(F)/det(H+G)) # Popovici
    R0=sqrt(det(F)/det(H+G)) # F is 4x4, H is 8x8, |aN|=a^n|N| if N is nxn. ResLibCal is missing (8log(2))^-2 because it uses FWHM
    R0::ResFloat,M::ResMatrix=_calcR0M(t,R0,Cov)
    rot4=basisrotation4xy(ra) # basis rotation for the resolution matrix
    M=rot4'*M*rot4 # rotate like ResLib's ResMatS
    return (R0,M,F,G)
end
"""
    popoviciR0M(t::TripleAxis)
returns a `ResRM` object containing the source-normalized resolution prefactor, `R₀`,
and the Popovici approximation for the resolution matrix, `M`, describing a 4-D ellipsoid.
"""
popoviciR0M(t::TripleAxis)=ResRM(_innerpopoviciR0M(t)...)
"""
    popoviciRmM(t::TripleAxis)
returns a `ResRM` object containing the monitor-normalized resolution prefactor, `Rₘ`,
and the Popovici approximation for the resolution matrix, `M`, describing a 4-D ellipsoid.
"""
popoviciRmM(t::TripleAxis)=ResRM(convertpopoviciR02Rm(t,_innerpopoviciR0M(t)...)...)
#@vectorize_1arg TripleAxis popoviciR0M
#@vectorize_1arg TripleAxis popoviciRmM

_makeAB(T::DataType)=(zeros(T,6,8),zeros(T,4,6))                   # A:6x8, B:4x6
_makeC(T::DataType)=zeros(T,4,8)                                   # C:4x8
_makeFG(T::DataType)=(zeros(T,4,4),zeros(T,8,8))                   # F:4x4, G:8x8
_makeDST(T::DataType)=(zeros(T,8,13),zeros(T,13,13),zeros(T,4,13)) # D:8x13, S:13x13, T:4x13
"""Functions to make empty matricies for Cooper-Nathans and/or Popovici resolution approximation"""
_makeAB,_makeC,_makeFG,_makeDST

function _fillAB!(t::TripleAxis,A::ResMatrix,B::ResMatrix)
    # fill matrix A
    A[1,1]=getki(t)*cot(t.a[2]/2)/2
    A[1,2]=-A[1,1]
    A[2,2]=getki(t)
    A[3,4]=getki(t)
    A[4,5]=getkf(t)*cot(t.a[6]/2)/2
    A[4,6]=-A[4,5]
    A[5,5]=getkf(t)
    A[6,7]=getkf(t)
    # fill matrix B
    B[1,1]= cos(getphi(t))
    B[1,2]= sin(getphi(t))
    B[1,4]=-cos(getphi(t)-t.a[4])
    B[1,5]=-sin(getphi(t)-t.a[4])
    B[2,1]=-B[1,2]
    B[2,2]= B[1,1]
    B[2,4]=-B[1,5]
    B[2,5]= B[1,4]
    B[3,3]= 1
    B[3,6]=-1
    B[4,1]= 2*ħ²2mₙ*getki(t)
    B[4,4]=-2*ħ²2mₙ*getkf(t)
end
function _fillC!(t::TripleAxis,C::ResMatrix)
    C[1,1]= 1/2
    C[1,2]= 1/2
    C[2,3]= 1/(2sin(t.a[2]/2))
    C[2,4]=-C[2,3] # mistake in paper
    C[3,5]= 1/2
    C[3,6]= 1/2
    C[4,7]= 1/(2sin(t.a[6]/2))
    C[4,8]=-C[4,7]
end
function _fillFG!(t::TripleAxis,F::ResMatrix,G::ResMatrix)
    # unused "horizontal focussing analyzer" correction
    # actually intended for if the analyzer is multi-bladed, e.g. RITA-II or  MACS
    # instead use the Popovici approximation.
    # check the conversion factor if this ever gets implemented:
    # (t.horifoc>0) && (αs[3]/=sqrt(12))
    #
    # fill matrix F (4x4) F[1:5:16] selects diagonal
    # getηs returns horizontal and vertical mosaics as Gaussian widths, σ
    F[1:5:16]=1./[getηs(t.mono);getηs(t.ana)].^2
    # fill matrix G (8x8) G[1:9:64] selects diagonal
    αs=getαs(t); βs=getβs(t); # these are σ (t.horcol and t.vercol are FWHM)
    G[1:9:64]=1./[αs[1:2];βs[1:2];αs[3:4];βs[3:4]].^2
end
function _fillDST!(t::TripleAxis,D::ResMatrix,S::ResMatrix,T::ResMatrix)
    # fill matrix D
    D[1,1] =-1/t.arms[1]
    D[1,3] =-cos(t.a[2]/2)/t.arms[1]
    D[1,4] = sin(t.a[2]/2)/t.arms[1]
    D[3,2] = D[1,1]
    D[3,5] =-D[1,1]
    D[2,3] = cos(t.a[2]/2)/t.arms[2]
    D[2,4] = sin(t.a[2]/2)/t.arms[2]
    D[2,6] = sin(t.a[4]/2)/t.arms[2]
    D[2,7] = cos(t.a[4]/2)/t.arms[2]
    D[4,5] =-1/t.arms[2]
    D[4,8] =-D[4,5]
    D[5,6] = sin(t.a[4]/2)/t.arms[3]
    D[5,7] =-cos(t.a[4]/2)/t.arms[3]
    D[5,9] =-cos(t.a[6]/2)/t.arms[3]
    D[5,10]= sin(t.a[6]/2)/t.arms[3]
    D[7,8] =-1/t.arms[3]
    D[7,11]=-D[7,8]
    D[6,9] = cos(t.a[6]/2)/t.arms[4]
    D[6,10]= sin(t.a[6]/2)/t.arms[4]
    D[6,12]= 1/t.arms[4]
    D[8,11]=-D[6,12]
    D[8,13]= D[6,12]
    # Define matrix S, for which it is necessary to rotate the sample shape
    # since (here) t.sample.shape is defined with x along Q and we need here x
    # along the bisector of ki and kf (Popovici's x₂)
    rot=basisrotation3xy(t.a[4]/2-getphi(t)) # phi is the angle between Q and ki, so this rotates x_sample along x2 (bisecting 2θs)
    # S=inv(S) reasigns the label S without changing the pointer that was passed in.
    # Since we define S by elements of its inverse ("Sinv") we need to allocate a matrix here to enable saving inv(Sinv) into S.
    Sinv=zeros(S) # makes a zero-filled matrix of the same size and type as S; (similar(S) would create a matrix the same type and size as S, but would not zero the elements)
    Sinv[1:2,1:2]=getshape(t.source)          # 2x2
    Sinv[3:5,3:5]=getshape(t.mono)            # 3x3
    Sinv[6:8,6:8]=rot*getshape(t.sample)*rot.'# 3x3 # correctly rotates shape description from x_sample to x2
    Sinv[9:11,9:11]=getshape(t.ana)           # 3x3
    Sinv[12:13,12:13]=getshape(t.detector)    # 2x2
    @assert isInvertable(Sinv) "The shape matrix is not invertable. Are all sizes reasonable?"
    Sinv=inv(Sinv)
    S[:]=Sinv[:]
    # fill matrix T
    T[1,1] =-1/(2t.arms[1])  #mistake in paper
    T[1,3] = cos(t.a[2]/2)*(1/t.arms[2]-1/t.arms[1])/2
    T[1,4] = sin(t.a[2]/2)*(1/t.arms[2]+1/t.arms[1]-2/(getmonoHR(t)*sin(t.a[2]/2)))/2
    T[1,6] = sin(t.a[4]/2)/(2t.arms[2])
    T[1,7] = cos(t.a[4]/2)/(2t.arms[2])
    T[2,2] =-1/(2t.arms[1]*sin(t.a[2]/2))
    T[2,5] = (1/t.arms[1]+1/t.arms[2]-2sin(t.a[2]/2)/getmonoVR(t))/(2sin(t.a[2]/2))
    T[2,8] =-1/(2t.arms[2]*sin(t.a[2]/2))
    T[3,6] = sin(t.a[4]/2)/(2t.arms[3])
    T[3,7] =-cos(t.a[4]/2)/(2t.arms[3])
    T[3,9] = cos(t.a[6]/2)*(1/t.arms[4]-1/t.arms[3])/2
    T[3,10]= sin(t.a[6]/2)*(1/t.arms[3]+1/t.arms[4]-2/(getanaHR(t)*sin(t.a[6]/2)))/2
    T[3,12]= 1/(2t.arms[4])
    T[4,8] =-1/(2t.arms[3]*sin(t.a[6]/2))
    T[4,11]= (1/t.arms[3]+1/t.arms[4]-2sin(t.a[6]/2)/getanaVR(t))/(2sin(t.a[6]/2))
    T[4,13]=-1/(2t.arms[4]*sin(t.a[6]/2))
end


# 3-matrix form
function convertcoopernathansR02Rm(t::TripleAxis,R0::ResFloat,C::ResMatrix,F::ResMatrix,G::ResMatrix)
    c=C[1:2,1:4]
    f=F[1:2,1:2]
    g=G[1:4,1:4]
    Rmon=(1/getki(t))*pi*getki(t)^3*cot(t.a[2]/2)/(2sin(t.a[2]/2))*sqrt(det(f)/det(g+c'*f*c)) # Zheludev's eq 9
    R0/=Rmon
    return R0
end
# 0-matrix form
function convertcoopernathansR02Rm(t::TripleAxis,R0::ResFloat)
    C=_makeC(ResFloat); F,G=_makeFG(ResFloat)
    _fillC!(t,C);_fillFG!(t,F,G)
    convertcoopernathansR02Rm(t,R0,C,F,G)
end
# 4-matrix form
convertcoopernathansR02Rm(t::TripleAxis,R0::ResFloat,M::ResMatrix,C::ResMatrix,F::ResMatrix,G::ResMatrix)=(convertcoopernathansR02Rm(t,R0,C,F,G),M)
# 1-matrix form
convertcoopernathansR02Rm(t::TripleAxis,R0::ResFloat,M::ResMatrix)=(convertcoopernathansR02Rm(t,R0),M)
# 1-ResRM, 0-matrix form; returns a ResRM (plus Array versions)
convertcoopernathansR02Rm(t::TripleAxis,RM::ResRM)=ResRM(convertcoopernathansR02Rm(t,RM.R),RM.M)
convertcoopernathansR02Rm{R<:ResRM}(t::TripleAxis,RM::Array{R})=map(x->convertcoopernathansR02Rm(t,x),RM)
convertcoopernathansR02Rm{T<:TripleAxis,R<:ResRM}(t::Array{T},RM::Array{R})=(@assert compatible(t,RM);map(convertcoopernathansR02Rm,t,RM))
# 1-ResRM, 0-matrix form; modifies Rm.R in place (plus Array versions)
convertcoopernathansR02Rm!(t::TripleAxis,RM::ResRM)=(RM.R=convertcoopernathansR02Rm(t,RM.R))
convertcoopernathansR02Rm!{R<:ResRM}(t::TripleAxis,RM::Array{R})=map(x->convertcoopernathansR02Rm!(t,x),RM)
convertcoopernathansR02Rm!{T<:TripleAxis,R<:ResRM}(t::Array{T},RM::Array{R})=(@assert compatible(t,RM);map(convertcoopernathansR02Rm!,t,RM))

"""
Functions to convert from a source-normalized resolution prefactor to a monitor-normalized
resolution prefactor in the Cooper-Nathans approximation.
"""
convertcoopernathansR02Rm,convertcoopernathansR02Rm!

_maketsd(T::DataType)=(zeros(T,2,7),zeros(T,7,7),zeros(T,4,7)) # t:2x7, s:7x7; d:4x7
function _filltsd!(ta::TripleAxis,t::ResMatrix,s::ResMatrix,d::ResMatrix)
     # Pull out relevant parameters from the TripleAxis structure
    L0=ta.arms[1]
    Lm=(length(ta.arms)>4) ? ta.arms[5] : ta.arms[2]
    # fill matrix t
    t[1,1]=-1/(2L0)  #mistake in paper
    t[1,3]= cos(ta.a[2]/2)*(1/Lm-1/L0)/2
    t[1,4]= sin(ta.a[2]/2)*(1/L0+1/Lm-2/(getmonoHR(ta)*sin(ta.a[2]/2)))/2
    t[1,7]= 1/(2Lm)
    t[2,2]=-1/(2L0*sin(ta.a[2]/2))
    t[2,5]=(1/L0+1/Lm-2sin(ta.a[2]/2)/getmonoVR(ta))/(2sin(ta.a[2]/2))
    # fill matrix s
    #s[:]=0 # just to be sure
    sinv=zeros(s)
    sinv[1:2,1:2]=getshape(ta.source)
    sinv[3:5,3:5]=getshape(ta.mono)
    sinv[6:7,6:7]=getshape(ta.monitor)
    @assert isInvertable(sinv) "The pre-sample shape matrix is not invertable. Are source, monochromator and monitor sizes reasonable?"
    sinv=inv(sinv) # invert s^-1, this is s
    s[:]=sinv[:]
    # fill matrix d
    d[1,1]=-1/L0
    d[1,3]=-cos(ta.a[2]/2)/L0
    d[1,4]= sin(ta.a[2]/2)/L0
    d[3,2]= d[1,1]
    d[3,5]=-d[1,1]
    d[2,3]= cos(ta.a[2]/2)/Lm
    d[2,4]= sin(ta.a[2]/2)/Lm
    d[2,7]= 1/Lm
    d[4,5]=-1/Lm
end

# 5-matrix form
function convertpopoviciR02Rm(ta::TripleAxis,R0::ResFloat,F::ResMatrix,G::ResMatrix,t::ResMatrix,s::ResMatrix,d::ResMatrix)
    f=F[1:2,1:2]
    g=G[1:4,1:4]
    _filltsd!(ta,t,s,d)
    Rmon=(1/getki(ta))*pi*getki(ta)^3*cot(ta.a[2]/2)/(2sin(ta.a[2]/2))*sqrt(det(f)/det(inv(d*inv(s+t'*f*t)*d')+g)) # Zheludev's eq 10
    R0/=Rmon
    return R0::ResFloat
end
# 2-matrix form
function convertpopoviciR02Rm(ta::TripleAxis,R0::ResFloat,F::ResMatrix,G::ResMatrix)
    t,s,d=_maketsd(ResFloat)
    convertpopoviciR02Rm(ta,R0,F,G,t,s,d)::ResFloat
end
# 0-matrix form
function convertpopoviciR02Rm(t::TripleAxis,R0::ResFloat)
    F,G=_makeFG(ResFloat);_fillFG!(t,F,G);
    convertpopoviciR02Rm(t,R0,F,G)::ResFloat
end
# 6-matrix form
convertpopoviciR02Rm(ta::TripleAxis,R0::ResFloat,M::ResMatrix,F::ResMatrix,G::ResMatrix,t::ResMatrix,s::ResMatrix,d::ResMatrix)=(convertpopoviciR02Rm(ta,R0,F,G,t,s,d)::ResFloat,M)
# 3-matrix form <- likely used variant
convertpopoviciR02Rm(t::TripleAxis,R0::ResFloat,M::ResMatrix,F::ResMatrix,G::ResMatrix)=(convertpopoviciR02Rm(t,R0,F,G)::ResFloat,M)
# 1-matrix form
convertpopoviciR02Rm(t::TripleAxis,R0::ResFloat,M::ResMatrix)=(convertpopoviciR02Rm(t,R0)::ResFloat,M)
# 1-ResRM, 0-matrix form; returns a ResRM (plus Array versions)
convertpopoviciR02Rm(t::TripleAxis,RM::ResRM)=ResRM(convertpopoviciR02Rm(t,RM.R),RM.M)
convertpopoviciR02Rm{R<:ResRM}(t::TripleAxis,RM::Array{R})=map(x->convertpopoviciR02Rm(t,x),RM)
convertpopoviciR02Rm{T<:TripleAxis,R<:ResRM}(t::Array{T},RM::Array{R})=(@assert compatible(t,RM);map(convertpopoviciR02Rm,t,RM))
# 1-ResRM, 0-matrix form; modifies RM.R in place (plus Array versions)
convertpopoviciR02Rm!(t::TripleAxis,RM::ResRM)=(RM.R=convertpopoviciR02Rm(t,RM.R))
convertpopoviciR02Rm!{R<:ResRM}(t::TripleAxis,RM::Array{R})=map(x->convertpopoviciR02Rm!(t,x),RM)
convertpopoviciR02Rm!{T<:TripleAxis,R<:ResRM}(t::Array{T},RM::Array{R})=(@assert compatible(t,RM);map(convertpopoviciR02Rm!,t,RM))
"""
Functions to convert from a source-normalized resolution prefactor to a monitor-normalized
resolution prefactor in the Popovici approximation.
"""
convertpopoviciR02Rm,convertpopoviciR02Rm!

# overloading for Scatterd objects:
coopernathansR0M(a::Scatterd)=coopernathansR0M(a.instrument)
coopernathansRmM(a::Scatterd)=coopernathansRmM(a.instrument)
popoviciR0M(a::Scatterd)=popoviciR0M(a.instrument)
popoviciRmM(a::Scatterd)=popoviciRmM(a.instrument)
