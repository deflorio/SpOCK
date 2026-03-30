!#######################################################################
! MSIS® (NRL-SOF-014-1) SOFTWARE
! NRLMSIS® empirical atmospheric model software. Use is governed by the
! Open Source Academic Research License Agreement contained in the file
! nrlmsis2.1_license.txt, which is part of this software package. BY
! USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
! CONDITIONS OF THE LICENSE.  
!#######################################################################

!!! ===========================================================================
!!! NRLMSIS 2.1:
!!! Neutral atmosphere empirical model from the surface to lower exosphere
!!! ===========================================================================
!**************************************************************************************************
! MSIS_CONSTANTS Module: Contains constants and hardwired parameters
!**************************************************************************************************
module msis_constants

  implicit none

  ! Floating Point Precision
#ifdef DBLE
  integer, parameter         :: rp = 8
#else
  integer, parameter         :: rp = 4
#endif

  ! Missing density value
  real(kind=rp),parameter    :: dmissing = 9.999e-38_rp

  ! Trigonometric constants
  real(kind=rp), parameter   :: pi = 3.1415926535897932384626433832795_rp
  real(kind=rp), parameter   :: deg2rad = pi / 180.0_rp
  real(kind=rp), parameter   :: doy2rad = 2.0_rp*pi / 365.0_rp
  real(kind=rp), parameter   :: lst2rad = pi / 12.0_rp
  !real(kind=rp), parameter   :: tanh1 = 0.761594155955765485_rp  ! tanh(1.0)
  real(kind=rp), parameter   :: tanh1 = tanh(1.0_rp)

  ! Thermodynamic constants
  ! Boltzmann constant (CODATA 2018) (J/kg)
  real(kind=rp), parameter   :: kB = 1.380649e-23_rp
  ! Avogadro constant (CODATA 2018)
  real(kind=rp), parameter   :: NA = 6.02214076e23_rp
  ! Reference gravity (CIMO Guide 2014) (m/s^2) (specified separately in alt2gph, in msis_utils.F90)
  real(kind=rp), parameter   :: g0 = 9.80665_rp
  ! Species molecular masses (kg/molecule) (CIPM 2007)
  real(kind=rp), parameter   :: specmass(1:10) = (/  0.0_rp,                          & ! Mass density (dummy value)
                                                    28.0134_rp,                       & ! N2
                                                    31.9988_rp,                       & ! O2
                                                    31.9988_rp/2.0_rp,                & ! O
                                                     4.0_rp,                          & ! He
                                                     1.0_rp,                          & ! H
                                                    39.948_rp,                        & ! Ar
                                                    28.0134_rp/2.0_rp,                & ! N
                                                    31.9988_rp/2.0_rp,                & ! Anomalous O
                                                    (28.0134_rp+31.9988_rp)/2.0_rp /) & ! NO
                                                    / (1.0e3_rp * NA)                   ! Convert from g/mol to kg/molecule
  ! Dry air mean mass in fully mixed atmosphere (CIPM 2007) (includes CO2 and other trace species that are not yet in MSIS)
  real(kind=rp), parameter   :: Mbar = 28.96546_rp / (1.0e3_rp * NA)   ! kg/molecule
  ! Dry air log volume mixing ratios (CIPM 2007)
  real(kind=rp), parameter   :: lnvmr(1:10) = log( (/ 1.0_rp,        & ! Mass density (dummy value)
                                                      0.780848_rp,   & ! N2
                                                      0.209390_rp,   & ! O2
                                                      1.0_rp,        & ! O (dummy value)
                                                      0.0000052_rp,  & ! He
                                                      1.0_rp,        & ! H (dummy value)
                                                      0.009332_rp,   & ! Ar
                                                      1.0_rp,        & ! N (dummy value)
                                                      1.0_rp,        & ! Anomalous O (dummy value)
                                                      1.0_rp /) )      ! NO (dummy value)
  ! Natural log of global average surface pressure (Pa)
  !real(kind=rp), parameter   :: lnP0 = 11.5080482 !+ 0.00759597 After calibration with MERRA2
  real(kind=rp), parameter   :: lnP0 = 11.515614
  ! Derived constants
  real(kind=rp), parameter   :: g0divkB = g0/kB * 1.0e3_rp  ! K/(kg km)
  real(kind=rp), parameter   :: Mbarg0divkB = Mbar*g0/kB * 1.0e3_rp   ! K/km
  ! References:
  ! CODATA Internationally recommended 2018 values of the fundamental physical constants.
  !   https://pml.nist.gov/cuu/Constants/; https://pml.nist.gov/cuu/pdf/wallet_2018.pdf
  ! Picard, A., Davis, R. S., Glaeser, M., and Fujii, K. (2007). Revised formula for the density of
  !   air (CIPM 2007). Metrologia 45, 149–155. doi:10.1088/0026-1394/45/2/004
  ! World Meteorological Organization (2014). WMO guide to meteorological instruments and methods of observation
  !   (the CIMO Guide). Part I, Chapter 12. https://www.wmo.int/pages/prog/www/IMOP/CIMO-Guide.html

  ! Vertical profile parameters
  integer, parameter         :: nspec = 11  !Number of species including temperature
  integer, parameter         :: nd = 27     !Number of temperature profile nodes
  integer, parameter         :: p = 4       !Spline order
  integer, parameter         :: nl = nd - p !Last temperature profile level index
  integer, parameter         :: nls = 9     !Last parameter index for each species (excluding O, NO splines)
  real(kind=rp), parameter   :: bwalt = 122.5_rp ! Reference geopotential height for Bates Profile
  real(kind=rp), parameter   :: zetaF = 70.0_rp  ! Fully mixed below this, uses constant mixing ratios
  real(kind=rp), parameter   :: zetaB = bwalt    ! Bates Profile above this altitude
  real(kind=rp), parameter   :: zetaA = 85.0_rp  ! Default reference height for active minor species
  real(kind=rp), parameter   :: zetagamma = 100.0_rp  ! Reference height of tanh taper of chemical/dynamical correction scale height
  real(kind=rp), parameter   :: Hgamma = 1.0_rp/30.0_rp  ! Inverse scale height of tanh taper of chemical/dynamical correction scale height
  real(kind=rp), parameter   :: nodesTN(0:nd+2) = &  !Nodes for temperature profile splines
      (/ -15., -10.,  -5.,   0.,   5., 10., 15., 20.,  25.,  30.,  35.,  40., 45., 50., &
          55.,  60.,  65.,  70.,  75., 80., 85., 92.5, 102.5, 112.5, 122.5, 132.5, 142.5, &
          152.5, 162.5, 172.5/)
  integer, parameter         :: izfmx = 13       ! fully mixed below this spline index
  integer, parameter         :: izfx = 14        ! Spline index at zetaF
  integer, parameter         :: izax = 17        ! Spline index at zetaA
  integer, parameter         :: itex = nl        ! Index of Bates exospheric temperature
  integer, parameter         :: itgb0 = nl - 1   ! Index of Bates temperature gradient at lower boundary
  integer, parameter         :: itb0 = nl - 2    ! Index of Bates temperature at lower boundary
  ! O1 Spline parameters
  integer, parameter         :: ndO1 = 13
  integer, parameter         :: nsplO1 = ndO1-5     !Number of unconstrained spline parameters for O1 (there are 2 additional C1-constrained splines)
  real(kind=rp), parameter   :: nodesO1(0:ndO1) = & !Nodes for O1 splines (Domain 50-85 km)
      (/ 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 92.5, 102.5, 112.5/)
  real(kind=rp), parameter   :: zetarefO1 = zetaA   !Joining height for O1 splines, and reference height for O1 density
  ! NO Spline parameters
  integer, parameter         :: ndNO = 13
  integer, parameter         :: nsplNO = ndNO-5     !Number of unconstrained spline parameters for NO (there are 2 additional C1-constrained splines)
  real(kind=rp), parameter   :: nodesNO(0:ndNO) = & !Nodes for NO splines (Domain 70-122.5 km)
      (/ 47.5, 55., 62.5, 70., 77.5, 85., 92.5, 100., 107.5, 115., 122.5, 130., 137.5, 145./)
  real(kind=rp), parameter   :: zetarefNO = zetaB   !Joining height for NO splines, and reference height for NO density
  !C2 Continuity matrix for temperature; Last 3 splines are constrained (must be recomputed if nodes change)
  real(kind=rp), parameter   :: c2tn(3,3) = reshape((/1.0_rp, -10.0_rp,  33.333333333333336_rp, &
                                                      1.0_rp,   0.0_rp, -16.666666666666668_rp, &
                                                      1.0_rp,  10.0_rp,  33.333333333333336_rp/), &
                                                    (/3,3/))
  !C1 Continuity for O1; Last 2 splines are constrained (must be recomputed if nodes change)
  real(kind=rp), parameter   :: c1o1(2,2) = reshape((/ 1.75_rp,               -2.916666573405061_rp, &
                                                      -1.624999900076852_rp,  21.458332647194382_rp /), &
                                                    (/2,2/))
  real(kind=rp), parameter   :: c1o1adj(2) = (/0.257142857142857_rp, -0.102857142686844_rp/) !Weights for coefficents on 3rd to last spline; product to be subtracted from RHS of continuity equation
  !C1 Continuity for NO; Last 2 splines are constrained (must be recomputed if nodes change)
  real(kind=rp), parameter   :: c1NO(2,2) = reshape((/ 1.5_rp, -3.75_rp, &
                                                       0.0_rp,  15.0_rp /), &
                                                    (/2,2/))
  real(kind=rp), parameter   :: c1NOadj(2) = (/0.166666666666667_rp, -0.066666666666667_rp/) !Weights for coefficents on 3rd to last spline; product to be subtracted from RHS of continuity equation
  ! Anomalous Oxygen parameters (legacy profile from NRLMSISE-00)
  real(kind=rp),parameter    :: zetarefOA = zetaB   !Reference height for anomalous oxygen density
  real(kind=rp),parameter    :: TOA = 4000.         !Temperature of anomalous oxygen density (K)
  real(kind=rp),parameter    :: HOA = (kB * TOA) / ( (16.0_rp/(1.0e3_rp*NA)) * g0 ) * 1.0e-3_rp  !Hydrostatic scale height of anomalous oxygen density (km)
    
  ! Horizontal and time-dependent basis function (gfn) parameters
  integer, parameter      :: maxnbf = 512   ! Number of basis functions to be allocated
  integer, parameter      :: maxn = 6       ! Maximum latitude (Legendre) spectral degree
  integer, parameter      :: maxl = 3       ! Maximum local time (tidal) spectral order
  integer, parameter      :: maxm = 2       ! Maximum longitude (stationary planetary wave) order
  integer, parameter      :: maxs = 2       ! Maximimum day of year (intra-annual) Fourier order
  integer, parameter      :: amaxn = 6      ! Maximum Legendre degree used in time independent and intra-annual zonal mean terms
  integer, parameter      :: amaxs = 2      ! Maximum intra-annual order used in zonal mean terms
  integer, parameter      :: tmaxl = 3      ! Maximum tidal order used
  integer, parameter      :: tmaxn = 6      ! Maximum Legendre degree coupled with tides
  integer, parameter      :: tmaxs = 2      ! Maximum intra-annual order coupled with tides
  integer, parameter      :: pmaxm = 2      ! Maximum stationary planetary wave order used
  integer, parameter      :: pmaxn = 6      ! Maximum Legendre degree coupled with SPW
  integer, parameter      :: pmaxs = 2      ! Maximum intra-annual order coupled with SPW
  integer, parameter      :: nsfx = 5       ! Number of linear solar flux terms
  integer, parameter      :: nsfxmod = 5    ! Number of nonlinear modulating solar flux terms (legacy NRLMSISE-00 terms)
  integer, parameter      :: nmag = 54      ! Number of terms in NRLMSISE-00 legacy geomagnetic parameterization
  integer, parameter      :: nut = 12       ! Number of terms in NRLMSISE-00 legacy UT parameterization
  integer, parameter      :: ctimeind = 0             ! Starting index of time-independent terms
  integer, parameter      :: cintann = ctimeind + (amaxn+1)   ! Starting index of zonal mean intra-annual terms
  integer, parameter      :: ctide = cintann + ((amaxn+1)*2*amaxs)   ! Starting index of zonal mean intra-annual terms
  integer, parameter      :: cspw = ctide + (4*tmaxs+2)*(tmaxl*(tmaxn+1)-(tmaxl*(tmaxl+1))/2) ! Starting index of SPW terms
  integer, parameter      :: csfx = cspw + (4*pmaxs+2)*(pmaxm*(pmaxn+1)-(pmaxm*(pmaxm+1))/2)   ! Starting index of linear solar flux terms
  integer, parameter      :: cextra = csfx + nsfx     ! Starting index of time-independent terms
  integer, parameter      :: mbf = 383                ! Last index of linear terms
  integer, parameter      :: cnonlin = mbf + 1        ! Starting index of nonlinear terms
  integer, parameter      :: csfxmod = cnonlin        ! Starting index of modulating solar flux terms
  integer, parameter      :: cmag = csfxmod + nsfxmod ! Starting index of daily geomagnetic terms
  integer, parameter      :: cut = cmag + nmag        ! Starting index of UT terms
    
  ! Weights for calculation log pressure spline coefficients from temperature coefficients (must be recalcuated if nodes change)
  real(kind=rp), parameter   :: gwht(0:3) =  (/ 5.0_rp/24.0_rp, 55.0_rp/24.0_rp, 55.0_rp/24.0_rp, 5.0_rp/24.0_rp /)

  ! Constants needed for analytical integration by parts of hydrostatic piecewise effective mass profile
  real(kind=rp), parameter   :: wbeta(0:nl) =  (nodesTN(4:nd)  - nodesTN(0:nl)) / 4.0_rp !Weights for 1st spline integration
  real(kind=rp), parameter   :: wgamma(0:nl) = (nodesTN(5:nd+1)- nodesTN(0:nl)) / 5.0_rp !Weights for 2nd spline integration
  ! Non-zero bspline values at zetaB (5th and 6th order) (must be recalcuated if nodes change)
  real(kind=rp), parameter   :: S5zetaB(0:3) = (/0.041666666666667_rp, 0.458333333333333_rp, 0.458333333333333_rp, &
                                                 0.041666666666667_rp/)
  real(kind=rp), parameter   :: S6zetaB(0:4) = (/0.008771929824561_rp, 0.216228070175439_rp, 0.550000000000000_rp, &
                                                 0.216666666666667_rp, 0.008333333333333_rp/)
  !Weights for calculating temperature gradient at zetaA (must be recalcuated if nodes change)
  real(kind=rp), parameter   :: wghtAxdz(0:2) = (/-0.102857142857_rp, 0.0495238095238_rp, 0.053333333333_rp/)
  !Non-zero bspline values at zetaA (4th, 5th and 6th order) (must be recalcuated if nodes change)
  real(kind=rp), parameter   :: S4zetaA(0:2) = (/0.257142857142857_rp, 0.653968253968254_rp, 0.088888888888889_rp/)
  real(kind=rp), parameter   :: S5zetaA(0:3) = (/0.085714285714286_rp, 0.587590187590188_rp, 0.313020313020313_rp, &
                                                 0.013675213675214_rp/)
  real(kind=rp), parameter   :: S6zetaA(0:4) = (/0.023376623376623_rp, 0.378732378732379_rp, 0.500743700743701_rp, &
                                                 0.095538448479625_rp, 0.001608848667672_rp/)
  !Non-zero bspline values at zetaF (4th and 5th order) (must be recalcuated if nodes change)
  real(kind=rp), parameter   :: S4zetaF(0:2) = (/0.166666666666667_rp, 0.666666666666667_rp, 0.166666666666667_rp/)
  real(kind=rp), parameter   :: S5zetaF(0:3) = (/0.041666666666667_rp, 0.458333333333333_rp, 0.458333333333333_rp, &
                                                 0.041666666666667_rp/)
  !Non-zero bspline values at zeta=0 (5th order) (must be recalcuated if nodes change)
  real(kind=rp), parameter   :: S5zeta0(0:2) = (/0.458333333333333_rp, 0.458333333333333_rp, 0.041666666666667_rp/)

end module msis_constants
!**************************************************************************************************
! MSIS_UTILS Module: Contains the following auxiliary subroutines:
!  alt2gph:  Converts geodetic altitude to geopotential height
!  gph2alt:  Converts geopotential height to geodetic altitude
!  bspline:  Computes B-splines using input nodes and up to specified order
!  dilog:    Computes dilogarithm function (expansion truncated at order 3, error < 1E-5)
!**************************************************************************************************

module msis_utils

contains

  !==================================================================================================
  ! ALT2GPH: Altitude to Geopotential Height
  ! References:
  !   DMA Technical Report TR8350.2 (1987),
  !     http://earth-info.nga.mil/GandG/publications/historic/historic.html
  !   Featherstone, W. E., and S. J. Claessens (2008), Closed-form transformation between
  !     geodetic and ellipsoidal coordinates, Studia Geophysica et Geodaetica, 52, 1-18
  !   Jekeli, C. (2009), Potential theory and static gravity field of the Earth, in 
  !     Treatise on Geophysics, ed. T. Herring, vol 3, 11-42
  !   NIMA Technical Report TR8350.2 (2000, 3rd edition, Amendment1), 
  !     http://earth-info.nga.mil/GandG/publications/tr8350.2/tr8350_2.html
  !==================================================================================================
  real(8) function alt2gph(lat,alt)

    implicit none

    ! Input variables
    real(8), intent(in) :: lat    !Geodetic latitude (deg)
    real(8), intent(in) :: alt    !Geodetic altitude (km)

    real(8), parameter  :: deg2rad = 0.017453292519943295d0

    ! WGS84 Defining parameters
    real(8), parameter  :: a = 6378.1370d0 * 1d3 !Semi-major axis of reference ellipsoid (m)
    real(8), parameter  :: finv = 298.257223563d0 ! 1/f = Reciprocal of flattening
    real(8), parameter  :: w = 7292115d-11 !Angular velocity of Earth rotation (rad/s)
    real(8), parameter  :: GM = 398600.4418 * 1d9 !Gravitational constant x Earth mass (m^3/s^2)

    ! WGS84 Derived parameters
    real(8), parameter  :: asq = a*a
    real(8), parameter  :: wsq = w*w
    real(8), parameter  :: f = 1.0d0 / finv
    real(8), parameter  :: esq = 2*f - f*f
    real(8), parameter  :: e = sqrt(esq)  !Ellipsoid eccentricity
    real(8), parameter  :: Elin = a*e     !Linear eccentricity of ellipsoid
    real(8), parameter  :: Elinsq = Elin*Elin
    real(8), parameter  :: epr = e / (1-f)  !Second eccentricity
    real(8), parameter  :: q0 = ((1.0d0 + 3.0d0/(epr*epr))*atan(epr) - 3.0d0/epr)/2.0d0  !DMA Technical Report tr8350.2, Eq. 3-25
    real(8), parameter  :: U0 = -GM*atan(epr)/Elin - wsq*asq/3d0 !Theoretical potential of reference ellipsoid (m^2/s^2)
    real(8), parameter  :: g0 = 9.80665d0 !Standard gravity (m/s^2), CGPM 1901; WMO
    real(8), parameter  :: GMdivElin = GM / Elin
  
    ! Parameters for centrifugal potential taper
    real(8), parameter  :: x0sq = 2d7**2   !Axial distance squared at which tapering begins (m^2)
    real(8), parameter  :: Hsq = 1.2d7**2  !Relaxation scale length of taper (m^2)

    ! Working variables
    real(8)             :: altm, sinsqlat, v, xsq, zsq
    real(8)             :: rsqminElinsq, usq, cossqdelta, epru, atanepru, q, U, Vc

    ! Compute Cartesian and ellipsoidal coordinates
    altm = alt * 1000.0d0
    sinsqlat = sin(lat*deg2rad)**2
    v = a / sqrt(1-esq*sinsqlat)           !Radius of curvature of the reference ellipsoid, Featherstone eq. 4
    xsq = (v + altm)**2 * (1 - sinsqlat)   !Squared x-coordinate of geocentric system, Featherstone eq. 1
    zsq = (v*(1-esq) + altm)**2 * sinsqlat !Squared z-coordinate of geocentric system, Featherstone eq. 3
    rsqminElinsq = xsq + zsq - Elinsq
    usq = rsqminElinsq/2.0d0 + sqrt(rsqminElinsq**2 / 4.0d0 + Elinsq*zsq)  !Ellipsoidal distance coordinate, Featherstone eq. 19 
    cossqdelta = zsq / usq                 !Ellipsoidal polar angle, Featherstone eq. 21

    ! Compute gravitational potential
    epru = Elin / sqrt(usq)                !Second eccentricity at ellipsoidal coordinate u
    atanepru = atan(epru)
    q = ((1+3.0d0/(epru*epru))*atanepru - 3.0d0/epru)/2.0d0   !Jekeli, eq. 114
    U = -GMdivElin * atanepru - wsq * ( asq * q * (cossqdelta - 1/3.0d0) / q0 ) / 2.0d0   !Jekeli, eq. 113

    ! Compute centrifugal potential and adjust total potential
    if (xsq .le. x0sq) then
      Vc = (wsq/2.0d0) * xsq
    else
      Vc = (wsq/2.0d0) * (Hsq*tanh((xsq-x0sq)/Hsq) + x0sq) !Centrifugal potential taper
    endif
    U = U - Vc
  
    ! Compute geopotential height
    alt2gph = (U - U0) / g0 / 1000.0d0

    return

  end function alt2gph

  !==================================================================================================
  ! GPH2ALT: Geopotential Height to Altitude
  !==================================================================================================
  real(8) function gph2alt(theta,gph)

    implicit none

    real(8), intent(in)  :: theta
    real(8), intent(in)  :: gph

    integer, parameter   :: maxn = 10
    real(8), parameter   :: epsilon = 0.0005

    real(8)              :: x,dx,y,dydz
    integer              :: n

    x = gph
    n = 0
    dx = epsilon + epsilon
    do while ((abs(dx) .gt. epsilon) .and. (n .lt. 10))
      y = alt2gph(theta,x)
      dydz = (alt2gph(theta,x+dx) - y)/dx
      dx = (gph - y)/dydz
      x = x + dx
      n = n + 1
    end do

    gph2alt = x

  end function gph2alt

  !==================================================================================================
  ! BSPLINE: Returns array of nonzero b-spline values, for all orders up to specified order (max 6)
  !==================================================================================================
  subroutine bspline(x,nodes,nd,kmax,eta,S,i)

    use msis_constants, only:  rp

    implicit none

    ! Input variables
    real(kind=rp), intent(in)         :: x              !Location at which splines are to be evaluated
    real(kind=rp),dimension(0:),intent(in)  :: nodes    !Spline node locations
    integer, intent(in)               :: nd             !Number of spline nodes minus one (0:nd)
    integer, intent(in)               :: kmax           !Maximum order (up to 6 allowed) of evaluated splines
    real(kind=rp), intent(in)         :: eta(0:30,2:6)  !Precomputed weights for recursion (reciprocals of node differences)
    ! Ouput variables
    real(kind=rp), intent(out)        :: S(-5:0,2:6)    !b-spline values (spline index relative to i (-5:0), spline order (2:6))
    integer, intent(out)              :: i              !Index of last nonzero b-spline

    ! Working variables
    integer                           :: j, k, l
    integer                           :: low, high
    real(kind=rp)                     :: w(-4:0)        !Weights for recursion relation
 
    ! Initialize to zero
    S(:,:) = 0.0_rp

    ! Find index of last (rightmost) nonzero spline
    if (x .ge. nodes(nd)) then
      i = nd
      return
    endif
    if (x .le. nodes(0)) then
      i = -1
      return
    endif
    low = 0
    high = nd
    i = (low + high)/2
    do while (x .lt. nodes(i) .or. x .ge. nodes(i + 1))
        if (x .lt. nodes(i)) then
          high = i
        else
          low = i
        endif
        i = (low + high)/2
    end do

    ! Initialize with linear splines
    S(0,2) = (x - nodes(i)) * eta(i,2)
    if (i .gt. 0) S(-1,2) = 1 - S(0,2)
    if (i .ge. nd-1) S(0,2) = 0.0_rp   !Reset out-of-bounds spline to zero

    ! k = 3 (quadratic splines)
    w(:) = 0.0_rp
    w(0) = (x - nodes(i)) * eta(i,3)
    if (i .ne. 0) w(-1) = (x - nodes(i-1)) * eta(i-1,3)
    if (i .lt. (nd-2)) S(0,3) = w(0)*S(0,2)
    if ( ((i-1) .ge. 0) .and. ((i-1) .lt. (nd-2)) ) &
        S(-1,3) = w(-1) * S(-1,2) + (1.0_rp - w(0))*S(0,2)
    if ((i-2) .ge. 0) S(-2,3) = (1.0_rp - w(-1))*S(-1,2)
    
    ! k = 4 (cubic splines)
    do l = 0, -2, -1
      j = i + l
      if (j .lt. 0) exit  !Skip out-of-bounds splines
      w(l) = (x - nodes(j)) * eta(j,4)
    enddo
    if (i .lt. (nd-3)) S(0,4) = w(0)*S(0,3)
    do l = -1, -2, -1
        if ( ((i+l) .ge. 0) .and. ((i+l) .lt. (nd-3)) ) &
            S(l,4) = w(l)*S(l,3) + (1.0_rp - w(l+1))*S(l+1,3)
    enddo
    if ((i-3) .ge. 0) S(-3,4) = (1.0_rp - w(-2))*S(-2,3)
  
    ! k = 5
    do l = 0, -3, -1
      j = i + l
      if (j .lt. 0) exit  !Skip out-of-bounds splines
      w(l) = (x - nodes(j)) * eta(j,5)
    enddo
    if (i .lt. (nd-4)) S(0,5) = w(0)*S(0,4)
    do l = -1, -3, -1
        if ( ((i+l) .ge. 0) .and. ((i+l) .lt. (nd-4)) ) &
            S(l,5) = w(l)*S(l,4) + (1.0_rp - w(l+1))*S(l+1,4)
    enddo
    if ((i-4) .ge. 0) S(-4,5) = (1.0_rp - w(-3))*S(-3,4)
    if (kmax .eq. 5) return  !Exit if only 5th order spline is needed

    ! k = 6
    do l = 0, -4, -1
      j = i + l
      if (j .lt. 0) exit  !Skip out-of-bounds splines
      w(l) = (x - nodes(j)) * eta(j,6)
    enddo
    if (i .lt. (nd-5)) S(0,6) = w(0)*S(0,5)
    do l = -1, -4, -1
      if ( ((i+l) .ge. 0) .and. ((i+l) .lt. (nd-5)) ) &
          S(l,6) = w(l)*S(l,5) + (1.0_rp - w(l+1))*S(l+1,5)
    enddo
    if ((i-5) .ge. 0) S(-5,6) = (1.0_rp - w(-4))*S(-4,5)

    return

  end subroutine bspline

  !==================================================================================================
  ! DILOG: Calculate dilogarithm in the domain [0,1)
  ! Retains terms up to order 3 in the expansion, which results in relative errors less than 1E-5.
  ! Reference: 
  !   Ginsberg, E. S., and D. Zaborowski (1975), The Dilogarithm function of a real argument, 
  !   Commun. ACM, 18, 200–202.
  !==================================================================================================
  real(kind=rp) function dilog(x0)

    use msis_constants, only     : rp, pi

    implicit none

    real(kind=rp), intent(in)   :: x0
    real(kind=rp), parameter    :: pi2_6 = pi*pi / 6.0_rp
    real(kind=rp)               :: x, xx, x4, lnx

    x = x0
    if (x .gt. 0.5_rp) then
      lnx = log(x)
      x = 1.0_rp - x          !Reflect argument into [0,0.5] range
      xx = x*x
      x4 = 4.0_rp*x
      dilog = pi2_6 - lnx*log(x) &
              - (4.0_rp*xx*(23.0_rp/16.0_rp + x/36.0_rp + xx/576.0_rp + xx*x/3600.0_rp) &
                  + x4 + 3.0_rp*(1.0_rp - xx)*lnx) / (1.0_rp + x4 + xx)
    else
      xx = x*x
      x4 = 4.0_rp*x
      dilog = (4.0_rp*xx*(23.0_rp/16.0_rp + x/36.0_rp + xx/576.0_rp + xx*x/3600.0_rp) &
                + x4 + 3.0_rp*(1.0_rp - xx)*log(1.0_rp - x)) / (1.0_rp + x4 + xx)
    endif

    return

  end function dilog

end module msis_utils
!!! ===========================================================================
!!! NRLMSIS 2.1:
!!! Neutral atmosphere empirical model from the surface to lower exosphere
!!! ===========================================================================
!!!
!!! MSISINIT: Initialization of MSIS parameters, switches, and options.
!
!     PREREQUISITES:
!       MSIS binary parameter file (msis207.parm)
!
!     CALLING SEQUENCE:
!       CALL MSISINIT([OPTIONAL ARGUMENTS])
!  
!     OPTIONAL ARGUMENTS:
!       parmpath        File path pointing to the MSIS parameter file.
!                         Default: Null string (current directory)
!       parmfile        Name of MSIS parameter file.
!                         Default: 'msis21.parm'
!       iun             File unit number for reading parameter file.
!                         Default: 67
!       switch_gfn      Logical array of 512 swtiches for individual terms. For
!                         advanced users.
!                         Default values: True (all switches on)
!       switch_legacy   Floating point array (1:25) of legacy switches that
!                         control groups of terms:
!                            1 - F10.7
!                            2 - Time independent
!                            3 - Symmetrical annual
!                            4 - Symmetrical semiannual
!                            5 - Asymmetrical annual
!                            6 - Asymmetrical semiannual
!                            7 - Diurnal
!                            8 - Semidiurnal
!                            9 - Geomagnetic activity:
!                                  1.0 = Daily Ap mode
!                                 -1.0 = Storm-time ap mode
!                           10 - All UT/long effects
!                           11 - Longitudinal
!                           12 - UT and mixed UT/long
!                           13 - Mixed Ap/UT/long
!                           14 - Terdiurnal
!                           15-25 - Not used in NRLMSIS 2.07
!                         For all switches:
!                           0.0 = Off
!                           1.0 = On
!                           2.0 = Main effects off, cross terms on
!                         Default values: 1.0
!       lzalt_type      Logical flag for altitude input type:
!                         True = Geodetic altitude (km)
!                         False = Geopotential height (km)
!                         Default: True (Geodetic altitude)
!       lspec_select    Logical array (1:10) flagging which densities to 
!                         calculate.
!                         True = Calculate, False = Do not calculate
!                            1 - Mass density
!                            2 - N2
!                            3 - O2
!                            4 - O
!                            5 - He
!                            6 - H
!                            7 - Ar
!                            8 - N
!                            9 - Anomalous O
!                           10 - NO
!                         Default values: True
!       lmass_include   Logical array (1:10) flagging which species to include
!                         in mass density calculation. Same ordering as 
!                         lspec_select.
!                         Default values: True
!       lN2_msis00      Logical flag for retrieving NRLMSISE-00 upper
!                         thermospheric N2 variation. See paper for details.
!                           False: Thermospheric N2 determined entirely by
!                             temperature profile and the constant mixing ratio
!                             of N2 in the lower atmosphere. 
!                           True: Upper thermospheric N2 relaxes to NRLMSISE-00
!                             Values.
!                         Default: False
!
!     NOTES:
!       - The switch_legacy optional argument performs the same function as
!         TSELEC(SW) in NRLSMSISE-00, except that switches 15-25 are not used in
!         NRLMSIS 2.07. The change in the switch-setting call is illustrated as
!         follows, where SW is the 25-element array of switches:
!           NRLMSISE-00: CALL TSELEC(SW)
!           NRLMSIS 2.07: call msisinit(switch_legacy=SW)
!
!!! ===========================================================================

!**************************************************************************************************
! MSIS_INIT Module: Contains initialization subroutines, model options, and model parameters
!**************************************************************************************************
module msis_init

  use msis_constants, only    : rp, nspec, nl, maxnbf, mbf

  implicit none
  
  !Model flags
  logical       :: initflag = .false.           !Flags whether model has been initialized
  logical       :: haveparmspace = .false.      !Flags whether parameter space has been initialized and allocated
  logical       :: zaltflag = .true.            !true: height input is geometric, false: height input is geopotential
  logical       :: specflag(1:nspec-1) = .true. !Array flagging which species densities are required
  logical       :: massflag(1:nspec-1) = .true. !Array flagging which species should be included in mass density
  logical       :: N2Rflag = .false.            !Flag for retrieving NRLMSISE-00 thermospheric N2 variations
  logical       :: zsfx(0:mbf) = .false.        !Flags zonal mean terms to be modulated by F1 (MSISE-00 legacy multiplier)
  logical       :: tsfx(0:mbf) = .false.        !Flags tide terms to be modulated by F2 (MSISE-00 legacy multiplier)
  logical       :: psfx(0:mbf) = .false.        !Flags SPW terms to be modulated by F3 (MSISE-00 legacy multiplier)
  logical       :: smod(0:nl) = .false.         !Flags which temperature levels get solar flux modulation; loadparmset turns flags on based on parameter values
  logical       :: swg(0:maxnbf-1) = .true.     !Switch array for globe subroutine.
  real(kind=rp) :: masswgt(1:nspec-1)  = 0.0_rp !Weights for calculating mass density
  real(4)       :: swleg(1:25)=1.0, swc(1:25), sav(1:25) !Legacy switch arrays

  ! Model parameter arrays
  type basissubset
    sequence
    character(8)               :: name
    integer                    :: bl,nl
    real(kind=rp), allocatable :: beta(:,:)
    logical, allocatable       :: active(:,:)
    integer, allocatable       :: fitb(:,:)
  end type basissubset
  type (basissubset)     :: TN
  type (basissubset)     :: PR
  type (basissubset)     :: N2
  type (basissubset)     :: O2
  type (basissubset)     :: O1
  type (basissubset)     :: HE
  type (basissubset)     :: H1
  type (basissubset)     :: AR
  type (basissubset)     :: N1
  type (basissubset)     :: OA   !Anomalous O
  type (basissubset)     :: NO
  integer                :: nvertparm
  
  ! Reciprocal node difference arrays (constant values needed for B-spline calculations)
  real(kind=rp)          :: etaTN(0:30,2:6) = 0.0_rp
  real(kind=rp)          :: etaO1(0:30,2:6) = 0.0_rp
  real(kind=rp)          :: etaNO(0:30,2:6) = 0.0_rp

  ! C1 constraint terms for O and NO related to the tapered logistic correction
  real(kind=rp)          :: HRfactO1ref, dHRfactO1ref, HRfactNOref, dHRfactNOref

contains

  !==================================================================================================
  ! MSISINIT: Entry point for initializing model and loading parameters
  !==================================================================================================
  subroutine msisinit(parmpath,parmfile,iun,switch_gfn,switch_legacy, &
                      lzalt_type,lspec_select,lmass_include,lN2_msis00)

    use msis_constants, only : specmass, nspec, maxnbf 

    implicit none

    character(len=*), intent(in), optional    :: parmpath                 !Path to parameter file
    character(len=*), intent(in), optional    :: parmfile                 !Parameter file name
    integer, intent(in), optional             :: iun                      !File unit number for reading parameter file
    logical, intent(in), optional             :: switch_gfn(0:maxnbf-1)   !Switch array for globe subroutine.
    real(4), intent(in), optional             :: switch_legacy(1:25)      !Legacy switch array
    logical, intent(in), optional             :: lzalt_type               !true: height input is geometric, false: height input is geopotential
    logical, intent(in), optional             :: lspec_select(1:nspec-1)  !Array flagging which species densities are required
    logical, intent(in), optional             :: lmass_include(1:nspec-1) !Array flagging which species should be included in mass density
    logical, intent(in), optional             :: lN2_msis00               !Flag for retrieving NRLMSISE-00 thermospheric N2 variations

    character(len=128)                        :: parmpath1
    character(len=128)                        :: parmfile1
    integer                                   :: iun1

    ! Path to parameter file
    if (present(parmpath)) then
      parmpath1 = parmpath
    else
      parmpath1 = '../extlib/Atmosphere/NRLMSIS-2.1/FORTRAN/'
    endif

    ! Parameter file name
    if (present(parmfile)) then
      parmfile1 = parmfile
    else
      parmfile1 = 'msis21.parm'
    endif

    ! Initialize model parameter space
    if (.not. haveparmspace) call initparmspace()

    ! Load parameter set
    if (present(iun)) then
      iun1 = iun
    else
      iun1 = 67
    endif
    call loadparmset(trim(parmpath1)//trim(parmfile1),iun1)

    ! Set switches
    swg(:) = .true.
    swleg(:) = 1.0
    if (present(switch_gfn)) then
      swg = switch_gfn
    else
      if (present(switch_legacy)) then
        swleg = switch_legacy
        call tselec(swleg)
      endif
    endif

    ! Input altitude type flag
    if (present(lzalt_type)) then
      zaltflag = lzalt_type
    else
      zaltflag = .true.
    endif

    ! Species flags for number density and mass density
    if (present(lspec_select)) then
      specflag = lspec_select
    else
      specflag(:) = .true.
    endif
    if (specflag(1)) then
      if (present(lmass_include)) then
        massflag = lmass_include
      else
        massflag(:) = .true.
      endif
    else
      massflag(:) = .false.
    endif
    where(massflag) specflag = .true.
    masswgt(:) = 0.0_rp
    where(massflag) masswgt = 1.0_rp
    masswgt(1) = 0.0_rp
    masswgt = masswgt * specmass
    masswgt(10) = 0.0_rp

    ! Flag for retrieving NRLMSISE-00 thermospheric N2 variations
    if (present(lN2_msis00)) then
      N2Rflag = lN2_msis00
    else
      N2Rflag = .false.
    endif

    ! Set model initialization flag
    initflag = .true.

    return

  end subroutine msisinit

  !==================================================================================================
  ! INITPARMSPACE: Initialize and allocate the model parameter space
  !==================================================================================================
  subroutine initparmspace()

    use msis_constants, only : nl, nls, nodesTN, ndO1, nsplO1, nodesO1, nsplNO, ndNO, nodesNO, &
                               zetagamma, Hgamma, zetarefO1, zetarefNO, maxnbf, ctide, cspw

    implicit none

    integer             :: n, m, j, k
    real(kind=rp)       :: gammaterm0

    ! Vertical parameter counter (number of vertical parameters in the parmeter file)
    nvertparm = 0

    ! Model formulation parameter subsets
    call initsubset(TN,0,nl,        maxnbf,'TN')
    call initsubset(PR,0,nl,        maxnbf,'PR')
    call initsubset(N2,0,nls,       maxnbf,'N2')
    call initsubset(O2,0,nls,       maxnbf,'O2')
    call initsubset(O1,0,nls+nsplO1,maxnbf,'O1')
    call initsubset(HE,0,nls,       maxnbf,'HE')
    call initsubset(H1,0,nls,       maxnbf,'H1')
    call initsubset(AR,0,nls,       maxnbf,'AR')
    call initsubset(N1,0,nls,       maxnbf,'N1')
    call initsubset(OA,0,nls,       maxnbf,'OA')
    call initsubset(NO,0,nls+nsplNO,maxnbf,'NO')

    ! Add the surface pressure parameter to the vertical parameter counter
    nvertparm = nvertparm + 1

    ! Set solar flux modulation flags
    zsfx(:) = .false.
    tsfx(:) = .false.
    psfx(:) = .false.
    ! F1, solar flux modulation of the zonal mean asymmetric annual terms
    zsfx(9:10) = .true.    !Pl(1,0) annual terms
    zsfx(13:14) = .true.   !Pl(3,0) annual terms
    zsfx(17:18) = .true.   !Pl(5,0) annual terms
    ! F2, solar sflux modulation of the tides
    tsfx(ctide:cspw-1) = .true.
    ! F3, solar sflux modulation of stationary planetary wave 1
    psfx(cspw:cspw+59) = .true. 

    ! Calculate reciprocal node difference arrays
    do k = 2, 6
      do j = 0, nl
        etaTN(j,k) = 1.0_rp / (nodesTN(j+k-1) - nodesTN(j))
      enddo
    enddo
    do k = 2, 4
      do j = 0, ndO1-k+1
        etaO1(j,k) = 1.0_rp / (nodesO1(j+k-1) - nodesO1(j))
      enddo
      do j = 0, ndNO-k+1
        etaNO(j,k) = 1.0_rp / (nodesNO(j+k-1) - nodesNO(j))
      enddo
    enddo

    ! Calculate C1 constraint terms for O and NO related to the tapered logistic correction
    gammaterm0 = tanh((zetarefO1 - zetagamma)*Hgamma)
    HRfactO1ref = 0.5_rp * (1.0_rp + gammaterm0)
    dHRfactO1ref = (1.0_rp - (zetarefO1 - zetagamma)*(1.0_rp - gammaterm0)*Hgamma) / HRfactO1ref
    gammaterm0 = tanh((zetarefNO - zetagamma)*Hgamma)
    HRfactNOref = 0.5_rp * (1.0_rp + gammaterm0)
    dHRfactNOref = (1.0_rp - (zetarefNO - zetagamma)*(1.0_rp - gammaterm0)*Hgamma) / HRfactNOref

    ! Set parameter space initialization flag
    haveparmspace = .true.

    return

  contains

      !--------------------------------------------------------------------------------------------------
      ! INITSUBSET: Initialize and allocate a parameter subset
      !--------------------------------------------------------------------------------------------------
      subroutine initsubset(subset,bl,nl,maxnbf,name)

        implicit none

        type (basissubset), intent(inout) :: subset
        integer, intent(in)               :: bl
        integer, intent(in)               :: nl
        integer, intent(in)               :: maxnbf
        character(2), intent(in)          :: name

        integer                           :: iz

        ! Allocate and initialize subset structure
        subset%name = name
        subset%bl = bl
        subset%nl = nl
        allocate(subset%beta(0:maxnbf-1,bl:nl), &
                 subset%active(0:maxnbf-1,bl:nl), &
                 subset%fitb(0:maxnbf-1,bl:nl))
        subset%beta = 0.0_rp
        subset%active = .false.
        subset%fitb = 0
        
        ! Increment vertical parameter counter except for pressure
        if (name .ne. 'PR') nvertparm = nvertparm + nl - bl + 1

        return

      end subroutine initsubset

  end subroutine initparmspace

  !==================================================================================================
  ! LOADPARMSET: Read in a parameter file
  !==================================================================================================
  subroutine loadparmset(name,iun)

    use msis_constants, only      : maxnbf, csfxmod

    implicit none

    character(len=*), intent(in) :: name
    integer, intent(in)          :: iun

    integer                      :: i0, i1
    logical                      :: havefile
    real(8), allocatable         :: parmin(:,:)

    ! Check if file exists
    inquire(file=trim(name),exist=havefile)
    if (havefile) then
       open(unit=iun,file=trim(name),status='old',access='stream',convert='little_endian')
    else
       print *,"MSIS parameter set ",trim(name)," not found. Stopping."
       stop
    endif

    ! Read in parameter values into temporary double-precision array
    allocate(parmin(0:maxnbf-1,0:nvertparm-1))
    read(iun) parmin
    close(iun)

    ! Transfer parameters to structures
    i0 = 0
    i1 = TN%nl - TN%bl
    TN%beta = parmin(:,i0:i1)
    i0 = i1 + 1
    i1 = i0
    PR%beta(:,0) = parmin(:,i0)
    i0 = i1 + 1
    i1 = i0 + N2%nl - N2%bl
    N2%beta = parmin(:,i0:i1)
    i0 = i1 + 1
    i1 = i0 + O2%nl - O2%bl
    O2%beta = parmin(:,i0:i1)
    i0 = i1 + 1
    i1 = i0 + O1%nl - O1%bl
    O1%beta = parmin(:,i0:i1)
    i0 = i1 + 1
    i1 = i0 + HE%nl - HE%bl
    HE%beta = parmin(:,i0:i1)
    i0 = i1 + 1
    i1 = i0 + H1%nl - H1%bl
    H1%beta = parmin(:,i0:i1)
    i0 = i1 + 1
    i1 = i0 + AR%nl - AR%bl
    AR%beta = parmin(:,i0:i1)
    i0 = i1 + 1
    i1 = i0 + N1%nl - N1%bl
    N1%beta = parmin(:,i0:i1)
    i0 = i1 + 1
    i1 = i0 + OA%nl - OA%bl
    OA%beta = parmin(:,i0:i1)
    i0 = i1 + 1
    i1 = i0 + NO%nl - NO%bl
    NO%beta = parmin(:,i0:i1)

    !Set solar flux modulation flags; if on for a given vertical parameter, then sfluxmod is called by tfnparm
    smod(:) = .false.
    where((Tn%beta(csfxmod+0,:) .ne. 0) .or. &
          (Tn%beta(csfxmod+1,:) .ne. 0) .or. &
          (Tn%beta(csfxmod+2,:) .ne. 0)) smod = .true.

    ! Compute log pressure spline coefficients from temperature spline coeffcients
    call pressparm()

    return

  end subroutine loadparmset

  !==================================================================================================
  ! PRESSPARM: Compute log pressure spline coefficients from temperature spline coeffcients
  !==================================================================================================
  subroutine pressparm()

    use msis_constants, only    : Mbarg0divkB, izfmx, mbf, gwht

    implicit none

    integer                    :: j, b, iz
    real(kind=rp)              :: lnz

    !Integrate pressure on nodes up to the last fully mixed level
    do j = 0, mbf
        lnz = 0.0
        do b = 0, 3
            lnz = lnz + TN%beta(j,b)*gwht(b)*Mbarg0divkB
        enddo
        PR%beta(j,1) = -lnz
        do iz = 1, izfmx
            lnz = 0.0
            do b = 0, 3
                lnz = lnz + TN%beta(j,iz+b)*gwht(b)*Mbarg0divkB
            enddo
            PR%beta(j,iz+1) = PR%beta(j,iz) - lnz
        enddo
    enddo

    return

  end subroutine pressparm

  !==================================================================================================
  ! TSELEC: Legacy switches and mapping to new switches
  !==================================================================================================
  subroutine tselec(sv)
  
    use msis_constants, only  : nsfx, nsfxmod, nut, cspw, csfx, csfxmod, cmag, cut

    implicit none

    real(4), intent(in)  :: sv(1:25)

    integer              :: i
    
    !Set cross-terms flags
    do i = 1, 25
      sav(i) = sv(i)
      swleg(i) = amod(sv(i), 2.0)
      if(abs(sv(i)) .eq. 1.0 .or. abs(sv(i)) .eq. 2.0) then
        swc(i) = 1.0
      else
        swc(i) = 0.0
      endif
    enddo
    
    !Main effects
    swg(0)                           = .true.                !Global term must be on
    swg(csfx:csfx+nsfx-1)            = (swleg(1) .eq. 1.0)   !Solar flux
    swg(310)                         = (swleg(1) .eq. 1.0)   !Solar flux (truncated quadratic F10.7a function)
    swg(1:6)                         = (swleg(2) .eq. 1.0)   !Time independent
    swg(304:305)                     = (swleg(2) .eq. 1.0)   !Time independent (extra, F10.7a modulated terms)
    swg(311:312)                     = (swleg(2) .eq. 1.0)   !Time independent (extra, truncated quadratic F10.7a modulated terms)
    swg(313:314)                     = (swleg(2) .eq. 1.0)   !Time independent (extra, dF10.7 modulated terms)
    swg((/7,8,11,12,15,16,19,20/))   = (swleg(3) .eq. 1.0)   !Symmetric annual
    swg(306:307)                     = (swleg(3) .eq. 1.0)   !Global AO (extra, solar-flux modulated terms)
    swg((/21,22,25,26,29,30,33,34/)) = (swleg(4) .eq. 1.0)   !Symmetric semiannual
    swg(308:309)                     = (swleg(4) .eq. 1.0)   !Global SAO (extra, solar-flux modulated terms)
    swg((/9,10,13,14,17,18/))        = (swleg(5) .eq. 1.0)   !Asymmetric annual
    swg((/23,24,27,28,31,32/))       = (swleg(6) .eq. 1.0)   !Asymmetric semiannual
    swg(35:94)                       = (swleg(7) .eq. 1.0)   !Diurnal
    swg(300:303)                     = (swleg(7) .eq. 1.0)   !Solar zenith angle
    swg(95:144)                      = (swleg(8) .eq. 1.0)   !Semidiurnal
    swg(145:184)                     = (swleg(14) .eq. 1.0)  !Terdiurnal
    swg(cmag:cmag+1)                 = .false.               !Geomagnetic activity mode master switch
    if((swleg(9) .gt. 0) .or. (swleg(13) .eq. 1)) swg(cmag:cmag+1) = (/.true.,.true./)  !Daily mode master switch
    if(swleg(9) .lt. 0)                           swg(cmag:cmag+1) = (/.false.,.true./) !Storm-time mode master switch
    swg(cmag+2:cmag+12)              = (swleg(9) .eq. 1.0)   !Daily geomagnetic activity terms
    swg(cmag+28:cmag+40)             = (swleg(9) .eq. -1.0)  !Storm-time geomagnetic activity terms
    swg(cspw:csfx-1)                 = ((swleg(11) .eq. 1.0) .and. (swleg(10) .eq. 1.0))  !Longitudinal
    swg(cut:cut+nut-1)               = ((swleg(12) .eq. 1.0) .and. (swleg(10) .eq. 1.0))  !UT/Lon
    swg(cmag+13:cmag+25)             = ((swleg(13) .eq. 1.0) .and. (swleg(10) .eq. 1.0))  !Mixed UT/Lon/Geomag (Daily mode terms)
    swg(cmag+41:cmag+53)             = ((swleg(13) .eq. 1.0) .and. (swleg(10) .eq. 1.0))  !Mixed UT/Lon/Geomag (Storm-time mode terms)

    !Cross terms
    swg(csfxmod:csfxmod+nsfxmod-1)   = (swc(1) .eq. 1.0)    !Solar activity modulation
    if (swc(1) .eq. 0) then                                 !Solar activity modulation
      swg(302:303) = .false.                                   !Solar zenith angle
      swg(304:305) = .false.                                   !Time independent
      swg(306:307) = .false.                                   !Global AO
      swg(308:309) = .false.                                   !Global SAO
      swg(311:314) = .false.                                   !Time independent
      swg(447) = .false.                                       !UT/Lon
      swg(454) = .false.                                       !UT/Lon
    endif
    if (swc(2) .eq. 0) then                                 !Time independent (latitude terms) (in MSISE-00, SWC(2) is not used - latitude modulations are always on)
      swg(9:20) = .false.                                      !AO
      swg(23:34) = .false.                                     !SAO
      swg(35:184) = .false.                                    !All tides
      swg(185:294) = .false.                                   !All SPW
      swg(392:414) = .false.                                   !Daily geomagnetic activity
      swg(420:442) = .false.                                   !Storm-time geomagnetic activity
      swg(449:453) = .false.                                   !UT/Lon
    endif
    if (swc(3) .eq. 0) then                                 !Symmetric annual
      swg(201:204) = .false.                                   !SPW1 (2,1)
      swg(209:212) = .false.                                   !SPW1 (4,1)
      swg(217:220) = .false.                                   !SPW1 (6,1)
      swg(255:258) = .false.                                   !SPW2 (2,2)
      swg(263:266) = .false.                                   !SPW2 (4,2)
      swg(271:274) = .false.                                   !SPW2 (6,2)
      swg(306:307) = .false.                                   !Global AO solar flux modulation
    endif
    if (swc(4) .eq. 0) then                                 !Symmetric semiannual
      swg(225:228) = .false.                                   !SPW1 (2,1)
      swg(233:236) = .false.                                   !SPW1 (4,1)
      swg(241:244) = .false.                                   !SPW1 (6,1)
      swg(275:278) = .false.                                   !SPW2 (2,2)
      swg(283:286) = .false.                                   !SPW2 (4,2)
      swg(291:294) = .false.                                   !SPW2 (6,2)
      swg(308:309) = .false.                                   !Global SAO solar flux modulation
    endif
    if (swc(5) .eq. 0) then                                 !Asymmetric annual
      swg(47:50) = .false.                                     !Diurnal (1,1)
      swg(51:54) = .false.                                     !Diurnal (2,1) !In MSISE-00, swc(5) is applied to all annual modulated tides
      swg(55:58) = .false.                                     !Diurnal (3,1)
      swg(59:62) = .false.                                     !Diurnal (4,1)
      swg(63:66) = .false.                                     !Diurnal (5,1)
      swg(67:70) = .false.                                     !Diurnal (6,1)
      swg(105:108) = .false.                                   !Semidiurnal (2,2)
      swg(109:112) = .false.                                   !Semidiurnal (3,2)
      swg(113:116) = .false.                                   !Semidiurnal (4,2)
      swg(117:120) = .false.                                   !Semidiurnal (5,2)
      swg(121:124) = .false.                                   !Semidiurnal (6,2)
      swg(153:156) = .false.                                   !Terdiurnal (3,3)
      swg(157:160) = .false.                                   !Terdiurnal (4,3)
      swg(161:164) = .false.                                   !Terdiurnal (5,3)
      swg(165:168) = .false.                                   !Terdiurnal (6,3)
      swg(197:200) = .false.                                   !SPW1 (1,1)
      swg(205:208) = .false.                                   !SPW1 (3,1)
      swg(213:216) = .false.                                   !SPW1 (5,1)
      swg(259:262) = .false.                                   !SPW2 (3,2)
      swg(267:270) = .false.                                   !SPW2 (5,2)
      swg(394:397) = .false.                                   !Geomag (Daily mode terms)
      swg(407:410) = .false.                                   !Mixed UT/Lon/Geomag (Daily mode terms)
      swg(422:425) = .false.                                   !Geomag (Storm-time mode terms)
      swg(435:438) = .false.                                   !Mixed UT/Lon/Geomag (Storm-time mode terms)
      swg(446)     = .false.                                   !UT/Lon
    endif
    if (swc(6) .eq. 0) then                                 !Asymmetric semiannual
      swg(221:224) = .false.                                   !SPW1 (1,1)
      swg(229:232) = .false.                                   !SPW1 (3,1)
      swg(237:240) = .false.                                   !SPW1 (5,1)
      swg(279:282) = .false.                                   !SPW2 (3,2)
      swg(287:290) = .false.                                   !SPW2 (5,2)
    endif
    if (swc(7) .eq. 0) then                                 !Diurnal
      swg(398:401) = .false.                                   !Geomag (Daily mode terms)
      swg(426:429) = .false.                                   !Geomag (Storm-time mode terms)
    endif
    if (swc(11) .eq. 0) then                                !Longitude
      swg(402:410) = .false.                                   !Mixed UT/Lon/Geomag (Daily mode terms)
      swg(430:438) = .false.                                   !Mixed UT/Lon/Geomag (Storm-time mode terms)
      swg(452:453) = .false.                                   !UT/Lon
    endif
    if (swc(12) .eq. 0) then                                !UT/Lon
      swg(411:414) = .false.                                   !Mixed UT/Lon/Geomag (Daily mode terms)
      swg(439:440) = .false.                                   !Mixed UT/Lon/Geomag (Storm-time mode terms)
    endif
    
  end subroutine tselec

  !==================================================================================================
  ! TRETRV: Legacy routine for retrieving switch settings
  !==================================================================================================
  subroutine tretrv(svv)

    implicit none

    real(4), intent(out) :: svv(1:25)

    integer              :: i

    do i = 1, 25
      svv(i) = sav(i)
    enddo
  
  end subroutine tretrv

end module msis_init
!!! ===========================================================================
!!! NRLMSIS 2.1:
!!! Neutral atmosphere empirical model from the surface to lower exosphere
!!! ===========================================================================
!!!
!!! MSISCALC: Interface with re-ordered input arguments and output arrays.
!
!     PREREQUISITES:
!       Must first run MSISINIT to load parameters and set switches. The 
!       MSISCALC subroutine checks for initialization and does a default
!       initialization if necessary. This self-initialization will be removed
!       in future versions.
!
!     CALLING SEQUENCE:
!       CALL MSISCALC(DAY, UTSEC, Z, LAT, LON, SFLUXAVG, SFLUX, AP, TN, DN, [TEX])
!  
!     INPUT VARIABLES:
!       DAY       Day of year (1.0 to 365.0 or 366.0)
!       UTSEC     Universal time (seconds)
!       Z         Geodetic altitude (km) (default) or Geopotential height (km)
!       LAT       Geodetic latitude (deg)
!       LON       Geodetic longitude (deg)
!       SFLUXAVG  81 day average, centered on input time, of F10.7 solar
!                 activity index
!       SFLUX     Daily F10.7 for previous day
!       AP        Geomagnetic activity index array:
!                   (1) Daily Ap
!                   (2) 3 hr ap index for current time
!                   (3) 3 hr ap index for 3 hrs before current time
!                   (4) 3 hr ap index for 6 hrs before current time
!                   (5) 3 hr ap index for 9 hrs before current time
!                   (6) Average of eight 3 hr ap indices from 12 to 33 hrs
!                       prior to current time
!                   (7) Average of eight 3 hr ap indices from 36 to 57 hrs
!                       prior to current time
!                 AP(2:7) are only used when switch_legacy(9) = -1.0 in MSISINIT
!
!     NOTES ON INPUT VARIABLES: 
!       - The day-of-year dependence of the model only uses the DAY argument. If
!         a continuous day-of-year dependence is desired, this argument should
!         include the fractional day (e.g., DAY = <day of year> + UTSEC/86400.0
!       - If lzalt_type = .true. (default) in the MSISINIT call, then Z is
!         treated as geodetic altitude.
!         If lzalt_type = .false., then Z is treated as geopotential height.
!       - F107 and F107A values are the 10.7 cm radio flux at the Sun-Earth
!         distance, not the radio flux at 1 AU. 
!
!     OUTPUT VARIABLES:
!       TN     Temperature at altitude (K)
!       DN(1)  Total mass density (kg/m3)
!       DN(2)  N2 number density (m-3)
!       DN(3)  O2 number density (m-3)
!       DN(4)  O number density (m-3)
!       DN(5)  He number density (m-3)
!       DN(6)  H number density (m-3)
!       DN(7)  Ar number density (m-3)
!       DN(8)  N number density (m-3)
!       DN(9)  Anomalous oxygen number density (m-3)
!       DN(10) NO number density (m-3)
!       TEX    Exospheric temperature (K) (optional argument)
!
!     NOTES ON OUTPUT VARIABLES: 
!       - Missing density values are returned as 9.999e-38
!       - Species included in mass density calculation are set in MSISINIT
!
!!! =========================================================================

!**************************************************************************************************
! MSIS_GFN Module: Contains subroutines to calculate global (horizontal and time-dependent) model 
!                  basis functions
!**************************************************************************************************
module msis_gfn

  use msis_constants, only : rp, maxn
  use msis_init, only      : TN,PR,N2,O2,O1,HE,H1,AR,N1,OA,NO, swg

  implicit none
  
  real(kind=rp)                :: plg(0:maxn,0:maxn)
  real(kind=rp)                :: cdoy(2), sdoy(2)
  real(kind=rp)                :: clst(3), slst(3)
  real(kind=rp)                :: clon(2), slon(2)
  real(kind=rp)                :: sfluxavgref = 150.0 ! Reference F10.7 value (=150 in NRLMSISE-00)
  real(kind=rp)                :: sfluxavg_quad_cutoff = 150.0 ! Cutoff F10.7 for truncated quadratic F10.7a function
  real(kind=rp)                :: lastlat = -999.9
  real(kind=rp)                :: lastdoy = -999.9
  real(kind=rp)                :: lastlst = -999.9
  real(kind=rp)                :: lastlon = -999.9
 
contains

  !==================================================================================================
  ! GLOBE: Calculate horizontal and time-dependent basis functions
  !        (Same purpose as NRLMSISE-00 "GLOBE7" subroutine)
  !==================================================================================================
  subroutine globe(doy,utsec,lat,lon,sfluxavg,sflux,ap,bf)

    use msis_constants, only    : deg2rad, doy2rad, lst2rad, &
                                  maxnbf, mbf, maxn, amaxn, amaxs, tmaxl, tmaxn, tmaxs, pmaxm, pmaxn, pmaxs, &
                                  nsfx, nsfxmod, ctimeind, cintann, ctide, cspw, csfx, cextra, cnonlin, csfxmod, cmag, cut
    implicit none

    real(kind=rp), intent(in)  :: doy                       ! Day of year
    real(kind=rp), intent(in)  :: utsec                     ! Universal time in seconds
    real(kind=rp), intent(in)  :: lat                       ! Latitude
    real(kind=rp), intent(in)  :: lon                       ! Longitdue
    real(kind=rp), intent(in)  :: sfluxavg                  ! 81-day average F10.7
    real(kind=rp), intent(in)  :: sflux                     ! Daily F10.7
    real(kind=rp), intent(in)  :: ap(1:7)                   ! Ap geomagnetic activity index history array
    real(kind=rp), intent(out) :: bf(0:maxnbf-1)            ! Output array of basis function terms

    real(kind=rp)              :: lst
    real(kind=rp)              :: slat, clat, clat2, clat4, slat2
    real(kind=rp)              :: cosdoy, sindoy
    real(kind=rp)              :: coslon, sinlon
    real(kind=rp)              :: pl
    real(kind=rp)              :: coslst, sinlst
    real(kind=rp)              :: dfa, df
    real(kind=rp)              :: theta
    real(kind=rp)              :: sza
    integer                    :: n, m, l, s, c

    ! Associated Legendre polynomials
    if (lat .ne. lastlat) then
      clat = sin(lat*deg2rad)  ! clat <=> sin, Legendre polyomial defined in colat
      slat = cos(lat*deg2rad)  ! slat <=> cos, Legendre polyomial defined in colat
      clat2 = clat*clat
      clat4 = clat2*clat2
      slat2 = slat*slat

      plg(0,0) = 1.0_rp
      plg(1,0) = clat
      plg(2,0) = 0.5_rp * (3.0_rp * clat2 - 1.0_rp)
      plg(3,0) = 0.5_rp * (5.0_rp * clat * clat2 - 3.0_rp * clat)
      plg(4,0) = (35.0_rp * clat4 - 30.0_rp * clat2 + 3.0_rp)/8.0_rp
      plg(5,0) = (63.0_rp * clat2 * clat2 * clat - 70.0_rp * clat2 * clat + 15.0_rp * clat)/8.0_rp
      plg(6,0) = (11.0_rp * clat * plg(5, 0) - 5.0_rp * plg(4, 0))/6.0_rp

      plg(1,1) = slat
      plg(2,1) = 3.0_rp * clat * slat
      plg(3,1) = 1.5_rp * (5.0_rp * clat2 - 1.0_rp) * slat
      plg(4,1) = 2.5_rp * (7.0_rp * clat2 * clat - 3.0_rp * clat) * slat
      plg(5,1) = 1.875_rp * (21.0_rp * clat4 - 14.0_rp * clat2 + 1.0_rp) * slat
      plg(6,1) = (11.0_rp * clat * plg(5, 1) - 6.0_rp * plg(4, 1))/5.0_rp

      plg(2,2) = 3.0_rp * slat2
      plg(3,2) = 15.0_rp * slat2 * clat
      plg(4,2) = 7.5_rp * (7.0_rp * clat2 - 1.0_rp) * slat2
      plg(5,2) = 3.0_rp * clat * plg(4, 2) - 2.0_rp * plg(3, 2)
      plg(6,2) = (11.0_rp * clat * plg(5, 2) - 7.0_rp * plg(4, 2))/4.0_rp

      plg(3,3) = 15.0_rp * slat2 * slat
      plg(4,3) = 105.0_rp * slat2 * slat * clat
      plg(5,3) = (9.0_rp * clat * plg(4, 3) - 7.0_rp * plg(3, 3))/2.0_rp
      plg(6,3) = (11.0_rp * clat * plg(5, 3) - 8.0_rp * plg(4, 3))/3.0_rp

      lastlat = lat
    endif

    ! Fourier harmonics of day of year
    if (doy .ne. lastdoy) then
      cdoy(1) = cos(doy2rad*doy)
      sdoy(1) = sin(doy2rad*doy)
      cdoy(2) = cos(doy2rad*doy*2.0_rp)
      sdoy(2) = sin(doy2rad*doy*2.0_rp)
      lastdoy = doy
    endif

    ! Fourier harmonics of local time
    lst = mod(utsec/3600.0_rp + lon/15.0_rp + 24.0_rp, 24.0_rp)
    if (lst .ne. lastlst) then
      clst(1) = cos(lst2rad*lst)
      slst(1) = sin(lst2rad*lst)
      clst(2) = cos(lst2rad*lst*2.0_rp)
      slst(2) = sin(lst2rad*lst*2.0_rp)
      clst(3) = cos(lst2rad*lst*3.0_rp)
      slst(3) = sin(lst2rad*lst*3.0_rp)
      lastlst = lst
    endif

    ! Fourier harmonics of longitude
    if (lon .ne. lastlon) then
      clon(1) = cos(deg2rad*lon)
      slon(1) = sin(deg2rad*lon)
      clon(2) = cos(deg2rad*lon*2.0_rp)
      slon(2) = sin(deg2rad*lon*2.0_rp)
      lastlon = lon
    endif

    !---------------------------------------------
    ! Coupled Linear Terms
    !---------------------------------------------

    ! Reset basis functions
    bf(:) = 0.0_rp

    ! Time-independent (pure latitude dependence)
    c = ctimeind
    do n = 0, amaxn
      bf(c) = plg(n,0)
      c = c + 1
    enddo

    ! Intra-annual (annual and semiannual)
    if (c .ne. cintann) stop 'problem with basis definitions'
    do s = 1, amaxs
      cosdoy = cdoy(s)
      sindoy = sdoy(s)
      do n = 0, amaxn
        pl = plg(n,0)
        bf(c) = pl*cosdoy
        bf(c+1) = pl*sindoy
        c = c + 2
      enddo
    enddo

    ! Migrating Tides (local time dependence)
    if (c .ne. ctide) stop 'problem with basis definitions'
    do l = 1, tmaxl
      coslst = clst(l)
      sinlst = slst(l)
      do n = l, tmaxn
        pl = plg(n,l)
        bf(c) = pl*coslst
        bf(c+1) = pl*sinlst
        c = c + 2
      enddo
      ! Intra-annual modulation of tides
      do s = 1, tmaxs
        cosdoy = cdoy(s)
        sindoy = sdoy(s)
        do n = l, tmaxn
          pl = plg(n,l)
          bf(c) = pl*coslst*cosdoy
          bf(c+1) = pl*sinlst*cosdoy
          bf(c+2) = pl*coslst*sindoy
          bf(c+3) = pl*sinlst*sindoy
          c = c + 4
        enddo
      enddo
    enddo

    ! Stationary Planetary Waves (longitude dependence)
    if (c .ne. cspw) stop 'problem with basis definitions'
    do m = 1, pmaxm
      coslon = clon(m)
      sinlon = slon(m)
      do n = m, pmaxn
        pl = plg(n,m)
        bf(c) = pl*coslon
        bf(c+1) = pl*sinlon
        c = c + 2
      enddo
      ! Intra-annual modulation of SPWs
      do s = 1, pmaxs
        cosdoy = cdoy(s)
        sindoy = sdoy(s)
        do n = m, pmaxn
          pl = plg(n,m)
          bf(c) = pl*coslon*cosdoy
          bf(c+1) = pl*sinlon*cosdoy
          bf(c+2) = pl*coslon*sindoy
          bf(c+3) = pl*sinlon*sindoy
          c = c + 4
        enddo
      enddo
    enddo
    
    ! Linear solar flux terms
    if (c .ne. csfx) stop 'problem with basis definitions'
    dfa = sfluxavg - sfluxavgref
    df = sflux - sfluxavg
    bf(c) = dfa
    bf(c+1) = dfa*dfa
    bf(c+2) = df
    bf(c+3) = df*df
    bf(c+4) = df*dfa
    c = c + nsfx

    ! Additional linear terms
    if (c .ne. cextra) stop 'problem with basis definitions'
    sza = solzen(doy,lst,lat,lon)
    bf(c)    = -0.5_rp*tanh((sza-98.0_rp)/6.0_rp)  !Solar zenith angle logistic function for O, H (transition width 3 deg, transition sza for horizon at ~65 km altitude)
    bf(c+1)  = -0.5_rp*tanh((sza-101.5_rp)/20.0_rp) !Solar zenith angle logistic function for NO (transition width 10 deg, transition sza for horizon at ~130 km altitude)
    bf(c+2)  = dfa*bf(c)                        !Solar flux modulation of logistic sza term
    bf(c+3)  = dfa*bf(c+1)                      !Solar flux modulation of logistic sza term
    bf(c+4)  = dfa*plg(2,0)                     !Solar flux modulation of P(2,0) term
    bf(c+5)  = dfa*plg(4,0)                     !Solar flux modulation of P(4,0) term
    bf(c+6)  = dfa*plg(0,0)*cdoy(1)             !Solar flux modulation of global AO
    bf(c+7)  = dfa*plg(0,0)*sdoy(1)             !Solar flux modulation of global AO
    bf(c+8) = dfa*plg(0,0)*cdoy(2)              !Solar flux modulation of global SAO
    bf(c+9) = dfa*plg(0,0)*sdoy(2)              !Solar flux modulation of global SAO
    if (sfluxavg .le. sfluxavg_quad_cutoff) then !Quadratic F10.7a function with cutoff of quadratic term (for robust extrapolation)
      bf(c+10) = dfa*dfa
    else
      bf(c+10) = (sfluxavg_quad_cutoff-sfluxavgref)*(2.0_rp*dfa - (sfluxavg_quad_cutoff-sfluxavgref))
    endif
    bf(c+11)  = bf(c+10)*plg(2,0)               !P(2,0) modulation of truncated quadratic F10.7a term
    bf(c+12)  = bf(c+10)*plg(4,0)               !P(4,0) modulation of truncated quadratic F10.7a term
    bf(c+13)  = df*plg(2,0)                     !P(2,0) modulation of df --> (F10.7 - F10.7a)
    bf(c+14)  = df*plg(4,0)                     !P(4,0) modulation of df --> (F10.7 - F10.7a)
    
    !---------------------------------------------
    ! Nonlinear Terms
    !---------------------------------------------

    c = cnonlin

    ! Solar flux modulation terms
    if (c .ne. csfxmod) stop 'problem with basis definitions'
    bf(c) = dfa  
    bf(c+1) = dfa*dfa
    bf(c+2) = df 
    bf(c+3) = df*df
    bf(c+4) = df*dfa
    c = c + nsfxmod

    ! Terms needed for legacy geomagnetic activity dependence
    if (c .ne. cmag) stop 'problem with basis set'
    bf(c:c+6) = ap - 4.0
    bf(c+8) =   doy2rad*doy
    bf(c+9) =   lst2rad*lst
    bf(c+10) =  deg2rad*lon
    bf(c+11) =  lst2rad*utsec/3600.0
    bf(c+12) =  abs(lat)
    c = c + 13
    do m = 0,1
      do n = 0,amaxn
        bf(c) = plg(n,m)
        c = c + 1
      enddo
    enddo

    ! Terms needed for legacy UT dependence
    c = cut
    bf(c) =   lst2rad*utsec/3600.0
    bf(c+1) = doy2rad*doy
    bf(c+2) = dfa
    bf(c+3) = deg2rad*lon
    bf(c+4) = plg(1,0)
    bf(c+5) = plg(3,0)
    bf(c+6) = plg(5,0)
    bf(c+7) = plg(3,2)
    bf(c+8) = plg(5,2)

    !---------------------------------------------
    ! Apply Switches
    !---------------------------------------------
    where(.not. swg(0:mbf)) bf(0:mbf) = 0.0_rp
    
    return

  end subroutine globe

  !==================================================================================================
  ! SOLZEN: Calculate solar zenith angle (adapted from IRI subroutine)
  !==================================================================================================
  real(kind=rp) function solzen(ddd,lst,lat,lon)

    use msis_constants, only    : pi, deg2rad

    implicit none

    real(kind=rp), intent(in)  :: ddd
    real(kind=rp), intent(in)  :: lst
    real(kind=rp), intent(in)  :: lat
    real(kind=rp), intent(in)  :: lon

    real(kind=rp)              :: wlon,dec
    real(kind=rp)              :: teqnx,tf,teqt
    real(kind=rp)              :: rlat,phi,cosx
    real(kind=rp), parameter   :: humr = pi/12.0_rp
    real(kind=rp), parameter   :: dumr = pi/182.5_rp
    real(kind=rp), parameter   :: p(5) = (/0.017203534,0.034407068,0.051610602,0.068814136,0.103221204/)

    wlon = 360.0 - lon
    teqnx = ddd + (lst + wlon / 15.0_rp) / 24.0_rp + 0.9369_rp
    teqnx = ddd + 0.9369_rp

    ! Solar declination
    dec = 23.256_rp * sin(p(1) * (teqnx - 82.242_rp)) + 0.381_rp * sin(p(2)*(teqnx - 44.855_rp))  &
         + 0.167_rp * sin(p(3) * (teqnx - 23.355_rp)) - 0.013_rp * sin(p(4)*(teqnx + 11.97_rp)) &
         + 0.011_rp * sin(p(5) * (teqnx - 10.410_rp)) + 0.339137_rp
    dec = dec * deg2rad

    ! Equation of time
    tf = teqnx - 0.5_rp
    teqt = -7.38_rp * sin(p(1) * (tf -  4.0_rp)) - 9.87_rp * sin(p(2) * (tf +  9.0_rp)) &
          + 0.27_rp * sin(p(3) * (tf - 53.0_rp)) -  0.2_rp * cos(p(4) * (tf - 17.0_rp))

    phi = humr * (lst - 12.0_rp) + teqt * deg2rad / 4.0_rp
    rlat = lat * deg2rad

    ! Cosine of solar zenith angle
    cosx = sin(rlat) * sin(dec) + cos(rlat) * cos(dec) * cos(phi)
    if (abs(cosx) .gt. 1.0_rp) cosx = sign(1.0_rp,cosx)

    solzen = acos(cosx) / deg2rad

    return

  end function solzen

  !==================================================================================================
  ! SFLUXMOD: Legacy nonlinear modulation of intra-annual, tide, and SPW terms
  !==================================================================================================
  real(kind=rp) function sfluxmod(iz,gf,parmset,dffact)

    use msis_constants, only       : maxnbf, mbf, csfx, csfxmod
    use msis_init, only            : basissubset, zsfx, tsfx, psfx

    implicit none

    integer, intent(in)           :: iz
    real(kind=rp), intent(in)     :: gf(0:maxnbf-1)
    type(basissubset), intent(in) :: parmset
    real(kind=rp), intent(in)     :: dffact  !Turns on or adjusts the delta-F terms added to F1 and F2 (eqns. A22b and A22c in Hedin (1987)).

    real(kind=rp)                 :: f1, f2, f3, sum
    integer                       :: j

    ! Intra-annual modulation factor
    if (swg(csfxmod)) then
      f1 = parmset%beta(csfxmod,iz) * gf(csfxmod) &
           + (parmset%beta(csfx+2,iz) * gf(csfxmod+2) + parmset%beta(csfx+3,iz) * gf(csfxmod+3) ) * dffact
    else
      f1 = 0.0_rp
    endif

    ! Migrating tide (local time) modulation factor
    if (swg(csfxmod+1)) then
      f2 = parmset%beta(csfxmod+1,iz) * gf(csfxmod) &
           + (parmset%beta(csfx+2,iz) * gf(csfxmod+2) + parmset%beta(csfx+3,iz) * gf(csfxmod+3) ) * dffact
    else
      f2 = 0.0_rp
    endif

    ! SPW (longitude) modulation factor
    if (swg(csfxmod+2)) then
      f3 = parmset%beta(csfxmod+2,iz) * gf(csfxmod)
    else
      f3 = 0.0_rp
    endif

    sum = 0.0
    do j = 0, mbf
      ! Apply intra-annual modulation
      if (zsfx(j)) then
        sum = sum + parmset%beta(j,iz)*gf(j)*f1
        cycle
      endif
      ! Apply migrating tide modulation
      if (tsfx(j)) then
        sum = sum + parmset%beta(j,iz)*gf(j)*f2
        cycle
      endif
      ! Apply SPW modulation
      if (psfx(j)) then
        sum = sum + parmset%beta(j,iz)*gf(j)*f3
        cycle
      endif
    enddo

    sfluxmod = sum

    return

  end function sfluxmod

  !==================================================================================================
  ! GEOMAG: Legacy nonlinear ap dependence (daily ap mode and ap history mode), including mixed 
  !         ap/UT/Longitude terms.
  ! Master switch control is as follows:
  !   swg(cmag) .nor. swg(cmag+1)   Do nothing: Return zero
  !   swg(cmag) .and. swg(cmag+1)   Daily Ap mode
  !   swg(cmag) .neqv. swg(cmag+1)  3-hour ap history mode
  !==================================================================================================
  real(kind=rp) function geomag(p0,bf,plg)

    use msis_constants, only    : nmag, cmag

    implicit none

    real(kind=rp), intent(in)  :: p0(0:nmag-1)
    real(kind=rp), intent(in)  :: bf(0:12)
    real(kind=rp), intent(in)  :: plg(0:6,0:1)

    logical                    :: swg1(0:nmag-1) !Copy of switches
    real(kind=rp)              :: p(0:nmag-1)    !Copy of parameters used to apply switches
    real(kind=rp)              :: delA, gbeta, ex, sumex, G(1:6)
    integer(4)                 :: i

    ! Return zero if both master switches are off    
    if (.not. (swg(cmag) .or. swg(cmag+1))) then
      geomag = 0.0_rp
      return
    endif

    ! Copy parameters
    p = p0
    swg1 = swg(cmag:cmag+nmag-1)

    ! Calculate function
    if (swg1(0) .eqv. swg1(1)) then
      ! Daily Ap mode
      if (p(1) .eq. 0) then  !If k00s is zero, then cannot compute function
        geomag = 0.0_rp
        return
      endif
      where(.not. swg1(2:25)) p(2:25) = 0.0_rp !Apply switches
      p(8) = p0(8) !Need doy phase term
      delA = G0fn(bf(0),p(0),p(1))
      geomag = ( p(2)*plg(0,0) +  p(3)*plg(2,0) +  p(4)*plg(4,0)                     &  ! time independent
        + (p(5)*plg(1,0) + p(6)*plg(3,0) + p(7)*plg(5,0)) * cos(bf(8) - p(8))        &  ! doy modulation
        + (p(9)*plg(1,1) + p(10)*plg(3,1) + p(11)*plg(5,1)) * cos(bf(9) - p(12))     &  ! local time modulation
        + (1.0_rp + p(13)*plg(1,0)) *                                                &
          (p(14)*plg(2,1) + p(15)*plg(4,1) + p(16)*plg(6,1)) * cos(bf(10) - p(17))   &  ! longitude effect
        + (p(18)*plg(1,1) + p(19)*plg(3,1) + p(20)*plg(5,1)) * cos(bf(10) - p(21)) * &
          cos(bf(8) - p(8))                                                          &  ! longitude with doy modulaiton
        + (p(22)*plg(1,0) + p(23)*plg(3,0) + p(24)*plg(5,0)) * cos(bf(11) - p(25)) ) &  ! universal time
        *delA
    else
      ! 3-hour ap history mode
      if (p(28) .eq. 0) then  !If beta00 is zero, then cannot compute function
        geomag = 0.0
        return
      endif
      where(.not. swg1(30:)) p(30:) = 0.0  !Apply switches
      p(36) = p0(36) !Need doy phase term
      gbeta = p(28)/(1 + p(29)*(45.0_rp - bf(12)))
      ex = exp(-10800.0_rp*gbeta)
      sumex = 1 + (1 - ex**19.0_rp) * ex**(0.5_rp) / (1 - ex)
      do i = 1, 6
        G(i) = G0fn(bf(i),p(26),p(27))
      enddo
      delA = ( G(1)                                                                 &
                    + ( G(2)*ex + G(3)*ex*ex + G(4)*ex**3.0_rp                      &
                       +(G(5)*ex**4.0_rp + G(6)*ex**12.0_rp)*(1-ex**8.0_rp)/(1-ex) ) ) / sumex
      geomag = ( p(30)*plg(0,0) +  p(31)*plg(2,0) +  p(32)*plg(4,0)                  &  ! time independent
        + (p(33)*plg(1,0) + p(34)*plg(3,0) + p(35)*plg(5,0)) * cos(bf(8) - p(36))    &  ! doy modulation
        + (p(37)*plg(1,1) + p(38)*plg(3,1) + p(39)*plg(5,1)) * cos(bf(9) - p(40))    &  ! local time modulation
        + (1.0_rp + p(41)*plg(1,0)) *                                                &
          (p(42)*plg(2,1) + p(43)*plg(4,1) + p(44)*plg(6,1)) * cos(bf(10) - p(45))   &  ! longitude effect
        + (p(46)*plg(1,1) + p(47)*plg(3,1) + p(48)*plg(5,1)) * cos(bf(10) - p(49)) * &
          cos(bf(8) - p(36))                                                         &  ! longitude with doy modulaiton
        + (p(50)*plg(1,0) + p(51)*plg(3,0) + p(52)*plg(5,0)) * cos(bf(11) - p(53)) ) &  ! universal time
        *delA
    endif

    return

    contains

      real(kind=rp) function G0fn(a,k00r,k00s)
          real(kind=rp),intent(in)  :: a, k00r, k00s
          G0fn = a + (k00r - 1.0_rp) * (a + (exp(-a*k00s) - 1.0_rp)/k00s)
          return
      end function G0fn

  end function geomag

  !==================================================================================================
  ! UTDEP: Legacy nonlinear UT dependence
  !==================================================================================================
  real(kind=rp) function utdep(p0,bf)

    use msis_constants, only    : nut, cut

    implicit none

    real(kind=rp), intent(in)  :: p0(0:nut-1)
    real(kind=rp), intent(in)  :: bf(0:8)

    real(kind=rp)              :: p(0:nut-1)    !Copy of parameters used to apply switches
    logical                    :: swg1(0:nut-1) !Copy of switches

    !Copy parameters
    p = p0
    swg1 = swg(cut:cut+nut-1)
    where(.not. swg1(3:nut-1)) p(3:nut-1) = 0.0  !Apply switches

    ! Calculate function
    utdep = cos(bf(0)-p(0)) *                          &
            (1 + p(3)*bf(4)*cos(bf(1)-p(1)))  *        &
            (1 + p(4)*bf(2)) * (1 + p(5)*bf(4)) *      &
            (p(6)*bf(4) + p(7)*bf(5) + p(8)*bf(6)) +   &
            cos(bf(0)-p(2)+2*bf(3)) * (p(9)*bf(7) + p(10)*bf(8)) * (1 + p(11)*bf(2))

    return

  end function utdep

end module msis_gfn
!**************************************************************************************************
! MSIS_TFN Module: Contains vertical temperature profile parameters and subroutines, including 
!                  temperature integration terms.
!**************************************************************************************************
module msis_tfn  

  use msis_constants, only     : rp, nl
    
  type tnparm
    sequence
    real(kind=rp)             :: cf(0:nl)    ! Spline coefficients
    real(kind=rp)             :: tzetaF      ! Tn at zetaF
    real(kind=rp)             :: tzetaA      ! Tn at zetaA (reference altitude for O1, H1)
    real(kind=rp)             :: dlntdzA     ! log-temperature gradient at zetaA (km^-1)
    real(kind=rp)             :: lndtotF     ! ln total number density at zetaF (m^-3)
    real(kind=rp)             :: tex
    real(kind=rp)             :: tgb0
    real(kind=rp)             :: tb0
    real(kind=rp)             :: sigma
    real(kind=rp)             :: sigmasq
    real(kind=rp)             :: b           ! b = 1-tb0/tex
    real(kind=rp)             :: beta(0:nl)  ! 1st integration coefficients on k=5 splines 
    real(kind=rp)             :: gamma(0:nl) ! 2nd integration coefficients on k=6 splines 
    real(kind=rp)             :: cVs         ! 1st integration constant (spline portion)
    real(kind=rp)             :: cVb         ! 1st integration constant (Bates portion)
    real(kind=rp)             :: cWs         ! 2nd integration constant (spline portion)
    real(kind=rp)             :: cWb         ! 2nd integration constant (Bates portion)
    real(kind=rp)             :: VzetaF      ! 1st indefinite integral at zetaF
    real(kind=rp)             :: VzetaA      ! 1st indefinite integral at zetaA
    real(kind=rp)             :: WzetaA      ! 2nd indefinite integral at zetaA
    real(kind=rp)             :: Vzeta0      ! 1st indefinite integral at zeta=0 (needed for pressure calculation)
  end type tnparm

  contains
    
  !==================================================================================================
  ! TFNPARM: Compute the vertical temperature and species-independent profile parameters
  !==================================================================================================
  subroutine tfnparm(gf,tpro)

    use msis_constants, only   : kB, lnP0, Mbarg0divkB, zetaB, zetaA, izfx, izax, itex, itgb0, itb0, c2tn, &
                                 maxnbf, mbf, nmag, nut, cmag, cut, &
                                 wbeta, wgamma, S5zetaB, S6zetaB, wghtAxdz, S4zetaA, S5zetaA, S6zetaA, &
                                 S4zetaF, S5zetaF, S5zeta0
    use msis_init, only        : smod, TN, PR
    use msis_gfn, only         : sfluxmod, geomag, utdep
    use msis_utils, only       : dilog

    implicit none

    real(kind=rp), intent(in) :: gf(0:maxnbf-1) ! Array of horizontal and temporal basis function terms   
    type(tnparm), intent(out) :: tpro           ! Output structure containing temperature vertical profile parameters
    
    integer(4)                :: ix
    real(kind=rp)             :: bc(3)
        
    ! Unconstrained spline coefficients
    do ix = 0, itb0-1
      tpro%cf(ix) = dot_product(TN%beta(0:mbf,ix),gf(0:mbf))
    enddo
    do ix = 0, itb0-1
      if (smod(ix)) then
        tpro%cf(ix) = tpro%cf(ix) + sfluxmod(ix,gf,TN,1.0_rp/TN%beta(0,ix))    !sfluxmod adds F10.7 modulation of tides
      endif
    enddo

    ! Exospheric temperature
    tpro%tex = dot_product(TN%beta(0:mbf,itex),gf(0:mbf))
    tpro%tex = tpro%tex + sfluxmod(itex,gf,TN,1.0_rp/TN%beta(0,itex))         
    tpro%tex = tpro%tex + geomag(TN%beta(cmag:cmag+nmag-1,itex),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
    tpro%tex = tpro%tex + utdep(TN%beta(cut:cut+nut-1,itex),gf(cut:cut+8))

    ! Temperature gradient at zetaB (122.5 km)
    tpro%tgb0 = dot_product(TN%beta(0:mbf,itgb0),gf(0:mbf))
    if (smod(itgb0)) tpro%tgb0 = tpro%tgb0 + sfluxmod(itgb0,gf,TN,1.0_rp/TN%beta(0,itgb0))         
    tpro%tgb0 = tpro%tgb0 + geomag(TN%beta(cmag:cmag+nmag-1,itgb0),gf(cmag:cmag+12),gf(cmag+13:cmag+26))

    ! Temperature at zetaB (122.5 km)
    tpro%tb0 = dot_product(TN%beta(0:mbf,itb0),gf(0:mbf))
    if (smod(itb0)) tpro%tb0 = tpro%tb0 + sfluxmod(itb0,gf,TN,1.0_rp/TN%beta(0,itb0))         
    tpro%tb0 = tpro%tb0 + geomag(TN%beta(cmag:cmag+nmag-1,itb0),gf(cmag:cmag+12),gf(cmag+13:cmag+26))

    ! Shape factor
    tpro%sigma = tpro%tgb0/(tpro%tex-tpro%tb0)

    ! Constrain top three spline coefficients for C2 continuity
    bc(1) = 1.0_rp/tpro%tb0
    bc(2) = -tpro%tgb0/(tpro%tb0*tpro%tb0)
    bc(3) = -bc(2)*(tpro%sigma + 2.0_rp*tpro%tgb0/tpro%tb0)    
    tpro%cf(itb0:itex) = matmul(bc, c2tn)

    ! Reference temperature at zetaF (70 km)
    tpro%tzetaF = 1.0_rp / dot_product(tpro%cf(izFx:izFx+2),S4zetaF)

    ! Reference temperature and gradient at zetaA (85 km)
    tpro%tzetaA = 1.0_rp / dot_product(tpro%cf(izAx:izAx+2),S4zetaA)
    tpro%dlntdzA = -dot_product(tpro%cf(izAx:izAx+2),wghtAxdz) * tpro%tzetaA
        
    ! Calculate spline coefficients for first and second 1/T integrals 
    tpro%beta(0) = tpro%cf(0)*wbeta(0)
    do ix = 1, nl
      tpro%beta(ix) = tpro%beta(ix-1) + tpro%cf(ix)*wbeta(ix)
    enddo
    tpro%gamma(0) = tpro%beta(0)*wgamma(0)
    do ix = 1, nl
      tpro%gamma(ix) = tpro%gamma(ix-1) + tpro%beta(ix)*wgamma(ix)
    enddo
        
    ! Integration terms and constants
    tpro%b = 1 - tpro%tb0 / tpro%tex
    tpro%sigmasq = tpro%sigma * tpro%sigma
    tpro%cVS = -dot_product(tpro%beta(itb0-1:itb0+2),S5zetaB)
    tpro%cWS = -dot_product(tpro%gamma(itb0-2:itb0+2),S6zetaB)
    tpro%cVB = -log(1-tpro%b) / (tpro%sigma * tpro%tex)
    tpro%cWB = -dilog(tpro%b) / (tpro%sigmasq * tpro%tex)
    tpro%VzetaF = dot_product(tpro%beta(izfx-1:izfx+2),S5zetaF) + tpro%cVS
    tpro%VzetaA = dot_product(tpro%beta(izax-1:izax+2),S5zetaA) + tpro%cVS
    tpro%WzetaA = dot_product(tpro%gamma(izax-2:izax+2),S6zetaA) + tpro%cVS*(zetaA-zetaB) + tpro%cWS
    tpro%Vzeta0 = dot_product(tpro%beta(0:2),S5zeta0) + tpro%cVS

    ! Compute total number density at zetaF
    tpro%lndtotF = lnP0 - Mbarg0divkB*(tpro%VzetaF - tpro%Vzeta0) - log(kB*tpro%TzetaF)

    return
    
  end subroutine tfnparm
    
  !==================================================================================================
  ! TFNX: Compute the temperature at specified geopotential height
  !==================================================================================================
  real(kind=rp) function tfnx(z,iz,wght,tpro)

    use msis_constants, only   : zetaB

    implicit none

    real(kind=rp), intent(in) :: z            ! Geopotential height
    integer, intent(in)       :: iz           ! Bspline reference index
    real(kind=rp), intent(in) :: wght(-3:0)   ! Bspline weights
    type(tnparm), intent(in)  :: tpro         ! Structure containing temperature vertical profile parameters

    integer                   :: i, j
  
    if (z .lt. zetaB) then 
      ! Spline region
      i = max(iz-3,0)
      if (iz .lt. 3) then
        j = -iz
      else
        j = -3
      endif
      tfnx = 1.0_rp / dot_product(tpro%cf(i:iz),wght(j:0))
    else
      ! Bates profile region
      tfnx = tpro%tex - (tpro%tex - tpro%tb0)*exp(-tpro%sigma * (z - zetaB))
    endif

    return

  end function tfnx
    
end module msis_tfn
!**************************************************************************************************
! MSIS_DFN Module: Contains vertical species density profile parameters and subroutines
!**************************************************************************************************
module msis_dfn

  use msis_constants, only : rp, nl, nsplO1, nsplNO
  use msis_utils, only     : bspline, dilog

  type dnparm
    sequence
    real(kind=rp)         :: lnPhiF            ! (Except O, H) Natural log of mixing ratio at zetaF (70 km), before chemical and dynamical corrections are applied (ln m^-3) (global term only)
    real(kind=rp)         :: lndref            ! Natural log of number density at reference height
    real(kind=rp)         :: zetaM             ! "Turbopause Height": Height of midpoint of effective mass transition (km)
    real(kind=rp)         :: HML               ! Scale height of lower portion of effective mass profile (km)
    real(kind=rp)         :: HMU               ! Scale height of upper portion of effective mass profile (km)
    real(kind=rp)         :: C                 ! Chapman term coefficient
    real(kind=rp)         :: zetaC             ! Chapman term reference height (km)
    real(kind=rp)         :: HC                ! Chapman term scale height (km)
    real(kind=rp)         :: R                 ! Chemical/dynamical term coefficient
    real(kind=rp)         :: zetaR             ! Chemical/dynamical term reference height (km)
    real(kind=rp)         :: HR                ! Chemical/dynamical term scale height (km)
    real(kind=rp)         :: cf(0:nsplO1+1)    ! Merged spline coefficients (for chemistry-dominated region of O1, NO, and (eventually), H, N)
    real(kind=rp)         :: zref              ! Reference height for hydrostatic integral and ideal gas terms
    real(kind=rp)         :: Mi(0:4)           ! Effective mass at nodes of piecewise mass profile (derived from zetaM, HML, HMU)
    real(kind=rp)         :: zetaMi(0:4)       ! Height of nodes of piecewise mass profile (derived from zetaM, HML, HMU)
    real(kind=rp)         :: aMi(0:4) = 0.0_rp ! Slopes of piecewise mass profile segments (derived from zetaM, HML, HMU)
    real(kind=rp)         :: WMi(0:4) = 0.0_rp ! 2nd indefinite integral of 1/T at mass profile nodes
    real(kind=rp)         :: XMi(0:4) = 0.0_rp ! Cumulative adjustment to M/T integral due to changing effective mass
    real(kind=rp)         :: Izref             ! Indefinite hydrostatic integral at reference height
    real(kind=rp)         :: Tref              ! Temperature at reference height (for ideal gas law term)
    real(kind=rp)         :: zmin              ! Minimum height of profile (missing values below)
    real(kind=rp)         :: zhyd              ! Hydrostatic terms needed above this height
    integer(kind=rp)      :: ispec             ! Species index
  end type dnparm

  contains

  !==================================================================================================
  ! DFNPARM: Compute the species density profile parameters
  !==================================================================================================
  subroutine dfnparm(ispec,gf,tpro,dpro)

    use msis_constants, only   : tanh1, specmass, lnvmr, Mbar, g0divkB, &
                                 nd, zetaF, zetaB, zetaA, nodesTN, &
                                 nodesO1, zetarefO1, c1o1, c1o1adj, &
                                 nodesNO, zetarefNO, c1NO, c1NOadj, &
                                 zetarefOA, &
                                 maxnbf, mbf, nmag, nut, cmag, cut
    use msis_init, only        : etaTN, TN,PR,N2,O2,O1,HE,H1,AR,N1,OA,NO, N2Rflag, &
                                 HRfactO1ref, dHRfactO1ref, HRfactNOref, dHRfactNOref
    use msis_gfn, only         : sfluxmod, geomag, utdep
    use msis_tfn, only         : tnparm

    implicit none

    integer, intent(in)       :: ispec          ! Species index
    real(kind=rp), intent(in) :: gf(0:maxnbf-1) ! Array of horizontal and temporal basis function terms   
    type(tnparm), intent(in)  :: tpro           ! Structure containing temperature vertical profile parameters
    type(dnparm), intent(out) :: dpro           ! Output structure containing density vertical profile parameters
        
    integer                   :: izf, i, i1, iz
    real(kind=rp)             :: Cterm, Rterm0, Rterm
    real(kind=rp)             :: bc(2)
    real(kind=rp)             :: hbetaL,hbetaU
    real(kind=rp)             :: delM, delz
    real(kind=rp)             :: Wi           ! 2nd indefinite integral at a piecewise mass profile node
    real(kind=rp)             :: Si(-5:0,2:6) ! Array of b-spline values at a mass profile node
    real(kind=rp)             :: Mzref        ! Effective mass at reference altitude

    dpro%ispec = ispec

    select case(ispec)

    ! Molecular Nitrogen ----------------------
    case(2)
      ! Mixing ratio and reference number density
      dpro%lnPhiF = lnvmr(ispec)
      dpro%lndref = tpro%lndtotF + dpro%lnPhiF
      dpro%zref = zetaF
      dpro%zmin = -1.0_rp
      dpro%zhyd = zetaF
      ! Effective mass
      dpro%zetaM = dot_product(N2%beta(0:mbf,1),gf(0:mbf))
      dpro%HML   = N2%beta(0,2)
      dpro%HMU   = N2%beta(0,3)
      ! Photochemical correction
      dpro%R     = 0.0_rp
      if (N2Rflag) dpro%R = dot_product(N2%beta(0:mbf,7),gf(0:mbf))
      dpro%zetaR = N2%beta(0,8)
      dpro%HR    = N2%beta(0,9)

    ! Molecular Oxygen ------------------------
    case(3)
      ! Mixing ratio and reference number density
      dpro%lnPhiF = lnvmr(ispec)
      dpro%lndref = tpro%lndtotF + dpro%lnPhiF
      dpro%zref = zetaF
      dpro%zmin = -1.0_rp
      dpro%zhyd = zetaF
      ! Effective mass
      dpro%zetaM = O2%beta(0,1)
      dpro%HML   = O2%beta(0,2)
      dpro%HMU   = O2%beta(0,3)
      ! Photochemical correction
      dpro%R     = dot_product(O2%beta(0:mbf,7),gf(0:mbf))
      dpro%R     = dpro%R + geomag(O2%beta(cmag:cmag+nmag-1,7),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      dpro%zetaR = O2%beta(0,8)
      dpro%HR    = O2%beta(0,9)

    ! Atomic Oxygen --------------------------
    case(4)
      ! Reference number density
      dpro%lnPhiF = 0.0_rp
      dpro%lndref = dot_product(O1%beta(0:mbf,0),gf(0:mbf))
      dpro%zref = zetarefO1
      dpro%zmin = nodesO1(3)
      dpro%zhyd = zetarefO1
      ! Effective mass
      dpro%zetaM = O1%beta(0,1)
      dpro%HML   = O1%beta(0,2)
      dpro%HMU   = O1%beta(0,3)
      ! Chapman correction
      dpro%C     = dot_product(O1%beta(0:mbf,4),gf(0:mbf))
      dpro%zetaC = O1%beta(0,5)
      dpro%HC    = O1%beta(0,6)
      ! Dynamical correction
      dpro%R     = dot_product(O1%beta(0:mbf,7),gf(0:mbf))
      dpro%R     = dpro%R + sfluxmod(7,gf,O1,0.0_rp)         
      dpro%R     = dpro%R + geomag(O1%beta(cmag:cmag+nmag-1,7),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      dpro%R     = dpro%R + utdep(O1%beta(cut:cut+nut-1,7),gf(cut:cut+8))
      dpro%zetaR = O1%beta(0,8)
      dpro%HR    = O1%beta(0,9)
      ! Unconstrained splines
      do izf = 0, nsplO1-1
        dpro%cf(izf) = dot_product(O1%beta(0:mbf,izf+10),gf(0:mbf))
      enddo
      ! Constrained splines calculated after case statement

    ! Helium ----------------------
    case(5)
      ! Mixing ratio and reference number density
      dpro%lnPhiF = lnvmr(ispec)
      dpro%lndref = tpro%lndtotF + dpro%lnPhiF
      dpro%zref = zetaF
      dpro%zmin = -1.0_rp
      dpro%zhyd = zetaF
      ! Effective mass
      dpro%zetaM = HE%beta(0,1)
      dpro%HML   = HE%beta(0,2)
      dpro%HMU   = HE%beta(0,3)
      ! Dynamical correction
      dpro%R     = dot_product(HE%beta(0:mbf,7),gf(0:mbf))
      dpro%R     = dpro%R + sfluxmod(7,gf,HE,1.0_rp)         
      dpro%R     = dpro%R + geomag(HE%beta(cmag:cmag+nmag-1,7),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      dpro%R     = dpro%R + utdep(HE%beta(cut:cut+nut-1,7),gf(cut:cut+8))
      dpro%zetaR = HE%beta(0,8)
      dpro%HR    = HE%beta(0,9)

    ! Atomic Hydrogen ----------------------
    case(6)
      ! Reference number density
      dpro%lnPhiF = 0.0_rp
      dpro%lndref = dot_product(H1%beta(0:mbf,0),gf(0:mbf))
      dpro%zref = zetaA
      dpro%zmin = 75.0_rp
      dpro%zhyd = zetaF
      ! Effective mass
      dpro%zetaM = H1%beta(0,1)
      dpro%HML   = H1%beta(0,2)
      dpro%HMU   = H1%beta(0,3)
      ! Chapman correction
      dpro%C     = dot_product(H1%beta(0:mbf,4),gf(0:mbf))
      dpro%zetaC = dot_product(H1%beta(0:mbf,5),gf(0:mbf))
      dpro%HC    = H1%beta(0,6)
      ! Dynamical correction
      dpro%R     = dot_product(H1%beta(0:mbf,7),gf(0:mbf))
      dpro%R     = dpro%R + sfluxmod(7,gf,H1,0.0_rp)        
      dpro%R     = dpro%R + geomag(H1%beta(cmag:cmag+nmag-1,7),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      dpro%R     = dpro%R + utdep(H1%beta(cut:cut+nut-1,7),gf(cut:cut+8))
      dpro%zetaR = H1%beta(0,8)
      dpro%HR    = H1%beta(0,9)

    ! Argon ----------------------
    case(7)
      ! Mixing ratio and reference number density
      dpro%lnPhiF = lnvmr(ispec)
      dpro%lndref = tpro%lndtotF + dpro%lnPhiF
      dpro%zref = zetaF
      dpro%zmin = -1.0_rp
      dpro%zhyd = zetaF
      ! Effective mass
      dpro%zetaM = AR%beta(0,1)
      dpro%HML   = AR%beta(0,2)
      dpro%HMU   = AR%beta(0,3)
      ! Dynamical correction
      dpro%R     = dot_product(AR%beta(0:mbf,7),gf(0:mbf))
      dpro%R     = dpro%R + geomag(AR%beta(cmag:cmag+nmag-1,7),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      dpro%R     = dpro%R + utdep(AR%beta(cut:cut+nut-1,7),gf(cut:cut+8))
      dpro%zetaR = AR%beta(0,8)
      dpro%HR    = AR%beta(0,9)

    ! Atomic Nitrogen ----------------------
    case(8)
      ! Reference number density
      dpro%lnPhiF = 0.0_rp
      dpro%lndref = dot_product(N1%beta(0:mbf,0),gf(0:mbf))
      dpro%lndref = dpro%lndref + sfluxmod(0,gf,N1,0.0_rp)         
      dpro%lndref = dpro%lndref + geomag(N1%beta(cmag:cmag+nmag-1,0),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      dpro%lndref = dpro%lndref + utdep(N1%beta(cut:cut+nut-1,0),gf(cut:cut+8))
      dpro%zref = zetaB
      dpro%zmin = 90.0_rp
      dpro%zhyd = zetaF
      ! Effective mass
      dpro%zetaM = N1%beta(0,1)
      dpro%HML   = N1%beta(0,2)
      dpro%HMU   = N1%beta(0,3)
      ! Chapman correction
      dpro%C     = N1%beta(0,4)
      dpro%zetaC = N1%beta(0,5)
      dpro%HC    = N1%beta(0,6)
      ! Dynamical correction
      dpro%R     = dot_product(N1%beta(0:mbf,7),gf(0:mbf))
      dpro%zetaR = N1%beta(0,8)
      dpro%HR    = N1%beta(0,9)

    ! Anomalous Oxygen ----------------------
    case(9)
      dpro%lndref = dot_product(OA%beta(0:mbf,0),gf(0:mbf))
      dpro%lndref = dpro%lndref + geomag(OA%beta(cmag:cmag+nmag-1,0),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      dpro%zref = zetarefOA
      dpro%zmin = 120.0_rp
      dpro%zhyd = 0.0_rp
      dpro%C     = OA%beta(0,4)
      dpro%zetaC = OA%beta(0,5)
      dpro%HC    = OA%beta(0,6)
      return !No further parameters needed for legacy anomalous oxygen profile

    ! Nitic Oxide ----------------------
    ! Added geomag dependence 2/18/21
    case(10)
      ! Skip if parameters are not defined
      if (NO%beta(0,0) .eq. 0.0_rp) then
          dpro%lndref = 0.0_rp
          return
      endif
      ! Reference number density
      dpro%lnPhiF = 0.0_rp
      dpro%lndref = dot_product(NO%beta(0:mbf,0),gf(0:mbf))
      dpro%lndref = dpro%lndref + geomag(NO%beta(cmag:cmag+nmag-1,0),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      dpro%zref = zetarefNO
      !dpro%zmin = nodesNO(3)
      dpro%zmin = 72.5  !JTE 1/18/22 Cut off profile below 72.5 km, due to possible spline artefacts at edge of domain (70 km)
      dpro%zhyd = zetarefNO
      ! Effective mass
      dpro%zetaM = dot_product(NO%beta(0:mbf,1),gf(0:mbf))
      dpro%HML   = dot_product(NO%beta(0:mbf,2),gf(0:mbf))
      dpro%HMU   = dot_product(NO%beta(0:mbf,3),gf(0:mbf))
      ! Chapman correction
      dpro%C     = dot_product(NO%beta(0:mbf,4),gf(0:mbf))
      dpro%C     = dpro%C + geomag(NO%beta(cmag:cmag+nmag-1,4),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      dpro%zetaC = dot_product(NO%beta(0:mbf,5),gf(0:mbf))
      dpro%HC    = dot_product(NO%beta(0:mbf,6),gf(0:mbf))
      ! Dynamical correction
      dpro%R     = dot_product(NO%beta(0:mbf,7),gf(0:mbf))
      dpro%zetaR = dot_product(NO%beta(0:mbf,8),gf(0:mbf))
      dpro%HR    = dot_product(NO%beta(0:mbf,9),gf(0:mbf))
      ! Unconstrained splines
      do izf = 0,nsplNO-1
          dpro%cf(izf) = dot_product(NO%beta(0:mbf,izf+10),gf(0:mbf))
          dpro%cf(izf) = dpro%cf(izf) + geomag(NO%beta(cmag:cmag+nmag-1,izf+10),gf(cmag:cmag+12),gf(cmag+13:cmag+26))
      enddo
      ! Constrained splines calculated after case statement

! Failsafe -----   ---------------------------
    case default
      stop 'Species not yet implemented'

    endselect
        
    ! Compute piecewise mass profile values and integration terms
    dpro%zetaMi(0) = dpro%zetaM - 2.0_rp*dpro%HML
    dpro%zetaMi(1) = dpro%zetaM - dpro%HML
    dpro%zetaMi(2) = dpro%zetaM
    dpro%zetaMi(3) = dpro%zetaM + dpro%HMU
    dpro%zetaMi(4) = dpro%zetaM + 2.0_rp*dpro%HMU
    dpro%Mi(0) = Mbar
    dpro%Mi(4) = specmass(ispec)
    dpro%Mi(2) = (dpro%Mi(0) + dpro%Mi(4)) / 2.0_rp
    delM = tanh1 * (dpro%Mi(4) - dpro%Mi(0)) / 2.0_rp
    dpro%Mi(1) = dpro%Mi(2) - delM
    dpro%Mi(3) = dpro%Mi(2) + delM
    do i = 0, 3
      dpro%aMi(i) = (dpro%Mi(i+1) - dpro%Mi(i)) / (dpro%zetaMi(i+1) - dpro%zetaMi(i))
    enddo
    do i = 0, 4
      delz = dpro%zetaMi(i) - zetaB
      if (dpro%zetaMi(i) .lt. zetaB) then
        call bspline(dpro%zetaMi(i),nodesTN,nd+2,6,etaTN,Si,iz)
        dpro%WMi(i) = dot_product(tpro%gamma(iz-5:iz),Si(:,6)) + tpro%cVS*delz + tpro%cWS
      else
        dpro%WMi(i) = (0.5_rp*delz*delz + dilog(tpro%b*exp(-tpro%sigma*delz))/tpro%sigmasq)/tpro%tex &
                      + tpro%cVB*delz + tpro%cWB
      endif
    end do
    dpro%XMi(0) = -dpro%aMi(0) * dpro%WMi(0)
    do i = 1, 3
      dpro%XMi(i) = dpro%XMi(i-1) - dpro%WMi(i) * (dpro%aMi(i) - dpro%aMi(i-1))
    end do
    dpro%XMi(4) = dpro%XMi(3) + dpro%WMi(4) * dpro%aMi(3)

    ! Calculate hydrostatic integral at reference height, and copy temperature
    if (dpro%zref .eq. zetaF) then
      Mzref = Mbar
      dpro%Tref = tpro%TzetaF
      dpro%Izref = Mbar * tpro%VzetaF
    else if (dpro%zref .eq. zetaB) then
      Mzref = pwmp(dpro%zref,dpro%zetaMi,dpro%Mi,dpro%aMi)
      dpro%Tref = tpro%Tb0
      dpro%Izref = 0.0_rp
      if ((zetaB .gt. dpro%zetaMi(0)) .and. (zetaB .lt. dpro%zetaMi(4))) then
        i = 0
        do i1 = 1, 3
          if (zetaB .lt. dpro%zetaMi(i1)) then
            exit
          else
            i = i1
          endif
        enddo
        dpro%Izref = dpro%Izref -  dpro%XMi(i)
      else
        dpro%Izref = dpro%Izref - dpro%XMi(4)                
      endif
    else if (dpro%zref .eq. zetaA) then
      Mzref = pwmp(dpro%zref,dpro%zetaMi,dpro%Mi,dpro%aMi)
      dpro%Tref = tpro%TzetaA
      dpro%Izref = Mzref * tpro%VzetaA
      if ((zetaA .gt. dpro%zetaMi(0)) .and. (zetaA .lt. dpro%zetaMi(4))) then
        i = 0
        do i1 = 1, 3
          if (zetaA .lt. dpro%zetaMi(i1)) then
            exit
          else
            i = i1
          endif
        enddo
        dpro%Izref = dpro%Izref - (dpro%aMi(i)*tpro%WzetaA + dpro%XMi(i))
      else
        dpro%Izref = dpro%Izref - dpro%XMi(4)                
      endif
    else 
      stop 'Integrals at reference height not available'
    endif

    ! C1 constraint for O1 at 85 km
    if (ispec .eq. 4) then
      Cterm = dpro%C*exp(-(dpro%zref-dpro%zetaC)/dpro%HC)
      Rterm0 = tanh((dpro%zref-dpro%zetaR)/(HRfactO1ref*dpro%HR))
      Rterm = dpro%R*(1+Rterm0)
      bc(1) = dpro%lndref - Cterm + Rterm - dpro%cf(7)*c1o1adj(1)      !Reference density, Chapman term, logistic term, and subtraction of last unconstrained spline(7)
      bc(2) = -Mzref*g0divkB/tpro%tzetaA &   !Gradient of hydrostatic term
              -tpro%dlntdzA &  !Gradient of ideal gas law term
              +Cterm/dpro%HC & !Gradient of Chapman term
              +Rterm*(1-Rterm0)/dpro%HR*dHrfactO1ref  & !Gradient of tapered logistic term
              -dpro%cf(7)*c1o1adj(2)  !Subtraction of gradient of last unconstrained spline(7)
      ! Compute coefficients for constrained splines
      dpro%cf(8:9) = matmul(bc,c1o1)
    endif
        
    ! C1 constraint for NO at 122.5 km
    if (ispec .eq. 10) then
      Cterm = dpro%C*exp(-(dpro%zref - dpro%zetaC)/dpro%HC)
      Rterm0 = tanh((dpro%zref-dpro%zetaR)/(HRfactNOref*dpro%HR))
      Rterm = dpro%R*(1+Rterm0)
      bc(1) = dpro%lndref - Cterm + Rterm - dpro%cf(7)*c1noadj(1)      !Reference density, Chapman term, logistic term, and subtraction of last unconstrained spline(7)
      bc(2) = -Mzref*g0divkB/tpro%tb0 &   !Gradient of hydrostatic term
              -tpro%tgb0/tpro%tb0 &  !Gradient of ideal gas law term
              +Cterm/dpro%HC & !Gradient of Chapman term
              +Rterm*(1-Rterm0)/dpro%HR*dHrfactNOref  & !Gradient of tapered logistic term
              -dpro%cf(7)*c1noadj(2)  !Subtraction of gradient of last unconstrained spline(7)
      ! Compute coefficients for constrained splines
      dpro%cf(8:9) = matmul(bc,c1no)
    endif
        
  return

  end subroutine dfnparm

  !==================================================================================================
  ! DFNX: Compute a species density at specified geopotential height
  !==================================================================================================
  real(kind=rp) function dfnx(z,tnz,lndtotz,Vz,Wz,HRfact,tpro,dpro)

    use msis_constants, only   : dmissing, g0divkB, ndO1, nodesO1, ndNO, nodesNO, HOA
    use msis_init, only        : etaO1, etaNO
    use msis_tfn, only         : tnparm

    implicit none

    real(kind=rp), intent(in) :: z            ! Geopotential height
    real(kind=rp), intent(in) :: tnz, lndtotz ! Temperature, total number density at input z 
    real(kind=rp), intent(in) :: Vz, Wz       ! First and second indefinite integrals of 1/T at z
    real(kind=rp), intent(in) :: HRfact       ! Reduction factor for chemical/dynamical correction scale height below zetaF
    type(tnparm), intent(in)  :: tpro         ! Structure containing temperature vertical profile parameters
    type(dnparm), intent(in)  :: dpro         ! Structure containing density vertical profile parameters

    integer(4)                :: i, i1, iz
    real(kind=rp)             :: Mz
    real(kind=rp)             :: Sz(-5:0,2:6)
    real(kind=rp)             :: Ihyd         ! Hydrostatic definite integral
    real(kind=rp)             :: ccor         ! Chapman and logistical corrections

    ! Below minimum height of profile
    if (z .lt. dpro%zmin) then
      dfnx = dmissing
      return
    endif
        
    ! Anomalous Oxygen (legacy MSISE-00 formulation)
    if (dpro%ispec .eq. 9) then
      dfnx = dpro%lndref - (z - dpro%zref)/HOA - dpro%C*exp(-(z-dpro%zetaC)/dpro%HC)
      dfnx = exp(dfnx)
      return               !No further calculation needed for anomalous oxygen
    endif

    ! Nitric Oxide: Skip if parameters are not defined
    if (dpro%ispec .eq. 10) then
      if (dpro%lndref .eq. 0.0_rp) then
        dfnx = dmissing
        return
      endif
    endif
        
    ! Chapman and logistic corrections
    select case(dpro%ispec)
    case(2,3,5,7)             !For N2, O2, He, and Ar: logistic correction only
      ccor =   dpro%R*(1+tanh((z-dpro%zetaR)/(HRfact*dpro%HR)))
    case(4,6,8,10)            !For O, H, N, and NO: Chapman and logistic corrections
      ccor = - dpro%C*exp(-(z-dpro%zetaC)/dpro%HC) &
             + dpro%R*(1+tanh((z-dpro%zetaR)/(HRfact*dpro%HR)))
    endselect

    ! Below height where hydrostatic terms are needed
    if (z .lt. dpro%zhyd) then
      select case(dpro%ispec)
      case(2,3,5,7)           !For N2, O2, He, and Ar, apply mixing ratios and exit
        dfnx = exp(lndtotz + dpro%lnPhiF + ccor)
        return
      case(4)                 !For O, evaluate splines
        call bspline(z,nodesO1,ndO1,4,etaO1,Sz,iz)
        dfnx = exp(dot_product(dpro%cf(iz-3:iz),Sz(-3:0,4)))
        return
      case(10)                !For NO, evaluate splines
        call bspline(z,nodesNO,ndNO,4,etaNO,Sz,iz)
        dfnx = exp(dot_product(dpro%cf(iz-3:iz),Sz(-3:0,4)))
        return
      endselect
    endif
        
    ! Calculate hydrostatic term and apply to reference density
    Mz = pwmp(z,dpro%zetaMi,dpro%Mi,dpro%aMi)
    Ihyd = Mz * Vz - dpro%Izref
    if ((z .gt. dpro%zetaMi(0)) .and. (z .lt. dpro%zetaMi(4))) then
      i = 0
      do i1 = 1, 3
        if (z .lt. dpro%zetaMi(i1)) then
          exit
        else
          i = i1
        endif
      enddo
      Ihyd = Ihyd - (dpro%aMi(i)*Wz + dpro%XMi(i))
    else if (z .ge. dpro%zetaMi(4)) then
      Ihyd = Ihyd - dpro%XMi(4)                
    endif
    dfnx = dpro%lndref - Ihyd * g0divkB  + ccor
        
    ! Apply ideal gas law                
    dfnx = exp(dfnx) * dpro%Tref/tnz

    return

  end function dfnx

  !==================================================================================================
  ! PWMP: Piecewise effective mass profile interpolation
  !==================================================================================================
  real(kind=rp) function pwmp(z,zm,m,dmdz)

    use msis_constants, only   : rp

    real(kind=rp), intent(in) :: z
    real(kind=rp), intent(in) :: zm(0:4)
    real(kind=rp), intent(in) :: m(0:4)
    real(kind=rp), intent(in) :: dmdz(0:3)

    integer                   :: irng  !Index of piecwise interval
    integer                   :: inode

    ! Most probable case
    if (z .ge. zm(4)) then
      pwmp = m(4)
      return
    endif

    ! Second most probable case
    if (z .le. zm(0)) then
      pwmp = m(0)
      return
    endif

    ! None of the above
    do inode = 0,3
      if (z .lt. zm(inode+1)) then
        pwmp = m(inode) + dmdz(inode)*(z - zm(inode))
        return
      endif
    enddo

    ! If we are here this is a problem
    stop 'Error in pwmp'

  end function pwmp

end module msis_dfn
!**************************************************************************************************
! MSIS_CALC Module: Contains main MSIS entry point
!**************************************************************************************************
module msis_calc

contains

  !==================================================================================================
  ! MSISCALC: The main MSIS subroutine entry point
  !==================================================================================================
  subroutine msiscalc(day,utsec,z,lat,lon,sfluxavg,sflux,ap,tn,dn,tex)

    use msis_constants, only    : rp, dmissing, lnp0, Mbarg0divkB, kB, nspec, nodesTN, nd, zetaF, zetaB, &
                                  Hgamma, zetagamma, maxnbf
    use msis_init, only         : msisinit, initflag, zaltflag, specflag, massflag, masswgt, etaTN
    use msis_gfn, only          : globe
    use msis_tfn, only          : tnparm, tfnparm, tfnx
    use msis_dfn, only          : dnparm, dfnparm, dfnx
    use msis_utils, only        : alt2gph, bspline, dilog

    implicit none

    real(kind=rp), intent(in)  :: day
    real(kind=rp), intent(in)  :: utsec
    real(kind=rp), intent(in)  :: z
    real(kind=rp), intent(in)  :: lat
    real(kind=rp), intent(in)  :: lon
    real(kind=rp), intent(in)  :: sfluxavg,sflux,ap(1:7)
    real(kind=rp), intent(out) :: tn, dn(1:10)
    real(kind=rp), intent(out), optional :: tex
  
    real(kind=rp), save        :: lastday = -9999.0
    real(kind=rp), save        :: lastutsec = -9999.0
    real(kind=rp), save        :: lastlat = -9999.0
    real(kind=rp), save        :: lastlon = -9999.0
    real(kind=rp), save        :: lastz = -9999.0
    real(kind=rp), save        :: lastsflux = -9999.0
    real(kind=rp), save        :: lastsfluxavg = -9999.0
    real(kind=rp), save        :: lastap(1:7) = -9999.0
    real(kind=rp), save        :: gf(0:maxnbf-1)
    real(kind=rp), save        :: Sz(-5:0,2:6)
    integer, save              :: iz
    type(tnparm), save         :: tpro
    type(dnparm), save         :: dpro(1:nspec-1)

    real(8)                    :: zaltd, latd
    real(kind=rp)              :: zeta, lndtotz, Vz, Wz, HRfact, lnPz, delz
    integer                    :: i, j, kmax, ispec

    ! Check if model has been initialized; if not, perform default initialization
    if (.not. initflag) call msisinit()

    ! Calculate geopotential height, if necessary
    if(zaltflag) then
      zaltd = dble(z)
      latd = dble(lat)
      zeta = alt2gph(latd,zaltd)
    else
      zeta = z
    endif

    ! If only altitude changes then update the local spline weights
    if (zeta .lt. zetaB) then
      if (zeta .ne. lastz) then
        if (zeta .lt. zetaF) then
          kmax = 5
        else
          kmax = 6
        endif
        call bspline(zeta,nodesTN,nd+2,kmax,etaTN,Sz,iz)
        lastz = zeta
      endif
    endif

    ! If location, time, or solar/geomagnetic conditions change then recompute the profile parameters
    if ((day .ne. lastday)     .or. (utsec .ne. lastutsec)       .or. &
        (lat .ne. lastlat)     .or. (lon .ne. lastlon)           .or. &
        (sflux .ne. lastsflux) .or. (sfluxavg .ne. lastsfluxavg) .or. &
        any(ap .ne. lastap)) then
      call globe(day,utsec,lat,lon,sfluxavg,sflux,ap,gf)
      call tfnparm(gf,tpro)
      do ispec = 2, nspec-1
        if (specflag(ispec)) call dfnparm(ispec,gf,tpro,dpro(ispec))
      enddo
      lastday = day
      lastutsec = utsec
      lastlat = lat
      lastlon = lon
      lastsflux = sflux
      lastsfluxavg = sfluxavg
      lastap = ap
    endif

    ! Exospheric temperature
    if (present(tex)) then
      tex = tpro%tex
    endif

    ! Temperature at altitude
    tn = tfnx(zeta,iz,Sz(-3:0,4),tpro)

    ! Temperature integration terms at altitude, total number density
    delz = zeta - zetaB
    if (zeta .lt. zetaF) then
      i = max(iz-4,0)
      if (iz .lt. 4) then
        j = -iz
      else
        j = -4
      endif
      Vz = dot_product(tpro%beta(i:iz),Sz(j:0,5)) + tpro%cVS
      Wz = 0.0_rp
      lnPz = lnP0 - Mbarg0divkB*(Vz - tpro%Vzeta0)
      lndtotz = lnPz - log(kB*tn)
    else
      if (zeta .lt. zetaB) then
        Vz = dot_product(tpro%beta(iz-4:iz),Sz(-4:0,5)) + tpro%cVS
        Wz = dot_product(tpro%gamma(iz-5:iz),Sz(-5:0,6)) + tpro%cVS*delz + tpro%cWS
      else
        Vz = (delz + log(tn/tpro%tex)/tpro%sigma)/tpro%tex + tpro%cVB
        Wz = (0.5_rp*delz*delz + dilog(tpro%b*exp(-tpro%sigma*delz))/tpro%sigmasq)/tpro%tex &
              + tpro%cVB*delz + tpro%cWB
      endif
    endif
        
    ! Species number densities at altitude
    HRfact = 0.5_rp * (1.0_rp + tanh(Hgamma*(zeta - zetagamma)))  !Reduction factor for chemical/dynamical correction scale height below zetagamma
    do ispec = 2, nspec-1
      if (specflag(ispec)) then
        dn(ispec) = dfnx(zeta,tn,lndtotz,Vz,Wz,HRfact,tpro,dpro(ispec))
      else
        dn(ispec) = dmissing
      endif
    enddo

    ! Mass density
    if (specflag(1)) then
      dn(1) = dot_product(dn,masswgt)
    else
      dn(1) = dmissing
    endif

    return

  end subroutine msiscalc

end module msis_calc
!!! GTD8D: Legacy wrapper with input and output arguments used in NRLMSISE-00
!
!     PREREQUISITES:
!       Must first run MSISINIT to load parameters and set switches. The
!       MSISCALC subroutine (called by this wrapper) checks for initialization
!       and does a default initialization if necessary. This self-initialization
!       will be removed in future versions.
!
!     CALLING SEQUENCE:
!       CALL GTD8D(IYD, SEC, ALT, GLAT, GLONG, STL, F107A, F107, AP, MASS, D, T)
!  
!     INPUT VARIABLES:
!       IYD    Year and day as YYDDD (day of year from 1 to 365 (or 366))
!                (Year is ignored in current model)
!       SEC    Universal time (seconds)
!       ALT    Geodetic altitude (km)
!       GLAT   Geodetic latitude (deg)
!       GLONG  Geodetic longitude (deg)
!       STL    Local solar time (Ignored; calculated from SEC and GLONG)
!       F107A  81 day average, centered on input time, of F10.7 solar activity
!                index
!       F107   Daily F10.7 for previous day
!       AP     Geomagnetic activity index array:
!                (1) Daily Ap
!                (2) 3 hr ap index for current time
!                (3) 3 hr ap index for 3 hrs before current time
!                (4) 3 hr ap index for 6 hrs before current time
!                (5) 3 hr ap index for 9 hrs before current time
!                (6) Average of eight 3 hr ap indices from 12 to 33 hrs
!                    prior to current time
!                (7) Average of eight 3 hr ap indices from 36 to 57 hrs
!                    prior to current time
!              AP(2:7) are only used when switch_legacy(9) = -1.0 in MSISINIT
!       MASS   Mass number (Ignored in 2.0)
!
!     NOTES ON INPUT VARIABLES: 
!       - If lzalt_type = .false. in the MSISINIT call, then the ALT input
!         argument is treated as geopotential height.
!       - The STL input argument is ignored in NRLMSIS 2.0. Instead, local time
!         is computed from universal time and longitude.
!       - F107 and F107A values are the 10.7 cm radio flux at the Sun-Earth
!         distance, not the radio flux at 1 AU. 
!       - The MASS input argument is ignored in NRLMSIS 2.0; species to be 
!         calculated are set in MSISINIT.
!
!     OUTPUT VARIABLES:
!       D(1)  He number density (cm-3)
!       D(2)  O number density (cm-3)
!       D(3)  N2 number density (cm-3)
!       D(4)  O2 number density (cm-3)
!       D(5)  Ar number density (cm-3)
!       D(6)  Total mass density (g/cm3)
!       D(7)  H number density (cm-3)
!       D(8)  N number density (cm-3)
!       D(9)  Anomalous oxygen number density (cm-3)
!       D(10) NO number density (cm-3)
!       T(1)  Exospheric temperature (K)
!       T(2)  Temperature at altitude (K)
!
!     NOTES ON OUTPUT VARIABLES: 
!       - Missing density values are returned as 9.999e-38
!       - Species included in mass density calculation are set in MSISINIT
!
!!! =========================================================================

!==================================================================================================
! GTD8D: Legacy wrapper
!==================================================================================================
subroutine gtd8d(iyd,sec,alt,glat,glong,stl,f107a,f107,ap,mass,d,t)

  use msis_constants, only     : rp, dmissing
  use msis_init, only          : msisinit
  use msis_calc, only          : msiscalc

  implicit none

  ! MSIS Legacy subroutine arguments
  integer, intent(in)         :: iyd
  real(4), intent(in)         :: sec
  real(4), intent(in)         :: alt
  real(4), intent(in)         :: glat
  real(4), intent(in)         :: glong
  real(4), intent(in)         :: stl
  real(4), intent(in)         :: f107a
  real(4), intent(in)         :: f107
  real(4), intent(in)         :: ap(7)
  integer, intent(in)         :: mass
  real(4), intent(out)        :: d(10), t(2)

  ! MSIS 1.97 subroutine arguments
  real(kind=rp)               :: xday
  real(kind=rp)               :: xutsec
  real(kind=rp)               :: xalt
  real(kind=rp)               :: xlat
  real(kind=rp)               :: xlon
  real(kind=rp)               :: xsfluxavg, xsflux
  real(kind=rp)               :: xap(1:7)
  real(kind=rp)               :: xtn
  real(kind=rp)               :: xdn(1:10)
  real(kind=rp)               :: xtex

  ! Convert the legacy input arguments to the new interface values and precision
  xday = mod(iyd,1000)
  xutsec = sec
  xalt = alt
  xlat = glat
  xlon = glong
  xsfluxavg = f107a
  xsflux = f107
  xap = ap

  ! Call the new subroutine
  call msiscalc(xday,xutsec,xalt,xlat,xlon,xsfluxavg,xsflux,xap,xtn,xdn,tex=xtex)

  ! Convert the output arguments to the legacy format (mks to cgs, re-order species)
  t(1) = sngl(xtex)    ! Expospheric temperature
  t(2) = sngl(xtn)     ! Temperature at altitude
  where (xdn .ne. dmissing) xdn = xdn*1d-6
  if (xdn(1) .ne. dmissing) xdn(1) = xdn(1)*1e3_rp
  d(1) = sngl(xdn(5))  ! [He]
  d(2) = sngl(xdn(4))  ! [O]
  d(3) = sngl(xdn(2))  ! [N2]
  d(4) = sngl(xdn(3))  ! [O2]
  d(5) = sngl(xdn(7))  ! [Ar]
  d(6) = sngl(xdn(1))  ! Mass density
  d(7) = sngl(xdn(6))  ! [H]
  d(8) = sngl(xdn(8))  ! [N]
  d(9) = sngl(xdn(9))  ! [Anomalous O]
  d(10) = sngl(xdn(10))  ! [NO]

  return

end subroutine gtd8d












