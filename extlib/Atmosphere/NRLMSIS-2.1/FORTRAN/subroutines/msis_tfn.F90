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
