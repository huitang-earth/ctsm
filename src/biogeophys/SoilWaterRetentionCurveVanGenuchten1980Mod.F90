module SoilWaterRetentionCurveVanGenuchten1980Mod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Implementation of soil_water_retention_curve_type using the Clapp-Hornberg 1978
  ! parameterizations.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  public :: soil_water_retention_curve_vangenuchten_1980_type
  
  type, extends(soil_water_retention_curve_type) :: &
       soil_water_retention_curve_vangenuchten_1980_type
     private
   contains
     procedure :: soil_hk              ! compute hydraulic conductivity
     procedure :: soil_suction         ! compute soil suction potential
     procedure :: soil_suction_inverse ! compute relative saturation at which soil suction is equal to a target value
  end type soil_water_retention_curve_vangenuchten_1980_type

  interface soil_water_retention_curve_vangenuchten_1980_type
     ! initialize a new soil_water_retention_curve_vangenuchten_1980_type object
     module procedure constructor  
  end interface soil_water_retention_curve_vangenuchten_1980_type

contains

  !-----------------------------------------------------------------------
  type(soil_water_retention_curve_vangenuchten_1980_type) function constructor()
    !
    ! !DESCRIPTION:
    ! Creates an object of type soil_water_retention_curve_vangenuchten_1980_type.
    ! For now, this is simply a place-holder.
    !-----------------------------------------------------------------------

  end function constructor

  !-----------------------------------------------------------------------
  subroutine soil_hk(this, c, j, s, imped, soilstate_inst, hk, dhkds)
    !
    ! !DESCRIPTION:
    ! Compute hydraulic conductivity
    !
    ! !USES:
    use SoilStateType  , only : soilstate_type
    !
    ! !ARGUMENTS:
    class(soil_water_retention_curve_vangenuchten_1980_type), intent(in) :: this
    integer,  intent(in)             :: c        !column index
    integer,  intent(in)             :: j        !level index
    real(r8), intent(in)             :: s        !relative saturation, [0, 1]
    real(r8), intent(in)             :: imped    !ice impedance
    type(soilstate_type), intent(in) :: soilstate_inst
    real(r8), intent(out)            :: hk       !hydraulic conductivity [mm/s]
    real(r8), optional, intent(out)  :: dhkds    !d[hk]/ds   [mm/s]
    
                                                  ! l is the pore-connectivity parameter (−)
    
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'soil_hk'
    
    ! !LOCAL VARIABLES:
    real(8) :: n1, m1, alpha_van ! (-), pore-size-distribution parameter for Van Genuchten 1.07
    
    !-----------------------------------------------------------------------
    
    associate(& 
         hksat             =>    soilstate_inst%hksat_col(c,j)          , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
  !       bsw               =>    soilstate_inst%bsw_col(c,j)              & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
         )

    n1 = n_van 
    m1 = 1.0/n1   ! (1-1.0/n1) ?
    !alpha_van = alpha_van
    
    l=0.5 
    

    !watsat = spafhy_para%watsat  
    !watres = spafhy_para%watres     
    !vol_ice = 0.0  
    !eff_porosity = max(0.01, watsat - vol_ice)
    !satfrac = (vol_liq - watres) / (eff_porosity-watres)

    ! hydraulic conductivity (vanGenuchten - Mualem)
    hk = hksat * s**l  * (1.0 - (1.0 - s**(1/m1))**m1)**2.0    ! [m s-1]
  

    end associate 

  end subroutine soil_hk
  
  
  
  SUBROUTINE soil_water_retention_curve(vol_liq, spafhy_para, smp)
  ! Converts vol. water content to soil water potential (in MPa)
  ! Add restriction that smp can't drop too low?

    real(8), intent(in) :: vol_liq        ! v/v, volumetric of liq in soil bucket
    type(spafhy_para_type), intent(in)    :: spafhy_para ! parameters
    real(8), intent(out):: smp            ! soil suction, negative, MPa

  ! !LOCAL VARIABLES:

    real(8) :: vol_ice     ! v/v, volumetric ice in soil bucket 
    real(8) :: satfrac     ! parameter for Van Genuchten
    real(8) :: n1, m1, alpha_van, watsat, watres 
    real(8) :: eff_porosity! v/v, volume of ice
        
    n1 = spafhy_para%n_van 
    m1 = 1.0/n1  
    alpha_van = spafhy_para%alpha_van 
    watsat = spafhy_para%watsat  
    watres = spafhy_para%watres
    
    vol_ice = 0.0  
    eff_porosity = max(0.01, watsat - vol_ice)
    
    satfrac = (vol_liq-watres)/(eff_porosity-watres)
    smp = -(1.0/alpha_van)*(satfrac**(1.0/(m1-1.0)) - 1.0 )**m1 !kPa
    smp = smp * 0.001 !MPa

  END SUBROUTINE soil_water_retention_curve
  
  ! Convert matric potential to water head (m)
  ! hp = h = Ψp / ρwg
  
  subroutine soil_suction(this, c, j, s, soilstate_inst, smp, dsmpds)
    !j, 
    ! !DESCRIPTION:
    ! Compute soil suction potential
    !
    ! !USES:
    use SoilStateType  , only : soilstate_type
    !
    ! !ARGUMENTS:
    class(soil_water_retention_curve_vangenuchten_1980_type), intent(in) :: this
    integer,  intent(in)             :: c       !column index
    integer,  intent(in)             :: j        !level index
    real(r8), intent(in)             :: s        !relative saturation, [0, 1]
    type(soilstate_type), intent(in) :: soilstate_inst
    real(r8), intent(out)            :: smp      !soil suction, negative, [mm]
    real(r8), optional, intent(out)  :: dsmpds   !d[smp]/ds, [mm]
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'soil_suction'
    real(8) :: n1, m1, alpha_van ! (-), pore-size-distribution parameter for Van Genuchten 1.07
    
    
    !-----------------------------------------------------------------------
    
    associate(& 
   !      bsw               =>    soilstate_inst%bsw_col(c,j)            , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                       
         sucsat            =>    soilstate_inst%sucsat_col(c,j)           & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
         )

    n1 = n_van 
    m1 = 1.0/n1 

    !compute soil suction potential, negative
    smp = -(1.0/alpha_van)*(s**(1.0/(m1-1.0)) - 1.0 )**m1 !kPa
    smp = smp * 0.001 !MPa

    !compute derivative
    !if(present(dsmpds))then
    !   dsmpds=-bsw*smp/s
    !endif

    end associate 

  end subroutine soil_suction

  !-----------------------------------------------------------------------
  subroutine soil_suction_inverse(this, c, j, smp_target, soilstate_inst, s_target)
    !
    ! !DESCRIPTION:
    ! Compute relative saturation at which soil suction is equal to a target value.
    ! This is done by inverting the soil_suction equation to solve for s.
    !
    ! !USES:
    use SoilStateType  , only : soilstate_type
    !
    ! !ARGUMENTS:
    class(soil_water_retention_curve_vangenuchten_1980_type), intent(in) :: this
    integer,  intent(in)             :: c       !column index
    integer,  intent(in)             :: j        !level index
    type(soilstate_type), intent(in) :: soilstate_inst
    real(r8) , intent(in)  :: smp_target ! target soil suction, negative [mm]
    real(r8) , intent(out) :: s_target   ! relative saturation at which smp = smp_target [0,1]
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'soil_suction_inverse'
    !-----------------------------------------------------------------------
    
    associate(& 
         bsw               =>    soilstate_inst%bsw_col(c,j)            , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
         sucsat            =>    soilstate_inst%sucsat_col(c,j)           & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
         )

    s_target = (-smp_target/sucsat)**(-1._r8/bsw)

    end associate 

  end subroutine soil_suction_inverse

end module SoilWaterRetentionCurveVanGenuchten1980Mod


