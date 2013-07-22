    !********************************************************************************
    !
    !  MODULE:  trisoprop
    !
    !  PURPOSE: Contains functions to calculate the thermal conductiviy of the triso
    !           particles and fuel pebbles
    !
    !  FUNCTIONS:
    !  k_graphite
    !  k_SiC
    !  k_densePyC
    !  k_PyC
    !  k_UO2 --------------------------plus two more to fix/delete
    !  k_trisolayers
    !
    !******************************************************************************** 
    MODULE trisoprop
    
        IMPLICIT NONE

    CONTAINS
            
            !================================================================================
            !  FUNCTION: THERMAL CONDUCTIVITY OF GRAPHITE [W/m-K]
            !================================================================================
            ! INPUT     ::  Temperature [Kelvin] (Must be between 723 K and 1575 K)
            ! OUTPUT    ::  Thermal conductivity of graphite [W/m-K]
            ! REFERENCE ::  SINAP report?
            !================================================================================
            REAL(8) FUNCTION k_graphite(T)
                IMPLICIT NONE
                REAL(8), INTENT(IN)     :: T        ! Temperature [Kelvin]
                REAL(8)                 :: DOSIS
                REAL(8)                 :: A
                REAL(8)                 :: B
                
                DOSIS=3.3e13*1360.0*24.0*3600.0/1.0e21
                A=(-0.3906e-4*T+0.06829)/(DOSIS+1.931e-4*T+0.105)
                B=1.228e-4*T+0.042
                
                k_graphite=(A+B)*1.2768E2       ! Thermal conductivity of graphite [W/m-K]
                
            END FUNCTION k_graphite

            !================================================================================
            !  FUNCTION: THERMAL CONDUCTIVITY OF TRISO SiC [W/m-K]
            !================================================================================
            ! INPUT     ::  Temperature [Kelvin]
            ! OUTPUT    ::  Thermal conductivity of SiC in TRISO [W/m-K]
            ! REFERENCE ::  ? Improved Prediction of the Temperature Feedback in TRISO-Fueled 
            !               Reactos: Appendix D 53 (70), INL 2009
            !================================================================================
            REAL(8) FUNCTION k_SiC(T)
                IMPLICIT NONE
                REAL(8), INTENT(IN)     :: T        ! Temperature [Kelvin]
                REAL(8)                 :: k_unirr  ! Thermal conductivity of unirradiated SiC
                REAL(8)                 :: DNE      ! Neutron fluence in 10^25 cm^-2 DNE units
                
                k_unirr=17885.0/(T+273.15)+2.0               ! [W/m-K]
                DNE=3.3e13*1360.0*24*3600/1.0e21
                
                k_SiC=k_unirr*EXP(-0.1277*DNE)      ! Thermal conductivity of  [W/m-K]
                
            END FUNCTION k_SiC
            
            !================================================================================
            !  FUNCTION: THERMAL CONDUCTIVITY OF DENSE PyC [W/m-K]
            !================================================================================
            ! INPUT     ::  Temperature [Kelvin]
            ! OUTPUT    ::  Thermal conductivity of dense PyC [W/m-K]
            ! REFERENCE ::  Improved Prediction of the Temperature Feedback in TRISO-Fueled 
            !               Reactos: Appendix D 53 (70), INL 2009
            !================================================================================
            REAL FUNCTION k_densePyC(T)
                IMPLICIT NONE
                REAL(8), INTENT(IN)     :: T        ! Temperature [Kelvin]
                REAL(8)                 :: P        ! Porosity
                REAL(8)                 :: ALPHA    ! 1.5 for spheres, average of 2.2
                REAL(8)                 :: FM       ! Porosity correction factor
                REAL(8)                 :: k_unirr  ! Thermal conductivity of unirradiated PyC
                REAL(8)                 :: DNE      ! Neutron fluence in 10^25 cm^-2 DNE units
                
                P=(1.93-1.9)/1.93                   ! Density of 1.9g/cc
                ALPHA=1.5
                FM=(1.0-P)/(1.0+(ALPHA-1.0)*P)
                k_unirr=2.443*(T+273.15)**(-0.574)*FM        ! [W/cm-K]
                DNE=3.3e13*1360.0*24*3600/1.0e21    ! for DNE > 3.03 [10^25 cm^-2 DNE units]
                
                k_densePyC=k_unirr*(1.0-0.3662*(1.0-EXP(-1.1028*DNE))-0.03554*DNE)/FM ! [W/cm-K]
                k_densePyC=k_densePyC*1.0e2         ! Thermal conductivity of dense PyC [W/m-K]
                
            END FUNCTION k_densePyC
            
            !================================================================================
            !  FUNCTION: THERMAL CONDUCTIVITY OF PyC [W/m-K]
            !================================================================================
            ! INPUT     ::  Temperature [Kelvin]
            ! OUTPUT    ::  Thermal conductivity of PyC [W/m-K]
            ! REFERENCE ::  Improved Prediction of the Temperature Feedback in TRISO-Fueled 
            !               Reactos: Appendix D 53 (70), INL 2009
            !================================================================================
            REAL(8) FUNCTION k_PyC(T)
                IMPLICIT NONE
                REAL(8), INTENT(IN)     :: T        ! Temperature [Kelvin]
                REAL(8)                 :: P        ! Porosity
                REAL(8)                 :: ALPHA    ! 1.5 for spheres, average of 2.2
                REAL(8)                 :: FM       ! Porosity correction factor
                REAL(8)                 :: k_unirr  ! Thermal conductivity of unirradiated PyC
                REAL(8)                 :: DNE      ! Neutron fluence in 10^25 cm^-2 DNE units
                
                P=(1.93-1.05)/1.93                  ! Density of 1.05g/cc
                ALPHA=2.2
                FM=(1.0-P)/(1.0+(ALPHA-1.0)*P) 
                k_unirr=2.443*(T+273.15)**(-0.574)*FM        ! [W/cm-K]
                DNE=3.3e13*1360.0*24*3600/1.0e21    ! for DNE > 3.03 [10^25 cm^-2 DNE units]
                
                k_PyC=k_unirr*(1.0-0.3662*(1.0-EXP(-1.1028*DNE))-0.03554*DNE)  ! [W/cm-K]
                k_PyC=k_PyC*1.0e2                   ! Thermal conductivity of PyC [W/m-K]
                
            END FUNCTION k_PyC
            
            !================================================================================
            !  FUNCTION: THERMAL CONDUCTIVITY OF UO2 [W/m-K]
            !================================================================================
            ! INPUT     ::  Temperature [Kelvin]
            ! OUTPUT    ::  Thermal conductivity of UO2 [W/m-K]
            ! REFERENCE ::  
            !================================================================================
            !REAL FUNCTION k_UO2(T)
            !    IMPLICIT NONE
            !    REAL(8), INTENT(IN)     :: T        ! Temperature [Kelvin]
            !    
            !    k_UO2=2.5                           ! Thermal conductivity of UO2 [W/m-K]
            !    
            !END FUNCTION k_UO2
        
            !================================================================================
            !  FUNCTION: THERMAL CONDUCTIVITY OF UO2 [W/m-K]
            !================================================================================
            ! INPUT     ::  Temperature [Kelvin]
            ! OUTPUT    ::  Thermal conductivity of UO2 [W/m-K]
            ! REFERENCE ::  
            !================================================================================
            REAL(8) FUNCTION k_UO2(TT)
                IMPLICIT NONE
                REAL(8), INTENT(IN)     :: TT        ! Temperature [Celcius]
                REAL(8)                 :: T         ! Temperature [Kelvin]
                
                k_UO2=1.2768E2                         ! Thermal conductivity of  [W/m-K]
                
                    T=TT+273.15
                    IF(T-273.15D0) 10,20,30
                10	T=273.15D0
	                GOTO 20
                30	IF(T-3088.71D0) 20,20,40
                40	T=3088.71D0
                20	k_UO2= .1288788790D+02&
	                -.2113481610D-01*T+.1830378640D-04*T*T&
	                -.8458677310D-08*T**3+.2008269200D-11*T**4&
	                -.1877097800D-15*T**5
                RETURN

            END FUNCTION k_UO2

            !================================================================================
            !  FUNCTION: THERMAL CONDUCTIVITY OF  [W/m-K]
            !================================================================================
            ! INPUT     ::  Temperature [Kelvin]
            ! OUTPUT    ::  Thermal conductivity of  [W/m-K]
            ! REFERENCE ::  
            !================================================================================
            !REAL FUNCTION k_(T)
            !    IMPLICIT NONE
            !    REAL, INTENT(IN)    :: T           ! Temperature [Kelvin]
            !    
            !    k_=1.2768E2       ! Thermal conductivity of  [W/m-K]
            !    
            !    !IMPLICIT REAL(8) (A-H,O-Z)
            !    !T=TT
            !    !TC=38.24D0/(T+402.55D0)+4.788D-13*TT*TT*TT
            !    !TCNUO2=TC*100.0D0
            !    !RETURN
            !    
            !END FUNCTION k_
            
            !================================================================================
            !  FUNCTION: THERMAL CONDUCTIVITY OF CHOSEN TRISO LAYER [W/m-K]
            !================================================================================
            ! INPUT     ::  Temperature [Kelvin]
            ! OUTPUT    ::  Thermal conductivity of chosen TRISO layer [W/m-K]
            ! REFERENCE ::  
            !================================================================================
            REAL(8) FUNCTION k_TRISOlayer(n,T)
                IMPLICIT NONE
                INTEGER, INTENT(IN) :: n            ! TRISO layer selection
                REAL(8), INTENT(IN) :: T            ! Temperature [Celcius]
        
                IF(N==1)THEN
                k_TRISOlayer=k_UO2(T)
                ELSEIF(N==2)THEN
                k_TRISOlayer=k_PyC(T)
                ELSEIF(N==3)THEN
                k_TRISOlayer=k_densePyC(T)
                ELSEIF(N==4)THEN
                k_TRISOlayer=k_SiC(T)
                ELSEIF(N==5)THEN
                k_TRISOlayer=k_densePyC(T)
                ELSEIF(N==6)THEN
                k_TRISOlayer=k_graphite(T)
                ELSEIF(N==7)THEN
                k_TRISOlayer=k_graphite(T)
                ENDIF
            
            END FUNCTION k_TRISOlayer
            
    END MODULE trisoprop