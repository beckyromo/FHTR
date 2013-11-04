    !********************************************************************************
    !
    !  MODULE:  flibeprop
    !
    !  PURPOSE: Contains functions to calculate the thermophysical properties of flibe.
    !
    !  FUNCTIONS:
    !  flibe_cp             heat capacity of flibe [J/kg-C]
    !  flibe_k              thermal conductivity of flibe [W/m-C]
    !  flibe_rho            density of flibe [kg/m^3]
    !  flibe_mu             viscosity of flibe [Pa-s]
    !  flibe_enthaply       calculates enthalpy of flibe for given temperature [J/kg]
    !  flibe_temperature    calculates temperature of flibe for given enthalpy [Celcius]
    !
    !******************************************************************************** 
    MODULE flibeprop
            
        USE global, ONLY: T_REF
    
        IMPLICIT NONE


    CONTAINS
    
            !================================================================================
            !  FUNCTION: HEAT CAPACITY OF FLIBE [J/kg-K]
            !================================================================================
            ! INPUT     ::  Temperature [Kelvin]
            ! OUTPUT    ::  Heat Capacity of Flibe [J/kg-K]
            ! REFERENCE ::  Benes and Konings (2012)
            !================================================================================
            REAL(8) FUNCTION flibe_cp(T)
                IMPLICIT NONE
                REAL(8), INTENT(IN) :: T            ! Temperature [Kelvin]  
                flibe_cp=2380.0                     ! Heat Capacity of Flibe [J/kg-C]
            END FUNCTION flibe_cp
            
            !================================================================================
            !  FUNCTION: THERMAL CONDUCTIVITY OF FLIBE [W/m-K]
            !================================================================================
            ! INPUT     ::  Temperature [Kelvin]
            ! OUTPUT    ::  Thermal conductivity of flibe [W/m-K]
            ! REFERENCE ::  Williams et al. (2001), Cantor et al. (1968), Gierszewski et al. (1980)
            !================================================================================
            REAL(8) FUNCTION flibe_k(T)
                IMPLICIT NONE
                REAL(8), INTENT(IN) :: T            ! Temperature [Kelvin]             
                flibe_k=1.0                         ! Thermal conductivity of flibe [W/m-K]
            END FUNCTION flibe_k
            
            !================================================================================
            !  FUNCTION: DENSITY OF FLIBE [kg/m^3]
            !================================================================================
            ! INPUT     ::  Temperature [Kelvin]
            ! OUTPUT    ::  Density of flibe [kg/m^3] 
            ! REFERENCE ::  Compilation
            !================================================================================
            REAL(8) FUNCTION flibe_rho(T)
                IMPLICIT NONE
                REAL(8), INTENT(IN) :: T            ! Temperature [Celcius]
                !flibe_rho=2422.2-0.4859*T           ! Density of flibe [kg/m^3]
                flibe_rho=(2.28-4.9E-4*T)*1.0E3
            END FUNCTION flibe_rho
            
            !================================================================================
            !  FUNCTION: VISCOSITY OF FLIBE [Pa-s]
            !================================================================================
            ! INPUT     ::  Temperature [Kelvin]
            ! OUTPUT    ::  Viscosity of flibe [Pa-s]
            ! REFERENCE ::  Williams et al. (2006), Benes et al. (2012)
            !================================================================================
            REAL(8) FUNCTION flibe_mu(T)
                IMPLICIT NONE
                REAL(8), INTENT(IN) :: T                    ! Temperature [Celcius]              
                flibe_mu=0.000116*EXP(3755.0/(T+273.15))    ! Viscosity of flibe [Pa-s]
            END FUNCTION
            
            !================================================================================
            !  FUNCTION: CONVERT FLIBE TEMPERATURE TO ENTHALPY
            !================================================================================
            ! INPUT     ::  Temperature [Kelvin]
            ! OUTPUT    ::  Enthalpy of flibe []
            ! REFERENCE ::  
            !================================================================================
            REAL(8) FUNCTION flibe_enthalpy(T)
                IMPLICIT NONE
                REAL(8), INTENT(IN) :: T            ! Temperature [Celcius]
                flibe_enthalpy=flibe_cp(T)*(T-T_REF)
            END FUNCTION         
            
            !================================================================================
            !  FUNCTION: CONVERT FLIBE ENTHALPY TO TEMPERATURE
            !================================================================================
            ! INPUT     ::  Enthalpy []
            ! OUTPUT    ::  Temperature of flibe []
            ! REFERENCE ::  
            !================================================================================
            REAL(8) FUNCTION flibe_temperature(h)
                IMPLICIT NONE
                REAL(8), INTENT(IN) :: h            ! Enthalpy
                flibe_temperature=h/flibe_cp(T_REF)+T_REF
            END FUNCTION         
    
    
    END MODULE flibeprop