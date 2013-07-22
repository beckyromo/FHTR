    !********************************************************************************
    !
    !  MODULE:  pebblebed
    !
    !  PURPOSE: Contains functions for pebble bed core calculations
    !
    !  FUNCTIONS:
    !  porosity_PB
    !  HTC_PB
    !  FF_PB
    !
    !******************************************************************************** 
    MODULE pebblebed
            
        USE global
        USE flibeprop, ONLY: flibe_cp, flibe_k, flibe_rho, flibe_mu
        
        IMPLICIT NONE


    CONTAINS
    
            !================================================================================
            !  FUNCTION: POROSITY OF PEBBLE PED
            !================================================================================
            ! INPUT     ::  
            ! OUTPUT    ::  Porosity of the pebble bed
            ! REFERENCE ::  
            !================================================================================
            REAL(8) FUNCTION porosity_PB(R)
                IMPLICIT NONE
                REAL(8), INTENT(IN)     :: R     ! RADIUS OF AREA WHERE PEBBLES ARE [m]   
                
                IF (R>1.0.AND.R<=1.886) THEN
                    porosity_PB=1.0-2.0/3.0*R**(-3.0)/(2.0/R-1.0)**0.5
                ELSEIF (R>1.886.AND.R<2.033) THEN
                    porosity_PB=1.8578-0.6649*R
                ELSEIF (R>=2.033) THEN
                    porosity_PB=0.151/(R-1.0)+0.36
                ENDIF
                
                porosity_PB=1.0-0.6806 !REVISE duty ratio 0.32 !! SET POROSITY BY SINAP DESIGN !!
                
            END FUNCTION porosity_PB
            
            
            !================================================================================
            !  FUNCTION: CALCULATE HEAT TRANSFER COEFFICIENT IN PEBBLE BED
            !================================================================================
            ! INPUTS
            !   W       :: MASS FLOW RATE [kg/s]
            !   DC      :: CHANNEL DAIMETER [m]
            !   DP      :: PEBBLE DAIMETER [m]
            !   T       :: COOLANT TEMPERATURE [Celcius]
            !   N       :: LOGICAL TRIP
            !       N=1 Wakao Correlation 
            !       N=2 Kunii and Levenspiel Correlation
            ! OUTPUT    
            !   HTC_PB  ::  Heat transfer coefficient of the pebble bed [W/m^2-K]            
            ! REFERENCE
            !
            !================================================================================
            REAL(8) FUNCTION HTC_PB(W,DC,DP,T,N)
                IMPLICIT NONE
                REAL(8), INTENT(IN)     :: W    ! Mass flow rate [kg/s]
                REAL(8), INTENT(IN)     :: DC   ! Channel diameter [m]
                REAL(8), INTENT(IN)     :: DP   ! Pebble diameter [m]
                REAL(8), INTENT(IN)     :: T    ! Temperature [Celcius]
                INTEGER, INTENT(IN)     :: N    ! Correlation Selector
                REAL(8)                 :: G    ! Mass flux = mass flow rate / area
                REAL(8)                 :: Re   ! Reynolds number
                REAL(8)                 :: Pr   ! Prandtl number
                REAL(8)                 :: Nu   ! Nusselt number
                
                ! AVERAGE FLUID PROPERTIES
                G=W/(PI*(DC/2.0)**2.0)          ! [kg/m^2-s]
                write(*,*) G, "Mass flow rate"
                Re=G*DP/flibe_mu(T)
                write(*,*) flibe_mu(T), "mu"
                write(*,*)
                Pr=flibe_cp(T)*flibe_mu(T)/flibe_k(T)

                ! HEAT TRANSFER CORRELATION
                SELECT CASE(N)
                CASE(1)
                    ! Wakao Kaviany's heat transfer handbook Ud=SUPERFICIAL VELOCITY Re:15-8500
                    Nu=2.0+1.1*Pr**(1.0/3.0)*Re**0.6 
                    IF (Re<=15.0.OR.Re>=8500) THEN
                        WRITE(*,*) Re,"Re overflow - Wakao"
                    ENDIF
                CASE(2)
                    !Kunii and Levenspiel
                    Nu=2.0+1.8*Pr**(1.0/3.0)*Re**0.5 
                END SELECT

                HTC_PB=Nu*flibe_k(T)/DP
                
            END FUNCTION HTC_PB
            
            
            !================================================================================
            !  FUNCTION: CALCULATE FRICTION FACTOR IN PEBBLE BED
            !================================================================================
            ! INPUTS
            !   W       :: MASS FLOW RATE [kg/s]
            !   DC      :: CHANNEL DAIMETER [m]
            !   DP      :: PEBBLE DAIMETER [m]
            !   T       :: COOLANT TEMPERATURE [Celcius]
            !   N       :: LOGICAL TRIP
            !       N=1 Ergu's Law 
            !       N=2 
            ! OUTPUT    
            !   FF_PB  ::  friction factor of the pebble bed [none]           
            ! REFERENCE
            !
            !================================================================================
            REAL(8) FUNCTION FF_PB(W,DC,DP,T,N)
                IMPLICIT NONE
                REAL(8), INTENT(IN)     :: W    ! Mass flow rate [kg/s]
                REAL(8), INTENT(IN)     :: DC   ! Channel diameter [m]
                REAL(8), INTENT(IN)     :: DP   ! Pebble diameter [m]
                REAL(8), INTENT(IN)     :: T    ! Temperature [?]
                INTEGER, INTENT(IN)     :: N    ! Correlation Selector
                REAL(8)                 :: E    ! Porosity
                REAL(8)                 :: G    ! Mass flux = mass flow rate / area
                REAL(8)                 :: Re   ! Reynolds number
                REAL(8)                 :: Pr   ! Prandtl number
                REAL(8)                 :: A
                REAL(8)                 :: B

                ! AVERAGE FLUID PROPERTIES
                E=porosity_pb(DC/DP)
                G=W/(PI*(DC/2.0)**2.0)      ! [kg/m^2-s]
                Re=G*DP/flibe_mu(T)
                ! Pr=flibe_cp(T)*flibe_mu(T)/flibe_k(T)

                ! FRICTION FACTOR CORRELATION
                SELECT CASE(N)
                CASE(1)
                    ! Ergun's law
                    A=180.0
                    B=1.8
                    FF_PB=(1.0-E)/E**(3.0)*(A*(1-E)/Re+B) 
                CASE(2)
                END SELECT

            END FUNCTION FF_PB
            
    END MODULE