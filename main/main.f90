    !****************************************************************************
    !
    !  PROGRAM: main
    !
    !  PURPOSE: Entry point for the main application.
    !
    !**************************************************************************** 
    PROGRAM main    
    
        USE pebblebed, ONLY: pebblebedLSSSloop,pebblebedcore
            
        IMPLICIT NONE
        
        real(8)             :: POWER                ! Reactor power [W]
        !real(8)             :: DECAY_HEAT           ! Decay heat of reactor [W]
        !real(8)             :: POWER_NORM           ! Norminal Fisison Power [W]
        real(8)             :: W_core               ! Mass flow rate in the core [kg/s]
        real(8)             :: T_in                 ! Inlet temperature of the core [Celcius]        
        real(8)             :: T_fuel_limit         ! LSSS temoerature limit for fuel (max fuel temp in hot channel) 
        real(8)             :: T_coolant_limit      ! LSSS temperature limit for coolant (max coolant temp in hot channel)
        real(8)             :: T_out_limit          ! LSSS temperature limit for T_out (max T_out temp in average channel)
        real(8)             :: T_out                ! Maximum coolant temperature [Celcius]
        real(8)             :: T_coolant_max        ! Maximum core temperature  [Celcius]
        real(8)             :: T_core_max           ! Outlet temperature of the core [Celcius]     
        
        
        !================================================================================
        ! INPUTS
        !================================================================================
        POWER=20.0E6_8              ! [W]
        W_core=76.14_8              ! [kg/s]
        T_in=470.0_8                ! [Celcius]
        T_fuel_limit=1300.0_8       ! [Celcius]
        T_coolant_limit=1200.0_8    ! [Celcius]
        T_out_limit=720.0_8         ! [Celcius]
        
        
        !================================================================================
        ! Start Calculations
        !================================================================================
        write(*,*)
        write(*,*) 'Begin Program'

        open(10,FILE="output.txt")

        open(20,FILE="LSSS.txt")

        !================================================================================
        ! Calculates LSSS for pebblebed core
        !================================================================================
        call pebblebedLSSSloop(W_core,T_in,T_fuel_limit,T_coolant_limit,T_out_limit,T_out,T_coolant_max,T_core_max)
        
        write(10,*) "SINGLE POWER LEVEL OUTPUT OF ", POWER
        call pebblebedcore(POWER,W_core,T_in,T_out,T_coolant_max,T_core_max,1)

        
                                    
        close(10)
        close(20)
        write(*,*) 'Program Complete'
               
        
        !================================================================================
        ! Wait for user
        !================================================================================
        ! read(*,*)


    CONTAINS


    END PROGRAM main
