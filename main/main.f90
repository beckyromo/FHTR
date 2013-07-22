    !****************************************************************************
    !
    !  PROGRAM: main
    !
    !  PURPOSE: Entry point for the main application.
    !
    !**************************************************************************** 
    PROGRAM main
 
        USE trisoprop, ONLY: k_graphite, k_SiC, K_densePyC, k_PyC, k_UO2, k_TRISOlayer
            
        IMPLICIT NONE
        
        integer             :: I                ! Loop counter
        real(8)             :: TC               ! Temperature [Celcius]
        real(8)             :: POWER            ! Reactor power [W]
        real(8)             :: DECAY_HEAT       ! Decay heat of reactor [W]
        real(8)             :: POWER_NORM       ! Norminal Fisison Power [W]
        real(8)             :: W_core           ! Mass flow rate in the core [kg/s]
        real(8)             :: T_in             ! Inlet temperature of the core [Celcius]
        real(8)             :: T_out            ! Maximum coolant temperature [Celcius]
        real(8)             :: T_coolant_max    ! Maximum core temperature  [Celcius]
        real(8)             :: T_core_max       ! Outlet temperature of the core [Celcius]
        
        ! INPUTS
        POWER=20.0E6_8    ! [W]
        W_core=76.14_8     ! [kg/s]
        T_in=470.0_8      ! [Celcius]
        
        
        !================================================================================
        ! Start Calculations
        !================================================================================       
        write(*,*) 'Begin Program'
        
        open(10,FILE="output.txt")
        open(20,FILE="LSSS.txt")
        write(20,*) "                           LSSS Results                           "
        write(20,*) "   POWER      W      T_IN    T_OUT    TC_MAX   TF_MAX"
        
        
        !================================================================================
        ! Output Thermal Conductivities from 700 to 1500 Celcius
        !================================================================================
        !write(10,*) "  THERMAL CONDUCTIVITY OF TRISO PARTICLE MATERIALS   "
        !write(10,*) "       T       GF      SiC     DPyC      PyC      UO2"
        !do I=700,1500,50
        !    TC=I
        !    write(10,'(F9.2,F9.3,F9.3,F9.3,F9.3,F9.3)') &
        !            TC, k_graphite(TC), k_SiC(TC), k_densePyC(TC), k_PyC(TC), k_UO2(TC)
        !end do
        
        
        !================================================================================
        ! Calculate core temperatures for average channel
        !================================================================================
        write(10,*)
        write(10,*) "      TEMPERATURES IN AVERAGE CHANNEL       "
        write(10,*) "      TC       TW       TG      TCL     TCLT"
        
        call core(POWER,W_core,T_in,T_out,T_coolant_max,T_core_max,1)
        write(20,'(F9.2,F9.3,F9.3,F9.3,F9.3,F9.3)') &
                    POWER/1.0E6,W_core,T_in,T_out,T_coolant_max,T_core_max       
        
        !================================================================================
        ! Calculate core temperatures for hot channel
        !================================================================================
        write(10,*)
        write(10,*) "        TEMPERATURES IN HOT CHANNEL         "
        write(10,*) "      TC       TW       TG      TCL     TCLT"
        
        call core(POWER,W_core,T_in,T_out,T_coolant_max,T_core_max,2)
        
        


        
        close(10)
        close(20)
        write(*,*) 'Program Complete'
               
        
        !================================================================================
        ! Wait for user
        !================================================================================
        read(*,*)


    CONTAINS


    END PROGRAM main
