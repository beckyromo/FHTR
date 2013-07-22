SUBROUTINE core(POWER,W_core,T_in,t_out,T_coolant_max,T_core_max,channel)
    
        USE global
        USE flibeprop, ONLY: flibe_cp, flibe_k, flibe_rho, flibe_mu, flibe_enthalpy, flibe_temperature
        USE trisoprop, ONLY: k_graphite, k_SiC, K_densePyC, k_PyC, k_UO2, k_TRISOlayer
        USE pipes, ONLY: FF_PIPE, dP_PIPE
        USE pebblebed, ONLY: FF_PB, HTC_PB
        
        IMPLICIT NONE

        !================================================================================
        ! Declare Variables
        !--> Declare/Initialize, Declare arrays, initialize array inputs
        !================================================================================
        
        integer,intent(in)  :: channel          ! 1 for average, 2 for hot channel
        
        integer             :: I                ! Loop counter
        integer             :: m                ! Loop counter
        integer             :: N_core           ! Number of nodes(CVs) the core is split into
        ! INPUTS
        real(8)             :: AFP(21)
        real(8)             :: LAMBDA(6)
        real(8)             :: BETA(6)
        ! UO2,PyC,DPyC,SiC,DPyC,TRISO+Graphite,Graphite
        real(8)             :: RPF(7)           ! Fuel radii [m]
        real(8)             :: KPF(7)           ! Fuel thermal conductivity [W/m-C]
        real(8)             :: RHOPF(7)         ! Fuel density [kg/m^3]
        real(8)             :: CPPF(7)          ! Fuel heat capacity [J/kg-C]
        real(8)             :: D_pebble         ! Diameter of pebble [m]                     
        real(8)             :: DD1              ! Face to face distance 1 of hexagonal core
        real(8)             :: DD2              ! Face to face distance 2 of hexagonal core
        real(8)             :: H_core           ! Height of the core [m]
        real(8)             :: FH               ! HCF - enthalpy rise hot channel factor (HCF)
        real(8)             :: FDTW             ! HCF - film/wall temperature rise
        real(8)             :: FDTF             ! HCF - fuel temperature rise
        real(8)             :: FCORE            ! HCF - 
        real(8)             :: FFUEL            ! HCF - 
        real(8)             :: FFDF             ! HCF - channel flow disparity factor
        real(8)             :: FKZ              ! Radial power peaking factor
        real(8)             :: FKTRISO          ! TRISO power peaking factor          
        real(8),intent(in)  :: POWER            ! Reactor power [W]
        real(8)             :: DECAY_HEAT       ! Decay heat of reactor [W]
        real(8)             :: POWER_NORM       ! Norminal Fisison Power [W]
        real(8),intent(inout)  :: W_core           ! Mass flow rate in the core [kg/s]
        real(8),intent(in)  :: T_in             ! Inlet temperature of the core [Celcius]
        
        real(8)             :: NTRISO           ! Number of TRISO particles per pebble
        real(8)             :: VolPF1           ! Volume of first layer of TRISO (UO2) [m^3]
        real(8)             :: VolPF5           ! Volume of first 5 layers of TRISO [m^3]
        real(8)             :: VolPF6           ! Volume of fuel matrix (TRISOs and graphite) of pebble [m^3]
        real(8)             :: VolPF7           ! Volume of pebble (fuel matrix plus graphite cladding) [m^3]
        real(8)             :: A_core           ! Area of the core [m^2]
        real(8)             :: D_core           ! Hydraulic diameter of the core [m]
        
        real(8)             :: FF_core          ! Friction factor in core 
        real(8)             :: HTC_core         ! Heat transfer coefficient in core [W/m^2-C]
        real(8)             :: A_HTC            ! Area of core with convection heat transfer [m^2]
        real(8)             :: ACTUAL_POWER     ! Actual Power [W]
        real(8)             :: dP_core          ! Core pressure drop [Pa] 
        real(8),allocatable :: enthalpy_core(:) ! Enthalpy in each node [J?/kg]
        real(8),allocatable :: T_w_core(:)      ! Pebble surface temperature at each node  [Celcius]
        real(8),allocatable :: T_CL_core(:)     ! Centerline temperature at each node  [Celcius]
        real(8),allocatable :: T_CL_TRISO(:)    ! Centerline temperature of centerline TRISO at each node [C]
        real(8),intent(out) :: T_coolant_max    ! Maximum coolant temperature [Celcius]
        real(8),intent(out) :: T_core_max       ! Maximum core temperature  [Celcius]
        real(8),intent(out) :: T_out            ! Outlet temperature of the core [Celcius]
        real(8)             :: TC               ! Temperature of coolant in core [Celcius]
        real(8)             :: TG               ! Temperature at inner edge of graphite cladding [Celcius]
        real(8)             :: TW               ! Temperature at the pebble surface [Celcius]
        real(8)             :: T_temp           ! Temperature place holder [Celcius]
        real(8)             :: NPPN             ! Number of pebbles in each node(CV)
        real(8)             :: PPP              ! Power per pebble [W]
        real(8)             :: PPT              ! Power per TRISO particle [W]
        real(8)             :: TOL              ! While loop error tolerance
        real(8)             :: U                ! Thermal resistance

        
        
        !================================================================================
        ! Initialize Variables / Inputs
        !================================================================================
        
        ! POWER DISTRIBUTION FROM SINAP REPORT
        DATA AFP/1.3522E-002,4.7941E-002,4.6349E-002,4.7030E-002,4.8748E-002,5.0615E-002,5.2291E-002,&
                 5.3749E-002,5.4762E-002,5.5298E-002,5.5540E-002,5.5100E-002,5.4286E-002,5.3050E-002,&
                 5.1321E-002,4.9228E-002,4.7067E-002,4.4800E-002,4.2871E-002,4.2380E-002,3.4045E-002/
        
        ! LAMBDA AND BETA FOR POINT KINETICS
        DATA LAMBDA/1.271E-2,3.174E-2,1.16E-1,3.11E-1,1.4E0,3.87E0/
        DATA BETA/2.3794E-4, 1.5831E-3,1.4183E-3, 2.8550E-3,8.3273E-4,3.0192E-4/
        
        ! ENGINERRING HOT CHANNEL FACTORS -------------------------------------
        FH=1.173_8
        FDTW=1.275_8
        FDTF=1.123_8
        FCORE=0.965_8
        FFUEL=0.940_8
        FFDF=0.8_8
        FKZ=1.4_8
        FKTRISO=1.0_8

        ! CORE
        !POWER=20.0E6_8    ! [W]
        DECAY_HEAT=0.0_8  ! [W]
        POWER_NORM=1.0_8  ! [W]
        !W_core=84.6_8     ! [kg/s]
        !T_in=600.0_8      ! [Celcius]
        DD1=134.69E-2_8   ! Face to face distance 1 [m]
        DD2=139.0E-2_8    ! Face to face distance 2 [m]
        H_core=138.6_8    ! Height of the core [m]
        D_pebble=6.0E-2_8 ! Diameter of fuel pebble [m]
        N_core=21         ! Number of nodes(CVs) in the core

        
        ! FUEL: UO2,PyC,DPyC,SiC,DPyC,TRISO+Graphite,Graphite (
        DATA	RPF/0.25E-3,0.345E-3,0.385E-3,0.42E-3,0.46E-3,25.0E-3,30.0E-3/  ! Radius [m]
        DATA	KPF/2.5,0.5,4.0,13.9,4.0,30.0,30.0/                             ! Thermal conductivity [W/m-C]
        DATA	RHOPF/10.5E3,1.05E3,1.9E3,3.18E3,1.9E3,1.73E3,1.73E3/           ! Density [kg/m^3]
        DATA	CPPF/332,1.5,1.5,1.5,1.5,710,710/                               ! Heat capacity [J/kg-K]
        
        ! Calculated inputs
        VolPF5=4.0/3.0*PI*RPF(5)**3.0
        VolPF6=4.0/3.0*PI*RPF(6)**3.0
        VolPF7=4.0/3.0*PI*RPF(7)**3.0
        VolPF1=4.0/3.0*PI*RPF(1)**3.0
        NTRISO=VolPF6*7.5/100.0/VolPF5
        A_core=DD1**2-((DD1*2.0**0.5-DD2)/2.0**0.5)**2*2.0-0.35*0.35
        D_core=(4.0*A_core/PI)**0.5                                     ! HYDRAULIC RADIUS OF THE CORE [m]
        ACTUAL_POWER=POWER_NORM*POWER+DECAY_HEAT        

        
        ! Set criteria for average (1) or hot (2) channel
        IF (channel==1) THEN
            FH=1.0_8
            FKZ=1.0_8
            FDTW=1.0_8
            FDTF=1.0_8
            FFDF=1.0_8
        ELSE IF (channel==2) THEN
            W_core=W_core*FFDF
        END IF
            
                
       
        !================================================================================       
        ! Core Loop Calculations - Average channel
        !================================================================================         

        ! Initialize variables for loop
        dP_core=0.0
        T_coolant_max=0.0
        T_core_max=0.0        
        ! Allocate arrays
        allocate(enthalpy_core(N_core)) ! Array of enthalpies in each node(CV) of the core
        allocate(T_w_core(N_core))      ! Array of surface tempertures in each node(CV) of the core
        allocate(T_CL_core(N_core))     ! Array of central temperature in each node(CV) of the core
        allocate(T_CL_TRISO(N_core))    ! Array of central temperature in TRISO in each node(CV) of the core
        
        ! Loop through each core CV node
        DO I=1,N_core
            
            IF (I==1) THEN
                enthalpy_core(I)=flibe_enthalpy(T_in)+ACTUAL_POWER*AFP(I)*FH*FKZ/W_core
                TC=( flibe_temperature(enthalpy_core(I)) + T_in ) / 2.0       ! Average coolant temperature in node
            ELSE
                enthalpy_core(I)=enthalpy_core(I-1)+ACTUAL_POWER*AFP(I)*FH*FKZ/W_core
                TC=( flibe_temperature(enthalpy_core(I)) + flibe_temperature(enthalpy_core(I-1)) ) / 2.0
            END IF   
            
            ! Calculate pressure drop in node
            IF (channel==1) THEN
                FF_core=FF_PB(W_core,D_core,D_pebble,TC,1)
                dP_core=dP_core+0.5*(H_core/N_core)*(W_core/A_core)**2.0/(flibe_rho(TC)*D_pebble) &
                        + flibe_rho(TC)*GRAVITY*H_core/N_core
            END IF
            
            ! Calculate heat transfer coefficient in node
            HTC_core=HTC_PB(W_core,D_core,D_pebble,TC,1)
            
            ! Find wall temperature
            IF (I==1) THEN
                NPPN=140.0                      ! Number of pebbles
            ELSE IF (I==21) THEN
                NPPN=416.0                      ! Number of pebbles
            ELSE
                NPPN=556.0                      ! Number of pebbles
            END  IF
            A_HTC=NPPN*PI*D_pebble**2.0
            TW=ACTUAL_POWER*AFP(I)*FKZ/(HTC_core*A_HTC) + TC
            TW=(TW-TC)*FDTW+TC
            T_w_core(I)=TW
            
            ! Check if temperature is maximum
            IF (TW>T_coolant_max) THEN
                T_coolant_max=TW
            END IF
            
            ! Calculate power in pebble and centerline TRISO particle
            PPP=ACTUAL_POWER*AFP(I)*FKZ/NPPN        ! Power per pebble
            PPT=PPP/NTRISO                      ! Power per TRISO particle

            ! Calculate thermal conductivties and find temperatures
            TOL=1.0e-5           
                        
            ! GRAPHITE CLADDING
            T_temp=TW
            TG=0.0
            do while (abs(T_temp-TG)>TOL)
                T_temp=TG
                KPF(7)=k_graphite( (TW+T_temp)/2.0 )
                U=1.0/(4.0*PI*KPF(7)) * (1.0/RPF(6) - 1.0/RPF(7))
                TG=PPP*U + TW
            end do
            
            ! FUEL MATRIX (TRISOs IN GRAPHITE MATRIX)
            T_temp=TG
            T_CL_core(I)=0.0
            do while (abs(T_temp-T_CL_core(I))>TOL)
                T_temp=T_CL_core(I)
                KPF(6)=k_graphite( (TG+T_temp)/2.0 )
                T_CL_core(I)=PPP*(RPF(6)**2.0/6.0/KPF(6)/VolPF6) + TG
            end do
            
            ! CENTERLINE TRISO PARTICLE
            T_temp=T_CL_core(I)
            U=0.0
            do m=2,5
                U=1.0/(4.0*PI*k_trisolayer(m,T_temp)) * (1.0/RPF(m-1) - 1.0/RPF(m)) + U
            end do
            T_CL_TRISO(I)=PPT*(RPF(1)**2.0/6.0/k_trisolayer(1,T_temp)/VolPF1+U) + T_CL_core(I)
            T_CL_TRISO(I)=(T_CL_TRISO(I)-TW)*FDTF+TW
            
            ! Check if temperature is max
            if (T_CL_TRISO(I)>T_core_max) then
                T_core_max=T_CL_TRISO(I)
            end if
            
            write(*,*) PPP, PPT, TC, TW, TG, T_CL_core(I), T_CL_TRISO(I)
            write(10,'(F9.3,F9.3,F9.3,F9.3,F9.3)') TC, TW, TG, T_CL_core(I), T_CL_TRISO(I)
        
        end do
        
        
        ! W_core=W_core/FFDF
        T_out=flibe_temperature(enthalpy_core(N_core))    
        
        
        !================================================================================
        ! Print LSSS
        !================================================================================  
        write(*,*)
        write(*,*) "                           LSSS Results                           "
        write(*,*) "   POWER      W      T_IN    T_OUT    TC_MAX   TF_MAX"
        write(*,'(F9.2,F9.3,F9.3,F9.3,F9.3,F9.3)') &
                    POWER/1.0E6,W_core,T_in,T_out,T_coolant_max,T_core_max
        write(*,*)
        
        write(10,*)
        write(10,*) "                           LSSS Results                           "
        write(10,*) "   POWER      W      T_IN    T_OUT    TC_MAX   TF_MAX"
        write(10,'(F9.2,F9.3,F9.3,F9.3,F9.3,F9.3)') &
                    POWER/1.0E6,W_core,T_in,T_out,T_coolant_max,T_core_max
        write(10,*)
        write(10,*) "=================================================================="
        write(10,*)
                
        
        !================================================================================
        ! Deallocate all arrays
        !================================================================================
        deallocate(enthalpy_core)  
        deallocate(T_w_core)       
        deallocate(T_CL_core)      
        deallocate(T_CL_TRISO)
        
    END SUBROUTINE core