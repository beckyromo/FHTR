!********************************************************************************
!
!  MODULE:  prismatic
!
!  PURPOSE: Contains functions for prismatic core calculations and subroutines 
!           to do LSSS calculations for a prismatic core
!
!  FUNCTIONS:
!  HTC
!  FF
!
!  SUBROUTINES:
!  prismaticcore        calculates core temperatures for given power, mass flow 
!                       rate, and inlet temperature for either average or hot 
!                       channel, prints core temperature results and outputs the 
!                       outlet temperature, maximum coolant temperature, and 
!                       maximum core/fuel temperature
!  prismaticLSSSloop    calls subroutine pebblebedcore iterating on power to 
!                       print LSSS for given fuel, coolant, and outlet 
!                       temperature limits
!
!******************************************************************************** 
MODULE prismatic
            
    USE global
    USE flibeprop, ONLY: flibe_cp, flibe_k, flibe_rho, flibe_mu
        
    IMPLICIT NONE


CONTAINS
    
    
    !================================================================================
    !  SUBROUTINE: prismaticcore
    !================================================================================
    ! INPUTS
    !   POWER           :: Power [W]
    !   W_core          :: Mass flow rate [kg/s]
    !   T_in            :: Inlet temperature [Celcius]
    !   channel         :: Selects for average (1) or hot (2) channel calculation
    ! OUTPUTS    
    !   T_out           :: Outlet temperature [Celcius]
    !   T_coolant_max   :: Maximum coolant temperature [Celcius]
    !   T_core_max      :: Maximum core/fuel temperature [Celcius]
    ! REFERENCE
    !   Yao's Code, which uses power distribution data
    !================================================================================
    SUBROUTINE prismaticcore(POWER,W_core,T_in,T_out,T_coolant_max,T_core_max,channel)
    
        USE global
        USE flibeprop, ONLY: flibe_cp, flibe_k, flibe_rho, flibe_mu, flibe_enthalpy, flibe_temperature
        USE trisoprop, ONLY: k_graphite, k_SiC, K_densePyC, k_PyC, k_UO2, k_TRISOlayer
        USE pipes, ONLY: FF_PIPE, dP_PIPE
        
        IMPLICIT NONE

        !================================================================================
        ! Declare Variables
        !--> Declare/Initialize, Declare arrays, initialize array inputs
        !================================================================================
        
        integer,intent(in)  :: channel          ! 1 for average, 2 for hot channel
        
        integer             :: I                ! Loop counter
        integer             :: m                ! Loop counter
        integer             :: N_core           ! Number of nodes(CVs) the core is split into
        real(8)             :: AFP(21)
        
        ! UO2,PyC,DPyC,SiC,DPyC,TRISO+Graphite,Graphite
        real(8)             :: RPF(7)           ! Fuel radii [m]
        real(8)             :: DPF(7)           ! Fuel diameter [m]
        real(8)             :: KPF(7)           ! Fuel thermal conductivity [W/m-C]
        real(8)             :: RHOPF(7)         ! Fuel density [kg/m^3]
        real(8)             :: CPPF(7)          ! Fuel heat capacity [J/kg-C]

        real(8)             :: DD1              ! Core: Face to face distance 1 of hexagonal core
        real(8)             :: DD2              ! Core: Face to face distance 2 of hexagonal core
        real(8)             :: H_core           ! Core: height [m]
        real(8)             :: A_core           ! Core: area of the core [m^2]
        real(8)             :: H_block          ! Height of fuel assemply block [m]
        real(8)             :: D_core           ! Core: hydraulic diameter of the core [m]
        real(8)             :: D_fuel           ! Diameter of fuel channel [m]
        real(8)             :: D_cool           ! Dimaeter of coolant channel [m]       
        real(8)             :: d_cell           ! Unit cell: flat to flat distance (2x short radius) [m] d=sqrt(3)*t
        real(8)             :: t_cell           ! Unit cell: side length of hexagon
        real(8)             :: A_cell           ! Unit cell: area
        real(8)             :: A_fuel           ! Unit cell: fuel area
        real(8)             :: A_cool           ! Unit cell: coolant area
        real(8)             :: A_grapghite      ! Unit cell: graphite block area
         
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
        real(8),intent(in)  :: W_core           ! Mass flow rate in the core [kg/s]
        real(8)             :: W                ! Mass flow rate in the core [kg/s]
        real(8),intent(in)  :: T_in             ! Inlet temperature of the core [Celcius]
        
        real(8)             :: NTRISO           ! Number of TRISO particles per pebble
        real(8)             :: NPPN             ! Number of pebbles in each node(CV)
        real(8)             :: PPP              ! Power per pebble [W]
        real(8)             :: PPT              ! Power per TRISO particle [W]
        
        real(8)             :: VolPF1           ! Volume of first layer of TRISO (UO2) [m^3]
        real(8)             :: VolPF5           ! Volume of first 5 layers of TRISO [m^3]
        real(8)             :: VolPF6           ! Volume of fuel matrix (TRISOs and graphite) of pebble [m^3]
        real(8)             :: VolPF7           ! Volume of pebble (fuel matrix plus graphite cladding) [m^3]
       
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
        real(8)             :: TOL              ! While loop error tolerance
        real(8)             :: U                ! Thermal resistance

        
        
        !================================================================================
        ! Initialize Variables / Inputs
        !================================================================================
        
        ! POWER DISTRIBUTION FROM SINAP REPORT
        DATA AFP/1.3522E-002,4.7941E-002,4.6349E-002,4.7030E-002,4.8748E-002,5.0615E-002,5.2291E-002,&
                 5.3749E-002,5.4762E-002,5.5298E-002,5.5540E-002,5.5100E-002,5.4286E-002,5.3050E-002,&
                 5.1321E-002,4.9228E-002,4.7067E-002,4.4800E-002,4.2871E-002,4.2380E-002,3.4045E-002/
        
        ! ENGINERRING HOT CHANNEL FACTORS -------------------------------------
        FFDF=0.8_8
        FH=1.173_8      ! 1.173 is only MITR value, 1.172 also uses some HTC-10 pebblebed
        FDTW=1.275_8    ! 1.275 is only MITR value, 1.168 also uses some HTC-10 pebblebed
        FDTF=1.123_8    ! 1.123 is only MITR value, 1.172 also uses some HTC-10 pebblebed
        ! See ICONE paper of Yao's for MITR values, see Yao's report for HTC-10 values
        FCORE=0.965_8       ! NOT USED
        FFUEL=0.940_8       ! NOT USED
        FKZ=1.4_8
        FKTRISO=1.0_8       ! NOT USED

        ! CORE
        DECAY_HEAT=0.0_8        ! [W]
        POWER_NORM=1.0_8        ! [W]
        DD1=134.69E-2_8         ! Face to face distance 1 [m]
        DD2=139.0E-2_8          ! Face to face distance 2 [m]
        H_block=0.793           ! Height of a block [m]
        H_core=H_block*10.0_8   ! Height of the core [m]
        N_core=21               ! Number of nodes(CVs) in the core
        D_fuel=12.7E-3_8        ! Diameter of fuel channel [m]
        D_cool=9.53E-3_8        ! Diameter of coolant channel [m] or revised value of 1.4 cm
        
        ! Unit Cell
        !d_cell=0.24_8                          ! Flat to flat distance (2x short radius) [m] d=sqrt(3)*t
        t_cell=18.8E-3_8                        ! Unit cell side length of hexagon [m]
        A_cell=3.0*sqrt(3.0)/2.0*t_cell**2.0    ! sqrt(3)/2.0*d_cell**2.0
        A_fuel=PI*D_fuel**2.0/4.0*6.0/3.0
        A_cool=PI*D_cool**2.0/4.0
        A_grapghite=A_cell-A_fuel-A_cool
    
        
        ! FUEL: UO2,PyC,DPyC,SiC,DPyC,TRISO+Graphite,Graphite (
        DATA	DPF/0.350E-3,0.550E-3,0.620E-3,0.690E-3,0.770E-3,25.0E-3,30.0E-3/   ! Diameter [m]
        DPF(6)=D_fuel
        RPF=DPF/2.0                                                                 ! Radius [m]
        DATA	KPF/2.5,0.5,4.0,13.9,4.0,30.0,30.0/                                 ! Thermal conductivity [W/m-C] NOTE OVERWRITTEN
        !DATA	RHOPF/10.4E3,1.0E3,1.85E3,3.2E3,1.8E3,1.74E3,1.74E3/                 ! Density [kg/m^3]
        
        ! Calculated inputs
        VolPF1=4.0/3.0*PI*RPF(1)**3.0   ! Fuel kernel volume
        VolPF5=4.0/3.0*PI*RPF(5)**3.0   ! TRISO volume
        VolPF6=A_fuel*H_core/N_core     ! Fuel matrix volume
        VolPF7=4.0/3.0*PI*RPF(7)**3.0
        NTRISO=VolPF6*7.5/100.0/VolPF5
        A_core=DD1**2-((DD1*2.0**0.5-DD2)/2.0**0.5)**2*2.0-0.35*0.35
        D_core=9.2 !(4.0*A_core/PI)**0.5                                     ! HYDRAULIC RADIUS OF THE CORE [m]
        A_core=PI*D_core**2.0/4.0
        ACTUAL_POWER=POWER_NORM*POWER+DECAY_HEAT        

        
        ! Set criteria for average (1) or hot (2) channel
        IF (channel==1) THEN
            W=W_core
            FH=1.0_8
            FKZ=1.0_8
            FDTW=1.0_8
            FDTF=1.0_8
            FFDF=1.0_8
        ELSE IF (channel==2) THEN
            W=W_core*FFDF
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
                enthalpy_core(I)=flibe_enthalpy(T_in)+ACTUAL_POWER*AFP(I)*FH*FKZ/W
                TC=( flibe_temperature(enthalpy_core(I)) + T_in ) / 2.0       ! Average coolant temperature in node
            ELSE
                enthalpy_core(I)=enthalpy_core(I-1)+ACTUAL_POWER*AFP(I)*FH*FKZ/W
                TC=( flibe_temperature(enthalpy_core(I)) + flibe_temperature(enthalpy_core(I-1)) ) / 2.0
            END IF   
            
            ! Calculate pressure drop in node
            IF (channel==1) THEN
                FF_core=FF(W,D_cool,TC,1)
                dP_core=dP_core+0.5*(H_core/N_core)*(W/A_core)**2.0/(flibe_rho(TC)*D_cool) &
                        + flibe_rho(TC)*GRAVITY*H_core/N_core
            END IF
            
            ! Calculate heat transfer coefficient in node
            HTC_core=HTC(W,D_cool,TC,1)
            
            ! Find wall temperature
            A_HTC=PI*D_cool*H_core/N_core
            TW=ACTUAL_POWER*AFP(I)*FKZ/(HTC_core*A_HTC) + TC
            TW=(TW-TC)*FDTW+TC
            T_w_core(I)=TW
            
            ! Check if temperature is maximum
            IF (TW>T_coolant_max) THEN
                T_coolant_max=TW
            END IF
            
            ! Calculate power in pebble and centerline TRISO particle
            !PPP=ACTUAL_POWER*AFP(I)*FKZ/NPPN        ! Power per pebble
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
            
            ! Apply hot channel factors
            T_CL_TRISO(I)=(T_CL_TRISO(I)-TW)*FDTF+TW
            T_CL_core(I)=(T_CL_core(I)-TW)*FDTF+TW
            
            ! Check if temperature is max
            if (T_CL_core(I)>T_core_max) then
                T_core_max=T_CL_core(I)
            end if
            
            ! Print results to file output.txt
            if (I==1) then
                write(10,*)
                if (channel==1) then
                    write(10,*) "      TEMPERATURES IN AVERAGE CHANNEL       "
                else if (channel==2) then
                    write(10,*) "        TEMPERATURES IN HOT CHANNEL         "
                end if
                write(10,*) "      TC       TW       TG      TCL     TCLT"
            end if
            write(10,'(F9.3,F9.3,F9.3,F9.3,F9.3)') TC, TW, TG, T_CL_core(I), T_CL_TRISO(I)
            ! write to command line        
            !write(*,*) PPP, PPT, TC, TW, TG, T_CL_core(I), T_CL_TRISO(I)
            
        
        end do
        
        
        T_out=flibe_temperature(enthalpy_core(N_core))    
        
        
        !================================================================================
        ! Print LSSS
        !================================================================================  
        !write(*,*)
        !write(*,*) "                           LSSS Results                           "
        !write(*,*) "   POWER      W      T_IN    T_OUT    TC_MAX   TF_MAX"
        !write(*,'(F9.2,F9.3,F9.3,F9.3,F9.3,F9.3)') &
        !            POWER/1.0E6,W,T_in,T_out,T_coolant_max,T_core_max
        !write(*,*)
        
        write(10,*)
        write(10,*) "                           LSSS Results                           "
        write(10,*) "   POWER      W      T_IN    T_OUT    TC_MAX   TF_MAX"
        write(10,'(F9.2,F9.3,F9.3,F9.3,F9.3,F9.3)') &
                    POWER/1.0E6,W,T_in,T_out,T_coolant_max,T_core_max
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
        
    END SUBROUTINE prismaticcore
    
    
    
    
    
    
   ! !================================================================================
   ! !  SUBROUTINE: prismaticLSSSloop  [W/m-K]
   ! !================================================================================
   ! ! INPUTS    
   ! ! OUTPUTS   
   ! ! REFERENCE 
   ! !================================================================================
   ! SUBROUTINE prismaticLSSSloop(W_core,T_in,T_fuel_limit,T_coolant_limit,T_out_limit,T_out,T_coolant_max,T_core_max)
   !     
   !     IMPLICIT NONE
   !     
   !     integer             :: I                    ! Loop counter
   !     logical             :: fuelflag             ! Flag if max fuel temperature has been met in hot channel
   !     logical             :: coolantflag          ! Flag if max coolant temperature has been met in hot channel
   !     logical             :: Toutflag             ! Flag if max coolant outlet temperature has been met in average channel
   !     real(8)             :: T_fuel_limit         ! LSSS temoerature limit for fuel (max fuel temp in hot channel) 
   !     real(8)             :: T_coolant_limit      ! LSSS temperature limit for coolant (max coolant temp in hot channel)
   !     real(8)             :: T_out_limit          ! LSSS temperature limit for T_out (max T_out temp in average channel)
   !     real(8)             :: TC                   ! Temperature [Celcius]
   !     real(8)             :: POWER                ! Reactor power [W]
   !     real(8)             :: DECAY_HEAT           ! Decay heat of reactor [W]
   !     real(8)             :: POWER_NORM           ! Norminal Fisison Power [W]
   !     real(8)             :: W_core               ! Mass flow rate in the core [kg/s]
   !     real(8)             :: T_in                 ! Inlet temperature of the core [Celcius]
   !     real(8)             :: T_out                ! Maximum coolant temperature [Celcius]
   !     real(8)             :: T_out_avg            ! Holder for T_out in average channel
   !     real(8)             :: T_coolant_max        ! Maximum core temperature  [Celcius]
   !     real(8)             :: T_core_max           ! Outlet temperature of the core [Celcius]
   !     real(8)             :: step                 ! Step size for power do loop
   !     real(8)             :: T_in_coolantrange    
   !     real(8)             :: POWER_coolantrange   
   !     real(8)             :: T_in_fuelrange       
   !     real(8)             :: POWER_fuelrange
   !
   !         
   !     write(20,*) "                                               LSSS Results"
   !     write(20,*)
   !     write(20,'(A54)',ADVANCE='NO') "                           AVG CHANNEL                             "    
   !     write(20,'(A54)')              "                           HOT CHANNEL                             "
   !     write(20,'(A54)',ADVANCE='NO') " =================================================================="    
   !     write(20,'(A54)')              " =================================================================="
   !     write(20,'(A54)',ADVANCE='NO') "   POWER      W      T_IN    T_OUT    TC_MAX   TF_MAX"    
   !     write(20,'(A54)')              "   POWER      W      T_IN    T_OUT    TC_MAX   TF_MAX"
   !     write(20,'(A54)',ADVANCE='NO') " ------------------------------------------------------------------"    
   !     write(20,'(A54)')              " ------------------------------------------------------------------"
   !         
   !     !================================================================================
   !     ! Find power level for LSSS limits for given T_in and mass flow rate 
   !     !================================================================================
   !     ! Reset flags to .false.
   !     fuelflag=.false. 
   !     coolantflag=.false.
   !     Toutflag=.false.
   !     write(*,*)
   !     write(*,'(A,F5.2,A,F7.2,A)') "For primary mass flow rate of ", W_core, " and  inlet temperature of ", T_in, ":"
   !     write(*,'(A)') "======================================================================"
   !     write(*,*)
   !     do POWER=0.0E6,50.0E6,0.01E6
   !         !================================================================================
   !         ! Calculate core temperatures for average channel
   !                     call pebblebedcore(POWER,W_core,T_in,T_out,T_coolant_max,T_core_max,1)
   !         write(20,'(F9.2,F9.3,F9.3,F9.3,F9.3,F9.3)', ADVANCE='NO') &
   !                     POWER/1.0E6,W_core,T_in,T_out,T_coolant_max,T_core_max
   !         ! Check for exceedance of LSSS T_out temperature limit
   !         if (T_out>=T_out_limit .AND. Toutflag==.false.) then
   !             write(*,'(A,F7.2,A)') "The maximum outlet average channel temperature exceeds ", T_out_limit, " at:"
   !             write(*,'(F5.2,A,F7.2,A,F7.2)') POWER/1.0E6, " MW    T_in = ", T_in, "   T_out = ", T_out
   !             write(*,*) 
   !             Toutflag=.true.
   !         end if
   !         !================================================================================
   !         ! Calculate core temperatures for hot channel
   !                     call pebblebedcore(POWER,W_core,T_in,T_out,T_coolant_max,T_core_max,2)
   !         write(20,'(F9.2,F9.3,F9.3,F9.3,F9.3,F9.3)') &
   !                     POWER/1.0E6,W_core,T_in,T_out,T_coolant_max,T_core_max
   !         ! Check for exceedance of LSSS coolant temperature limit   
   !         if (T_coolant_max>=T_coolant_limit .AND. coolantflag==.false.) then
   !             write(*,'(A,F7.2,A)') "The maximum  coolant hot channel   temperature exceeds ", T_coolant_limit, " at:"
   !             write(*,'(F5.2,A,F7.2,A,F7.2,A,F7.2)') POWER/1.0E6, " MW    T_in = ", T_in, "   T_out = ", T_out, "   TCMAX = ", T_coolant_max
   !             write(*,*) 
   !             coolantflag=.true.
   !         end if
			!! Check for exceedance of LSSS fuel temperature limit
   !         if (T_core_max>=T_fuel_limit .AND. fuelflag==.false.) then
   !             write(*,'(A,F7.2,A)') "The maximum    fuel hot channel    temperature exceeds ", T_fuel_limit, " at:"
   !             write(*,'(F5.2,A,F7.2,A,F7.2,A,F7.2)') POWER/1.0E6, " MW    T_in = ", T_in, "   T_out = ", T_out, "   TFMAX = ", T_core_max
   !             write(*,*) 
   !             fuelflag=.true.
   !         end if
   !         
   !     end do
   !
   !     
   !     
   !     !================================================================================
   !     ! Find power level for LSSS limits for given T_out and mass flow rate :: Find T_in and POWER range
   !     !================================================================================
   !     ! Reset flags to .false.
   !     fuelflag=.false. 
   !     coolantflag=.false.
   !     Toutflag=.false.
   !     write(*,*)
   !     write(*,'(A,F5.2,A,F7.2,A)') "For primary mass flow rate of ", W_core, " and outlet temperature of ", T_out_limit, ":"
   !     write(*,'(A)') "======================================================================"
   !     do T_in=470.0_8,T_out_limit,1.0_8 
   !         ! Check if all LSSS limits met
   !         if (Toutflag==.true. .AND. coolantflag==.true. .AND. fuelflag==.true.) then
   !             write(*,*) "LSSS Limits Range Found"
   !             write(*,*)
   !             EXIT
   !         else          
   !             ! Loop through power level to find T_out=T_out_limit
   !             ! Reset T_out flags to .false.
   !             Toutflag=.false.
   !             do POWER=0.0E6,50.0E6,1.0E6
   !                 if (Toutflag==.false.) then
   !                     ! Calculate core temperatures for average channel
   !                     call pebblebedcore(POWER,W_core,T_in,T_out,T_coolant_max,T_core_max,1)
   !                     T_out_avg=T_out
   !                     ! Check for exceedance of LSSS T_out temperature limit
   !                     if (T_out>=T_out_limit .AND. Toutflag==.false.) then
   !                         Toutflag=.true.
   !                         ! Calculate core temperatures for hot channel
   !                         call pebblebedcore(POWER,W_core,T_in,T_out,T_coolant_max,T_core_max,2)
   !                         if (T_coolant_max <= T_coolant_limit .AND. coolantflag==.false.) then
   !                             !write(*,*) T_in, T_out_avg
   !                             T_in_coolantrange=T_in-5.0
   !                             POWER_coolantrange=POWER
   !                             !write(*,'(A,F7.2,A,F5.2,A)') & 
   !                             !"For coolant LSSS start T_in at ", T_in_coolantrange, " and the power level at ", POWER/1.0E6, " MW" 
   !                             coolantflag=.true.
   !                         end if
   !                         if (T_core_max <= T_fuel_limit .AND. fuelflag==.false.) then
   !                             !write(*,*) T_in,  T_out_avg
   !                             T_in_fuelrange=T_in-7.0
   !                             POWER_fuelrange=POWER
   !                             !write(*,'(A,F7.2,A,F5.2,A)') & 
   !                             !"For fuel LSSS start T_in at ", T_in_fuelrange, " and the power level at ", POWER/1.0E6, " MW" 
   !                             !write(*,*) T_core_max
   !                             fuelflag=.true.
   !                         end if
   !                     end if
   !                 end if
   !             end do
   !         end if
   !     end do
   !     !================================================================================
   !     ! Find power level for LSSS limits for given T_out and mass flow rate :: use range found for T_in and POWER
   !     !================================================================================
   !     ! Reset flags to .false.
   !     fuelflag=.false. 
   !     coolantflag=.false.
   !     Toutflag=.false.
   !     !write(*,*)
   !     !write(*,'(A,F5.2,A,F7.2,A)') "For primary mass flow rate of ", W_core, " and outlet temperature of ", T_out_limit, ":"
   !     do T_in=T_in_coolantrange,T_out_limit,0.01_8 
   !         ! Check if all LSSS limits met
   !         if (Toutflag==.true. .AND. coolantflag==.true.) then
   !             !write(*,*) "LSSS Coolant Limit Found"
   !             EXIT
   !         else          
   !             ! Loop through power level to find T_out=T_out_limit
   !             ! Reset T_out flags to .false.
   !             Toutflag=.false.
   !             do POWER=POWER_coolantrange,50.0E6,0.01E6_8
   !                 if (Toutflag==.false.) then
   !                     ! Calculate core temperatures for average channel
   !                     call pebblebedcore(POWER,W_core,T_in,T_out,T_coolant_max,T_core_max,1)
   !                     T_out_avg=T_out
   !                     ! Check for exceedance of LSSS T_out temperature limit
   !                     if (T_out>=T_out_limit .AND. Toutflag==.false.) then
   !                         Toutflag=.true.
   !                         ! Calculate core temperatures for hot channel
   !                         call pebblebedcore(POWER,W_core,T_in,T_out,T_coolant_max,T_core_max,2)
   !                         if (T_coolant_max <= T_coolant_limit .AND. coolantflag==.false.) then
   !                             write(*,'(A,F7.2,A)') "The maximum   coolant hot channel  temperature exceeds ", T_coolant_limit, " at:"
   !                             write(*,'(F5.2,A,F7.2,A,F7.2,A,F7.2)') POWER/1.0E6, " MW    T_in = ", T_in, "   T_out = ", T_out_avg, "   TCMAX = ", T_coolant_max
   !                             write(*,*)
   !                             coolantflag=.true.
   !                         end if
   !                     end if
   !                 end if
   !             end do
   !         end if
   !     end do
   !     !================================================================================
   !     ! Find power level for LSSS limits for given T_out and mass flow rate 
   !     !================================================================================
   !     ! Reset flags to .false.
   !     fuelflag=.false. 
   !     !coolantflag=.false.
   !     Toutflag=.false.
   !     !write(*,*)
   !     !write(*,'(A,F5.2,A,F7.2,A)') "For primary mass flow rate of ", W_core, " and outlet temperature of ", T_out_limit, ":"
   !     do T_in=T_in_fuelrange,T_out_limit,0.01_8
   !         ! Check if all LSSS limits met
   !         if (Toutflag==.true. .AND. fuelflag==.true.) then
   !             !write(*,*) "LSSS Fuel Limit Found"
   !             EXIT
   !         else          
   !             ! Loop through power level to find T_out=T_out_limit
   !             ! Reset flags to .false.
   !             Toutflag=.false.
   !             do POWER=POWER_fuelrange,50.0E6,0.01E6
   !                 if (Toutflag==.false.) then
   !                     ! Calculate core temperatures for average channel
   !                     call pebblebedcore(POWER,W_core,T_in,T_out,T_coolant_max,T_core_max,1)
   !                     T_out_avg=T_out
   !                     ! Check for exceedance of LSSS T_out temperature limit
   !                     if (T_out>=T_out_limit .AND. Toutflag==.false.) then
   !                         Toutflag=.true.
   !                         ! Calculate core temperatures for hot channel
   !                         call pebblebedcore(POWER,W_core,T_in,T_out,T_coolant_max,T_core_max,2)
   !                         if (T_core_max <= T_fuel_limit .AND. fuelflag==.false.) then
   !                             write(*,'(A,F7.2,A)') "The maximum    fuel hot channel    temperature exceeds ", T_fuel_limit, " at:"
   !                             write(*,'(F5.2,A,F7.2,A,F7.2,A,F7.2)') POWER/1.0E6, " MW    T_in = ", T_in, "   T_out = ", T_out_avg, "   TFMAX = ", T_core_max 
   !                             write(*,*)
   !                             fuelflag=.true.
   !                         end if
   !                     end if
   !                 end if
   !             end do
   !         end if
   !     end do
   !     
   ! END SUBROUTINE prismaticLSSSloop
    
    
    
    
    
            
            
    !================================================================================
    !  FUNCTION: CALCULATE HEAT TRANSFER COEFFICIENT
    !================================================================================
    ! INPUTS
    !   W       :: MASS FLOW RATE [kg/s]
    !   DC      :: CHANNEL DAIMETER [m]
    !   DP      :: PEBBLE DAIMETER [m]
    !   T       :: COOLANT TEMPERATURE [Celcius]
    !   N       :: LOGICAL TRIP
    !       N=1 Wakao Correlation 
    !       N=2 Kunii and Levenspiel Correlation
    !       n=3 Dittus-Boelter
    ! OUTPUT    
    !   HTC  ::  Heat transfer coefficient of the pebble bed [W/m^2-K]            
    ! REFERENCE
    !
    !================================================================================
    REAL(8) FUNCTION HTC(W,DC,T,N)
        IMPLICIT NONE
        REAL(8), INTENT(IN)     :: W    ! Mass flow rate [kg/s]
        REAL(8), INTENT(IN)     :: DC   ! Channel diameter [m]
        REAL(8), INTENT(IN)     :: T    ! Temperature [Celcius]
        INTEGER, INTENT(IN)     :: N    ! Correlation Selector
        REAL(8)                 :: G    ! Mass flux = mass flow rate / area
        REAL(8)                 :: Re   ! Reynolds number
        REAL(8)                 :: Pr   ! Prandtl number
        REAL(8)                 :: Nu   ! Nusselt number
                
        ! AVERAGE FLUID PROPERTIES
        G=W/(PI*(DC/2.0)**2.0)          ! [kg/m^2-s]
        Re=G*DC/flibe_mu(T)
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
            ! Kunii and Levenspiel
            Nu=2.0+1.8*Pr**(1.0/3.0)*Re**0.5 
        CASE(3)
            ! Dittus-Boelter (heated)
            if (Re<10000) then 
                write(*,*) "Reynold's number out of range for selected correlation."
            else if (Pr<0.7 .OR. Pr>100) then
                write(*,*) "Prandtl number out of range for selected correlation."
            end if
            Nu=0.023*Re**0.8*Pr**0.4
        END SELECT

        HTC=Nu*flibe_k(T)/DC
                
    END FUNCTION HTC
            
            
    !================================================================================
    !  FUNCTION: CALCULATE FRICTION FACTOR IN PEBBLE BED
    !================================================================================
    ! INPUTS
    !   W       :: MASS FLOW RATE [kg/s]
    !   DC      :: CHANNEL DAIMETER [m]
    !   T       :: COOLANT TEMPERATURE [Celcius]
    !   N       :: LOGICAL TRIP
    !       N=1 pipe flow
    !       N=2 churchill
    ! OUTPUT    
    !   FF  ::  friction factor [none]           
    ! REFERENCE
    !
    !================================================================================
    REAL(8) FUNCTION FF(W,DC,T,N)
        IMPLICIT NONE
        REAL(8), INTENT(IN)     :: W    ! Mass flow rate [kg/s]
        REAL(8), INTENT(IN)     :: DC   ! Channel diameter [m]
        REAL(8), INTENT(IN)     :: T    ! Temperature [?]
        INTEGER, INTENT(IN)     :: N    ! Correlation Selector
        REAL(8)                 :: G    ! Mass flux = mass flow rate / area
        REAL(8)                 :: Re   ! Reynolds number
        REAL(8)                 :: A
        REAL(8)                 :: B
        REAL(8)                 :: e    ! Surface roughness
        

        ! AVERAGE FLUID PROPERTIES
        G=W/(PI*(DC/2.0)**2.0)      ! [kg/m^2-s]
        Re=G*DC/flibe_mu(T)

        ! FRICTION FACTOR CORRELATION
        SELECT CASE(N)
        CASE(1)
            ! Laminar and turblent pipe flow
            if (Re<=4000) then
                FF=64.0/Re
            else if (Re>4000 .AND. Re<30000) then ! Blasius 
                FF=0.316*Re**(-0.25)
            else if (Re>30000 .AND. Re<1000000) then ! McAdams
                FF=0.184*Re**(-0.2)
            else
                write(*,*) "Reynold's number out of range for selected correlation."
            end if             
        CASE(2)
            ! Churchill
            if (Re<2100) then
                FF=64.0/Re
            else
                A=(-2.457*log( (7.0/Re)**0.9 + 0.27*e/DC ) )**16.0
                B=(37530.0/Re)**16.0
                e=30*10**(-6.0) ! [m] surface roughness of HastN from weighted average of Ni and Mo
                FF=8.0*( (8.0/Re)**12.0 + 1.0/((A+B)**1.5) )**(1.0/12.0)
            end if
        END SELECT

    END FUNCTION FF
            
END MODULE prismatic