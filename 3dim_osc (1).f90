program springparticles
!use mtmod 
!use rng
implicit none

real :: k, r_0, mass_p1, mass_p2, t, dt, r, E_kin_p1,E_kin_p2,E_kin_tot=0.0,E_pot=0.0,E_tot=0.0
        
real,dimension(3) :: position_p1, position_p2, velocity_p1, velocity_p2, acceleration_p1,acceleration_p2, &
                    lin_mom_p1,lin_mom_p2,ang_mom_p1,ang_mom_p2,center,total_lin_mom,total_ang_mom
integer :: n

!!!!! Get user input !!!!!

print *,"Please give in the following order x and y coordinate for: Position(particle 1), position(particle 2)"
read *,position_p1(1),position_p1(2),position_p1(3),position_p2(1),position_p2(2),position_p2(3)

!call sgrnd(1223)
velocity_p1(1) =rand(0)
velocity_p1(2) =rand(0)
velocity_p1(3) =rand(0)
velocity_p2(1) =rand(0)
velocity_p2(2) =rand(0)
velocity_p2(3) =rand(0)

!!!!! Initialize !!!!!

mass_p1 = 1.0
mass_p2 = 1.0
r_0 = 1.0
k = 1.0
t= 0
dt = 0.01

r       = sqrt((position_p1(1)-position_p2(1))**2+(position_p1(2)-position_p2(2))**2+(position_p1(3)-position_p2(3))**2)
center  = position_p1 + 0.5*(position_p2-position_p1)

!!!!! Prepare a file !!!!!

call execute_command_line('rm distance_2dim.xyz')
open(unit=1,file='distance_2dim.xyz')
    write(1,'(F0.5,x,F0.5,x,"0")',advance="no") t
    write(1,'(F0.5,x,F0.5,x,"0")',advance="no") r
    write(1,'(F0.5,x,F0.5,x,"0")',advance="no") r-r_0! This gives my perfect cosine 
close(unit=1,status='keep')


call execute_command_line('rm Trajectory.xyz')      ! This is my VMD readable file
open(unit=1,file='Trajectory.xyz')
    write(1,*) "2"
    write(1,*) "Time: ",t,n
    write(1,'("N",x,F0.5,x,F0.5,x,F0.5)') position_p1(1),position_p1(2),position_p1(3)
    write(1,'("O",x,F0.5,x,F0.5,x,F0.5)') position_p2(1),position_p2(2),position_p2(3)
close(unit=1,status='keep')

 ! Calculate momenta                        ! For simplicity the total momenta are arrays where both values are the same, and only one is relevant
    lin_mom_p1      =  (velocity_p1)*mass_p1 
    lin_mom_p2      =  (velocity_p2)*mass_p2 
    total_lin_mom   = lin_mom_p1 + lin_mom_p2
    
    ang_mom_p1      = ((position_p1(1)-center(1))*velocity_p1(2)-(position_p1(2)-center(2))*velocity_p1(1))*mass_p1
    ang_mom_p2      =  ((position_p2(1)-center(1))*velocity_p2(2)-(position_p2(2)-center(2))*velocity_p2(1))*mass_p2
    total_ang_mom   = ang_mom_p1 + ang_mom_p2

!print*,total_ang_mom(1),total_lin_mom(1)

! Calculate energies
    E_kin_p1    = 0.5*mass_p1* sqrt(velocity_p1(1)**2+velocity_p1(2)**2+velocity_p1(3)**2)**2
    E_kin_p2    = 0.5*mass_p2* sqrt(velocity_p2(1)**2+velocity_p2(2)**2+velocity_p2(3)**2)**2
    E_kin_tot   = (E_kin_p1 + E_kin_p2) 
    E_pot       = (r-1)**2
    E_tot       = E_kin_tot + E_pot
    print *,E_kin_p1,E_kin_p2,E_kin_tot,E_pot,E_tot
    
!!!!! Time evolution starts !!!!!

do n = 1,1000

    t = t + dt
    
    ! Calculate acceleration
    acceleration_p1 = -k*2*(r-1)* ((position_p1 - position_p2)/r) /mass_p1
    acceleration_p2 = -k*2*(r-1)* ((position_p2 - position_p1)/r) /mass_p2


    ! Calculate new velocity
    velocity_p1 = velocity_p1+acceleration_p1*dt
    velocity_p2 = velocity_p2+acceleration_p2*dt
    
    ! Calculate new position
    position_p1 = position_p1+velocity_p1*dt
    position_p2 = position_p2+velocity_p2*dt

    ! Calculate distance

    r       = sqrt((position_p1(1)-position_p2(1))**2+(position_p1(2)-position_p2(2))**2+(position_p1(3)-position_p2(3))**2)
    center  = position_p1 + 0.5*(position_p2-position_p1)
    
    ! WRITE FILE
     open(unit=1,file='distance_2dim.xyz',status='old',access='append')
         write(1,'(F0.5,x,F0.5,x,"0")',advance="no") t
         write(1,'(F0.5,x,F0.5,x,"0")',advance="no") r ! 
         write(1,'(F0.5,x,F0.5,x,"0")',advance="no") r-r_0 ! 
     close(unit=1,status='keep')
 
     open(unit=1,file='Trajectory.xyz',status='old',access='append')
         write(1,*) "2"
         write(1,*) "Time: ",t,n
         write(1,'("N",x,F0.5,x,F0.5,x,F0.5)') position_p1(1),position_p1(2),position_p1(3)
         write(1,'("O",x,F0.5,x,F0.5,x,F0.5)') position_p2(1),position_p2(2),position_p2(3)
     close(unit=1,status='keep')

 ! Calculate momenta

    lin_mom_p1      =  velocity_p1*mass_p1 
    lin_mom_p2      =  velocity_p2*mass_p2 
    total_lin_mom   = lin_mom_p1 + lin_mom_p2
    
    ang_mom_p1      = ((position_p1(1)-center(1))*velocity_p1(2)-(position_p1(2)-center(2))*velocity_p1(1))*mass_p1
    ang_mom_p2      =  ((position_p2(1)-center(1))*velocity_p2(2)-(position_p2(2)-center(2))*velocity_p2(1))*mass_p2
    total_ang_mom   = ang_mom_p1 + ang_mom_p2
!print*,total_ang_mom(1),total_lin_mom(1)

! Calculate energies
    E_kin_p1    = 0.5*mass_p1* sqrt((velocity_p1(1))**2+(velocity_p1(2))**2+velocity_p1(3)**2)**2
    E_kin_p2    = 0.5*mass_p2* sqrt((velocity_p2(1))**2+(velocity_p2(2))**2+velocity_p2(3)**2)**2
    E_kin_tot   = (E_kin_p1 + E_kin_p2)
    E_pot       = (r-1)**2 
    E_tot       = E_pot + E_kin_tot 

print *,E_kin_p1,E_kin_p2,E_kin_tot,E_pot,E_tot

end do

end program springparticles

