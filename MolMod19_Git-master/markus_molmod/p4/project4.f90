module defs
  !-----------------------Global Variables--------------------------!
  !Argon parameters

  integer, parameter:: S =50000          	!Time steps
  integer, parameter:: skip = 500          	!Don't record all configurations
  real*8, parameter :: sigma = 3.4d-10
  real*8, parameter :: epsilon1 = 1.65d-21
  real*8, parameter :: rho = 0.001633d0 	!density
  real*8, parameter :: m = 39.948d0	      	!Argon mass
  real*8, parameter :: dt=1d-8			!Timestep
  real*8, parameter :: gamma1 = 0.5d0 		!Langevin parameters
  real*8, parameter :: Lside = 40*3.4d-10	!!Cube length
   real*8, parameter :: Kb=1.38064852d-23	!Boltmann constant
  !character(len=5) :: fname


  !-----------------------Global Variables--------------------------!

end module defs



program test
  integer :: N, k
  real*8 :: t, mdstime,mcstime, scells, temp


  N=20; temp=10

  !call mcsimulation(N,mcstime,temp)
  call  mdsimulation(N,mdstime,temp)
  write(*,*)" N ", "  Temp    ", "  mdstime  ", "  mcstime "
  write(*,'(i4,f8.3,f12.6, f12.6)') N, temp, mdstime, mcstime


end program test

subroutine mcsimulation(N,stime,temp)
use defs

!!!!!!!!!!!!!!!!!!!!!!Algorithm

  !Establish an initial configuration. consider only the positions of the particles. The energy of the initial configuration is not important.
  !Choose a particle at random and make a trial change in its position.
  !Compute ΔE, the change in the potential energy of the system due to the trial move.
  !If ΔE is less than or equal to zero, accept the new configuration and go to step 8.
  !  If ΔE is positive, compute the quantity w = e-βΔE.
  !  Generate a uniform random number r in the unit interval [0,1].
  !  If r ≤ w, accept the trial move; otherwise retain the previous microstate.
  !Determine the value of the desired physical quantities.
  !Repeat steps (2) through (8) to obtain a sufficient number of microstates.
  !Accumulate data for various averages.



   integer :: pindex, i, j, k, N, step, dskip, ix
  real*8, dimension(:,:), allocatable  :: q, v, F
  real*8, dimension(:), allocatable :: x1, x2, r1, r2, r3, ri,rf
  real*8, intent(out) :: stime
  real*8, parameter :: cutoff = 2.5*sigma
  real*8 :: Ekin, Epot, diff, start, finish, temp,kT
  real*8 :: w,rn, dy, dz, ntry, nacp,aveep, aR
  real*8 :: particle, newEpot, deltaEpot, totEp, pairEp
  real*8 :: diffusion, rmsquare,  oldcrd(3), newcrd(3), dx(3)
  character(len=30) :: fname2, fn, ft
  !-----------------------Global Variables--------------------------!


  allocate(q(N, 3), v(N, 3), F(N, 3))
  allocate(x1(3), x2(3), r1(N), r2(N), r3(N), ri(N), rf(N))

  n_q =0;  totep = 0.0;  kT= temp;
  dskip =int(skip)

  write (ft,'(F6.2)')temp
  write (fn,'(i4.4)')N
  fname2 =trim(ft)//'_'//trim(fn)
  fname2 =trim(fname2)


  call random_number(q)   ! uniform random number in [0,1)
  q= q*Lside-Lside;

  do ix=1, N
     ri(ix)= sqrt(q(ix,1)**2 + q(ix,2)**2 + q(ix,3)**2)
  end do

  call resetPE()

  call cpu_time(start)
  open(unit = 1,file = 'mcoutput'//trim(fname2)//'.xyz', action = 'write')
  open(unit = 2,file = 'mcenergy'//trim(fname2)//'.dat', action = 'write')
  !open(unit = 4,file = 'mcdiffusion'//trim(fname2)//'.dat', action = 'write')

  do step=1, S

     !Choose a particle at random and make a trial change in its position (dx, dy, dz)

     call random_number(particle)
     pindex= mod(int(real(N)*particle), N) +1; !pick index of position in {0-N}
     call random_number(dx)


     do k=1,3

        newcrd(k)=q(pindex,k) + dx(k)*Lside -Lside

        ! periodic boundary condition
        if(newcrd(k) < 0)    newcrd(k) = newcrd(k)+ Lside
        if(newcrd(k) > Lside) newcrd(k) = newcrd(k) - Lside
        q(pindex,k)=newcrd(k)
     end do

     oldcrd(1:3) = q(pindex,1:3)

     !Compute ΔE, the change in the potential energy of the system due to the trial move.
     !change the energy measurement
     deltaEpot =0.0

     do i =1, N
        if(i == pindex) cycle
        x1(1:3)= q(i,1:3)

        call calc_force (newcrd, x1, Epot)
        deltaEpot = deltaEpot + Epot

        call calc_force (oldcrd, x1, Epot)
        deltaEpot = deltaEpot - Epot
     end do


     !If ΔE is less than or equal to zero,
     !accept the new configuration skip the negative if condition

     if(deltaEpot>0)then
        ! If ΔE is positive, compute the quantity w = e-βΔE.
        w = exp(-deltaEpot/kT);

        ! Generate a uniform random number r in the unit interval [0,1].
        call random_number(rn)

        ! If r ≤ w, accept the trial move; otherwise retain the previous microstate.
        !reverse position if conditions are not met

        if(rn>w)then
           do k=1,3
              q(pindex,k)= q(pindex,k) -r1(k);
           end do
        else
           nacp = nacp + 1.0
        end if
     end if

     !Adjust potential if criterion met
     !calculate everything here
     q(pindex,k)= q(pindex,k) + r1(k);
     totep = totep + deltaEpot;
     ntry = ntry + 1.0
     aveep = aveep + totep

     if (modulo(step,skip)==0) then
        aR = Real(nacp)/real(ntry) !number of success vs number of tries
        !average energy and acceptance ratio
        !average step and potential energy per prticle ==totep/N

        write(2, *) step/real(N), totep/real(N), aR
     end if

     !produce the xyz file
     if (modulo(step,skip)==0) then

        write(1,*) N
        write(1,*)

        do i=1, N
           write(1, *) "Ar",q(i,1),q(i,2) ,q(i,3)
        end do

     end if

     if(mod(step, dskip) == 0) then
        !use less steps to get a proper line
        diff =0.0

        do ix=1, N
           rf(ix)= sqrt(q(ix,1)**2 + q(ix,2)**2 + q(ix,3)**2)
        end do

        !populate diffusion
        !plot diffAverage vs time
        !diff = diffusion(N, ri, rf,step);
        !write(4,*) step, diff

          do ix=1, N
           ri(ix)=  rf(ix)
        end do

     end if

  end do

  close(1); close(2);   close(4)
!!!!!!!!!!!!!!!!!!!!!!call radial

  call cpu_time(finish)
  stime = finish - start
  call radial(N, Lside, q, 'mc',temp);

  deallocate( q, v, F, x1, x2, r1, r2, r3, ri,rf)


contains
  subroutine calc_force (dz, dy, Epot)           ! total potential energy calculated

    implicit none
    !!N, q, totEp, aveEp, ntry, nacp
    integer :: i, j, N
    real*8 :: dx(3), dy(3), dz(3), rr
    real*8 :: Epot, totEp, aveEp

    dx(:)= dz(:)- dy(:)
    dx(1:3) = Lside*0.5 - abs(Lside*0.5 - abs(dx(1:3)))
    rr = sqrt(sum(dx(1:3) ** 2))

    Epot= 0.0! calculate pairppotential
    if(rr < cutoff) then
       Epot = Epot+ 4 * epsilon1 *((sigma/rr)**12-(sigma/rr)**6)
    endif

    totEp = totEp + Epot
  end subroutine calc_force

  subroutine resetPE()
    implicit none
    integer :: i, j
    real*8 :: x1(3), x2(3), Epot
    totep = 0.0

    do i = 1, N-1
       x1(1:3) = q(i,1:3)
       do j = i + 1, N
          x2(1:3) = q(j,1:3)
          call calc_force(x1, x2, Epot)
          totep = totep + Epot
       end do
    end do
    aveep = 0.0; ntry = 0.0; nacp = 0.0
  end subroutine resetPE

end subroutine mcsimulation




subroutine mdsimulation(N,stime,temp)
  use defs
  implicit none

  type node
     integer :: i
     type(node), pointer :: next => null()
  end type node

   !-----------------------Global Variables--------------------------!
  real*8, intent(out) :: stime
  integer :: subcells, N, cell, i, j, step,ix, dskip
  integer, dimension(:), allocatable :: head
  real*8, dimension(:),   allocatable  :: ri,rf
  real*8, dimension(:,:), allocatable :: q, v, F, v0, q_abs,  q0         !q(x,y,z), v(vx, vy, vz), Fx, Fy, Fz
  type(node), dimension(:), allocatable :: linkedList
  !-----------------------Global Variables--------------------------!
  real*8, parameter :: cutoff = 2.5*sigma
  real*8 :: start, finish,Ekin, Epot, diff, manf,temp
  real*8 :: diffusion, rmsquare,kT
  character(len=30) :: fname, fn, ft

  subcells= floor(Lside/cutoff)

  !allocate(molecule(N))
  allocate(q(N,3), v(N,3), F(N,3), q0(N,3), q_abs(N,3),ri(N), rf(N))
  allocate(head(subcells**3))
  allocate(linkedList(N))

  !call init_random_seed()
  call random_seed()
  call random_number(q)   !random values for positions [0,1]
  q = q*Lside -Lside
  q0 = q
  q_abs = q
  v = 0.0d0               !Set initial velocities to zero.

  do i=1, N
     ri(i)= sqrt(q(i,1)**2 + q(i,2)**2 + q(i,3)**2)
  end do

  kT= Kb *temp;
  dskip =int(skip)
  write (ft,'(F6.2)')temp
  write (fn,'(i4.4)')N
  fname =trim(ft)//'_'//trim(fn)
  fname =trim(fname)

  call cpu_time(start)

  open(unit = 1,file = 'mdoutput'//trim(fname)//'.xyz', action = 'write')
  open(unit = 2,file = 'mdenergy'//trim(fname)//'.dat', action = 'write')
  open(unit = 4,file = 'mddiffusion'//trim(fname)//'.dat', action = 'write')

  do step = 1, S
     Ekin = 0.0d0
     Epot = 0.0d0
     if(mod(step, skip) == 0) then       !This to avoid saving all configurations.
        write(1,*) N
        write(1,*)
     end if
     head = -1                   !head = -1 means empty cell.

     do i = 1, N
        do j = 1, 3
           !Check that a molecule doesn't move over L in one time step.
           if(abs(q(i,j)) > 2*Lside) then
              write(*,*) 'Use smaller time-step', step, q(i,j) !This just in case a molecule moves over L distance in a single time step.
              stop
           end if
           !Periodic boundary conditions
           if(q(i,j) < 0) then
              q(i,j) = q(i,j) + Lside
           else if(q(i,j) >= Lside) then
              q(i,j) = q(i,j) - Lside
           end if
        end do
        !Check in which subcell molecule i is.
        cell = (floor(q(i,1)/cutoff))*subcells**2 + (floor(q(i,2)/cutoff))*subcells + (floor(q(i,3)/cutoff)) + 1

        !Construct linked list
        linkedList(i)%i = head(cell)
        head(cell) = i


        Ekin = Ekin + 0.5d0*m*norm2(v(i,:))**2
        if(mod(step, skip) == 0) then
           write(1,*) 'Ar', q(i,1), q(i,2), q(i,3)
        end if
     end do

     call calcForces

     if(mod(step, skip) == 0) then
        write(2,*) step, Epot, Ekin, Ekin+Epot

     end if


     if(mod(step, dskip) == 0) then
        !use less steps to get a proper line
        diff =0.0

        do ix=1, N
           rf(ix)= sqrt(q(ix,1)**2 + q(ix,2)**2 + q(ix,3)**2)
        end do

        !populate diffusion
        !plot diffAverage vs time
        diff = diffusion(N, q_abs, q0);
        write(4,*) step, diff

     end if

     !Velocity Verlet
     v = v + F/(2*m)*dt        !v(t + dt/2)
     q_abs = q_abs + v*dt
     q = q + v*dt              !q(t + dt)
     call calcForces
     v = v + F/(2*m)*dt        !v(t + dt)
  end do
  close(4);    close(2) ;  close(1)


  call cpu_time(finish)
  call radial(N, Lside, q, 'md',temp);
  stime =finish-start;

  deallocate(q, v, F, ri, rf, head)
contains

  subroutine calcForces
    F = 0.0d0
    !Get random impulses.
    call random_seed()
    !call init_random_seed()
    call random_number(F)
    F = F - 0.5d0
    F = F*dsqrt(2.0d0*gamma1*m*temp/step)

    !Account for viscious forces
    F = F - gamma1*v*m

    call getInterPartForce
  end subroutine calcForces

  subroutine getInterPartForce
    integer i, j, k, l, m, n, cell1, cell2, x, y, z, particle1, particle2
    integer, dimension(3) :: shift
    real*8, dimension(3) :: dq
    real*8 :: r

    do i = 1, subcells
       do j = 1, subcells
          do k = 1, subcells
             cell1 = ((i-1)*subcells**2) + (j-1)*subcells + k
             do l = -1, 1
                do m = -1, 1
                   do n = -1, 1
                      x = i+l
                      y = j+m
                      z = k+n
                      shift = 0
                      if(x < 1) then              !Periodic boundary for subcells at edges
                         x = subcells
                         shift(1) = Lside
                      end if
                      if(y < 1) then
                         y = subcells
                         shift(2) = Lside
                      end if
                      if(z < 1) then
                         z = subcells
                         shift(3) = Lside
                      end if
                      if(x > subcells) then
                         x = 1
                         shift(1) = -Lside
                      end if
                      if(y > subcells) then
                         y = 1
                         shift(2) = -Lside
                      end if
                      if(z > subcells) then
                         z = 1
                         shift(3) = -Lside
                      end if
                      cell2 = ((x-1)*subcells**2) + (y-1)*subcells + z
                      particle1 = head(cell1)
                      do while(particle1 /= -1)
                         particle2 = head(cell2)
                         do while(particle2 /= -1)
                            if(particle1 < particle2) then                         !Don't double count pair potential.
                               dq = q(particle1,:) - (q(particle2,:)+shift)
                               r = sqrt(dq(1)**2 + dq(2)**2 + dq(3)**2)
                               if(r < cutoff) then
                                  F(particle1,:) = F(particle1,:) &
                                       - 4*dq(:)*epsilon1*(6*sigma**6 *r**(-8)-12*sigma**(12)*r**(-14))
                                  F(particle2,:) = F(particle2,:) - F(particle1,:)
                                  Epot = Epot + 4 * epsilon1 *((sigma/r)**12-(sigma/r)**6)
                               end if
                            end if
                            particle2 = linkedList(particle2)%i
                         end do
                         particle1 = linkedList(particle1)%i
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine getInterPartForce


end subroutine mdsimulation


real*8 function rmsquare(N, rs)
  !Root means square

  real*8 :: avr, rms
  integer, intent(in):: N
  real*8, intent(in) :: rs(N)

  avr =0; rms=0;

  do i=1,N
     avr = avr + rs(i);
  end do


  avr = avr/N;

  do i=1, N
     rms= rms +(rs(i)-avr)*(rs(i)-avr);
  end do
  rmsquare= rms/N;

end function rmsquare


real*8 function diffusion(N, q_abs, q0)
  !calculate diffusion
  real*8 :: sumr, avg
  integer, intent(in):: N
  real*8, intent(in) :: q_abs(N,3), q0(N,3);
  integer ::i

  sumr=0;
  do i = 1, N
    sumr = sumr + norm2(q_abs(i,:)-q0(i,:))**2;
  end do
  !sumr = sum((ri-rf)**2);
  diffusion = sumr/N; ! return the average
end function diffusion


!============================================================================
real*8 function DiffCoef(vf,vi,N,dt)
  !----------------------------------------------------------------------------
  ! Evaluates the diffusion coefficient from the Green-Kubo formula:
  !             1  inf  N
  ! D = lim    --- Int <Sum v(t)v(0)> dt
  !     t->inf 3N   0   i=1
  !----------------------------------------------------------------------------

  implicit none
  real*8 :: sumr, avg, diff, dt, fact
  integer, intent(in):: N
  real*8, intent(in) :: vi(N,3), vf(N,3);
  integer ::i,steps
  logical init /.true./

  if (init) then
     Diff = 0.d0
     fact = dt/(3.d0*N)
     init = .false.
  end if

  do i = 1,N
     Diff = Diff + fact*(vf(i,1)*vi(i,1) +vf(i,2)*vi(i,2) + vf(i,3)*vi(i,3))
  end do

  DiffCoef = Diff
end function DiffCoef


subroutine radial(N, Lside, r, mcd,temp);
use defs, only: rho, sigma
  integer, parameter :: bins=200;
   real*8, parameter :: cutoff = 2.5*sigma
  real*8, parameter :: pi = 3.141592653
  real*8 ::aa, volume, g(2*bins), dh(3), r1, Lside, temp;
  integer :: i,j, k, ka;
  integer, intent(in):: N
  real*8, intent(in) :: r(N,3);
  character(len=2) :: mcd!indicate if mc or md
 character(len=30) :: fname, fn, ft



  aa=0.4*Lside/bins;
  volume = 4*pi/3.0;
  g=0;

  write (ft,'(F6.2)')temp
  write (fn,'(i4.4)')N
  fname = trim(mcd)//'radial'//trim(ft)//'_'//trim(fn)//'.dat'
  fname = trim(fname)


  do i=1,N-1
     do j=1+i,N

        do k=1, 3
           dh(k) = r(i,k) - r(j,k);
           if (dh(k)>Lside/2)then
              dh(k)= dh(k)-Lside;
           else if (dh(k)<-Lside/2)then
              dh(k)=dh(k)+Lside;
           end if
        end do
        r1 = sqrt((dh(3)*dh(3))+(dh(1)*dh(1))+(dh(2)*dh(2)));

        if(r1 < 0.4*Lside) then
           ka = int(r1/aa)
           g(ka) = g(ka)+ 2.0;
        end if

     end do
  end do


  !write to file normalized radial functions
  open(unit = 44,file = fname, action = 'write')

  do i=1,bins
     g(i)=g(i)/( volume*((i+1)**3-i**3)*aa**3*rho);
     write(44,*) (real(i)+0.5)*aa, g(i)/real(N);
  end do

  close(44);
end subroutine radial



SUBROUTINE init_random_seed()
  !borrowed : https://gcc.gnu.org/onlinedocs/gcc-4.6.2/gfortran/RANDOM_005fSEED.html
  integer :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed

  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))

  CALL SYSTEM_CLOCK(COUNT=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  DEALLOCATE(seed)
END SUBROUTINE init_random_seed
