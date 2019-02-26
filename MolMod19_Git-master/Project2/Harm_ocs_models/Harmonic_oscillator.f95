Program harmonic_oscillator

  IMPLICIT NONE
  
  ! parameters of simulation, there are two particles interacting through harmonic bond
  ! simulation is carried out with a set timestep length for a set number of timesteps
  
  INTEGER(4),PARAMETER :: number_particles = 2
  INTEGER(4),PARAMETER :: number_timesteps = 1000
  REAL(8),PARAMETER :: timestep = 0.01
  REAL(8),PARAMETER :: bond_length = 1D0
  REAL(8),PARAMETER :: bond_potential = 1D0
  
  INTEGER(4) :: I
  
  ! the x and y axis positions, velocities and forces 
  
  REAL(8) :: xpos(number_particles),ypos(number_particles)
  REAL(8) :: xvel(number_particles),yvel(number_particles)
  REAL(8) :: xforce(number_particles),yforce(number_particles)
  
  ! initialize position of both particles with distance between set to stretch bond
  
  xpos(1) = -7.5D-1
  xpos(2) =  7.55D-1
  ypos = 0D0
  
  ! initialize forces (I chose zero, you can start them random)
  
  xvel = 0D0
  yvel = 0D0
  
  CALL calculate_forces(xforce,yforce,xpos,ypos)
  
  DO I = 1,number_timesteps
  
  ! performs velocity verlet algorithm, step 1 (see notes)
  
    xvel = xvel + xforce*5D-1*timestep
    yvel = yvel + yforce*5D-1*timestep
    
  ! performs velocity verlet algorithm, step 2 (see notes)

    xpos = xpos + xvel*timestep
    ypos = ypos + yvel*timestep

  ! calculates forces
  
    CALL calculate_forces(xforce,yforce,xpos,ypos)
    
  ! performs velocity verlet algorithm, step 3 (see notes)

    xvel = xvel + xforce*5D-1*timestep
    yvel = yvel + yforce*5D-1*timestep
        
  ENDDO    
  
  ! subroutine that calculates the forces as a function of the positions of the particles
  ! enters the particles positions and outputs the forces on the particles
  
  CONTAINS
  
  SUBROUTINE calculate_forces(xf,yf,xp,yp)
  
    IMPLICIT NONE
    REAL(8),INTENT(IN) :: xp(number_particles),yp(number_particles)
    REAL(8),INTENT(OUT):: xf(number_particles),yf(number_particles)
    REAL(8) :: distance,xvec,yvec,xhat,yhat,scalar_force
    INTEGER(4) :: I
    
    
    DO I = 1,number_particles,2
    
    ! calculates the directional vector and distance between particles
    
      xvec = xp(I+1) - xp(I)
      yvec = yp(I+1) - yp(I)
      
      distance = DSQRT(xvec*xvec + yvec*yvec)
      
      xhat = xvec/distance
      yhat = yvec/distance
      
    ! calculates scalar value of force due to harmonic bond of length bond_length and 
    ! bond potential constant bond_potential
        
      scalar_force = 2D0*bond_potential*(distance - bond_length)
      
    ! calculate the vector components of the force
    
      xf(I)   =  xhat*scalar_force  
      yf(I)   =  yhat*scalar_force  
      xf(I+1) = -xhat*scalar_force  
      yf(I+1) = -yhat*scalar_force  
      
    ENDDO    
  
  END SUBROUTINE calculate_forces
  
END PROGRAM harmonic_oscillator