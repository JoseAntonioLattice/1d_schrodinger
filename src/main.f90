program main

  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

  integer(i4), parameter :: n = 2000
  real(dp), parameter :: m = 1.0_dp
  real(dp), parameter :: a = 1.0_dp
  real(dp), parameter :: dx = a / (n + 1)
  real(dp), parameter :: lambda = 1/(2*m*dx**2)
  
  real(dp), dimension(0:N+1) :: x
  
  real(dp), dimension(0:n+1) :: Psi
  real(dp), dimension(n,n)   :: Hamiltonian
  real(dp), dimension(n)     :: V_potential

  ! Energy
  real(dp), dimension(n) :: Ei, Er

  !Lapack parameters
  integer, parameter :: lwork = 13000
  real(dp), dimension(lwork) :: work
  real(dp), dimension(n,n)   :: vl, vr
  integer :: info
  
  integer(i4) :: i
  
  ! Infinite Square Well boundary conditions
  Psi(0)   = 0.0_dp
  Psi(N+1) = 0.0_dp

  ! Pontential inside the Infinite Well
  V_potential = 0.0_dp

  ! spatial grid
  x(0)   = 0.0_dp
  x(n+1) = a
  x(1:n) = [(x(0) + dx*i, i = 1, N)]

  ! Constructing the Hamiltonian Matrix
  Hamiltonian = 0.0_dp
  do i = 1, n
     Hamiltonian(i,i) = (2.0_dp + V_potential(i)) * lambda
     if( i /= n ) Hamiltonian(i+1,i) = -1.0_dp * lambda
     if( i /= 1 ) Hamiltonian(i-1,i) = -1.0_dp * lambda
  end do

  !do i = 1, n
  !   write(*, '(5(f10.3,2x))') hamiltonian(:,i)
  !end do
  
  call dgeev('N','V', n, Hamiltonian, n, Er, Ei,vl, n, vr, n, work, lwork, info)
  print*, info, work(1)
  !print*, Er

  open(unit = 100, file = 'data.dat')
  do i = 1, n
     write(100,*) i, Er(i), vr(i,:)
  end do
  
end program main
