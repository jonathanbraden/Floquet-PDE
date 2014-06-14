module Hamiltonian

  use fftw3
  use Model

  implicit none

  integer, parameter :: n_Hamiltonian_pieces = 2  ! number of terms in split Hamiltonian

  real(dl), dimension(1:n_Hamiltonian_pieces) :: tglobal   ! stores current value of the time that the fields and momenta are being stored at
  real(dl) :: k2

  real(C_DOUBLE), pointer :: laplace(:)
  complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)
  type(C_PTR) :: planf, planb

#define SPECTRAL 1
!#define DISCRETE_2 1
!#define DISCRETE_4 1

#define PERIODICBC
!#define DIRICHLETBC 1

 contains

   subroutine set_k2eff(keff)
     real(dl) :: keff
     k2 = keff
   end subroutine set_k2eff

   subroutine Hamiltonian_Split(dt, term_index, fld, fldp, tglob)
     real(dl) :: dt
     integer :: term_index  ! labels which term of the split hamiltonian to use
     real(dl), dimension(imin:imax) :: fld, fldp
     real(dl), dimension(1:2) :: tglob

     select case(term_index)
     case(1)
        tglob(1) = tglob(1) + dt
        call Hamiltonian_field(dt, fld, fldp)
     case(2)
        tglob(2) = tglob(2) + dt
        call Hamiltonian_momentum(dt, fld, fldp, tglob(1))
     case default
        print*, "Undefined Hamiltonian term"
        stop
     end select
   end subroutine Hamiltonian_Split

   subroutine Hamiltonian_field(dt, fld, fldp)
     real(dl) :: dt
     real(dl), dimension(imin:imax) :: fld, fldp
     integer :: i

     fld(1:nlat) = fld(1:nlat) + fldp(1:nlat)*dt
   end subroutine Hamiltonian_field

! 
! Just need to fill in this subroutine
!
! Should use a macro to define my laplacian, and place the potential somewhere
! Also, i want to scan of k^2 values, so this needs to be permitted somehow
!
   subroutine Hamiltonian_momentum(dt, fld, fldp, time)
     real(dl) :: dt
     real(dl), dimension(imin:imax) :: fld, fldp
     real(dl) :: time

     real(dl) :: lap
     integer :: i

!     call boundary_conditions()

#ifdef SPECTRAL
     laplace(1:nlat) = fld(1:nlat)
     call laplacian_1d(nlat, laplace, Fk, twopi/length, planf, planb)
     do i=1,nlat
        fldp(i) = fldp(i) -  &
             dt * ( (k2 + effective_mass(i*dx,time) )*fld(i) - laplace(i) )
     enddo
#endif
#ifdef DISCRETE_2
     call boundary_conditions()
     do i=1,nlat
        lap = (fld(i+1)+fld(i-1)-2.*fld(i)) / dx**2
        fldp(i) = fldp(i) - dt * ( (k2 + effective_mass(i*dx,time) )*fld(i) - lap )
     enddo
#endif
#ifdef DISCRETE_4
     call boundary_conditions()
     do i=1,nlat
        lap = - (fld(i+2) + fld(i-2)) + 16.*(fld(i+1)+fld(i-1)) - 30.*fld(i)
        lap = lap / 12./dx**2
        fldp(i) = fldp(i) - dt * ( (k2+effective_mass(i*dx,time) )*fld(i) -lap )
     enddo
#endif
   end subroutine Hamiltonian_momentum

  subroutine boundary_conditions()
    integer :: j

! subroutine to implement periodic boundary conditions
#ifdef PERIODICBC
    do j=1,pad
       fld(1-j) = fld(nlat + 1 - j)
       fld(nlat+j) = fld(j)
       fldp(1-j) = fldp(nlat + 1 - j)
       fldp(nlat+j) = fldp(j)
    enddo
#endif
#ifdef DIRICHLETBC
    fld(0) = 0.
    fld(nlat+1) = 0.
#endif

  end subroutine boundary_conditions

   subroutine init_arrays()
     call allocate_1d_array(nlat, laplace, Fk)
     planf = fftw_plan_dft_r2c_1d(nlat, laplace, Fk, FFTW_PATIENT+FFTW_DESTROY_INPUT)
     planb = fftw_plan_dft_c2r_1d(nlat, Fk, laplace, FFTW_PATIENT+FFTW_DESTROY_INPUT)
   end subroutine init_arrays

end module Hamiltonian
