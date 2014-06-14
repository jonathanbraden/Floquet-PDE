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

!#define GL_ORDER4 1
#define GL_ORDER6 1
!#define GL_ORDER10 1
 
#ifdef GL_ORDER6 
#define BVECTOR (/ 5.0/18.0, 4.0/9.0, 5.0/18.0 /) 
#endif

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
   
     real(dl) :: r

     select case(term_index)
     case(1)
        tglob(1) = tglob(1) + dt
        call Hamiltonian_field(dt, fld, fldp)
     case(2)
        tglob(2) = tglob(2) + dt
        r = r_evolution(tglob(1))
        call Hamiltonian_momentum(dt, fld, fldp, r)
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
   subroutine Hamiltonian_momentum(dt, fld, fldp, rcur)
     real(dl) :: dt
     real(dl), dimension(imin:imax) :: fld, fldp
     real(dl), intent(in) :: rcur
     real(dl) :: lap
     integer :: i

! Evolve the collective coordinate forward in time here

#ifdef SPECTRAL
     laplace(1:nlat) = fld(1:nlat)
     call laplacian_1d(nlat, laplace, Fk, twopi/length, planf, planb)
     do i=1,nlat
        fldp(i) = fldp(i) -  &
             dt * ( (k2 + effective_mass_rad(i*dx,rcur) )*fld(i) - laplace(i) )
     enddo
#endif
#ifdef DISCRETE_2
     call boundary_conditions()
     do i=1,nlat
        lap = (fld(i+1)+fld(i-1)-2.*fld(i)) / dx**2
        fldp(i) = fldp(i) - dt*( (k2+effective_mass_rad(i*dx,rvar) )*fld(i) - lap )
     enddo
#endif
#ifdef DISCRETE_4
     call boundary_conditions()
     do i=1,nlat
        lap = - (fld(i+2) + fld(i-2)) + 16.*(fld(i+1) + fld(i-1)) - 30.*fld(i)
        lap = lap / (12.*dx**2)
        fldp(i) = fldp(i) - dt*( (k2+effective_mass_rad(i*dx,rvar) )*fld(i) -lap )
     enddo
#endif
   end subroutine Hamiltonian_momentum

   subroutine gl(y,dt)
     real(dl), intent(inout) :: y(2)
     real(dl), intent(in) :: dt

#ifdef GL_ORDER4
     integer, parameter :: s = 2
     real*8, parameter :: a(s,s) = reshape((/ 0.25, 0.25 - 0.5/sqrt(3.0), 0.25 + 0.5/sqrt(3.0), 0.25 /), [s,s])
     real*8, parameter :: b(s) = (/ 0.5, 0.5 /)
#endif
#ifdef GL_ORDER6
     integer, parameter :: s = 3
     real(dl), parameter :: a(s,s) = reshape( (/ &
       5./36.0, 2.0/9.0 - 1./sqrt(15.0), 5.0/36.0 - 0.5/sqrt(15.0), &
       5.0/36.0 + sqrt(15.0)/24.0, 2.0/9.0, 5.0/36.0-sqrt(15.0)/24.0, &
       5.0/36.0 + 0.5/sqrt(15.0), 2.0/9.0+1.0/sqrt(15.0), 5.0/36.0 &
       /), [s,s]) 
     real(dl), parameter :: b(s) = BVECTOR
#endif
#ifdef GL_ORDER10
     integer, parameter :: s = 5
     real(dl), parameter :: a(s,s) = reshape( (/ & 
           0.5923172126404727187856601017997934066D-1, -1.9570364359076037492643214050884060018D-2, &
           1.1254400818642955552716244215090748773D-2, -0.5593793660812184876817721964475928216D-2, &
           1.5881129678659985393652424705934162371D-3,  1.2815100567004528349616684832951382219D-1, &
           1.1965716762484161701032287870890954823D-1, -2.4592114619642200389318251686004016630D-2, &
           1.0318280670683357408953945056355839486D-2, -2.7689943987696030442826307588795957613D-3, &
           1.1377628800422460252874127381536557686D-1,  2.6000465168064151859240589518757397939D-1, &
           1.4222222222222222222222222222222222222D-1, -2.0690316430958284571760137769754882933D-2, &
           4.6871545238699412283907465445931044619D-3,  1.2123243692686414680141465111883827708D-1, &
           2.2899605457899987661169181236146325697D-1,  3.0903655906408664483376269613044846112D-1, &
           1.1965716762484161701032287870890954823D-1, -0.9687563141950739739034827969555140871D-2, &
           1.1687532956022854521776677788936526508D-1,  2.4490812891049541889746347938229502468D-1, &
           2.7319004362580148889172820022935369566D-1,  2.5888469960875927151328897146870315648D-1, &
           0.5923172126404727187856601017997934066D-1 /) , [s,s])
     real(dl), parameter :: b(s) = (/ &
           1.1846344252809454375713202035995868132D-1,  2.3931433524968323402064575741781909646D-1, &
           2.8444444444444444444444444444444444444D-1,  2.3931433524968323402064575741781909646D-1, &
           1.1846344252809454375713202035995868132D-1 /)
#endif
     real(dl), dimension(2,s) :: g
     integer :: i,k
     integer, parameter :: niter=16

     g=0.
     do k=1,niter
        g = matmul(g,a)
        do i=1,s
           call derivs_bg( y+g(:,i)*dt, g(:,i) )
        enddo
     enddo
     y=y+matmul(g,b)*dt
   end subroutine gl

  subroutine boundary_conditions(fld)
    real(dl), dimension(imin:imax),intent(inout) :: fld
    integer :: j

! subroutine to implement periodic boundary conditions
#ifdef PERIODICBC
    do j=1,pad
       fld(1-j) = fld(nlat + 1 - j)
       fld(nlat+j) = fld(j)
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
