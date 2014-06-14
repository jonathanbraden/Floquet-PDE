!
! To do : make the time-stepping in here adaptive
!

module Hamiltonian
  use fftw3
  use Model
  implicit none

  real(dl) :: k2

  real(C_DOUBLE), pointer :: laplace(:)
  complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)
  type(C_PTR) :: planf, planb

  integer, parameter :: nvar = 2*nlat+2
  real(dl), dimension(nvar) :: yvec

  contains
    subroutine step(dt, nsteps, fld, fldp, rcur)
      real(dl), intent(in) :: dt
      integer, intent(in) :: nsteps
      real(dl), dimension(nlat) :: fld, fldp
      real(dl), dimension(2) :: rcur
      integer :: i

      yvec(1:nlat) = fld(1:nlat)
      yvec(nlat+1:2*nlat) = fldp(1:nlat)
      yvec(2*nlat+1) = rcur(1)
      yvec(2*nlat+2) = rcur(2)
      do i=1,nsteps
         call gl10(yvec, dt)
         write(40,*) i*dt, yvec(2*nlat+1:2*nlat+2)
      enddo
      fld(1:nlat) = yvec(1:nlat)
      fldp(1:nlat) = yvec(nlat+1:2*nlat)
      rcur(1) = yvec(2*nlat+1)
      rcur(2) = yvec(2*nlat+2)
    end subroutine step

    subroutine set_k2eff(keff)
      real(dl) :: keff
      k2=keff
    end subroutine set_k2eff

    subroutine init_arrays()
      call allocate_1d_array(nlat, laplace, Fk)
      planf = fftw_plan_dft_r2c_1d(nlat, laplace, Fk, FFTW_PATIENT+FFTW_DESTROY_INPUT)
      planb = fftw_plan_dft_c2r_1d(nlat, Fk, laplace, FFTW_PATIENT+FFTW_DESTROY_INPUT)
    end subroutine init_arrays

!
! 10th order Gauss-Legendre integrator
! Used for pieces of Hamiltonian without exact solutions
! To do : implement entire evolution using this
! To do : pack into smaller vector as required (just chop a, so probably not worth it)
!
!
! Hmm, need to sort out the packing of the vector for this to work (c.f. matmul's that appear)
!
    subroutine gl10( y, dt )
      real*8 :: y(nvar)
      real*8 :: dt

      integer, parameter :: s = 5
      real*8 :: g(nvar, s)
      real*8 :: g_bg(2*nlat,s)
      real*8 :: g_fluc(2*nlat,s)

      ! Butcher tableau for 8th order Gauss-Legendre method
      real*8, parameter :: a(s,s) = reshape( (/ &
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
      real, parameter :: b(s) = (/ &
           1.1846344252809454375713202035995868132D-1,  2.3931433524968323402064575741781909646D-1, &
           2.8444444444444444444444444444444444444D-1,  2.3931433524968323402064575741781909646D-1, &
           1.1846344252809454375713202035995868132D-1 /)

      integer :: i,k

      g = 0.
      do k=1,16
#ifdef USEBLAS
         call DGEMM('N','N', nvar, s, s, 1., g, nvar, a, s, 0., g,nvar)
#else
         g = matmul(g,a)
#endif
         do i=1,s
            call derivs( y+g(:,i)*dt , g(:,i) )
         enddo
      enddo
#ifdef USEBLAS
      call DGEMV('N','N', nvar,s, dt,g,nvar, b,1 ,1.,y,1)
#else
      y = y + matmul(g,b)*dt
#endif
    end subroutine gl10
!
! Evolution vector in phase space for GL integration
!
    subroutine derivs(yc, yp)
      real*8, dimension(1:nvar) :: yc, yp
      integer :: i
! dphi/dt
      yp(1:nlat) = yc(nlat+1:2*nlat)
! dpi/dt
      laplace(:) = yc(1:nlat)
      call laplacian_1d(nlat, laplace, Fk, twopi/length, planf, planb)
      do i=1,nlat
         yp(nlat+i) = laplace(i) - ( k2 + effective_mass_rad(i*dx,yc(2*nlat+1:2*nlat+2)) )*yc(i)
      enddo

      yp(2*nlat+1) = yc(2*nlat+2)
      yp(2*nlat+2) = -vprime(yc(2*nlat+1),yc(2*nlat+2)) ! / gamma**3 / 2./ tension
    end subroutine derivs

  end module Hamiltonian
