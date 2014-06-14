program odetest

implicit none

integer, parameter :: fields = 2
integer :: l
real*8, parameter :: rinit=4.

real*8 :: tcur, dt1, dt2, dr
real*8, parameter :: period = 0.5*3.14159/6.**0.5*exp(2.*rinit)!498.58635008
real*8, parameter :: dt = period / 2**15
real*8, parameter :: tend = 5.*period
real*8, dimension(2) :: y
y = (/ rinit, 0.0 /)

open(unit=30,file="collective_background.dat")

dr = 1./40.
do while (tcur < tend)
   dt1 = (1./20.)/( 1.e-7+abs(y(2)) ) 
   dt2 = dr/( 1.e-7+abs(vprime(y(1))) )
   dt1 = min(dt1,dt2,1./10.)

   write(30,'(4(E13.7,2x))') tcur/period, 2.**0.5*y, (potential(y(1))+0.5*2.**0.5/0.75*y(2)**2)
   tcur = tcur + dt1
   if ( tend - tcur < 0 ) then
      print*,"last step ",dt1, tend - ( tcur - dt1 )
      dt1 = tend - ( tcur - dt1 )
      tcur = tend 
   endif
   call gl10(y,dt1)
enddo
write(30,'(4(E13.7,2x))') tcur/period, 2.**0.5*y, (potential(y(1))+05.*2.**05./0.75*y(2)**2)

do l=-50,300
   dr=0.01
   write(25,'(4(E13.7,2x))') l*dr, vprime(l*dr)
enddo

contains

!
! Implicit Gauss-Legendre integrators
!
subroutine evalf(y, yprime)
  real*8, dimension(fields) :: y, yprime
  real*8 :: csh, cot, th, gamma
  
  csh = 1./(sinh(2.*y(1)))
  cot = 1./(tanh(2.*y(1)))
  th = tanh(2.*y(1))

  gamma = (1.-y(2)**2)**(-0.5)

  yprime(1) = y(2)
  yprime(2) = -vprime(y(1)) !/gamma**3
end subroutine evalf

real function potential(x)
  real*8 :: x
  real*8 :: csh, cot, gamma

  csh = 1./(sinh(2.*x))
  cot = 1./(tanh(2.*x))
  potential=-8./3. + 8.*x + 12.*cot - (24.*x+8.)*cot**2 + 16.*x*cot**3
  potential=2.**0.5*potential
end function potential

real function vprime(r)
  real*8 :: r
  real*8 :: th, cot, csh

  cot = 1./tanh(2.*r)
  th = tanh(2.*r)

  if (abs(r)<1.e-7) then
     vprime = 64.*r/5. - 128.*r**2/5. + 512.*r**3/105. + 512*r**4/21. - 2048.*r**5/175.
! - 4096*x**6/225. + 212992*r**7/17325 + 8192*r**8/693 - 141197312*x**9/14189175.
  else
     vprime = 32. - (96.*r+32.)*cot + (96.*r-48.)*cot**2 + (96.*r+48.)*cot**3 - 96.*r*cot**4
  endif
!  vprime = 2.**0.5/th**4
  vprime = 0.75*vprime*0.5
end function vprime

subroutine gl4(y, dt)
  integer, parameter :: s=2, n=2

  real*8 :: y(n), g(n,s), dt
  integer :: i,k

  ! Butcher tableu for 4th-order integrator
  real*8, parameter :: a(s,s) = reshape((/ 0.25, 0.25 - 0.5/sqrt(3.0), 0.25 + 0.5/sqrt(3.0), 0.25 /), [s,s])
  real*8, parameter :: b(s) = (/ 0.5, 0.5 /)

  g = 0.
  do k=1,16
     g = matmul(g,a)
     do i=1,s
        call evalf(y + g(:,i)*dt, g(:,i))
     enddo
  enddo

  y = y + matmul(g,b)*dt
end subroutine gl4

!
! Higher order integrators are created exactly the same way, except s is increased and the corresponding Butcher tableau is used
!
subroutine gl6(y, dt)
  integer, parameter :: s=3

  real*8 :: y(fields), g(fields,s), dt
  integer :: i,k

  ! Butcher tableu for 4th-order integrator
  real*8, parameter :: a(s,s) = reshape( (/ &
       5./36.0, 2.0/9.0 - 1./sqrt(15.0), 5.0/36.0 - 0.5/sqrt(15.0), &
       5.0/36.0 + sqrt(15.0)/24.0, 2.0/9.0, 5.0/36.0-sqrt(15.0)/24.0, &
       5.0/36.0 + 0.5/sqrt(15.0), 2.0/9.0+1.0/sqrt(15.0), 5.0/36.0 &
       /), [s,s])
  real*8, parameter :: b(s) = (/ 5.0/18.0, 4.0/9.0, 5.0/18.0 /)

  g = 0.
  do k=1,16
     g = matmul(g,a)
     do i=1,s
        call evalf(y + g(:,i)*dt, g(:,i))
     enddo
  enddo

  y = y + matmul(g,b)*dt
end subroutine gl6

    subroutine gl10( y, dt )
      real*8 :: y(fields)
      real*8 :: dt

      integer, parameter :: s = 5
      real*8 :: g(fields, s)

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
         g = matmul(g,a)
         do i=1,s
            call evalf( y+g(:,i)*dt , g(:,i) )
         enddo
      enddo
      y = y + matmul(g,b)*dt
    end subroutine gl10

end program odetest
