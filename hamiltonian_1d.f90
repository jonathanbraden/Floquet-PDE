module Hamiltonian
  integer :: n_Hamiltonian_pieces = 2  ! number of terms in split Hamiltonian

  real(dl), dimension(n_Hamiltonian_pieces) :: tglobal   ! stores current value of the time that the fields and momenta are being stored at

 contains

   subroutine Hamiltonian_Split(dt, term_index)
     real(dl) :: dt
     integer :: term_index  ! labels which term of the split hamiltonian to use

     select case(term_index)
     case(1)
        tglobal(1) = tglobal(1) + dt
        call Hamiltonian_field(dt)
     case(2)
        tglobal(2) = tglobal(2) + dt
        call Hamiltonian_momentum(dt)
     case default
        print*, "Undefined Hamiltonian term"
        stop
     end select
   end subroutine Hamiltonian_Split

   subroutine Hamiltonian_field(dt)
     real(dl) :: dt
     integer :: i

     do i=1,nlat
        fld(i) = fld(i) + fldp(i)*dt
     enddo
   end subroutine Hamiltonian_field

! 
! Just need to fill in this subroutine
!
! Should use a macro to define my laplacian, and place the potential somewhere
! Also, i want to scan of k^2 values, so this needs to be permitted somehow
!
   subroutine Hamiltonian_momentum(dt)
     real(dl) :: dt
     real(dl) :: lap

     call boundary_conditions()

     do i=1,nlat
        lap = fld(i+1) + fld(i-1) - 2.*fld(i)
        lap = lap / dx**2
        fldp(i) = fldp(i) - dt * ( k2 + lap + vprime(fld(i)) )
     enddo
   end subroutine Hamiltonian_momentum

! Second derivative of the potential as a function of field phi
!
! Fix up to allow multiple fields
   real(dl) function vprime(phi)
     real(dl) :: phi

     vprime = 
   end function effective_mass

end module Hamiltonian
