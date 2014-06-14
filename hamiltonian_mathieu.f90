module Hamiltonian

  use Model

  implicit none

  integer, parameter :: n_Hamiltonian_pieces = 2  ! number of terms in split Hamiltonian

  real(dl), dimension(1:n_Hamiltonian_pieces) :: tglobal   ! stores current value of the time that the fields and momenta are being stored at
  real(dl) :: k2

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

     fld = fld + fldp*dt
   end subroutine Hamiltonian_field

! 
! Just need to fill in this subroutine
!
! Should use a macro to define my laplacian, and place the potential somewhere
! Also, i want to scan of k^2 values, so this needs to be permitted somehow
!
   subroutine Hamiltonian_momentum(dt)
     real(dl) :: dt
     integer :: i

     call boundary_conditions()

     fldp = fldp - dt*( k2 + sin(tglobal(1)))*fld
   end subroutine Hamiltonian_momentum

end module Hamiltonian
