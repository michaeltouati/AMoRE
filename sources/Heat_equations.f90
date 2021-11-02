!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!! Angular Momentum Model Of Relativistic Electron beam (AMoRE) code !!
!!                                                                   !!
!! Copyright © 2015 Michaël J TOUATI                                 !!
!!                                                                   !!
!! This file is part of AMoRE.                                       !!
!!                                                                   !!
!! AMoRE is free software: you can redistribute it and/or modify     !!
!! it under the terms of the GNU General Public License as published !!
!! by the Free Software Foundation, either version 3 of the License, !!
!! or (at your option) any later version.                            !!
!!                                                                   !!
!! AMoRE is distributed in the hope that it will be useful,          !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with AMoRE. If not, see <https://www.gnu.org/licenses/>.    !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Initial commit written by Michaël J TOUATI - Oct. 2015
module heat_equations

  use acuracy
  use constants
  use transport_coef
  use input

  implicit none
  private :: electron_thermal_conduction
  public  :: hydrodynamic_quantities

  contains

  subroutine electron_thermal_conduction(Z, ni, Te_old, Ti_old, dt, Te, Ti)   
    ! input  : Z      = table containing the atomic number of the material in ()
    !          ni     = table containing the ion density in the material in (/cm^3)
    !          Te_old = table containing the electron temperature in the material at time t(n)
    !                   in (K)
    !          Ti_old = table containing the ion temperature in the material at time t(n) in (K)
    !          dt     = time step in (fs)
    !          Te     = table containing the electron temperature in the material at time t(n+1) 
    !                   in (K) before electron thermal diffusion
    !          Ti     = table containing the ion temperature in the material at time t(n+1) in (K)
    !                   before electron thermal diffusion
    ! output : Te = table containing the electron temperature in the material at time t(n+1) 
    !               in (K) after electron thermal diffusion
    !          Ti = table containing the ion temperature in the material at time t(n+1) in (K)
    !               after electron thermal diffusion
    implicit none
    real(PR), intent(in)                                :: dt
    real(PR), dimension(0:N_x+1,0:N_z+1), intent(in)    :: Z, ni
    real(PR), dimension(0:N_x+1,0:N_z+1), intent(in)    :: Te_old, Ti_old
    real(PR), dimension(0:N_x+1,0:N_z+1), intent(inout) :: Te, Ti
    integer                                             :: i,k
    real(PR)                                            :: kapa_iph_k, kapa_imh_k 
    real(PR)                                            :: kapa_i_kph, kapa_i_kmh
    real(PR)                                            :: div_q
    !$omp  PARALLEL DO DEFAULT(NONE) &
    !$omp& SHARED(N_z, N_x, dt, d_z, d_x, Z, ni, Te_old, Ti_old, Te) &
    !$omp& PRIVATE(i, k, kapa_iph_k, kapa_imh_k,kapa_i_kph, kapa_i_kmh,div_q) &
    !$omp& COLLAPSE(2)
    do k=1,N_z,1
      do i=1,N_x,1         
        ! Diffusion coefficients          
        kapa_iph_k = 2._PR * cond(Z(i+1,k), ni(i+1,k), Te_old(i+1,k), Ti_old(i+1,k)) &
        &* cond(Z(i,k), ni(i,k), Te_old(i,k), Ti_old(i,k))&
        / ( cond(Z(i+1,k), ni(i+1,k), Te_old(i+1,k), Ti_old(i+1,k)) &
        &+ cond(Z(i,k), ni(i,k), Te_old(i,k), Ti_old(i,k)) )          
        kapa_imh_k = 2._PR * cond(Z(i-1,k), ni(i-1,k), Te_old(i-1,k), Ti_old(i-1,k)) &
        &* cond(Z(i,k), ni(i,k), Te_old(i,k), Ti_old(i,k))&
        / ( cond(Z(i-1,k), ni(i-1,k), Te_old(i-1,k), Ti_old(i-1,k)) &
        &+ cond(Z(i,k), ni(i,k), Te_old(i,k), Ti_old(i,k)) )         
        kapa_i_kph = 2._PR * cond(Z(i,k+1), ni(i,k+1), Te_old(i,k+1), Ti_old(i,k+1)) &
        &* cond(Z(i,k), ni(i,k), Te_old(i,k), Ti_old(i,k))&
        / ( cond(Z(i,k+1), ni(i,k+1), Te_old(i,k+1), Ti_old(i,k+1)) &
        &+ cond(Z(i,k), ni(i,k), Te_old(i,k), Ti_old(i,k)) )         
        kapa_i_kmh = 2._PR * cond(Z(i,k-1), ni(i,k-1), Te_old(i,k-1), Ti_old(i,k-1)) &
        &* cond(Z(i,k), ni(i,k), Te_old(i,k), Ti_old(i,k))&
        / ( cond(Z(i,k-1), ni(i,k-1), Te_old(i,k-1), Ti_old(i,k-1)) &
        &+ cond(Z(i,k), ni(i,k), Te_old(i,k), Ti_old(i,k)) )          
        ! Thermal flux divergence div q = div( -kappa grad(Te) )      
        div_q=-(&
        (kapa_iph_k/(d_x*microns)**2._PR)*(Te_old(i+1,k)-Te_old(i,k)) &
        &- (kapa_imh_k/(d_x*microns)**2._PR)*(Te_old(i,k)-Te_old(i-1,k))&
        +&
        (kapa_i_kph/(d_z*microns)**2._PR)*(Te_old(i,k+1)-Te_old(i,k)) &
        &- (kapa_i_kmh/(d_z*microns)**2._PR)*(Te_old(i,k)-Te_old(i,k-1))&
        )                   
        ! 2nd order explicit scheme for the electron heat equation     
        Te(i,k) = Te(i,k) - ( ( dt*fs / capacity(Z(i,k), ni(i,k), Te_old(i,k)) ) * div_q )           
      end do
    end do
    !$omp END PARALLEL DO
    ! boundary conditions
    do k=1,N_z,1
      Te(0,k)     = Te(1,k)
      Te(N_x+1,k) = Te(N_x,k)
    end do
    do i=1,N_x,1
      Te(i,0)     = Te(i,1)
      Te(i,N_z+1) = Te(i,N_z)
    end do
    if (biTemperature.eqv..FALSE.) then
      Ti = Te
    end if    
  end subroutine electron_thermal_conduction

  subroutine hydrodynamic_quantities(dt, norm, Sva, phi, &
                                     A, Z, ni, Te, Ti, &
                                     E_x, E_z, jrx, jrz, &
                                     n_b, jbx, jbz, &
                                     Pdepos, Plost_col, Plost_res)
    ! input  : dt     = time step in (fs)
    !          norm   = coefficient of normalization for the distribution function in (/cm^3/keV)
    !          Sva    = table containing the kinetic energy grid in (keV)
    !          phi    = table containing the two first angular moment components of 
    !                   the distribution function at time t(n) in ()
    !          A      = table containing the atomic weight of the material in ()
    !          Z      = table containing the atomic number of the material in ()
    !          ni     = table containing the ion density in the material in (/cm^3)
    !          Te     = table containing the electron temperature in the material 
    !                   at time t(n) in (K)
    !          Ti     = at time t(n) table containing the ion temperature in the material in (K)
    !          Ex(z)  = table containing the electric field component on x(z)-axis in 
    !                   (statVolt/cm)
    !          jrx(z) = table containing the return current density component on x(z)-axis
    !                   in (statAmpere/cm)
    !          jbx(z) = table containing the fast e- beam current density component on x(z)-axis
    !                   in (statAmpere/cm)
    ! output : n_b       = table containing the fast e- beam density in (/cm^3)
    !          Pdepos    = table containing the density of power deposited on the background
    !                      electron in (erg/cm^3/s)
    !          Plost_col = table containing the density of power lost by the fast e- beam
    !                      due to collisions in (/erg/cm^3/s)
    !          Plost_res = table containing the density of power lost by the fast e- beam
    !                      due to the resistive electric field in (/erg/cm^3/s)
    real(PR), intent(in)                                                                      :: dt, norm
    real(PR), dimension(1:N_x,1:N_z,1:N_eps), intent(in)                                      :: Sva
    real(PR), dimension(forward:backward,psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2), intent(in) :: phi
    real(PR), dimension(0:N_x+1,0:N_z+1), intent(in)                                          :: A, Z, ni
    real(PR), dimension(0:N_x+1,0:N_z+1), intent(inout)                                       :: Te, Ti
    real(PR), dimension(1:N_x,1:N_z), intent(in)                                              :: E_x, E_z
    real(PR), dimension(1:N_x,1:N_z), intent(in)                                              :: jbx, jbz
    real(PR), dimension(1:N_x,1:N_z), intent(in)                                              :: jrx, jrz
    real(PR), dimension(1:N_x,1:N_z), intent(out)                                             :: n_b, Pdepos
    real(PR), dimension(1:N_x,1:N_z), intent(out)                                             :: Plost_col
    real(PR), dimension(1:N_x,1:N_z), intent(out)                                             :: Plost_res
    ! locals
    integer                              :: k, i
    real(PR)                             :: Omega, Cve, Cvi, Cv
    real(PR), dimension(0:N_x+1,0:N_z+1) :: Te_old, Ti_old
    Te_old = Te
    Ti_old = Ti
    !$omp PARALLEL DO DEFAULT(NONE) &
    !$omp SHARED(N_z, N_x, N_eps, dt, norm, d_eps, Sva, Z, A, ni) &
    !$omp SHARED(phi, n_b, Plost_col, Plost_res, Pdepos, Te, Ti) &
    !$omp SHARED(backward, biTemperature, E_x, E_z, jbx, jbz, jrx, jrz) &
    !$omp PRIVATE(k, i, Omega, Cv, Cve, Cvi) COLLAPSE(2)
    do k=1,N_z,1
      do i=1,N_x,1
        n_b(i,k)       = sum(phi(forward,psi0,i,k,1:N_eps)) * norm * d_eps
        Plost_col(i,k) = sum(Sva(i,k,1:N_eps) * phi(forward,psi0,i,k,1:N_eps)) * norm * d_eps
        if (backward == 2) then
          n_b(i,k)       = n_b(i,k) + sum(phi(backward,psi0,i,k,1:N_eps)) * norm * d_eps
          Plost_col(i,k) = Plost_col(i,k) + sum(Sva(i,k,1:N_eps) * phi(backward,psi0,i,k,1:N_eps)) * norm * d_eps
        end if
        Plost_res(i,k) = -( E_x(i,k) * jbx(i,k) ) - ( E_z(i,k) * jbz(i,k) )
        Pdepos(i,k)    = Plost_col(i,k) + ( E_x(i,k) * jrx(i,k) ) + ( E_z(i,k) * jrz(i,k) )
        if (biTemperature) then
          Omega = Omega_ei(A(i,k), Z(i,k), ni(i,k), Te(i,k), Ti(i,k))
          Cve   = capacity(Z(i,k), ni(i,k), Te(i,k))
          Cvi   = ion_capacity(Z(i,k), ni(i,k), Ti(i,k))
          Te(i,k) = Te(i,k) + ((Pdepos(i,k)/Cve)*dt*fs)&
          &- ((Omega/Cve)*(Te(i,k)-Ti(i,k))*dt*fs)
          Ti(i,k) = Ti(i,k) &
          &+ ((Omega/Cvi)*(Te(i,k)-Ti(i,k))*dt*fs)
        else
          Cv      = capacity(Z(i,k), ni(i,k), Te(i,k)) + ion_capacity(Z(i,k), ni(i,k), Ti(i,k))
          Te(i,k) = Te(i,k) + ((Pdepos(i,k)/Cv)*dt*fs)   
          Ti(i,k) = Te(i,k)
        end if    
      end do
    end do
    !$omp END PARALLEL DO
    if (thermal_conduction) then
    call electron_thermal_conduction(Z, ni, Te_old, Ti_old, dt, Te, Ti)
    end if
  end subroutine hydrodynamic_quantities

end module heat_equations