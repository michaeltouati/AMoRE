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
module vfp

use acuracy
use constants
use physics_library
use fokker_planck_coef
use input

implicit none
private  :: anisotropy_vector_norm, closure_parameter 
private :: flux_hll, flux_limitation, flux_hll_2nd_order
private :: collisional_derivative
public  :: F_z, F_x, cfl_scheme, grid, currents
public  :: spatial_derivatives, energy_derivative
public  :: source_terms, update, boundary_cond

contains

!=====================================================================================
!                             Minerbo's closure parameter
!=====================================================================================

pure subroutine anisotropy_vector_norm(phi, omega)    
! input  : phi = array containing the zero and first order angular moments
!                (Psi_0 Psi_1x Psi_1z)
! output : omega = anisotropy vector norm ( |Psi_1| / Psi_0 )
implicit none
real(PR), intent(in), dimension(psi0:psi1z)  :: phi
real(PR), intent(out)                        :: omega       
if (phi(psi0) <= 0._PR) then
omega = 0._PR               ! chosen isotrope by default
else
omega = sqrt((phi(psi1x)**2._PR) + (phi(psi1z)**2._PR))
omega = omega / phi(psi0) 
end if
! ensure 0 <= omega <= 1 :
omega = max(omega, 0._PR)
omega = min(omega,1._PR)
end subroutine anisotropy_vector_norm

pure subroutine closure_parameter(phi,mu)  
! input  : phi = array containing the zero and first order angular moments
!                (Psi_0 Psi_1x Psi_1z)
! output : mu = Minerbo's closure parameter
implicit none
real(PR), intent(in), dimension(psi0:psi1z) :: phi
real(PR), intent(out)              :: mu
real(PR)                           :: omega    
call anisotropy_vector_norm(phi, omega)
mu = ((omega**2._PR) / 2._PR) * ( 1._PR + (omega**2._PR))
mu = max(mu, 0._PR)
mu = min(mu,1._PR)
end subroutine closure_parameter

!=====================================================================================
!                                    Space fluxes
!=====================================================================================

pure subroutine F_x(phi, Fx)
! input  : phi = array containing the zero and first order angular moments
!                (Psi_0 Psi_1x Psi_1z)
! output : F_x = array containing the first and second order angular moment
!                projections on the x-axis (Psi_1x Psi_2xx Psi_2xz)
!                according to the Minerbo's closure.
implicit none
real(PR), intent(in), dimension(psi0:psi1z)  :: phi
real(PR), intent(out), dimension(psi0:psi1z) :: Fx
real(PR)                                     :: mu, psi_2xx, psi_2xz, ax, az             
Fx(psi0) = phi(psi1x)
call closure_parameter(phi,mu)
psi_2xx =  (1._PR/3._PR) * phi(psi0) * ( 1._PR - mu ) 
if (sqrt((phi(psi1x)**2._PR)+(phi(psi1z)**2._PR)).ne.0._PR) then 
ax = phi(psi1x)**2._PR
ax = ax / ( (phi(psi1x)**2._PR)+(phi(psi1z)**2._PR))
psi_2xx= psi_2xx + mu * phi(psi0) * ax
end if
Fx(psi1x) = psi_2xx
psi_2xz =  0._PR 
if (sqrt((phi(psi1x)**2._PR)+(phi(psi1z)**2._PR)).ne.0._PR) then 
az = phi(psi1x)*phi(psi1z)
az = az / ( (phi(psi1x)**2._PR) + (phi(psi1z)**2._PR) )
psi_2xz = mu * phi(psi0) * az
end if
Fx(psi1z) = psi_2xz
end subroutine F_x

pure subroutine F_z(phi, Fz)
! input  : phi = array containing the zero and first order angular moments
!                (Psi_0 Psi_1x Psi_1z)
! output : F_x = array containing the first and second order angular moment
!                projections on the z-axis (Psi_1z Psi_2zx Psi_2zz)
!                according to the Minerbo's closure
implicit none
real(PR), intent(in), dimension(psi0:psi1z)  :: phi
real(PR), intent(out), dimension(psi0:psi1z) :: Fz
real(PR)                                     :: mu, psi_2zx, psi_2zz, ax, az             
Fz(psi0) = phi(psi1z)
call closure_parameter(phi,mu)
psi_2zx =  0._PR 
if (sqrt((phi(psi1x)**2._PR)+(phi(psi1z)**2._PR)).ne.0._PR) then 
ax = phi(psi1x)*phi(psi1z)
ax = ax / ( (phi(psi1x)**2._PR) + (phi(psi1z)**2._PR) )
psi_2zx = mu * phi(psi0) * ax
end if
Fz(psi1x) = psi_2zx
psi_2zz =  (1._PR/3._PR) * phi(psi0) * ( 1._PR - mu ) 
if (sqrt((phi(psi1x)**2._PR)+(phi(psi1z)**2._PR)).ne.0._PR) then 
az = phi(psi1z)**2._PR
az = az / ( (phi(psi1x)**2._PR)+(phi(psi1z)**2._PR))
psi_2zz= psi_2zz + mu * phi(psi0) * az
end if
Fz(psi1z) = psi_2zz
end subroutine F_z

!=====================================================================================
!                                     HLL scheme
!=====================================================================================

pure subroutine flux_hll(phi_l, phi_r, flux_l, flux_r, fluxhll, velocity)
! input  : phi_l    = array containing the values of the zero and first order 
!                     angular moments (Psi_0 Psi_1x Psi_1z) in the cell at left side 
!                     of the interface
!          phi_r    = array containing the values of the zero and first order 
!                     angular moments (Psi_0 Psi_1x Psi_1z) in the cell at right side
!                     of the interface
!		   flux_l   = array containing the values of the flux in the cell 
!                     at left side of the interface
!		   flux_r   = array containing the values of the flux in the cell 
!                     at right side of the interface
!          velocity = flux velocity at the interface
! output : flux_hll = first order HLL flux crossing the interface according to
!          Berthon C. et al., Adv. Appl. Math. Mech. 2, 3, p 259 (2010)
implicit none
real(PR), intent(in), dimension(psi0:psi1z)  :: phi_l, phi_r, flux_l, flux_r
real(PR), intent(in)                         :: velocity
real(PR), intent(out), dimension(psi0:psi1z) :: fluxhll
fluxhll = 0.5_PR*(flux_l + flux_r)
fluxhll = velocity*fluxhll + (0.5_PR*abs(velocity)*(phi_l - phi_r)) 
end subroutine flux_hll

pure subroutine flux_limitation(phi, dphi)
! input  : phi = array containing the zero and first order angular moments
!                (Psi_0 Psi_1x Psi_1z)
!          dphi = array containing the second order correction to get the value 
!                 of phi at +/- half the grid size phi+/- = phi +/- dphi
! output : dphi = array containing the modified second order correction 
!                 computed in order to preserve Psi_0 >= 0 and |Psi_1| <= Psi_0
real(PR), intent(in),    dimension(psi0:psi1z) :: phi
real(PR), intent(inout), dimension(psi0:psi1z) :: dphi
real(PR)                                       :: theta, a, b, c, delta
real(PR), dimension(4)                         :: root
if (hll_order == 1) then 
theta  = 0._PR        ! 1st order
else                     ! 2nd order
a   = (dphi(psi0)**2._PR) - (dphi(psi1z)**2._PR) - (dphi(psi1x)**2._PR)
b   = phi(psi0) * dphi(psi0) &
&- (phi(psi1z) * dphi(psi1z)) - (phi(psi1x) * dphi(psi1x))
c   = (phi(psi0)**2._PR) - (phi(psi1z)**2._PR) - (phi(psi1x)**2._PR)
! solve a theta^2 +/- b theta + c >= 0 :
delta = (b**2) - a * c   
if (c < 0._PR) then
theta = 0._PR ! we come back to the first order  
! because phi is out of the realizability domain (|Psi_1| > Psi_0)  
else if (delta < 0._PR .or. (a == 0._PR .and. b == 0._PR)) then
theta = 1._PR ! we do not need to limit the flux to preserve 
! Psi_0 > 0 and |Psi_1| < Psi_0
else if (a == 0._PR) then
root(1:2) = -0.5_PR * c / b
root(3:4) =  0.5_PR * c / b
theta = min(1._PR, minval(root, root > 0._PR))
else
root(1) = (-b - sqrt(delta)) / a
root(2) = (-b + sqrt(delta)) / a
root(3) = ( b - sqrt(delta)) / a
root(4) = ( b + sqrt(delta)) / a
theta = min(1._PR, minval(root, root > 0._PR))
end if
end if
dphi = theta * dphi
end subroutine flux_limitation

pure subroutine flux_hll_2nd_order(psi1mu, phi_l2, phi_l, phi_m, phi_r, phi_r2,&
flux_balance, velocity, alpha_l, alpha_r)
! input  : psi1mu   = direction index psi1x or psi1z
!          phi_l2   = array containing the values of the zero and first order 
!                     angular moments (Psi_0 Psi_1x Psi_1z) two cells at left side 
!                     of the considered cell m
!          phi_l    = array containing the values of the zero and first order 
!                     angular moments (Psi_0 Psi_1x Psi_1z) one cell at left side 
!                     of the considered cell m
!          phi_m    = array containing the values of the zero and first order 
!                     angular moments (Psi_0 Psi_1x Psi_1z) in the considered cell m
!          phi_r    = array containing the values of the zero and first order 
!                     angular moments (Psi_0 Psi_1x Psi_1z) one cell at right side 
!                     of the considered cell m
!          phi_r2   = array containing the values of the zero and first order 
!                     angular moments (Psi_0 Psi_1x Psi_1z) two cells at right side 
!                     of the considered cell m
!          velocity = flux velocity in the considered cell m
!          alpha_l  = correction factor to get the flux at the left interface
!                     of the considered cell 
!          alpha_r  = correction factor to get the flux at the right interface
!                     of the considered cell
! output : flux_balance = flux balance between the HLL flux at the right and left 
!                         interface
implicit none
integer, intent(in)                          :: psi1mu
real(PR), intent(in), dimension(psi0:psi1z)  :: phi_l2,phi_l, phi_m, phi_r, phi_r2
real(PR), intent(in)                         :: velocity, alpha_l, alpha_r
real(PR), intent(out), dimension(psi0:psi1z) :: flux_balance
real(PR), dimension(psi0:psi1z)              :: p_g, p_m, p_d
real(PR), dimension(psi0:psi1z)              :: phi_l_plus, phi_m_minus
real(PR), dimension(psi0:psi1z)              :: phi_m_plus, phi_r_minus 
real(PR), dimension(psi0:psi1z)              :: flux_l_plus, flux_m_minus
real(PR), dimension(psi0:psi1z)              :: flux_m_plus, flux_r_minus 
real(PR), dimension(psi0:psi1z)              :: fluxhll_l, fluxhll_r
real(PR), dimension(psi0:psi1z)              :: phi_l_plus_temp, flux_l_plus_temp
real(PR), dimension(psi0:psi1z)              :: phi_r_minus_temp, flux_r_minus_temp
p_m = max(0._PR,min(phi_r-phi_m, phi_m-phi_l)) &
+ min(0._PR,max(phi_r-phi_m, phi_m-phi_l))
p_m = 0.5_PR*p_m
call flux_limitation(phi_m, p_m)
p_g = max(0._PR,min(phi_m-phi_l, phi_l-phi_l2)) &
+ min(0._PR,max(phi_m-phi_l, phi_l-phi_l2))
p_g = 0.5_PR*p_g
call flux_limitation(phi_l, p_g)
p_d = max(0._PR,min(phi_r2-phi_r, phi_r-phi_m)) &
+ min(0._PR,max(phi_r2-phi_r, phi_r-phi_m))
p_d = 0.5_PR*p_d
call flux_limitation(phi_r, p_d)
phi_l_plus  = phi_l + p_g
phi_m_minus = phi_m - p_m
phi_m_plus  = phi_m + p_m
phi_r_minus = phi_r - p_d
if (psi1mu == psi1z) then
call F_z(phi_l_plus, flux_l_plus)
call F_z(phi_m_minus,flux_m_minus)
call F_z(phi_m_plus, flux_m_plus)
call F_z(phi_r_minus, flux_r_minus)
else if (psi1mu == psi1x) then
call F_x(phi_l_plus, flux_l_plus)
call F_x(phi_m_minus,flux_m_minus)
call F_x(phi_m_plus, flux_m_plus)
call F_x(phi_r_minus, flux_r_minus)
end if
phi_r_minus_temp  = alpha_r*phi_r_minus
flux_r_minus_temp = alpha_r*flux_r_minus
call flux_hll(phi_m_plus, phi_r_minus_temp, flux_m_plus,&
flux_r_minus_temp, fluxhll_r, velocity)
phi_l_plus_temp  = alpha_l*phi_l_plus
flux_l_plus_temp = alpha_l*flux_l_plus
call flux_hll(phi_l_plus_temp, phi_m_minus, flux_l_plus_temp,&
flux_m_minus, fluxhll_l, velocity)
flux_balance = fluxhll_r - fluxhll_l
end subroutine flux_hll_2nd_order

!=====================================================================================
!                                  CFL condition
!=====================================================================================

pure subroutine cfl_scheme(epsa, Sva, nua, Ex, Ez, By, dt)
! input  : epsa = array containing the kinetic energy grid in (keV)
!          Sva  = table containing the total stopping power
!                 of fast e- with kinetic energy eps(l) multiplied by their 
!                 corresponding velocity in (erg/s)
!          nua  = table containing the total angular collision rate of fast e-
!                 with kinetic energy eps(l) 
!          Ex   = table containing the electric field component on the x-axis
!          Ez   = table containing the electric field component on the z-axis
!          By   = table containing the magnetic field component on the y-axis
! output : dt   = step time computed according to an approximate CFL accounting
!                 for all numerical schemes used to solve the angular moment 
!                 equations. It is relaxed when using the implicit scheme for
!                 collisional effects. 
implicit none
real(PR), dimension(1:N_eps), intent(in)             :: epsa
real(PR), dimension(1:N_x,1:N_z,1:N_eps), intent(in) :: Sva, nua
real(PR), dimension(1:N_x,1:N_z), intent(in)         :: Ex, Ez, By
real(PR), intent(out)                                :: dt
real(PR)                                             :: Vmax, Vepsmax, Exx, Ezz, Byy
real(PR)                                             :: nnu, pvmin, pcmin
Vmax  = vit(L_eps)   * (fs / microns)
pvmin = mom(epsa(1)) * vit(epsa(1)) / keV
pcmin = mom(epsa(1)) * c           / keV
if (coll_implicit_scheme) then
Vepsmax = 0._PR
nnu     = 0._PR
else
Vepsmax = maxval(maxval(maxval(Sva,dim=3),dim=2)) * fs / keV
nnu     = maxval(maxval(maxval(nua,dim=3),dim=2)) * fs
end if
Exx = e*maxval(maxval(abs(Ex(:,:)),dim=2)) * (Vmax*microns/fs) * fs/keV
Ezz = e*maxval(maxval(abs(Ez(:,:)),dim=2)) * (Vmax*microns/fs) * fs/keV
Byy = e*maxval(maxval(abs(By(:,:)),dim=2)) * (Vmax*microns/fs) * fs/keV
dt   = cfl / &
( &
(real(hll_order)*Vmax * ((1._PR/d_x)+(1._PR/d_z))) &
+ (Vepsmax/d_eps) &
+ nnu&
+ (real(hll_order)*sqrt((Exx**2._PR)+(Ezz**2._PR))/d_eps)&
+ (sqrt((Exx**2._PR)+(Ezz**2._PR))/pvmin) &
+ (Byy/pcmin) &
)
end subroutine cfl_scheme

!=====================================================================================
!                                        loops
!=====================================================================================

pure subroutine grid(za,xa,epsa)
! output : za   = array containing the spatial z-axis grid in (microns)
!          xa   = array containing the spatial x-axis grid in (microns)
!          epsa = array containing the kinetic energy grid in (keV)
implicit none
real(PR), dimension(1:N_z), intent(out)               :: za
real(PR), dimension(1:N_x), intent(out)               :: xa
real(PR), dimension(1:N_eps), intent(out)             :: epsa
integer                                               :: i,k,l
do l=1,N_eps,1
epsa(l) = eps_min + 0.5_PR*d_eps + (real(l- 1) * d_eps)
end do
do i=1,N_x,1
xa(i) = -(L_x/2._PR) + 0.5_PR*d_x + (real(i-1)*d_x)
end do
do k=1,N_z,1
za(k) = 0.5_PR*d_z + (real(k-1)*d_z)
end do
end subroutine grid

subroutine currents(norma,epsa,phi,By,jbx,jbz,jrx,jrz)
! input  : norma = coefficient of normalization for the distribution function in (/cm^3/keV)
!          epsa  = array containing the kinetic energy grid in (keV)
!          phi   = table containing the two first angular moment components of 
!                  the distribution function in ()
!          By    = table containing the magnetic field in (gauss)
! output : jbx   = table containing the beam current density component on the 
!                  x-axis computed by summing over the first order angular moment over
!                  the kinetic energy grid in (statAmpere/cm)
!          jbz   = table containing the beam current density component on the 
!                  z-axis computed by summing the first order angular moment over
!                  the kinetic energy grid in (statAmpere/cm)
!          jrx   = table containing the return current component on the x-axis
!                  computed according to the Ampere's equation in (statAmpere/cm)
!          jrz   = table containing the return current component on the z-axis
!                  computed according to the Ampere's equation in (statAmpere/cm)
implicit none
real(PR), intent(in)                                                     :: norma
real(PR), dimension(1:N_eps), intent(in)                                 :: epsa
real(PR), dimension(forward:backward,psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2), intent(in) :: phi
real(PR), dimension(1:N_x,1:N_z), intent(in)                             :: By
real(PR), dimension(1:N_x,1:N_z), intent(out)                            :: jbx,jbz,jrx,jrz
integer                                                                  :: i,k
!$omp PARALLEL DO DEFAULT(SHARED) PRIVATE(k,i) COLLAPSE(2)
do k=1,N_z,1
do i=1,N_x,1
jbx(i,k) = sum(phi(forward,psi1x,i,k,1:N_eps)*vit(epsa(1:N_eps)))*(-e)*norma*d_eps
jbz(i,k) = sum(phi(forward,psi1z,i,k,1:N_eps)*vit(epsa(1:N_eps)))*(-e)*norma*d_eps
if (backward == 2) then
jbx(i,k) = jbx(i,k) + sum(phi(backward,psi1x,i,k,1:N_eps)*vit(epsa(1:N_eps)))*(-e)*norma*d_eps
jbz(i,k) = jbz(i,k) + sum(phi(backward,psi1z,i,k,1:N_eps)*vit(epsa(1:N_eps)))*(-e)*norma*d_eps
end if
if ((k.ne.1).and.(k.ne.N_z)) then
jrx(i,k) = - jbx(i,k) - ((c/(4._PR*pi))*((By(i,k+1)-By(i,k-1))/(2._PR*d_z*microns)))
else if (k == 1) then
jrx(i,k) = - jbx(i,k) - ((c/(4._PR*pi))*((By(i,k+1)-By(i,k  ))/(2._PR*d_z*microns)))
else
jrx(i,k) = - jbx(i,k) - ((c/(4._PR*pi))*((By(i,k  )-By(i,k-1))/(2._PR*d_z*microns)))
end if
if ((i.ne.1).and.(i.ne.N_x)) then
jrz(i,k) = - jbz(i,k) + ((c/(4._PR*pi))*((By(i+1,k)-By(i-1,k))/(2._PR*d_x*microns)))
else if (i == 1) then
jrz(i,k) = - jbz(i,k) + ((c/(4._PR*pi))*((By(i+1,k)-0._PR    )/(2._PR*d_x*microns)))
else
jrz(i,k) = - jbz(i,k) + ((c/(4._PR*pi))*((0._PR    -By(i-1,k))/(2._PR*d_x*microns)))
end if
end do
end do
!$omp END PARALLEL DO
end subroutine currents

subroutine spatial_derivatives(eps_a,phi_n,dF_dz,dF_dx)
! input  : eps_a = array containing the kinetic energy grid in (keV)
!          phi_n = table containing the two first angular moment components of 
!                  the distribution function in ()
! output : dF_dz = table containing the z-derivative components of the angular 
!                  moment equations computed according to the HLL scheme
!          dF_dx = table containing the x-derivative components of the angular 
!                  moment equations computed according to the HLL scheme
implicit none
real(PR), dimension(1:N_eps), intent(in)                                                  :: eps_a
real(PR), dimension(forward:backward,psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2), intent(in) :: phi_n
real(PR), dimension(forward:backward,psi0:psi1z,1:N_x,1:N_z,1:N_eps), intent(out)         :: dF_dz, dF_dx
! locals
integer                                                                                   :: m,l,k,i
real(PR), dimension(psi0:psi1z)                                                           :: phi_r2,phi_r,phi_m,phi_l,phi_l2
real(PR), dimension(psi0:psi1z)                                                           :: flux_balance
real(PR)                                                                                  :: velocity, alpha_l, alpha_r
alpha_l = 1._PR 
alpha_r = 1._PR
!$omp PARALLEL DO DEFAULT(SHARED) PRIVATE(m,l,k,i,phi_r2,phi_r,phi_m,phi_l,phi_l2,velocity,flux_balance) COLLAPSE(4)
do m=forward,backward,1
do l=1,N_eps,1
do k=1,N_z,1
do i=1,N_x,1
phi_l2 = phi_n(m,psi0:psi1z,i,k-2,l)
phi_l  = phi_n(m,psi0:psi1z,i,k-1,l)
phi_m  = phi_n(m,psi0:psi1z,i,k  ,l)
phi_r  = phi_n(m,psi0:psi1z,i,k+1,l)
phi_r2 = phi_n(m,psi0:psi1z,i,k+2,l)
velocity = vit(eps_a(l))*fs/microns
call flux_hll_2nd_order(psi1z, phi_l2, phi_l, phi_m, phi_r, phi_r2,&
flux_balance, velocity, alpha_l, alpha_r)
dF_dz(m,psi0:psi1z,i,k,l) = flux_balance(psi0:psi1z)/d_z
end do
end do
end do
end do
!$omp END PARALLEL DO
!$omp PARALLEL DO DEFAULT(SHARED) PRIVATE(m,l,k,i,phi_r2,phi_r,phi_m,phi_l,phi_l2,velocity,flux_balance) COLLAPSE(4)
do m=forward,backward,1
do l=1,N_eps,1
do k=1,N_z,1       
do i=1,N_x,1
phi_l2 = phi_n(m,psi0:psi1z,i-2,k,l)
phi_l  = phi_n(m,psi0:psi1z,i-1,k,l)
phi_m  = phi_n(m,psi0:psi1z,i  ,k,l)
phi_r  = phi_n(m,psi0:psi1z,i+1,k,l)
phi_r2 = phi_n(m,psi0:psi1z,i+2,k,l)
velocity = vit(eps_a(l))*fs/microns
call flux_hll_2nd_order(psi1x, phi_l2, phi_l, phi_m, phi_r, phi_r2,&
flux_balance, velocity, alpha_l, alpha_r)
dF_dx(m,psi0:psi1z,i,k,l)=flux_balance(psi0:psi1z)/d_x
end do
end do
end do
end do
!$omp END PARALLEL DO
end subroutine spatial_derivatives

pure subroutine collisional_derivative(phin, Sva_l, Sva_lp1, der)
! input  : phin  = array containg the two first angular moment components of 
!                  the distribution function in ()
!          Sva_l = array containing the total stopping power
!                  of fast e- with kinetic energy eps(l) multiplied by their 
!                  corresponding velocity in (erg/s)
!          Sva_lp1 = array containing the total stopping power
!                  of fast e- with kinetic energy eps(l+1) multiplied by their 
!                  corresponding velocity in (erg/s)
! output : der = array containing the kinetic energy derivative of the angular 
!                moment equations due to the collisional slowing down of fast e-.
!                It is computed according to the explicit downwind scheme
implicit none
real(PR), dimension(psi0:psi1z,-1:N_eps+2), intent(in) :: phin 
real(PR), dimension(1:N_eps), intent(in)               :: Sva_l, Sva_lp1
real(PR), dimension(psi0:psi1z,1:N_eps), intent(out)   :: der
! locals
integer                                                :: l
real(PR), dimension(psi0:psi1z)                        :: Fe_l, Fe_r
do l=1,N_eps,1
Fe_r(psi0:psi1z) = Sva_lp1(l) * (fs / keV) * phin(psi0:psi1z,l+1) 
Fe_l(psi0:psi1z) = Sva_l(l)   * (fs / keV) * phin(psi0:psi1z,l)
der(psi0:psi1z,l) = - (Fe_r(psi0:psi1z) - Fe_l(psi0:psi1z)) / d_eps
end do
end subroutine collisional_derivative

subroutine energy_derivative(E_x,E_z,eps_a, Sv_a,Sv_a_lp1,phi_n,dF_deps)
! input  : E_x(z) = table containing the electric field component on the x(z)-axis
!                   in (statVolt/cm)
!          eps_a  = array containing the kinetic energy grid in (keV)
!          Sv_a   = table containing the total stopping power
!                   of fast e- with kinetic energy eps(l) multiplied by their 
!                   corresponding velocity in (erg/s)
!          Sv_a_lp1 = Sva_lp1 = table containing the total stopping power
!                     of fast e- with kinetic energy eps(l+1) multiplied by their 
!                     corresponding velocity in (erg/s)
!          phi_n = table containing the two first angular moment components of 
!                  the distribution function in ()
! output : dF_deps = table containing the kinetic energy derivative components 
!                    of the angular moment equations computed according to 
!                    the HLL scheme
implicit none
real(PR), dimension(1:N_x,1:N_z), intent(in)                             :: E_x, E_z
real(PR), dimension(1:N_eps), intent(in)                                 :: eps_a
real(PR), dimension(1:N_x,1:N_z,1:N_eps), intent(in)                     :: Sv_a,Sv_a_lp1
real(PR), dimension(forward:backward,psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2), intent(in) :: phi_n
real(PR), dimension(forward:backward,psi0:psi1z,1:N_x,1:N_z,1:N_eps), intent(out)         :: dF_deps
! locals
integer                                                                  :: m,k,i,l
real(PR)                                                                 :: E_x_temp, E_z_temp 
real(PR), dimension(:), allocatable                                      :: Sv_a_temp,Sv_a_lp1_temp
real(PR), dimension(:,:), allocatable                                    :: phi_temp_eps
real(PR), dimension(:,:), allocatable                                    :: der_collisional
real(PR), dimension(psi0:psi1z)                                          :: phi_l2,phi_l,phi_m,phi_r,phi_r2
real(PR), dimension(psi0:psi1z)                                          :: flux_balance_x
real(PR), dimension(psi0:psi1z)                                          :: flux_balance_z
real(PR)                                                                 :: alpha_r, alpha_l
!$omp PARALLEL DO DEFAULT(SHARED) &
!$omp& PRIVATE(m,l,k,i,alpha_l,alpha_r,E_x_temp, E_z_temp,&
!$omp& phi_l2,phi_l,phi_m,phi_r,phi_r2,&
!$omp& flux_balance_x,flux_balance_z) COLLAPSE(4)
do m=forward,backward,1
do l=1,N_eps,1
do k=1,N_z,1
do i=1,N_x,1  
E_x_temp  = -e*E_x(i,k) * vit(eps_a(l)) * fs / keV
E_z_temp  = -e*E_z(i,k) * vit(eps_a(l)) * fs / keV
alpha_l = vit(eps_a(l)-d_eps) / vit(eps_a(l))
alpha_r = vit(eps_a(l)+d_eps) / vit(eps_a(l))
phi_l2 = phi_n(m,psi0:psi1z,i,k,l-2)
phi_l  = phi_n(m,psi0:psi1z,i,k,l-1)
phi_m  = phi_n(m,psi0:psi1z,i,k,l  )
phi_r  = phi_n(m,psi0:psi1z,i,k,l+1)
phi_r2 = phi_n(m,psi0:psi1z,i,k,l+2)
call flux_hll_2nd_order(psi1x, phi_l2, phi_l, phi_m, phi_r, phi_r2,&
flux_balance_x, E_x_temp, alpha_l, alpha_r)
call flux_hll_2nd_order(psi1z, phi_l2, phi_l, phi_m, phi_r, phi_r2,&
flux_balance_z, E_z_temp, alpha_l, alpha_r)
dF_deps(m,psi0:psi1z,i,k,l) = (flux_balance_x(psi0:psi1z) &
&+  flux_balance_z(psi0:psi1z))/d_eps
end do    
end do
end do
end do
!$omp END PARALLEL DO
if (coll_implicit_scheme.eqv..false.) then
!$omp PARALLEL DO DEFAULT(SHARED) &
!$omp& PRIVATE(m,k,i,phi_temp_eps,Sv_a_temp,Sv_a_lp1_temp,der_collisional) COLLAPSE(3)
do m=forward,backward,1
do k=1,N_z,1
do i=1,N_x,1 
allocate(phi_temp_eps(psi0:psi1z,-1:N_eps+2),der_collisional(psi0:psi1z,1:N_eps))
allocate(Sv_a_temp(1:N_eps),Sv_a_lp1_temp(1:N_eps))
phi_temp_eps(psi0:psi1z,-1:N_eps+2) = phi_n(m,psi0:psi1z,i,k,-1:N_eps+2)
Sv_a_temp(1:N_eps)     = Sv_a(i,k,1:N_eps)
Sv_a_lp1_temp(1:N_eps) = Sv_a_lp1(i,k,1:N_eps)
call collisional_derivative(phi_temp_eps(psi0:psi1z,-1:N_eps+2),&
Sv_a_temp(1:N_eps), Sv_a_lp1_temp(1:N_eps),&
der_collisional(psi0:psi1z,1:N_eps))
dF_deps(m,psi0:psi1z,i,k,1:N_eps) = dF_deps(m,psi0:psi1z,i,k,1:N_eps) &
&+ der_collisional(psi0:psi1z,1:N_eps)
deallocate(phi_temp_eps,der_collisional,Sv_a_temp,Sv_a_lp1_temp)
end do
end do
end do
end if
end subroutine energy_derivative

subroutine source_terms(E_x,E_z,B_y,eps_a,nu_a,phi_n,gama_E,gama_B,gama_nu)
! input  : E_x(z)= table containing the electric field component on the x(z)-axis
!                  in (statVolt/cm)
!          B_y   = table containing the magnetic field component on the y-axis in (gauss)
!          eps_a = array containing the kinetic energy grid in (keV)
!          nu_a  = table containing the total angular collision rate of fast e-
!                  with kinetic energy eps(l) in (/s)
!          phi_n = table containing the two first angular moment components of 
!                  the distribution function in ()
! output : gama_E  = table containing the source term of angular moment equations
!                     depending on the electric fields
!          gama_B  = table containing the source term of angular moment equations
!                     depending on the magnetic fields
!          gama_nu = table containing the source term of angular moment equations
!                    depending on the angular collision rate
implicit none
real(PR), dimension(1:N_x,1:N_z), intent(in)                             :: E_x, E_z, B_y
real(PR), dimension(1:N_eps), intent(in)                                 :: eps_a
real(PR), dimension(1:N_x,1:N_z,1:N_eps), intent(in)                     :: nu_a
real(PR), dimension(forward:backward,psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2), intent(in) :: phi_n
real(PR), dimension(forward:backward,psi0:psi1z,1:N_x,1:N_z,1:N_eps), intent(out)         :: gama_E, gama_B 
real(PR), dimension(forward:backward,psi0:psi1z,1:N_x,1:N_z,1:N_eps), intent(inout)       :: gama_nu
! locals
real(PR)                                                                 :: Exv, Ezv, pv, Byv, pc, nnu
real(PR), dimension(psi0:psi1z)                                          :: phi_temp, Fx, Fz
integer                                                                  :: m,l,k,i
!$omp PARALLEL DO DEFAULT(SHARED) PRIVATE(m,l,k,i,Ezv,Exv,pv,phi_temp,Fz,Fx,Byv,pc,nnu) COLLAPSE(4)
do m=forward,backward,1
do l=1,N_eps,1
do k=1,N_z,1      
do i=1,N_x,1                     
! 'Electric Force'
Ezv = -e*E_z(i,k)*vit(eps_a(l))* fs / keV
Exv = -e*E_x(i,k)*vit(eps_a(l))* fs / keV
pv = mom(eps_a(l))*vit(eps_a(l)) / keV
phi_temp(psi0:psi1z) = phi_n(m,psi0:psi1z,i,k,l)                  
call F_z(phi_temp, Fz)
call F_x(phi_temp, Fx)
gama_E(m,psi0,i,k,l)  = 0._PR
gama_E(m,psi1x,i,k,l) = (Exv/pv) * phi_n(m,psi0,i,k,l)  &
&- ( Fx(psi1x)*(Exv/pv) + Fx(psi1z)*(Ezv/pv) ) 
gama_E(m,psi1z,i,k,l) = (Ezv/pv) * phi_n(m,psi0,i,k,l)  &
&- ( Fz(psi1x)*(Exv/pv) + Fz(psi1z)*(Ezv/pv) )                  
! 'Magnetic Force'
Byv = -e*B_y(i,k)*vit(eps_a(l))* fs / keV
pc = mom(eps_a(l)) * c / keV
gama_B(m,psi0,i,k,l)  = 0._PR
gama_B(m,psi1x,i,k,l) = - (Byv/pc) * phi_n(m,psi1z,i,k,l)
gama_B(m,psi1z,i,k,l) =   (Byv/pc) * phi_n(m,psi1x,i,k,l)                 
if (coll_implicit_scheme) then
gama_nu(m,psi0:psi1z,i,k,l) = 0._PR
else
! 'Diffusion Strenght'
nnu = fs * nu_a(i,k,l)
gama_nu(m,psi0,i,k,l)  = 0._PR
gama_nu(m,psi1x,i,k,l) =  - nnu * phi_n(m,psi1x,i,k,l)
gama_nu(m,psi1z,i,k,l) =  - nnu * phi_n(m,psi1z,i,k,l)
end if
end do
end do
end do
end do
!$omp END PARALLEL DO
end subroutine source_terms

subroutine update(phi_n, d_t, dF_dx, dF_dz, dF_deps, gama_e, gama_b, gama_nu, Sv_a, Sv_a_lp1, nu_a, phi_np1)
! input  : phi_n    = table containing the two first angular moment components of 
!                     the distribution function at time t(n)
!          d_t      = time step at the considered iteration n
!          dF_dx    = table containing the x-derivative components of the angular 
!                     moment equations computed according to the HLL scheme
!          dF_dz    = table containing the z-derivative components of the angular 
!                     moment equations computed according to the HLL scheme
!          dF_deps  = table containing the kinetic energy derivative components 
!                     of the angular moment equations computed according to 
!                     the HLL scheme
!          gama_E  = table containing the source term of angular moment equations
!                     depending on the electric fields
!          gama_B  = table containing the source term of angular moment equations
!                     depending on the magnetic fields
!          gama_nu = table containing the source term of angular moment equations
!                     depending on the angular collision rate
!          Sv_a     = table containing the total stopping power
!                     of fast e- with kinetic energy eps(l) multiplied by their 
!                     corresponding velocity in (erg/s)
!          Sv_a_lp1 = table containing the total stopping power
!                     of fast e- with kinetic energy eps(l+1) multiplied by their 
!                     corresponding velocity in (erg/s)
!          nu_a     = table containing the total angular collision rate of fast e-
!                     with kinetic energy eps(l) in (/s)
! output : phi_np1  = table containing the two first angular moment components of 
!                     the distribution function at time t(n+1) = t(n) + d_t
implicit none
real(PR), dimension(forward:backward,psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2), intent(in)  :: phi_n
real(PR), intent(in)                                                      :: d_t
real(PR), dimension(forward:backward,psi0:psi1z,1:N_x,1:N_z,1:N_eps), intent(in)           :: dF_dx, dF_dz, dF_deps
real(PR), dimension(forward:backward,psi0:psi1z,1:N_x,1:N_z,1:N_eps), intent(in)           :: gama_e, gama_b, gama_nu
real(PR), dimension(1:N_x,1:N_z,1:N_eps), intent(in)                      :: Sv_a, Sv_a_lp1, nu_a
real(PR), dimension(forward:backward,psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2), intent(out) :: phi_np1
integer                                                                   :: m, l, k, i, psi_mu
!$omp PARALLEL DO DEFAULT(SHARED) PRIVATE(m,l,k,i,psi_mu) COLLAPSE(5)
do m= forward,backward,1
do l=1,N_eps,1
do k=1,N_z,1      
do i=1,N_x,1
do psi_mu=psi0,psi1z,1                   
phi_np1(m,psi_mu,i,k,l) = phi_n(m,psi_mu,i,k,l) &
& - ( ( dF_dx(m,psi_mu,i,k,l) + dF_dz(m,psi_mu,i,k,l) ) * d_t )&
& - ( dF_deps(m,psi_mu,i,k,l) * d_t) &
& + ( ( gama_e(m,psi_mu,i,k,l) &
&     + gama_b(m,psi_mu,i,k,l) &
&     + gama_nu(m,psi_mu,i,k,l) ) * d_t)                   
end do
end do
end do
end do
end do
!$omp END PARALLEL DO
if (coll_implicit_scheme) then
! implicit resolution of collisional terms (implicit downwind scheme)
!$omp PARALLEL DO DEFAULT(SHARED) PRIVATE(m,k,i) COLLAPSE(3)
do m=forward,backward,1
do k=1,N_z,1
do i=1,N_x,1
do l=N_eps,1,-1
phi_np1(m,psi0,i,k,l)  = ( phi_np1(m,psi0,i,k,l)  &
&+ (Sv_a_lp1(i,k,l)*(fs/keV)*d_t*phi_np1(m,psi0,i,k,l+1)/d_eps) ) &
&/ ( 1._PR + (Sv_a(i,k,l)*(fs/keV) * d_t / d_eps) )
phi_np1(m,psi1x,i,k,l) = ( phi_np1(m,psi1x,i,k,l) &
&+ (Sv_a_lp1(i,k,l)*(fs/keV)*d_t*phi_np1(m,psi1x,i,k,l+1)/d_eps) ) &
&/ ( 1._PR + (Sv_a(i,k,l)*(fs/keV) * d_t / d_eps) &
&+ (nu_a(i,k,l) * fs * d_t) )
phi_np1(m,psi1z,i,k,l) = ( phi_np1(m,psi1z,i,k,l) &
&+ (Sv_a_lp1(i,k,l)*(fs/keV)*d_t*phi_np1(m,psi1z,i,k,l+1)/d_eps) ) &
&/ ( 1._PR + (Sv_a(i,k,l)*(fs/keV) * d_t / d_eps) &
&+ (nu_a(i,k,l) * fs * d_t) )                     
end do             
end do
end do
end do
!$omp END PARALLEL DO
end if     
end subroutine update

!=========================================================================================
!                                   Boundary conditions 
!           and computation of injected and escaping fast e-beam kinetic energies
!=========================================================================================

subroutine boundary_cond(norma, time, dt, epsa, Ue, Usf, Usb, Usu, Usd, phin)
! input  : norma = coefficient of normalization for the distribution function in (/cm^3/keV)
!          time  = time in (fs)
!          dt    = time step in (fs)
!          epsa  = array containing the kinetic energy grid in (keV)
! in(out)put : Ue = time integrated fast e- beam energy injected in the simulation box
!                   at z = 0 in (J)
!              Usf = time integrated fast e- beam energy escaping from the simulation box
!                    at z = L_z in (J)
!              Usb = time integrated fast e- beam energy escaping from the simulation box
!                    at z = 0 in (J)
!              Usu = time integrated fast e- beam energy escaping from the simulation box
!                    at x = L_x/2 in (J)
!              Usd = time integrated fast e- beam energy escaping from the simulation box
!                    at x = -L_x/2 in (J)
!              phin  = table containing the two first angular moments of the distribution 
!                      function in ()
implicit none
real(PR), intent(in)                                                                         :: norma
real(PR), intent(in)                                                                         :: time, dt
real(PR), dimension(1:N_eps), intent(in)                                                     :: epsa
real(PR), intent(inout)                                                                      :: Ue, Usf, Usb, Usu, Usd
real(PR), dimension(forward:backward,psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2), intent(inout) :: phin
! locals
real(PR), dimension(psi0:psi1z)                                                              :: phitemp, Fx_l, Fx_r, Fz_l, Fz_r
integer                                                                                      :: m, i, k, l, psi_mu
real(PR),dimension(:,:,:),allocatable                                                        :: usu_temp_x, usd_temp_x
real(PR),dimension(:,:),allocatable                                                          :: ue_temp_z, usf_temp_z, usb_temp_z
!
allocate(usu_temp_x(forward:backward,1:N_z,1:N_eps), usd_temp_x(forward:backward,1:N_z,1:N_eps),&
ue_temp_z(1:N_x,1:N_eps), usf_temp_z(1:N_x,1:N_eps), usb_temp_z(1:N_x,1:N_eps))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Energy boundary condition (kinetic energy indices l = 0,-1 and l = N_eps+1,N_eps+2) : !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp PARALLEL DO DEFAULT(SHARED) PRIVATE(m,k,i,psi_mu) COLLAPSE(4)    
do m = forward,backward,1
do k =-1,N_z+2,1
do i =-1,N_x+2,1
do psi_mu=psi0,psi1z,1
phin(m,psi_mu,i,k,0)  = phin(m,psi_mu,i,k,1)
phin(m,psi_mu,i,k,-1) = phin(m,psi_mu,i,k,0) 
phin(m,psi_mu,i,k,N_eps+1) = 0._PR 
phin(m,psi_mu,i,k,N_eps+2) = phin(m,psi_mu,i,k,N_eps+1)
end do
end do
end do
end do
!$omp END PARALLEL DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Space boundary conditions (transverse indices i = 0,-1 and i = N_x+1, N_x+2) : !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp PARALLEL DO DEFAULT(SHARED) PRIVATE(m,l,k,phitemp,Fx_l,Fx_r) COLLAPSE(3)
do m=forward,backward,1    
do l=1,N_eps,1
do k=-1,(N_z+2),1            
if (phin(m,psi1x,1,k,l) < 0._PR) then
phin(m,psi0:psi1z,0,k,l) = phin(m,psi0:psi1z,1,k,l)
else 
phin(m,psi0:psi1z,0,k,l) = 0._PR
end if             
phin(m,psi0:psi1z,-1,k,l) = phin(m,psi0:psi1z,0,k,l)              
if (phin(m,psi1x,N_x,k,l) > 0._PR) then
phin(m,psi0:psi1z,N_x+1,k,l) = phin(m,psi0:psi1z,N_x,k,l)
else 
phin(m,psi0:psi1z,N_x+1,k,l) = 0._PR
end if            
phin(m,psi0:psi1z,N_x+2,k,l) = phin(m,psi0:psi1z,N_x+1,k,l)   
! Energy diagnostics :          
if (k.ne.(N_z+1).and.(k.ne.(N_z+2)).and.(k.ne.-1).and.(k.ne.0)) then
phitemp(psi0:psi1z)=phin(m,psi0:psi1z,N_x+1,k,l)
call F_x(phitemp,Fx_r)
usu_temp_x(m,k,l) = Fx_r(psi0)*vit(epsa(l))*norma*(epsa(l)*keV)*d_eps*&
                  (dt*fs*d_z*sim_box_thickness)*(microns**2._PR)/Joules
phitemp(psi0:psi1z)=phin(m,psi0:psi1z,0,k,l)
call F_x(phitemp,Fx_l)
usd_temp_x(m,k,l) = -Fx_l(psi0)*vit(epsa(l))*norma*(epsa(l)*keV)*d_eps*&
                   (dt*fs*d_z*sim_box_thickness)*(microns**2._PR)/Joules
end if         
end do
end do
end do
!$omp END PARALLEL DO
usu = usu + sum(usu_temp_x(forward:backward,1:N_z,1:N_eps))
usd = usd + sum(usd_temp_x(forward:backward,1:N_z,1:N_eps))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Space boundary conditions (longitudinal indices k = N_z + 1, N_z+2 and k = 1,0) : !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp PARALLEL DO DEFAULT(SHARED) PRIVATE(l,i,phitemp,Fz_l,Fz_r) COLLAPSE(2)
do l=1,N_eps,1
do i=-1,N_x+2,1
! absorbing boundary condition at the rear side of the target for fast electrons propagating forward :                 
if (phin(forward,psi1z,i,N_z,l) > 0._PR) then
phin(forward,psi0:psi1z,i,N_z+1,l) = phin(forward,psi0:psi1z,i,N_z,l)
else 
phin(forward,psi0:psi1z,i,N_z+1,l) = 0._PR
end if
phin(forward,psi0:psi1z,i,N_z+2,l) = phin(forward,psi0:psi1z,i,N_z+1,l)
! absorbing boundary condition at the irradiated side of the target for fast electrons propagating forward : :
if (irradiated_side_refluxing .eqv. .false.) then
if ((time > source_duration) .and. (phin(forward,psi1z,i,1,l) < 0._PR)) then
phin(forward,psi0:psi1z,i,0,l) = phin(forward,psi0:psi1z,i,1,l)
else if ((time > source_duration) .and. (phin(forward,psi1z,i,1,l) > 0._PR)) then
phin(forward,psi0:psi1z,i,0,l) = 0._PR
end if
phin(forward,psi0:psi1z,i,-1,l) = phin(forward,psi0:psi1z,i,0,l)
end if
if (backward.eq.2) then
! if rear side refluxing : initialization of refluxing electrons propagating backward assuming a specular reflection of fast electrons propagating forward
if (phin(forward,psi1z,i,N_z,l) > 0._PR) then
phin(backward,psi0,i,N_z+1,l)  =  phin(forward,psi0,i,N_z,l)
phin(backward,psi1x,i,N_z+1,l) =  phin(forward,psi1x,i,N_z,l)
phin(backward,psi1z,i,N_z+1,l) = -phin(forward,psi1z,i,N_z,l)
else 
phin(backward,psi0:psi1z,i,N_z+1,l) = 0._PR
end if
phin(backward,psi0:psi1z,i,N_z+2,l) = phin(backward,psi0:psi1z,i,N_z+1,l)           
! absorbing boundary condition at the irradiated side of the target
if (phin(backward,psi1z,i,1,l) < 0._PR) then
phin(backward,psi0:psi1z,i,0,l) = phin(backward,psi0:psi1z,i,1,l)
else 
phin(backward,psi0:psi1z,i,0,l) = 0._PR
end if
phin(backward,psi0:psi1z,i,-1,l) = phin(backward,psi0:psi1z,i,0,l)
! absorbing boundary condition at the rear side of the target
if (phin(backward,psi1z,i,N_z,l) > 0._PR) then
phin(backward,psi0:psi1z,i,N_z+1,l) = phin(backward,psi0:psi1z,i,N_z+1,l) + phin(backward,psi0:psi1z,i,N_z,l)
end if
phin(backward,psi0:psi1z,i,N_z+2,l) = phin(backward,psi0:psi1z,i,N_z+1,l)
end if
! Energy diagnostics :
if ((i.ne.0).and.(i.ne.(-1)).and.(i.ne.(N_x+2)).and.(i.ne.(N_x+1))) then
phitemp(psi0:psi1z) = phin(forward,psi0:psi1z,i,0,l)      
call F_z(phitemp, Fz_r) 
ue_temp_z(i,l)  = max(0._PR, Fz_r(psi0)*vit(epsa(l))*norma*(epsa(l)*keV)*d_eps*&
                      (dt*fs*d_x)*sim_box_thickness*microns/Joules)
if (backward.eq.1) then
phitemp(psi0:psi1z) = phin(forward,psi0:psi1z,i,N_z,l)
call F_z(phitemp, Fz_r) 
usf_temp_z(i,l) = max(0._PR, Fz_r(psi0)*vit(epsa(l))*norma*(epsa(l)*keV)*d_eps*&
                      (dt*fs*d_x*sim_box_thickness)*microns/Joules)

phitemp(psi0:psi1z) = phin(forward,psi0:psi1z,i,1,l)      
call F_z(phitemp, Fz_l)
usb_temp_z(i,l) = max(0._PR,-Fz_l(psi0)*vit(epsa(l))*norma*(epsa(l)*keV)*d_eps*&
                      (dt*fs*d_x)*sim_box_thickness*microns/Joules)
else if ( (backward.eq.2) .and. (irradiated_side_refluxing .eqv. .false.) ) then
usf_temp_z(i,l) = 0._PR
phitemp(psi0:psi1z) = phin(forward,psi0:psi1z,i,1,l)      
call F_z(phitemp, Fz_l)
usb_temp_z(i,l) = max(0._PR,-Fz_l(psi0)*vit(epsa(l))*norma*(epsa(l)*keV)*d_eps*&
                      (dt*fs*d_x)*sim_box_thickness*microns/Joules)
phitemp(psi0:psi1z) = phin(backward,psi0:psi1z,i,1,l)      
call F_z(phitemp, Fz_l)
usb_temp_z(i,l) = usb_temp_z(i,l) + max(0._PR,-Fz_l(psi0)*vit(epsa(l))*norma*(epsa(l)*keV)*d_eps*&
                      (dt*fs*d_x)*sim_box_thickness*microns/Joules)
else if ( (backward.eq.2) .and. (irradiated_side_refluxing .eqv. .true.) ) then
usf_temp_z(i,l) = 0._PR
usb_temp_z(i,l) = 0._PR
end if  
!irradiated side refluxing : initialization of refluxing electrons propagating forward assuming a specular reflection of fast electrons propagating forward
if ( (irradiated_side_refluxing .eqv. .true.) .and. (phin(backward,psi1z,i,1,l) < 0._PR) )then
phin(forward,psi0,i,0,l)  = phin(forward,psi0,i,0,l)  + phin(backward,psi0,i,1,l)
phin(forward,psi1x,i,0,l) = phin(forward,psi1x,i,0,l) + phin(backward,psi1x,i,1,l)
phin(forward,psi1z,i,0,l) = phin(forward,psi1z,i,0,l) - phin(backward,psi1z,i,1,l)
if (phin(forward,psi0,i,0,l) < 0.) phin(forward,:,i,0,l) = 0._PR
phin(forward,psi0:psi1z,i,-1,l) = phin(forward,psi0:psi1z,i,0,l)
end if
end if
end do
end do
!$omp END PARALLEL DO
ue  = ue  + sum(ue_temp_z(1:N_x, 1:N_eps))
usf = usf + sum(usf_temp_z(1:N_x,1:N_eps))
usb = usb + sum(usb_temp_z(1:N_x,1:N_eps))
deallocate(usu_temp_x, usd_temp_x, ue_temp_z, usf_temp_z, usb_temp_z)
end subroutine boundary_cond

end module vfp