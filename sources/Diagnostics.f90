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
module diagnostics

use acuracy
use constants
use physics_library
use fokker_planck_coef
use transport_coef
use input
use initialization
use vfp

implicit none
private :: diagnose_energy, diagnose_hydro,diagnose_fb_x, diagnose_fb_z
private :: store_fb, store_energy, store_hydro
private :: Kioni_cross_section, Kalpha_diagnostic
public  :: initialize_diagnostics, plasma_diagnostics
public  :: Kalpha_emission, write_results, energy_balance

contains

!========================================================================================= 
!                              Kalpha and Kbeta diagnostics                    
!=========================================================================================

elemental function Kioni_cross_section(E_ioniK, U, J)
! input  : E_ioniK = K-shell electron ionisation energy in (eV)
!          U       = kinetic energy of the fast e- projectile normalized by E_ioniK in ()
!          J       = m_e c^2 / E_ioniK in ()
! output : Kioni_cross_section = K-shell ionization cross section in (cm^2)
!          according to Hombourger C., Jour. Phys. B 31, 16, 3693 (1998)
implicit none
real(PR),intent(in) :: E_ioniK , U , J
real(PR)            :: Kioni_cross_section
real(PR)            :: E_Rydberg , R_Bohr 
real(PR)            :: Gr, D, C
! Fundamental Hydrogen state energy in (eV) :
E_Rydberg = 13.606_PR     
! Bohr radius in (cm) :
R_Bohr = 5.2918e-9_PR
! Grysinski coefficient
Gr = (1._PR+(2._PR*J))*((U+J)**2._PR)*(((1._PR+U)*(U+(2._PR*J))*((1._PR+J)**2._PR))**1.5_PR) &
  &/ ((U+(2._PR*J))*((1._PR+J)**2._PR)*((J**2._PR)*(1._PR+(2._PR*J))+(U*(U+(2._PR*J))*(1.+(J**2_PR))))**1.5_PR)
D = ( 3.125_PR + ( (-4.172_PR) / U ) + ( 1.877_PR  / (U**2_PR) ) ) * (log(U) / U)
C = 2.0305_PR + ( (-0.3160_PR) / U ) + ( 0.1545_PR / (U**2_PR) )
Kioni_cross_section = 2._PR * Pi * ((R_Bohr)**2_PR) * Gr * D * ((E_Rydberg / E_ioniK)**C)
end function Kioni_cross_section

pure subroutine Kalpha_diagnostic(Kalphatab, norma, dt, Z, ni,&
                                  epsa, phi, nK_Holes, nKalpha, nKbeta, ionizrate)
! input  : Kalphatab = table containing data provided in the file 'sources/data/Kalpha_tab.dat' 
!          norma     = coefficient of normalization for the distribution function in (/cm^3/keV) 
!          dt        = time step in (fs)
!          Z         = atomic number of the material at a given location in ()
!          ni        = ion density in the material at a given location in (/cm^3)
!          epsa      = array containing the kinetic energy grid in (keV)
!          phi       = array containing the 0th order angular moment at a given location
!                      in space, that is, the fast e- spectrum at this location in ()
!          nK_Holes  = K-shell holes density in (/cm^3) at time t(n) and a given location
!          nKalpha   = Time integrated density of Kalpha (Kalpha1+Kalpha2) photons emitted
!                      per steradian in (/cm^3/sr) at time t(n) and a given location
!          nKbeta    = Time integrated density of Kbeta photons emitted per steradian 
!                      in (/cm^3/sr) at time t(n) and a given location
! output : nK_Holes  = K-shell holes density at time t(n+1) in (/cm^3)
!          nKalpha   = Time integrated density of Kalpha (Kalpha1+Kalpha2) photons emitted
!                      per steradian in (/cm^3/sr) at time t(n+1) and a given location
!          nKbeta    = Time integrated density of Kbeta photons emitted per steradian 
!                      in (/cm^3/sr) at time t(n+1) and a given location
!          ionizrate = K-shell electron ionization rate in (/s) at time t(n) and a given 
!             location according to A. Thomas et al., New Jour. Phys. 15, 1, 015017 (2013)
implicit none
real(PR), dimension(1:79,1:7),intent(in) :: Kalphatab
real(PR), intent(in)                     :: norma, dt, Z, ni
real(PR), dimension(1:N_eps), intent(in) :: epsa, phi
real(PR), intent(inout)                  :: nK_Holes, nKalpha, nKbeta
real(PR), intent(out)                    :: ionizrate
real(PR)                                 :: U, J
real(PR)                                 :: a1, a2, a3, b1, b2, b3
real(PR)                                 :: sigmaKalpha, Sigma_K, tau_K, gama, Wk
real(PR)                                 :: Fkalpha, Fkalpha1, Fkalpha2, Fkbeta 
real(PR)                                 :: nK_Holes_old, nK_Holes_mean
real(PR)                                 :: dnKalpha1, dnKalpha2, dnKbeta
integer                                  :: l, Zint
! K-shell hole life time tau_K (s) 
! according to A. Thomas et al., New Jour. Phys. 15, 1, 015017 (2013) :
a1 = 0.0002725_PR
a2 = 0.09932_PR
a3 = 2.160_PR
tau_K = hbar / ( exp( - a1 * (Z**2._PR) + a2 * Z - a3 ) * 1.e-3_PR * keV )
! Kalpha1, Kalpha2 and Kbeta contributions in the total K-shell fluoresence yield 
! (other contributions are neglected) :
Zint = int(Z)
FKalpha  = ( 1._PR + Kalphatab(Zint,7) )**(-1._PR)    !Kalpha_tab(7) = I_Kbeta/I_Kalpha
Fkalpha1 = ( 1._PR + Kalphatab(Zint,5) )**(-1._PR)    !Kalpha_tab(5) = I_Kalpha2/I_Kalpha1
Fkalpha2 = 1._PR - Fkalpha1                              
FKbeta   = 1._PR - FKalpha                               
! K-shell electron ionization rate :
J   = 511.e3_PR / Kalphatab(Zint,2)
Sigma_K=0._PR
do l=1,N_eps,1
U   = epsa(l) * 1.e3_PR / Kalphatab(Zint,2)                 !Kalpha_tab(2) = E_ioni_K (eV)
sigmaKalpha   = Kioni_cross_section(Kalphatab(Zint,2),U,J)        
Sigma_K = Sigma_K + ((phi(l)*norma) * vit(epsa(l)) * sigmaKalpha * d_eps) 
end do
ionizrate = Sigma_K
! Computation of hole density at time t(n+1) and 
! the average value nK_holes_mean between t(n) and t(n+1)
! according to A. Thomas et al., New Jour. Phys. 15, 1, 015017 (2013) :
nK_Holes_old = nK_Holes
gama = ionizrate + (1._PR/tau_K)
nK_Holes = nK_Holes_old &
&+ ( nK_Holes_old - (2 * ni * ionizrate / gama) ) * ( exp(-gama*dt) - 1._PR )
nk_Holes_mean = ( nK_Holes - nK_Holes_old + ( 2._PR * ni * ionizrate * dt * fs ) ) &
&/ ( gama * dt *fs )
! K-shell Fluorescence Yield 
! according to Kahoul et al., Rad. Phys. Chem., 80, 3, 369 (2011) :
b1 = 0.985_PR
b2 = 30.896_PR
b3 = 3.847_PR
Wk = b1 * ((Z/b2)**b3) / ( 1._PR + ((Z/b2)**b3) )
! Density of Kalpha1, Kalpha2 and Kbeta photons emitted between t(n) and t(n+1) in (/cm^3):
dnKalpha1 = ( ((FKalpha*FKalpha1*Wk/tau_K) * nK_Holes_mean) ) * dt * fs
dnKalpha2 = ( ((FKalpha*FKalpha2*Wk/tau_K) * nK_Holes_mean) ) * dt * fs
dnKbeta   = ( ((FKbeta          *Wk/tau_K) * nK_Holes_mean) ) * dt * fs
! Time integrated density of Kalpha and Kbeta Photons emitted per steradian :
nKalpha = nKalpha +  ((1._PR/(4._PR * pi)) * (dnKalpha1 + dnKalpha2))
nKbeta  = nKbeta  +  ((1._PR/(4._PR * pi)) * dnKbeta)
end subroutine Kalpha_diagnostic

subroutine Kalpha_emission(Kalpha_tab, norm, d_t, Z, ni, eps_a, phi_n,&
n_K_Holes, n_Kalpha, n_Kbeta, ioniz_rate) 
! input  : Kalphatab = table containing data provided in the file 'sources/data/Kalpha_tab.dat' 
!          norm      = coefficient of normalization for the distribution function in (/cm^3/keV) 
!          d_t        = time step in (fs)
!          Z         = table containing the atomic number in each spatial cell in ()
!          ni        = table containing the ion density in each spatial cell in (/cm^3)
!          epsa      = array containing the kinetic energy grid in (keV)
!          phi       = array containing the 0th and 1st order angular moment components
!                      of the distribution function in ()
!          nK_Holes  = table containing the K-shell holes density in (/cm^3) at time t(n) 
!                      in each spatial cell
!          nKalpha   = table containing the time integrated density of Kalpha 
!                      (Kalpha1+Kalpha2) photons emitted per steradian in (/cm^3/sr)
!                       at time t(n) in each spatial cell
!          nKbeta    = table containing the time integrated density of Kbeta photons 
!                      emitted per steradian in (/cm^3/sr) at time t(n) in each spatial cell
! output : nK_Holes  = table containing the K-shell holes density at time t(n+1)
!                      in each spatial cell
!          nKalpha   = table containing the time integrated density of Kalpha (Kalpha1+Kalpha2) 
!                      photons emittedper steradian in (/cm^3/sr) at time t(n+1) in each 
!                      spatial cell
!          nKbeta    = table containing the time integrated density of Kbeta photons 
!                      emitted per steradian in (/cm^3/sr) at time t(n+1) in each spatial cell
!          ionizrate = K-shell electron ionization rate in (/s) in each spatial cell
implicit none
real(PR), dimension(1:79,1:7), intent(in)                                                  :: Kalpha_tab
real(PR), intent(in)                                                                       :: norm, d_t
real(PR), dimension(0:N_x+1,0:N_z+1), intent(in)                                           :: Z, ni
real(PR), dimension(1:N_eps), intent(in)                                                   :: eps_a
real(PR), dimension(forward:backward,psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2), intent(in)  :: phi_n
real(PR), dimension(1:N_x,1:N_z), intent(inout)                                            :: n_K_Holes, n_Kalpha
real(PR), dimension(1:N_x,1:N_z), intent(inout)                                            :: n_Kbeta, ioniz_rate
! locals
integer                                                                                    :: m,k,i
real(PR), dimension(:), allocatable                                                        :: phi_temp
!
!$omp PARALLEL DO DEFAULT(SHARED) PRIVATE(k,i,phi_temp) COLLAPSE(2)
do k=1,N_z,1
do i=1,N_x,1
allocate(phi_temp(1:N_eps))
phi_temp(1:N_eps) = phi_n(forward,psi0,i,k,1:N_eps)
if (backward == 2) phi_temp(1:N_eps) = phi_temp(1:N_eps) + phi_n(backward,psi0,i,k,1:N_eps)
call Kalpha_diagnostic(Kalpha_tab, norm, d_t,&
Z(i,k), ni(i,k), eps_a, phi_temp,&
n_K_Holes(i,k), n_Kalpha(i,k), n_Kbeta(i,k), ioniz_rate(i,k))
deallocate(phi_temp)
end do
end do
!$omp END PARALLEL DO
end subroutine Kalpha_emission

!=========================================================================================
!     Write some material properties in txt files at the beginning of the simulation
!=========================================================================================

subroutine plasma_diagnostics(eta_tab)
  ! input : eta_tab   = tabulated resistivity if provided by the user 
  implicit none
  real(PR), dimension(1:N_eta_tab,1:2), intent(in) :: eta_tab
  integer                                          :: i
  real(PR)                                         :: T, Zf, Cv, res, cond_t
  real(PR)                                         :: S, nu, energy, om
  real(PR)                                         :: Z, A, ni
  call get_A_Z_ni(A,Z,ni)
  call system('mkdir -p results/'//trim(simu))
  open (unit=1,file ='results/'//trim(simu)//'/zeffvsTe[eV].dat'                           ,form='formatted',status='unknown')
  open (unit=2,file ='results/'//trim(simu)//'/electron_capacity[SI]vsTe[eV].dat'          ,form='formatted',status='unknown')
  open (unit=3,file ='results/'//trim(simu)//'/ion_capacity[SI]vsTe[eV].dat'               ,form='formatted',status='unknown')
  open (unit=4,file ='results/'//trim(simu)//'/resistivity[SI]vsTe[eV]_Te_eq_Ti.dat'       ,form='formatted',status='unknown')
  open (unit=7,file ='results/'//trim(simu)//'/resistivity[SI]vsTe[eV]_Ti_eq_Tamb.dat'     ,form='formatted',status='unknown')
  open (unit=8,file ='results/'//trim(simu)//'/conductivity[SI]vsTe[eV]_Te_eq_Ti.dat'      ,form='formatted',status='unknown')
  open (unit=9,file ='results/'//trim(simu)//'/conductivity[SI]vsTe[eV]_Ti_eq_Tamb.dat'    ,form='formatted',status='unknown')
  open (unit=10,file='results/'//trim(simu)//'/G[SI]vsTe[eV]_Te_eq_Ti.dat'                 ,form='formatted',status='unknown')
  open (unit=11,file='results/'//trim(simu)//'/G[SI]vsTe[eV]_Ti_eq_Tamb.dat'               ,form='formatted',status='unknown')
  open (unit=12,file='results/'//trim(simu)//'/ang_coll_rate[s-1]vsEps[keV].dat'           ,form='formatted',status='unknown')
  open (unit=13,file='results/'//trim(simu)//'/stopping_power[keV_Microns-1]vsEps[keV].dat',form='formatted',status='unknown')
  do i=-2000,4000
    T= (10._PR**(1.e-3_PR*i)) * eV 
    Zf = zeff(Z,ni,T)
    write(1,*) T / eV, Zf
    Cv   = capacity(Z,ni,T) * 1.E-1_PR               ! Cv[CGS] * 1.E-1_PR = Cv[SI]
    write(2,*) T / eV, Cv 
    Cv   = ion_capacity(Z,ni,T) * 1.E-1_PR           ! Cv[CGS] * 1.E-1_PR = Cv[SI]
    write(3,*) T / eV, Cv 
    res   = resis(eta_tab, Z,ni,T,T) * 9.E9_PR       ! res[CGS]*9.E9_PR = res[SI]
    write(4,*) T / eV, res 
    res   = resis(eta_tab, Z,ni, T, Tamb) * 9.E9_PR  ! res[CGS]*9.E9_PR = res[SI]
    write(7,*) T / eV, res 
    cond_t   = cond(Z,ni,T,T)
    write(8,*) T / eV, cond_t * 1.E-5_PR            ! cond_t[CGS] * 1.E-5_PR = cond_t[SI]
    cond_t   = cond(Z,ni,T,Tamb)
    write(9,*) T / eV, cond_t * 1.E-5_PR            ! cond_t[CGS] * 1.E-5_PR = cond_t[SI]
    om = Omega_ei(A, Z, ni, T, T)
    write(10,*) T / eV, om * 1.E-1_PR               ! om[CGS] * 1.E-1_PR = om[SI]
    om = Omega_ei(A, Z, ni, T, Tamb)
    write(11,*) T / eV, om                          ! om[CGS] * 1.E-1_PR = om[SI]
  end do
  do i=1,100
    energy = (10_PR*(i**2)) 
    nu   = nu_tot(Z, ni, Tamb, Tamb, energy)
    write(12,*) energy , nu
    S   = S_tot(A, Z, ni, Tamb, Tamb, energy) * microns / keV
    write(13,*) energy , S
  end do
  close(1)
  close(2)
  close(3)
  close(4)
  close(7)
  close(8)
  close(9)
  close(10)
  close(11)
  close(12)
  close(13)
end subroutine plasma_diagnostics

!=========================================================================================
!     Write some electron properties in txt files at the beginning of the simulation
!=========================================================================================

subroutine electron_diagnostics(eps_tab,xa,epsa)
implicit none
! inputs
real(PR), dimension(1:N_eps_tab,1:2), intent(in) :: eps_tab
real(PR), dimension(1:N_x), intent(in)   :: xa
real(PR), dimension(1:N_eps), intent(in) :: epsa
! locals
integer                              :: i, k, l, n
real(PR)                             :: delta, omega0, a2, theta0, alpha, alpha_x, alpha_z
real(PR), dimension(1:N_eps)         :: dN_dE
real(PR), dimension(1:N_x)           :: dN_dx
integer, parameter                   :: N_theta   = 360
real(PR), parameter                  :: theta_max = 90., dtt = 0.25
integer                              :: N_tt
real(PR), dimension(:),allocatable   :: tt, dN_dtt
real(PR), dimension(1:N_theta)       :: theta
real(PR), dimension(1:N_x,1:N_theta) :: dN_dtheta
!
! Create and or Open Files where data will be stored 
open (unit=1,file ='results/'//trim(simu)//'/fast_electron_spectrum.dat'     ,form='formatted',status='unknown')
open (unit=2,file ='results/'//trim(simu)//'/fast_electron_angular_distr.dat',form='formatted',status='unknown')
open (unit=3,file ='results/'//trim(simu)//'/fast_electron_spatial_distr.dat',form='formatted',status='unknown')
open (unit=4,file ='results/'//trim(simu)//'/fast_electron_temporal_distr.dat',form='formatted',status='unknown')
! Fast Electron Kinetic Energy Spectrum
delta = epsa(2) - epsa(1)
do l=1,N_eps,1
dN_dE(l) = function_energy(eps_tab,epsa(l))
end do
dN_dE = dN_dE / (sum(dN_dE(1:N_eps))*delta)
do l=1,N_eps,1
write(1,*) epsa(l), dN_dE(l)
end do
! Fast Electron Angular Distribution
delta = 2.*theta_max / N_theta
if (sigma_theta.ne.0._PR) then
a2 = 1._PR / ((pi*sigma_theta/180._PR)**2._PR)
omega0 = (1._PR/tanh(a2)) - (1._PR /a2)
else
omega0 = 1._PR
end if
do i=1,N_x,1
theta0 = function_theta(xa(i))
alpha   = 3. * omega0               / ( 1. - ((0.5*(omega0**2.))*(1.+(omega0**2.))) )
alpha_z = 3. * omega0 * cos(theta0) / ( 1. - ((0.5*(omega0**2.))*(1.+(omega0**2.))) )
alpha_x = 3. * omega0 * sin(theta0) / ( 1. - ((0.5*(omega0**2.))*(1.+(omega0**2.))) )
do k=1,N_theta,1
theta(k) = - theta_max + (real(k)-1.)*delta
theta(k) = pi*theta(k)/180.
dN_dtheta(i,k) = alpha * exp( ( alpha_x * sin(theta(k)) ) + ( alpha_z * cos(theta(k)) ) ) / ( 4.*pi * sinh(alpha) )
theta(k) = 180.*theta(k)/pi
end do
dN_dtheta(i,1:N_theta) = dN_dtheta(i,1:N_theta) / ( sum(dN_dtheta(i,1:N_theta))*delta )
end do
do i=1,N_x,1
do k=1,N_theta,1
write(2,*) xa(i), theta(k), dN_dtheta(i,k)
end do
end do
! Fast Electron Spatial Distribution
delta = xa(2) - xa(1)
do i=1,N_x,1
dN_dx(i) = function_space(xa(i))
end do
dN_dx = dN_dx / (sum(dN_dx(1:N_x))*delta)
do i=1,N_x,1
write(3,*) xa(i), dN_dx(i)
end do
! Fast Electron Temporal Distribution
N_tt = int(L_t / dtt)
allocate(tt(1:N_tt),dN_dtt(1:N_tt))
do n=1,N_tt,1
tt(n) = (real(n)-1.)*dtt
dN_dtt(n) = function_time(1._PR, tt(n))
end do
dN_dtt = dN_dtt / sum(dN_dtt(:))
do n=1,N_tt,1
write(4,*) tt(n), dN_dtt(n)
end do
deallocate(tt,dN_dtt)
! close the files
close(1)
close(2)
close(3)
close(4)
end subroutine electron_diagnostics

!=========================================================================================
!                  Open files for storing the computation results and 
!  store the needed data for the computation of the emission of Kalpha and Kbeta photons
!=========================================================================================

subroutine initialize_diagnostics(Kalphatab)
! output : Kalphatab = table containing data provided in the file 'sources/data/Kalpha_tab.dat'
!          array 1 corresponds to the atomic number Z () 
!          array 2 to the K-shell electron ionisation energy E_ionization (eV) 
!          array 3 to the Kalpha1 photon energy E_Kalpha1 (eV) if emitted
!          array 4 to the Kalpha2 photon energy E_Kalpha2 (eV)  
!          array 5 to the ratio IKalpha2/IKalpha1 of single photon signal intensities
!          array 6 to the Kbeta photon energy E_Kbeta (eV) and
!          array 7 to the ratio IKbeta/IKalpha of single photon signal intensities.
!          Data from Hydrogen H (Z=1) to Gold Au (Z=79) are provided according to
!          Thomson A. et al., X-ray data booklet (2009)
implicit none
real(PR), dimension(1:79,1:7), intent(out) :: Kalphatab
integer                                    :: i
open (unit=10 ,file='results/'//trim(simu)//'/psi0_z[cm-3_keV-1].dat'     ,form='formatted',status='unknown')
open (unit=11 ,file='results/'//trim(simu)//'/psi1x_z[cm-3_keV-1].dat'    ,form='formatted',status='unknown')
open (unit=12,file='results/'//trim(simu)//'/psi1z_z[cm-3_keV-1].dat'     ,form='formatted',status='unknown')
open (unit=13,file='results/'//trim(simu)//'/psi0_x[cm-3_keV-1].dat'      ,form='formatted',status='unknown')
open (unit=14,file='results/'//trim(simu)//'/psi1x_x[cm-3_keV-1].dat'     ,form='formatted',status='unknown')
open (unit=15,file='results/'//trim(simu)//'/psi1z_x[cm-3_keV-1].dat'     ,form='formatted',status='unknown')
open (unit=20,file='results/'//trim(simu)//'/nb[cm-3].dat'                ,form='formatted',status='unknown')
open (unit=21,file='results/'//trim(simu)//'/jb_x[A_cm-2].dat'            ,form='formatted',status='unknown')
open (unit=22,file='results/'//trim(simu)//'/jb_z[A_cm-2].dat'            ,form='formatted',status='unknown')
open (unit=23,file='results/'//trim(simu)//'/E_x[V_m-1].dat'              ,form='formatted',status='unknown')
open (unit=24,file='results/'//trim(simu)//'/E_z[V_m-1].dat'              ,form='formatted',status='unknown')
open (unit=25,file='results/'//trim(simu)//'/B_y[Tesla].dat'              ,form='formatted',status='unknown')
open (unit=26,file='results/'//trim(simu)//'/We[erg_s-1_cm-3].dat'        ,form='formatted',status='unknown')
open (unit=27,file='results/'//trim(simu)//'/Te[eV].dat'                  ,form='formatted',status='unknown')
open (unit=28,file='results/'//trim(simu)//'/Wi[erg_s-1_cm-3].dat'        ,form='formatted',status='unknown')
open (unit=29,file='results/'//trim(simu)//'/Ti[eV].dat'                  ,form='formatted',status='unknown')
open (unit=30,file='results/'//trim(simu)//'/resis[Ohm.m].dat'            ,form='formatted',status='unknown')
open (unit=31,file='results/'//trim(simu)//'/Kappa_e[erg_m-1_K-1_s-1].dat',form='formatted',status='unknown')
open (unit=32,file='results/'//trim(simu)//'/K_shell_ioniz_rate_[s-1].dat',form='formatted',status='unknown')
open (unit=33,file='results/'//trim(simu)//'/n_Kalpha[cm-3].dat'          ,form='formatted',status='unknown')
open (unit=34,file='results/'//trim(simu)//'/n_Kbeta[cm-3].dat'           ,form='formatted',status='unknown')
open (unit=35,file='results/'//trim(simu)//'/ni[cm-3].dat'                ,form='formatted',status='unknown')
open (unit=40,file='results/'//trim(simu)//'/U_e[J].dat'                  ,form='formatted',status='unknown')
open (unit=41,file='results/'//trim(simu)//'/U_b[J].dat'                  ,form='formatted',status='unknown')
open (unit=42,file='results/'//trim(simu)//'/Ud_col[J].dat'               ,form='formatted',status='unknown')
open (unit=43,file='results/'//trim(simu)//'/Ud_res[J].dat'               ,form='formatted',status='unknown')
open (unit=44,file='results/'//trim(simu)//'/U_el[J].dat'                 ,form='formatted',status='unknown')
open (unit=45,file='results/'//trim(simu)//'/U_ma[J].dat'                 ,form='formatted',status='unknown')
open (unit=46,file='results/'//trim(simu)//'/U_sf[J].dat'                 ,form='formatted',status='unknown')
open (unit=47,file='results/'//trim(simu)//'/U_sb[J].dat'                 ,form='formatted',status='unknown')
open (unit=48,file='results/'//trim(simu)//'/U_su[J].dat'                 ,form='formatted',status='unknown')
open (unit=49,file='results/'//trim(simu)//'/U_sd[J].dat'                 ,form='formatted',status='unknown')
open (unit=50,file='results/'//trim(simu)//'/ne[cm-3].dat'                ,form='formatted',status='unknown')
open (unit=52,file='results/'//trim(simu)//'/je_x[A_cm-2].dat'            ,form='formatted',status='unknown')
open (unit=53,file='results/'//trim(simu)//'/je_z[A_cm-2].dat'            ,form='formatted',status='unknown')
! Open, read and store the file 'sources/data/Kalpha_tab.dat' in the table Kalphatab :
open (unit=51,file='sources/data/Kalpha_tab.dat'          ,form='formatted',status='unknown')      
read(51,*)
do i=1,79,1
read(51,*) Kalphatab(i,1),Kalphatab(i,2),Kalphatab(i,3),Kalphatab(i,4),&
Kalphatab(i,5),Kalphatab(i,6),Kalphatab(i,7)
end do
close(51)
end subroutine initialize_diagnostics

!=========================================================================================
!                Write the different energy contribution in txt files
!=========================================================================================

subroutine store_energy(N_file, t, u)
! input : N_file = file unit
!         t      = time in (fs)
!         u      =  energy in (J)
implicit none
integer, intent(in)  :: N_file
real(PR), intent(in) :: t, u
if (abs(u) <= zero) then
write(N_file,'(2E23.15)') t, 0._PR
else
write(N_file,'(2E23.15)') t, u
end if
end subroutine store_energy

subroutine diagnose_energy(t, Ue, Ub, Udcol, Udres, Uel, Uma, Usf, Usb, Usu, Usd)
! input : t     = time in (fs)
!         Ue    = time integrated fast e- beam energy injected in the simulation box at 
!                 time t in (J)
!         Ub    = instantaneous beam energy in the simulation box at time t in (J)
!         Udcol = time integrated energy lost by fast e- due to collisions at time t in (J) 
!         Udres = time integrated energy lost by fast e- due to electric fields at time t 
!                 in (J)
!         Uel   = instantaneous electric energy in the simulation box at time t in (J)
!         Uma   = instantaneous magnetic energy in the simulation box at time t in (J)
!         Usf   = time integrated fast e- beam energy escaping from the simulation box
!                 at z = L_z at time t in (J)
!         Usb   = time integrated fast e- beam energy escaping from the simulation box
!                 at z = 0 at time t in (J)
!         Usu   = time integrated fast e- beam energy escaping from the simulation box
!                 at x = +L_x/2 at time t in (J)
!         Usd   = time integrated fast e- beam energy escaping from the simulation box
!                 at x = -L_x/2 at time t in (J)
implicit none
real(PR), intent(in) :: t, Ue, Ub, Udcol, Udres, Uel, Uma, Usf, Usb, Usu, Usd 
call store_energy(40,t,Ue)
call store_energy(41,t,Ub)
call store_energy(42,t,Udcol)
call store_energy(43,t,Udres)
call store_energy(44,t,Uel)
call store_energy(45,t,Uma)
call store_energy(46,t,Usf)
call store_energy(47,t,Usb)
call store_energy(48,t,Usu)
call store_energy(49,t,Usd)
end subroutine diagnose_energy

!=========================================================================================
!                Compute the different energy contributions at each time step
!               and print in the console the energy balance every Delta_t_diag
!=========================================================================================

subroutine energy_balance(diag_condition, N_t, time, d_t, norm, eps_a, phi_n,&
  p_depos, p_lost_col, p_lost_res, E_x, E_z, B_y_n,&
  U_e, U_b, U_d, U_d_col, U_d_res, U_el, U_ma,&
  U_sf, U_sb, U_su, U_sd)
  ! input :      diag_condition = logical ensuring that the energy balance is printed
  !                               in the console every Delta_t_diag
  !              N_t            = time iteration
  !              time           = time in (fs)
  !              d_t            = time step in (fs)
  !              norm           = coefficient of normalization for the distribution function
  !                               in (/cm^3/keV)
  !              eps_a          = array containing the kinetic energy grid in (keV)
  !              phi_n          = table containing the 0th and 1st order angular moment 
  !                               components of the distribution function in ()
  !              p_depos        = table containing the density of power deposited on 
  !                               at t = time in (erg/cm^3/s)
  !              p_lost_col    = table containing the density of power lost by the fast e-
  !                              beam due to collisions in (erg/cm^3/s)
  !              p_lost_res    = table containing the density of power lost by the fast e-
  !                              beam due to the self-generated electric field 
  !                              in (erg/cm^3/s)
  !              E_x           = table containing the self-generated electric field 
  !                              component on the x-axis in (statVolt/cm)
  !              E_z           = table containing the self-generated electric field 
  !                              component on the z-axis in (statVolt/cm)
  !              B_y_n         = table containing the self-generated magnetic field 
  !                              component on the x-axis in (gauss)
  ! in(out)put : Ue    = injected fast e- beam energy at time t in (J)
  !              Ub    = instantaneous beam energy in the simulation box at time t in (J)
  !              Udcol = time integrated energy lost by fast e- due to collisions at t=time
  !                      in (J) 
  !              Udres = time integrated energy lost by fast e- due to electric fields at 
  !                      t = time in (J)
  !              Uel   = instantaneous electric energy in the simulation box at t=time in (J)
  !              Uma   = instantaneous magnetic energy in the simulation box at t=time in (J)
  !              Usf   = time integrated fast e- beam energy escaping from the simulation box
  !                      at z = L_z at time t in (J)
  !              Usb   = time integrated fast e- beam energy escaping from the simulation box
  !                      at z = 0 at t=time in (J)
  !              Usu   = time integrated fast e- beam energy escaping from the simulation box
  !                      at x = +L_x/2 at t=time in (J)
  !              Usd   = time integrated fast e- beam energy escaping from the simulation box
  !                      at x = -L_x/2 at t=time in (J)
  implicit none
  logical, intent(in)                                                      :: diag_condition
  integer, intent(in)                                                      :: N_t
  real(PR), intent(in)                                                     :: time, d_t, norm
  real(PR), dimension(1:N_eps), intent(in)                                 :: eps_a
  real(PR), dimension(forward:backward,psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2), intent(in) :: phi_n
  real(PR), dimension(1:N_x,1:N_z), intent(in)                             :: p_depos, p_lost_col, p_lost_res
  real(PR), dimension(1:N_x,1:N_z), intent(in)                             :: E_x, E_z, B_y_n
  real(PR), intent(inout)                                                  :: U_e, U_b, U_d
  real(PR), intent(inout)                                                  :: U_d_col, U_d_res, U_el, U_ma
  real(PR), intent(inout)                                                  :: U_sf, U_sb, U_su, U_sd
  integer                                                                  :: m,l
  real(PR)                                                                 :: U_s, error
  do m=forward,backward,1
    do l=1,N_eps,1
      U_b     = U_b + sum(phi_n(m,psi0,1:N_x,1:N_z,l)) * d_x * sim_box_thickness * d_z * (microns**2_PR) &
      * (eps_a(l)*keV) * d_eps * norm / Joules 
    end do   
  end do
  U_d     = U_d     + sum(p_depos)    * d_x *  sim_box_thickness*d_z*((microns)**2._PR)*d_t*fs/Joules
  U_d_col = U_d_col + sum(p_lost_col) * d_x *  sim_box_thickness*d_z*((microns)**2._PR)*d_t*fs/Joules
  U_d_res = U_d_res + sum(p_lost_res) * d_x *  sim_box_thickness*d_z*((microns)**2._PR)*d_t*fs/Joules
  U_el    = (sum(E_x**2._PR+E_z**2._PR)/ 2._PR)*d_x*sim_box_thickness*d_z*((microns)**2._PR)/Joules 
  U_ma    = (sum(B_y_n**2._PR)/ 2._PR) *d_x*sim_box_thickness*d_z*((microns)**2._PR)/Joules 
  U_s     = U_sf+U_sb+U_su+U_sd
  error   = 100._PR*((U_e-(U_b+U_d_col+U_d_res+U_s+U_el+U_ma))/U_e)
  ! Print in the console :
  if (diag_condition) then
    write(*,*)'=========================='
    write(*,'(A,1E14.7,A)')' t =', time, ' fs'
    write(*,'(A,1I12,A)')' ( iteration :',N_t,' )'
    write(*,*)'=========================='
    write(*,'(A,1E14.7,A)')' * Time integrated injected energy =', U_e, ' J'
    write(*,'(A,1E14.7,A)')' * Instantaneous beam energy       =', U_b, ' J'
    write(*,'(A,1E14.7,A)')' * Energy deposited :              =', U_d, ' J'
    write(*,'(A,1E14.7,A)')'    * due to collisions            =', U_d_col, ' J'
    write(*,'(A,1E14.7,A)')'    * due to Ohmic heating         =', U_d_res, ' J'
    write(*,'(A,1E14.7,A)')' * Instantaneous electric energy   =', U_el, ' J'
    write(*,'(A,1E14.7,A)')' * Instantaneous magnetic energy   =', U_ma, ' J'
    write(*,'(A,1E14.7,A)')' * Time integrated escaped energy  =', U_sf, ' J'
    write(*,'(A,1E14.7,A)')'    * forward                      =', U_sf, ' J'
    write(*,'(A,1E14.7,A)')'    * backward                     =', U_sb, ' J'
    write(*,'(A,1E14.7,A)')'    * upward                       =', U_su, ' J'
    write(*,'(A,1E14.7,A)')'    * downward                     =', U_sd, ' J'
    write(*,'(A,1F6.2,A)')' -> ENERGY CONSERVATION ERROR      =', error, ' %'
  end if
  ! write in txt files :
  call diagnose_energy(time, U_e, U_b, U_d_col, U_d_res, U_el, U_ma, U_sf, U_sb, U_su, U_sd)
end subroutine energy_balance

!=========================================================================================
!                   Write the two first order angular moment components
!                      and the main physical quantities in txt files
!=========================================================================================

subroutine store_fb(N_file, t, x, z, nrj, f0)
! input : N_file = file unit
!         t      = time in (fs)
!         z      = longitudinal location z in (microns)
!         nrj    = kinetic energy value in (keV)
!         x      = transverse location x in (microns) at which the fast e- beam density 
!                  has the maximum value for a given depth z
!         f0     = psi_0, Psi_1x or Psi_1z in (/cm^3/keV)
implicit none
integer, intent(in)  :: N_file
real(PR), intent(in) :: t, x, z, nrj
real(PR), intent(in) :: f0
if (abs(f0) <= zero) then
write(N_file,'(5E23.15)') t, z, nrj, x, 0._PR
else
write(N_file,'(5E23.15)') t, z, nrj, x, f0
end if
end subroutine store_fb

subroutine diagnose_fb_z(norma, t, x, z, nrj, f)
! input : norma = coefficient of normalization for the distribution function in (/cm^3/keV) 
!         t     = time in (fs)
!         x     = transverse location x in (microns) at which the fast e- beam density 
!                 has the maximum value for a given depth z
!         z     = longitudinal location z in (microns)
!         nrj   = kinetic energy value in (keV)
!         f     = array containing the 0th and 1st order angular moment components
!                 (psi_0 Psi_1x Psi_1z)^T
implicit none
real(PR), intent(in)                        :: t, x, z, nrj
real(PR), dimension(psi0:psi1z), intent(in) :: f
real(PR), intent(in)                        :: norma
real(PR)                                    :: f0
f0 = f(psi0)*norma  
call store_fb(10, t, x, z, nrj, f0)
f0 = f(psi1x)*norma
call store_fb(11, t, x, z, nrj, f0)
f0 = f(psi1z)*norma
call store_fb(12, t, x, z, nrj, f0)
end subroutine diagnose_fb_z

subroutine diagnose_fb_x(norma, t, x, z, nrj, f)
! input : norma = coefficient of normalization for the distribution function in (/cm^3/keV) 
!         t     = time in (fs)
!         x     = transverse location x in (microns) 
!         z     = longitudinal location z in (microns) at which the fast e- beam density 
!                 has the maximum value for a given transverse location x
!         nrj   = kinetic energy value in (keV)
!         f     = array containing the 0th and 1st order angular moment components
!                 (psi_0 Psi_1x Psi_1z)^T
implicit none
real(PR), intent(in)                 :: t, x, z, nrj
real(PR), dimension(1:3), intent(in) :: f
real(PR), intent(in)                 :: norma
real(PR)                             :: f0
f0 = f(psi0)*norma 
call store_fb(13, t, z, x, nrj, f0)
f0 = f(psi1x)*norma
call store_fb(14, t, z, x, nrj, f0)
f0 = f(psi1z)*norma
call store_fb(15, t, z, x, nrj, f0)
end subroutine diagnose_fb_x

subroutine store_hydro(N_file,t,x,z,f)
! input : t = time in (fs)
!         x = transverse position in (microns)
!         z = longitudinal posiition in (microns)
!         f = physical quantity to store
implicit none
integer, intent(in)  :: N_file
real(PR), intent(in) :: t,x,z,f
if (abs(f) <= zero) then
write(N_file,'(4E23.15)') t, x, z, 0._PR
else
write(N_file,'(4E23.15)') t, x, z, f
end if
end subroutine store_hydro

subroutine diagnose_hydro(eta_tab, A, Zn, ni, nb, jbx, jbz, jrx, jrz, Ex, Ez, By, &
pdepos, Te, Ti, nKalpha, nKbeta, ionizrate,&
t, x, z)
! input : eta_tab   = tabulated resistivity if provided by the user 
!                     (tabulated_resistivity = .true.)
!         A         = Atomic weight at a given location  in ()
!         Zn        = Atomic number at a given location in ()
!         ni        = ion density at a given location in (/cm^3)
!         nb        = beam density at time t at a given location in (/cm^3)
!         jbx       = beam current density component on x-axis at time t at a given location  
!                     in (statAmpere/cm^2)
!         jbz       = beam current density component on z-axis at time t at a given location 
!                     in (statAmpere/cm^2)
!         Ex        = Self-generated electric field on x-axis at time t at a given location 
!                     in (statvolt/cm)
!         Ez        = Self-generated electric field on z-axis at time t at a given location 
!                     in (statvolt/cm)
!         By        = Self-generated magnetic field on y-axis at time t at a given location 
!                     in (gauss)
!         pdepos    = density of power deposited on background electrons at time t 
!                     at a given location in (erg/cm^3/s)
!         Te        = background electron temperature at time t at a given location in (K)
!         Ti        = background ion temperature at time t at a given location in (K)
!         nKalpha   = Time integrated density of Kalpha photons (Kalpha1+Kalpha2)
!                     emitted at time t at a given location in (/cm^3/sr)
!         nKbeta    = Time integrated density of Kbeta photons emitted at time t at a
!                     given location in (/cm^3/sr)
!         ionizrate = K-shell electron ionization rate in (/s) at time t and a given 
!                     location in (/s)
!         t         = time in (fs)
!         x         = transverse position in (microns)
!         z         = longitudinal position in (microns)
implicit none
real(PR), dimension(1:N_eta_tab,1:2), intent(in) :: eta_tab
real(PR), intent(in)                             :: A, Zn, ni
real(PR), intent(in)                             :: nb, jbx, jbz, jrx, jrz, Ex, Ez, By
real(PR), intent(in)                             :: pdepos, Te, Ti
real(PR), intent(in)                             :: nKalpha, nKbeta, ionizrate, t, x, z
real(PR)                                         :: f
call store_hydro(20,t,x,z,nb)
call store_hydro(21,t,x,z,jbx/3.E9_PR) ! jbx(statAmpere/cm^2)/3.E9_PR = jbx(Ampere/cm^2)
call store_hydro(22,t,x,z,jbz/3.E9_PR) ! jbz(statAmpere/cm^2)/3.E9_PR = jbz(Ampere/cm^2)
call store_hydro(23,t,x,z,Ex*3.e4_PR)  ! Ex(statvolt/cm)*3.e4_PR = Ex(volt/m)
call store_hydro(24,t,x,z,Ez*3.e4_PR)  ! Ez(statvolt/cm)*3.e4_PR = Ez(volt/m)
call store_hydro(25,t,x,z,By/1.e4_PR)  ! By(gauss)/1.e4_PR = By(T)
call store_hydro(26,t,x,z,pdepos)
call store_hydro(27,t,x,z,Te/eV)       ! Te(K)/eV = Te(eV)
f = Omega_ei(A, Zn, ni, Te, Ti)*(Te - Ti)
call store_hydro(28,t,x,z,f)
call store_hydro(29,t,x,z,Ti/eV)       ! Ti(K)/eV = Ti(eV)
f = resis(eta_tab, Zn,ni,Te,tamb) * 9.E9_PR
call store_hydro(30,t,x,z,f)           ! eta(s)*9.E9_PR = eta(Ohm.m)
f = cond(Zn, ni, Te, Ti)
call store_hydro(31,t,x,z,f)
call store_hydro(32,t,x,z,ionizrate)
call store_hydro(33,t,x,z,nKalpha)
call store_hydro(34,t,x,z,nKbeta)
call store_hydro(35,t,x,z,ni)
f = zeff(Zn, ni, Te)*ni
call store_hydro(50,t,x,z,f)
call store_hydro(52,t,x,z,jrx/3.E9_PR) ! jex(statAmpere/cm^2)/3.E9_PR = jbx(Ampere/cm^2)
call store_hydro(53,t,x,z,jrz/3.E9_PR) ! jez(statAmpere/cm^2)/3.E9_PR = jbx(Ampere/cm^2)
end subroutine diagnose_hydro

subroutine write_results(eta_tab, A, Z, ni, n_b, j_b_x_n, j_b_z_n,&
 j_r_x_n, j_r_z_n, E_x, E_z, B_y_n, &
p_depos, Te, Ti,&
n_Kalpha, n_Kbeta, ioniz_rate,&
phi_n, norm, time, x_a, z_a, eps_a)
! input : see subroutines diagnose_hydro, diagnose_fb_x and diagnose_fb_z
implicit none
real(PR), dimension(1:N_eta_tab,1:2), intent(in)                                          :: eta_tab
real(PR), dimension(0:N_x+1,0:N_z+1), intent(in)                                          :: A, Z, ni, Te, Ti
real(PR), dimension(1:N_x,1:N_z), intent(in)                                              :: n_b, j_b_x_n, j_b_z_n
real(PR), dimension(1:N_x,1:N_z), intent(in)                                              :: j_r_x_n, j_r_z_n
real(PR), dimension(1:N_x,1:N_z), intent(in)                                              :: E_x, E_z, B_y_n
real(PR), dimension(1:N_x,1:N_z), intent(in)                                              :: n_Kalpha, n_Kbeta
real(PR), dimension(1:N_x,1:N_z), intent(in)                                              :: p_depos, ioniz_rate
real(PR), dimension(forward:backward,psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2), intent(in) :: phi_n
real(PR), intent(in)                                                                      :: norm, time
real(PR), dimension(1:N_x), intent(in)                                                    :: x_a
real(PR), dimension(1:N_z), intent(in)                                                    :: z_a
real(PR), dimension(1:N_eps), intent(in)                                                  :: eps_a
! locals
integer                                                                                   :: i,k,l
integer, dimension(1:N_z)                                                                 :: i_nbmax
integer, dimension(1:N_x)                                                                 :: k_nbmax
real(PR), dimension(psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2)                              :: phi_temp
!
phi_temp(psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2) = phi_n(forward,psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2)
if (backward == 2) then
phi_temp(psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2) = phi_temp(psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2)&
                                                  + phi_n(backward,psi0:psi1z,-1:N_x+2,-1:N_z+2,-1:N_eps+2)
end if
! Plasma and Beam Densities (they are printed into files every Delta_t_diag)
do i=1,N_x,1
do k=1,N_z,1
call diagnose_hydro(eta_tab, A(i,k),Z(i,k), ni(i,k), n_b(i,k), j_b_x_n(i,k), j_b_z_n(i,k),&
j_r_x_n(i,k), j_r_z_n(i,k), E_x(i,k), E_z(i,k), B_y_n(i,k), &
p_depos(i,k), Te(i,k), Ti(i,k),&
n_Kalpha(i,k), n_Kbeta(i,k), ioniz_rate(i,k),&
time, x_a(i), z_a(k))
end do
end do       
! Electron Beam Kinetics 
!(they are printed into files every Delta_t_diag steptimes and :
! _only at the position x where the electron beam density is maximum for the given
!  position z -> diagnose_fb_z
! _only at the position z where the electron beam density is maximum for the given
!  position x -> diagnose_fb_x)
do k=1,N_z,1
! search for the index at which the beam density is maximum :
i_nbmax(k) = 1
do i=2,N_x,1
if (n_b(i,k) >= n_b(i_nbmax(k),k)) i_nbmax(k)=i
end do          
do l=1,N_eps,1    
! write :      
call diagnose_fb_z(norm, time, &
x_a(i_nbmax(k)), z_a(k), eps_a(l),&
phi_temp(psi0:psi1z,i_nbmax(k),k,l))
end do     
end do
do i=1,N_x,1
! search for the index at which the beam density is maximum :           
k_nbmax(i) = 1
do k=2,N_z,1
if (n_b(i,k) >= n_b(i,k_nbmax(i))) k_nbmax(i)=k
end do        
! write :        
do l=1,N_eps,1          
call diagnose_fb_x(norm, time, &
x_a(i), z_a(k_nbmax(i)), eps_a(l),&
phi_temp(psi0:psi1z,i,k_nbmax(i),l))
end do     
end do 
end subroutine write_results

end module diagnostics
