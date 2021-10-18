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
module physics_library

use constants

implicit none
private :: lambda_De_Broglie_clas,  lambda_landau_clas_ei, lambda_landau_clas_ee
private :: lambda_Debye_e, lambda_De_Broglie_rel_ei, lambda_landau_rel_ei
private :: lambda_De_Broglie_rel_ee, lambda_landau_rel_ee, lambda_e
public  :: gama, mom, vit, beta
public  :: log_clas_ei, log_clas_ee, log_rel_ei, log_rel_ee 
public  :: zeff, Tfermi, omega_pe, Ioni_m, vTh_e, lambda_i, lambda_Debye
public  :: Chemical_potential, Fermi_integrale

contains

!=====================================================================================
!                                 Special relativity
!=====================================================================================

elemental function gama(nrj)
! input  : nrj   = electron energy in (keV)
! output : gama = Lorentz factor in ()
implicit none
real(PR),intent(in) :: nrj
real(PR)            :: gama
gama = 1._PR + (nrj / mec2)
end function gama

elemental function beta(nrj)
! input  : nrj  = electron energy in (keV)
! output : beta = v / c in ()  
implicit none
real(PR),intent(in) :: nrj
real(PR)            :: beta
beta =  sqrt( 1._PR - (1._PR / (gama(nrj)**2._PR)) )
end function beta

elemental function mom(nrj)
! input  : nrj = electron energy in (keV)
! output : mom = electron momentum gama*me*v in (g.cm/s)     
implicit none
real(PR),intent(in) :: nrj
real(PR)            :: mom
mom = me * c * sqrt((gama(nrj)**2._PR) - 1._PR)
end function mom

elemental function vit(nrj)
! input  : nrj = electron energy in (keV)
! output : v   = electron velocity in (cm/s)      
implicit none
real(PR),intent(in) :: nrj
real(PR)            :: vit
vit = c * beta(nrj)
end function vit

!=====================================================================================
!              Ionic effective charge, Fermi Temperature, Debye length,
!   interionic mean distance, Langmuir plasma pulsation and mean ionization energy
!=====================================================================================

elemental function zeff(Z, ni, Te)
! input  : Z  = atomic number in ()
!          ni = atomic density in (/cm^3)
!          Te = electron Temperature in (K)
! output : zeff = ionization state 
!                 according to More E., Adv. At. Mol. Phys. 21, 305 (1985)
implicit none
real(PR), intent(in) :: Z, ni, Te
real(PR), parameter  :: a1 = 0.003323_PR, a2 = 0.9718_PR, a3 = 9.26148E-5_PR, a4 = 3.10165_PR
real(PR), parameter  :: b0 = -1.763_PR, b1 = 1.43175_PR, b2 = 0.31546_PR
real(PR), parameter  :: c1 = -0.366667_PR, c2 = 0.983333_PR
real(PR)             :: T, T0, R, Tf, AZ, BZ, CZ, Q1, Q, XZ, zeff0
real(PR)             :: K, fe, zeff
! More formula
T    =  Te / 1.1604e4_PR 
T0   =  T / (real(Z)**(4._PR / 3._PR))
R    =  (ni * mu) / Z
Tf   =  T0 / (1._PR + T0)
AZ   =  a1 * (T0**a2) + a3 * (T0**a4)
BZ   = -exp(b0 + b1 * Tf + b2 * (Tf**7._PR))
CZ   =  c1 * Tf + c2
Q1   =  AZ*(R**BZ)
Q    =  ((R**CZ) + (Q1**CZ))**(1._PR / CZ)
XZ   =  14.3139_PR * (Q**0.6624_PR)
zeff0 =  Z * XZ / (1._PR + XZ + sqrt(1._PR + 2._PR * XZ))
! Modifications for metals in order to impose a smooth transition between 
! the ambiant Temperature (zeff = 3 for Al, 1 for Cu and 2 for Ta) and
! the Temperature at which the metal becomes liquid
if (Z == 13.) then
K=0.5_PR*( 1._PR + tanh( (Te - Tfermi(Z, ni) ) / (0.05_PR*Tfermi(Z, ni)) ) )
fe=K**(1/(zeff0**2._PR))
zeff = (1._PR-fe)*3._PR + fe*zeff0
elseif (Z == 29.) then
K=0.5_PR*( 1._PR + tanh((Te-(35._PR*Tfermi(Z, ni))) / (1.925_PR*Tfermi(Z, ni))) )
fe=K**(1/(zeff0**2._PR))
zeff = ( ( ((1._PR-fe)*1._PR)**2._PR) + ((fe*zeff0)**2._PR) )**0.5_PR
elseif (Z == 73.) then
K=0.5_PR*( 1._PR + tanh( (Te-(20._PR*Tfermi(Z, ni))) / (1.25_PR*Tfermi(Z, ni))) )
fe=K**(1/(zeff0**2._PR))
zeff = ( ( ((1._PR-fe)*2._PR)**2._PR) + ((fe*zeff0)**2._PR) )**0.5_PR
else
zeff=zeff0
end if
end function zeff

elemental function Tfermi(Z, ni)
! input  : Z  = Atomic number in ()
!          ni = ion density in (/cm^3)
! output : Tfermi = Fermi Temperature in (K)
implicit none
real(PR), intent(in) :: Z, ni
real(PR)             :: Zc, Tfermi
if (Z == 13._PR) then
Zc = 3._PR
else if (Z == 29._PR) then
Zc = 1._PR
else if (Z == 73._PR) then
Zc = 2._PR
else 
Zc = Z
end if
Tfermi = (1._PR / kB) * (hbar**2._PR) / (2._PR * me) 
Tfermi = Tfermi * ((3._PR * (pi**2._PR) * Zc * ni)**(2._PR/3._PR))
end function Tfermi

elemental function omega_pe(Z, ni, Te)
! input  : Z  = atomic number in ()
!          ni = atomic density in (/cm^3)
!          Te = electron Temperature in (K)
! output : omega_pe = Langmuir frequency in (rad/s)
implicit none
real(PR), intent(in) :: Z, ni, Te
real(PR)             :: omega_pe
omega_pe = sqrt(4 * pi * Zeff(Z,ni,Te) * ni * (e**2._PR) / me)
end function omega_pe

elemental function Ioni_m(Z, ni, Te)
! input  : Z  = atomic number in ()
!          ni = atomic density in (/cm^3)
!          Te = electron Temperature in (K)
! output : Ioni_m = mean excitation potential of bound electrons in an ion in (erg)
!                   according to More E., Adv. At. Mol. Phys. 21, 305 (1985)
implicit none
real(PR), intent(in) :: Z, ni, Te
real(PR)             :: q, Ioni_m
q = Zeff(Z, ni, Te) / Z
Ioni_m = 10._PR* 1.6022E-12_PR * Z 
Ioni_m = Ioni_m * exp(1.29_PR * q**(0.72_PR - 0.18_PR*q))/(1._PR- q)**(1._PR/2._PR)
end function Ioni_m

elemental function vTh_e(Z, ni, Te)
! input  : Z  = atomic number in ()
!          ni = atomic density in (/cm^3)
!          Te = electron Temperature in (K)
! output : vTh_e = Thermal velocity in (cm/s) with the degeneracy correction
!                   according to Lee and More, Phys. Fluids, 27(5), p. 1273 (1984) 
implicit none
real(PR), intent(in) :: Z, ni, Te
real(PR)             :: vTh_e
vTh_e = (3._PR* kb * (( (Te**2._PR) + (Tfermi(Z, ni)**2._PR) )**0.5_PR) / me )**0.5_PR
end function vTh_e

elemental function lambda_De_Broglie_clas(Z, ni, Te)
! input  : Z  = atomic number in ()
!          ni = atomic density in (/cm^3)
!          Te = electron Temperature in (K)
! output : lambda_De_Broglie_clas = non relativistic electron De Broglie
!                                   wavelength in (cm)
implicit none
real(PR), intent(in) :: Z, ni, Te
real(PR)             :: lambda_De_Broglie_clas
lambda_De_Broglie_clas = hbar / (2._PR * me * vTh_e(Z, ni, Te)) 
end function lambda_De_Broglie_clas

elemental function lambda_Debye_e(Z, ni, Te)
! input  : Z  = atomic number in ()
!          ni = atomic density in (/cm^3)
!          Te = electron Temperature in (K)
! output : lambda_Debye_e = elecron Debye shielding length in (cm) with degeneracy 
!                           correction according to 
!                           Lee and More, Phys. Fluids, 27(5), p. 1273 (1984)   
implicit none
real(PR), intent(in) :: Z, ni, Te
real(PR)             :: lambda_Debye_e
lambda_Debye_e = (4 * pi * zeff(Z, ni, Te) * ni * e**2._PR /&
& ( kB * ( (Te**2._PR) + ((Tfermi(Z, ni))**2._PR)  )**0.5_PR) )**(-0.5_PR)     
end function lambda_Debye_e

elemental function lambda_Debye(Z, ni, Te, Ti)
! input  : Z  = atomic number in ()
!          ni = atomic density in (/cm^3)
!          Te = electron Temperature in (K)
!          Ti = ion Temperature in (K)
! output : lambda_Debye = Debye shielding length in (cm) with degeneracy 
!                         accounting for ion contribution correction according to 
!                         Lee and More, Phys. Fluids, 27(5), p. 1273 (1984)      
implicit none
real(PR), intent(in) :: Z, ni, Te, Ti
real(PR)             :: lambda_Debye
lambda_Debye = ( (lambda_Debye_e(Z, ni, Te)**(-2._PR)) &
& +  (4 * pi * zeff(Z, ni, Te) * ni * e**2._PR /(kB*Ti)  ) )**(-0.5_PR) 
end function lambda_Debye

!=====================================================================================
!                    Classical electron-ion Coulomb logarithm
!=====================================================================================

elemental function lambda_i(ni)
! input  : ni = atomic density in (/cm^3)
! output : lambda_i = Mean atomic distance in (/cm^3)
implicit none
real(PR), intent(in) :: ni
real(PR)             :: lambda_i
lambda_i = (3 / (4 * pi *ni))**(1._PR/3._PR)
end function lambda_i

elemental function lambda_Landau_clas_ei(Z, ni, Te)
! input  : Z  = atomic number in ()
!          ni = atomic density in (/cm^3)
!          Te = electron Temperature in (K)
! output : lambda_Landau_clas_ei = classical Landau electron-ion mean interaction
!                                  distance in (cm)
implicit none
real(PR), intent(in) :: Z, ni, Te
real(PR)             :: lambda_Landau_clas_ei
lambda_Landau_clas_ei = zeff(Z, ni, Te)*(e**2._PR)/(0.5_PR*me*vTh_e(Z, ni, Te)**2._PR)
end function lambda_Landau_clas_ei

elemental function log_clas_ei(Z, ni, Te, Ti)
! input  : Z  = atomic number in ()
!          ni = atomic density in (/cm^3)
!          Te = electron Temperature in (K)
!          Ti = ion Temperature in (K)
! output : log_clas_ei = classical electron-ion Coulomb logarithm in () according to 
!                        Lee and More, Phys. Fluids, 27(5), p. 1273 (1984)
implicit none
real(PR), intent(in) :: Z, ni, Te, Ti
real(PR)             :: bmax, bmin, log_clas_ei
! maximal impact parameter
bmax = ( (lambda_Debye(Z, ni, Te, Ti)**2._PR)  &
&+   (lambda_i(ni)**2._PR)                  )**(0.5_PR)
! minimal impact parameter
bmin = ( (lambda_Landau_clas_ei(Z, ni, Te)**2._PR) &
&+   (lambda_De_Broglie_clas(Z, ni, Te)**2._PR) )**(0.5_PR)
log_clas_ei = ( (2._PR**2._PR) &
&+   ((0.5_PR * log(1 + ((bmax**2._PR) / (bmin**2._PR))))**2._PR) )**0.5_PR
end function log_clas_ei

!=====================================================================================
!                   Classical electron-electron Coulomb logarithm
!=====================================================================================

elemental function lambda_e(Z, ni, Te)
! input  : Z  = atomic number in ()
!          ni = atomic density in (/cm^3)
!          Te = electron Temperature in (K)
! output : lambda_i = Mean electronic distance in (/cm^3)
implicit none
real(PR), intent(in) :: Z, ni, Te
real(PR)             :: lambda_e
lambda_e = (3 / (4 * pi * Zeff(Z, ni, Te) * ni))**(1._PR/3._PR)
end function lambda_e

elemental function lambda_Landau_clas_ee(Z, ni, Te)
! input  : Z  = atomic number in ()
!          ni = atomic density in (/cm^3)
!          Te = electron Temperature in (K)
! output : lambda_Landau_clas_ee = classical Landau electron-electron mean interaction
!                                  distance in (cm)
implicit none
real(PR), intent(in) :: Z, ni, Te
real(PR)             :: lambda_Landau_clas_ee
lambda_Landau_clas_ee =e**2._PR / (0.5_PR * me * vTh_e(Z, ni, Te)**2._PR )
end function lambda_Landau_clas_ee

elemental function log_clas_ee(Z, ni, Te)
! input  : Z  = atomic number in ()
!          ni = atomic density in (/cm^3)
!          Te = electron Temperature in (K)
! output : log_clas_ee = classical electron-electron Coulomb logarithm in () according
!                        to Lee and More, Phys. Fluids, 27(5), p. 1273 (1984)
implicit none
real(PR), intent(in) :: Z, ni, Te
real(PR)             :: bmax, bmin, log_clas_ee
! maximal impact parameter
bmax = ( (lambda_Debye_e(Z, ni, Te)**2._PR)+(lambda_e(Z, ni, Te)**2._PR) )**(0.5_PR)
! minimal impact parameter
bmin = ( (lambda_Landau_clas_ee(Z, ni, Te)**2._PR) &
&+   (lambda_De_Broglie_clas(Z, ni, Te)**2._PR) )**0.5_PR
log_clas_ee = ( (2._PR**2._PR) &
&+ ((0.5_PR * log(1 + (bmax**2._PR / bmin**2._PR)))**2._PR) )**0.5_PR
end function log_clas_ee

!=====================================================================================
!                       Relativistic electron-ion coulomb logarithm
!=====================================================================================

elemental function lambda_De_Broglie_rel_ei(nrj)
! input  : nrj = electron energy in (keV)
! output : lambda_De_Broglie_rel_ei = relativistic electron De Broglie wavelength
!                                     in an ion potential in (cm)
implicit none
real(PR), intent(in) :: nrj
real(PR)             :: lambda_De_Broglie_rel_ei
lambda_De_Broglie_rel_ei = 2._PR * hbar / mom(nrj) 
end function lambda_De_Broglie_rel_ei

elemental function lambda_Landau_rel_ei(Z, nrj)
! input  : Z   = atomic number in ()
!          nrj = electron energy in (keV)
! output : lambda_Landau_rel_ei = relativistic Landau electron-ion mean interaction
!                                  distance in (cm)
implicit none
real(PR), intent(in) :: Z, nrj
real(PR)             :: lambda_Landau_rel_ei
lambda_Landau_rel_ei = Z * (e**2._PR) / (nrj * keV)
end function lambda_Landau_rel_ei

elemental function log_rel_ei(Z, ni, Te, Ti, nrj)
! input  : Z   = atomic number in ()
!          ni  = atomic density in (/cm^3)
!          Te  = electron Temperature in (K)
!          Ti  = ion Temperature in (K)
!          nrj =  in (keV)  
! output : log_rel_ei = relativistic electron-ion Coulomb logarithm in ()
implicit none
real(PR), intent(in) :: Z, ni, Te, Ti, nrj
real(PR)             :: bmax, bmin, log_rel_ei
! maximal impact parameter
bmax = max(lambda_Debye(Z, ni, Te, Ti), lambda_i(ni))
! minimal impact parameter
bmin = max(lambda_Landau_rel_ei(Z, nrj), lambda_De_Broglie_rel_ei(nrj))
log_rel_ei = 0.5_PR * log(1 + (bmax**2._PR / bmin**2._PR))
end function log_rel_ei

!=====================================================================================
!                     Relativistic electron-electron coulomb logarithm
!=====================================================================================

elemental function lambda_De_Broglie_rel_ee(nrj)
! input  : nrj = electron energy in (keV)
! output : lambda_De_Broglie_rel_ee = relativistic electron De Broglie
!                                     wavelength in (cm) when interacting with 
!                                     another electron   
implicit none
real(PR), intent(in) :: nrj
real(PR)             :: lambda_De_Broglie_rel_ee
lambda_De_Broglie_rel_ee = 2._PR * hbar / (me * c * sqrt((gama(nrj) - 1._PR) / 2._PR))  
end function lambda_De_Broglie_rel_ee

elemental function lambda_Landau_rel_ee(Z, nrj)
! input  : Z   = atomic number in ()
!          nrj = electron energy in (keV)
! output : lambda_Landau_rel_ei = relativistic Landau electron-ion mean interaction
!                                  distance in (cm)    
implicit none
real(PR), intent(in) :: Z, nrj
real(PR)             :: aa, lambda_Landau_rel_ee
aa = me * (c**2._PR) *( sqrt( (gama(nrj) + 1._PR) / 2._PR ) - 1._PR )
lambda_Landau_rel_ee = Z * (e**2._PR) / aa
end function lambda_Landau_rel_ee

elemental function log_rel_ee(Z, ni, Te, Ti, nrj)
! input  : Z   = atomic number in ()
!          ni  = atomic density in (/cm^3)
!          Te  = electron Temperature in (K)
!          Ti  = ion Temperature in (K)
!          nrj =  in (keV)  
! output : log_rel_ei = relativistic electron-ion Coulomb logarithm in ()    
implicit none
real(PR), intent(in) :: Z, ni, Te, Ti,nrj
real(PR)             :: bmax, bmin, log_rel_ee
! maximal impact parameter
bmax = max(lambda_Debye(Z, ni, Te, Ti), lambda_i(ni))
! mimimal impact parameter
bmin = max(lambda_Landau_rel_ee(Z, nrj), lambda_De_Broglie_rel_ee(nrj))
log_rel_ee = 0.5_PR * log(1 + (bmax**2._PR / bmin**2._PR))
end function log_rel_ee

!=====================================================================================
!                                   Fermi-Dirac Integrals
!=====================================================================================


elemental function Chemical_potential(Z, ni, Te)
! input  : Z   = atomic number in ()
!          ni  = atomic density in (/cm^3)
!          Te  = electron Temperature in (K)
! output : Chemical potential in (erg) according to 
!          ANTIA, H. M., Astrophysical Journal Supplement Series, 84, 101 (1993)
implicit none
real(PR), intent(in) :: Z, ni, Te
real(PR)             :: Chemical_potential
real(PR), parameter  :: a0=1.999266880833e4_PR, a1=5.702479099336e3_PR
real(PR), parameter  :: a2=6.610132843877e2_PR, a3=3.818838129486e1_PR, a4=1._PR
real(PR), parameter  :: b0= 1.771804140488e4_PR, b1=-2.014785161019e3_PR
real(PR), parameter  :: b2= 9.130355392717e1_PR, b3=-1.670718177489
real(PR), parameter  :: c0=-1.277060388085e-2_PR, c1= 7.187946804945e-2_PR
real(PR), parameter  :: c2=-4.262314235106e-1_PR, c3= 4.997559426872e-1_PR
real(PR), parameter  :: c4=-1.285579118012, c5=-3.930805454272e-1_PR, c6= 1._PR
real(PR), parameter  :: d0=-9.745794806288e-3_PR, d1= 5.485432756838e-2_PR
real(PR), parameter  :: d2=-3.299466243260e-1_PR, d3= 4.077841975923e-1_PR
real(PR), parameter  :: d4=-1.145531476975, d5=-6.067091689181e-2_PR
real(PR)             :: y, R_p, yy, R_m, R
y = (2._PR/3._PR) * zeff(Z, ni, Te) * ( (Tfermi(Z, ni)/Te)**(3._PR/2._PR) )
R_p= ( a0 + (a1*y) + (a2*(y**2._PR)) + (a3*(y**3._PR)) + (a4*(y**4._PR)) ) / &
& ( b0 + (b1*y) + (b2*(y**2._PR)) + (b3*(y**3._PR)) ) 
yy=y**(-2._PR/3._PR)
R_m= ( c0 + (c1*yy) + (c2*(yy**2._PR)) + (c3*(yy**3._PR)) + (c4*(yy**4._PR)) + (c5*(yy**5._PR)) + (c6*(yy**6._PR)) ) / &
& ( d0 + (d1*yy) + (d2*(yy**2._PR)) + (d3*(yy**3._PR)) + (d4*(yy**4._PR)) + (d5*(yy**5._PR)) )
if (y < 4) then
R= log(y*R_p)
else
R=R_m/yy
end if
Chemical_potential = R*kb*Te
end function Chemical_potential

!=====================================================================================
!                                  /infinity
! Fermi-Dirac integral : F_j(Te) = |   x^j / [ exp( x- (mu/kTe) )+ 1 ] dx
!                                  /0
!=====================================================================================

elemental function Fermi_integrale(j, Z, ni, Te)
! input  : j  = 1/2 or integer
!          Z  = atomic number in ()
!          ni = atomic density in (/cm^3)
!          Te = electron Temperature in (K)
! output : Fermi_integrale = Value of the Fermi-Dirac integral at Te F_j(Te) according
!          to Aymerich-Humet X. et al., Jour. of Applied Phys 54(5),2850 (1983)
implicit none
real(PR), intent(in) :: j, Z, ni, Te
real(PR)             :: muc, Fermi_integrale
real(PR)             :: x, a, b, c, d
muc = Chemical_potential(Z, ni, Te)
x  = muc/(kb*Te)
a = ( 1._PR + ((15._PR/4._PR)*(real(j+1,PR))) + (((real(j+1,PR))**2._PR)/40._PR) )**0.5_PR
b = 1.8_PR + (0.61_PR*real(j,PR))
c = 2._PR + ( ( 2._PR - sqrt(2._PR) ) * (2._PR**(-real(j,PR))) )
d = (b + x + ( (abs(x-b)**c) + (a**c) )**(1._PR/c) )**(real(j+1,PR))
Fermi_integrale = ( ( real(j+1,PR) * (2._PR**real(j+1,PR)) / d ) &
&+   (exp(-x)/gamma(real(j+1,PR)))           )**(-1._PR)
end function Fermi_integrale

end module physics_library