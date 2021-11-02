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
module fokker_planck_coef

  use acuracy
  use constants
  use physics_library

  implicit none

  private :: S_i, S_e_free, S_e_bound, S_p  
  private :: nu_i, nu_e_free, nu_e_bound, nu_p
  public  :: S_tot, nu_tot 

  contains

  !=====================================================================================
  !                                   Stopping powers
  !=====================================================================================

  elemental function S_i(A, Z, ni, Te, Ti, nrj)
    ! input  : A   = Atomic weight of the material in ()
    !          Z   = Atomic number of the material in ()
    !          ni  = Ion nuclei density in the material in (/cm^3)
    !          Te  = Electron temperature in the material in (K)
    !          Ti  = Ion temperature in the material in (K)
    !          nrj = Electron projectile kinetic energy in (keV)
    ! output : S_i = Relativistic electron stopping power contribution due to 
    !                collisions with ion nuclei in (erg/cm) 
    !                according to Mott, N. F., Proc. R. Soc. Lond. A 135, 429 (1932) 
    implicit none
    real(PR), intent(in) :: A, Z, ni, Te, Ti, nrj
    real(PR)             :: chi, L, S_i
    chi = 4 * pi * ni * (Z**2._PR) * (e**4._PR)  / (A * mu * (vit(nrj)**2._PR))
    L   = log_rel_ei(Z, ni, Te, Ti, nrj) 
    L   = L - (1._PR + (beta(nrj)**2._PR) / 2._PR)
    S_i = chi * L
  end function S_i

  elemental function S_e_free(Z, ni, Te, Ti, nrj)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion nuclei density in the material in (/cm^3)
    !          Te  = Electron temperature in the material in (K)
    !          Ti  = Ion temperature in the material in (K)
    !          nrj = Electron projectile kinetic energy in (keV)
    ! output : S_e_free = Relativistic electron stopping power contribution due to 
    !                     collisions with free electrons in (erg/cm)   
    !                     according to Möller, C., Ann. Phys. 14, 531 (1932)
    implicit none
    real(PR), intent(in) :: Z, ni, Te, Ti, nrj
    real(PR)             :: Zf, chi, g, L,  S_e_free    
    Zf  = Zeff(Z, ni, Te)
    chi = 4 * pi * Zf * ni * (e**4._PR)  / (me * (vit(nrj)**2._PR))
    g   = gama(nrj)
    L   = log_rel_ee(Z, ni, Te, Ti,nrj) 
    L   = L - log(2._PR) + 0.5_PR - (((2._PR*g-1._PR)/(2._PR*(g**2._PR)))*log(2._PR))
    L   = L + (1._PR / 16._PR) * ((g-1._PR)/g)**2._PR
    S_e_free = chi * L
  end function S_e_free

  elemental function S_e_bound(Z, ni, Te, nrj)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion nuclei density in the material in (/cm^3)
    !          Te  = Electron temperature in the material in (K)
    !          nrj = Electron projectile kinetic energy in (keV)
    ! output : S_e_bound = Relativistic electron stopping power contribution due to 
    !                      collisions with bound electrons in (erg/cm) 
    !                      according to Bethe, H. Z. f. Physik, 76, 293 (1932)
    !                               and Möller, C., Ann. Phys. 14, 531 (1932) 
    implicit none
    real(PR), intent(in) :: Z, ni, Te, nrj
    real(PR)             :: Zf, chi, g, L,  S_e_bound    
    Zf  = Zeff(Z, ni, Te)
    chi = 4 * pi * (Z - Zf) * ni * (e**4._PR)  / (me * (vit(nrj)**2._PR))
    g   = gama(nrj)
    L   = log(nrj * keV * sqrt((g + 1._PR)/2._PR )/ Ioni_m(Z, ni, Te)) 
    L   = L + (1/(2._PR*(g**2._PR))) - (((2._PR*g-1._PR)/(2._PR*(g**2._PR)))*log(2._PR))
    L   = L + (1._PR / 16._PR) * (((g-1._PR)/g)**2._PR)
    S_e_bound = chi * L
  end function S_e_bound

  elemental function S_p(Z, ni, Te, Ti, nrj)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion nuclei density in the material in (/cm^3)
    !          Te  = Electron temperature in the material in (K)
    !          Ti  = Ion temperature in the material in (K)
    !          nrj = Electron projectile kinetic energy in (keV)
    ! output : S_p = Relativistic electron stopping power contribution due to 
    !                collisions with screened free electrons (plasmons) in (erg/cm)  
    !                according to Pines D. and Bohm D., Phys. rev. 85, 2 (1952)
    !                but with a modified screening length    
    implicit none
    real(PR), intent(in) :: Z, ni, Te, Ti, nrj
    real(PR)             :: Zf, chi, bmax, L, S_p
    Zf   = Zeff(Z, ni, Te)
    chi  = 4 * pi * Zf * ni * (e**4._PR)  / (me * (vit(nrj)**2._PR))
    bmax = ( (lambda_Debye(Z, ni, Te, Ti)**2._PR)  &
    &+   (lambda_i(ni)**2._PR)                  )**(0.5_PR)
    L    = vit(nrj) / ( omega_pe(Z, ni, Te) * bmax)
    L    = log(L)
    S_p = chi * L
  end function S_p

  elemental function S_rad(Z, ni, Te, nrj)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion nuclei density in the material in (/cm^3)
    !          Te  = Electron temperature in the material in (K)
    !          nrj = Electron projectile kinetic energy in (keV)
    ! output : S_rad = Relativistic electron stopping power contribution due to 
    !                 Bremsstrahlung radiation losses in (erg/cm)    
    !                 according to Heitler W. and Sauter F., Nature 132, 892 (1933)
    implicit none
    real(PR), intent(in) :: Z, ni, Te, nrj
    real(PR)             :: Zf, chi, L, S_rad
    Zf  = Zeff(Z, ni, Te)
    chi = 4 * pi * (alpha/pi) * (real(Z)-Zf)*(real(Z)-Zf+1)*ni*(e**4._PR)*gama(nrj) 
    chi = chi / (me *(c**2._PR))
    L   = log(2._PR*gama(nrj)) - 1._PR / 3._PR
    S_rad = chi * L
  end function S_rad

  elemental function S_tot(A, Z, ni, Te, Ti, nrj)
    ! input  : A   = Atomic weight of the material in ()
    !          Z   = Atomic number of the material in ()
    !          ni  = Ion nuclei density in the material in (/cm^3)
    !          Te  = Electron temperature in the material in (K)
    !          Ti  = Ion temperature in the material in (K)
    !          nrj = Electron projectile kinetic energy in (keV)
    ! output : S_tot = Relativistic electron stopping power in (erg/cm)  
    implicit none
    real(PR), intent(in) :: A, Z, ni, Te, Ti, nrj
    real(PR) :: Si, Sef, Seb, Sp, Srad, S_tot     
    Si    = S_i(A,Z, ni, Te, Ti, nrj)
    Sef   = S_e_free(Z, ni, Te, Ti, nrj)
    Seb   = S_e_bound(Z, ni, Te, nrj)
    Sp    = S_p(Z, ni, Te, Ti, nrj)
    Srad  = S_rad(Z, ni, Te, nrj)
    S_tot = Si + Sef + Seb + Sp + Srad
  end function S_tot

  !=====================================================================================
  !                              Angular collision rates
  !=====================================================================================
  ! The following isotropization rates nu_alpha by colliding alpha particles 
  ! (ion nuclei, free electrons, bound electrons or screened free electrons) 
  ! are computed according to the formula : nu_alpha = ( m_alpha v / p^2 ) S_alpha
  ! where m_alpha is the targetted particle mass, 
  !       v and p are respectively the velocity and the momentum of the relativistic 
  !       electron projectile and
  !       S_alpha is the corresponding relativistic electron stopping power contribution
  !       due to collisions with alpha particles
  ! as demonstrated by Touati, M. et al., New. Jour. Phys. 16, 073014 (2014)

  elemental function nu_i(Z, ni, Te, Ti, nrj)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion nuclei density in the material in (/cm^3)
    !          Te  = Electron temperature in the material in (K)
    !          Ti  = Ion temperature in the material in (K)
    !          nrj = Electron projectile kinetic energy in (keV)
    ! output :nu_i = Relativistic electron isotropization rate contribution due to
    !                collisions with ion nuclei in (/s)     
    implicit none
    real(PR), intent(in) :: Z, ni, Te, Ti, nrj
    real(PR)             :: chi, g, L, nu_i
    g   = gama(nrj)
    chi = 4 * pi * ni * (Z**2._PR) * (e**4._PR)
    chi = chi / ((g**2._PR) * (me**2._PR) * (vit(nrj)**3._PR))
    L   = log_rel_ei(Z, ni, Te, Ti, nrj) 
    L   = L - ((1._PR + (beta(nrj)**2._PR)) / 2._PR)
    nu_i = chi * L 
  end function nu_i

  elemental function nu_e_free(Z, ni, Te, Ti, nrj)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion nuclei density in the material in (/cm^3)
    !          Te  = Electron temperature in the material in (K)
    !          Ti  = Ion temperature in the material in (K)
    !          nrj = Electron projectile kinetic energy in (keV)
    ! output : nu_e_free = Relativistic electron isotropization rate contribution due to 
    !                      collisions with free electrons in (/s)      
    implicit none
    real(PR), intent(in) :: Z, ni, Te, Ti, nrj
    real(PR)             :: Zf, chi, g, L,  nu_e_free
    Zf  = Zeff(Z, ni, Te)
    g   = gama(nrj)
    chi = 4 * pi * Zf * ni * (e**4._PR)  / (g**2._PR * (me**2._PR) * (vit(nrj)**3._PR))
    L   = log_rel_ee(Z, ni, Te, Ti,nrj) 
    L   = L - log(2._PR) + 0.5_PR - (((2._PR*g-1._PR)/(2._PR*(g**2._PR)))*log(2._PR))
    L   = L + (1._PR / 16._PR) * ((g-1._PR)/g)**2._PR
    nu_e_free = chi * L 
  end function nu_e_free

  elemental function nu_e_bound(Z, ni, Te, nrj)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion nuclei density in the material in (/cm^3)
    !          Te  = Electron temperature in the material in (K)
    !          nrj = Electron projectile kinetic energy in (keV)
    ! output : nu_e_bound = Relativistic electron isotropization rate contribution due to 
    !                      collisions with bound electrons in (/s)     
    implicit none
    real(PR), intent(in) :: Z, ni, Te, nrj
    real(PR)             :: Zf, chi, g, L,  nu_e_bound
    Zf  = Zeff(Z, ni, Te)
    g   = gama(nrj)
    chi = 4 * pi * (Z - Zf) * ni * (e**4._PR)  / (g**2._PR * (me**2._PR) * (vit(nrj)**3._PR))
    L   = log(nrj * keV * sqrt((g + 1._PR)/2._PR )/ Ioni_m(Z, ni, Te)) 
    L   = L + (1/(2._PR*(g**2._PR))) - (((2._PR*g-1._PR)/(2._PR*(g**2._PR)))*log(2._PR))
    L   = L + (1._PR / 16._PR) * ((g-1._PR)/g)**2._PR
    nu_e_bound = chi * L
  end function nu_e_bound

  elemental function nu_p(Z, ni, Te, Ti, nrj)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion nuclei density in the material in (/cm^3)
    !          Te  = Electron temperature in the material in (K)
    !          Ti  = Ion temperature in the material in (K)
    !          nrj = Electron projectile kinetic energy in (keV)
    ! output : nu_p = Relativistic electron isotropization rate contribution due to 
    !                collisions with screened free electrons (plasmons) in (/s)         
    implicit none
    real(PR), intent(in) :: Z, ni, Te, Ti, nrj
    real(PR)             :: Zf, g, chi, bmax, L, nu_p
    Zf   = Zeff(Z, ni, Te)
    g    = gama(nrj)
    chi  = 4 * pi * Zf * ni * (e**4._PR)  / (g**2._PR * (me**2._PR) * (vit(nrj)**3._PR))
    bmax = ( (lambda_Debye(Z, ni, Te, Ti)**2._PR)  &
    &+   (lambda_i(ni)**2._PR)                  )**(0.5_PR)
    L    = vit(nrj) / ( omega_pe(Z, ni, Te) * bmax)
    L    = log(L)
    nu_p = chi * L
  end function nu_p

  elemental function nu_tot(Z, ni, Te, Ti, nrj)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion nuclei density in the material in (/cm^3)
    !          Te  = Electron temperature in the material in (K)
    !          Ti  = Ion temperature in the material in (K)
    !          nrj = Electron projectile kinetic energy in (keV)
    ! output :nu_tot = Relativistic electron isotropization rate in (/s)
    implicit none
    real(PR), intent(in) :: Z, ni, Te, Ti, nrj
    real(PR) :: nui, nuef, nueb, nup, nu_tot 
    nui    = nu_i(Z, ni, Te, Ti, nrj)
    nuef   = nu_e_free(Z, ni, Te, Ti, nrj)
    nueb   = nu_e_bound(Z, ni, Te, nrj)
    nup    = nu_p(Z, ni, Te, Ti, nrj) 
    nu_tot = nui + nuef + nueb + nup
  end function nu_tot

end module fokker_planck_coef