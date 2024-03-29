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
module transport_coef

  use acuracy
  use constants
  use physics_library
  use input

  implicit none
  private :: nu_ee, nu_eph, nu_c, nu_ei, nu_res, nu_cond
  private :: eps_delta_T, gama_E, gama_Lorenz
  private :: A_alpha, A_beta, gama_ii, G_gama, S_11
  private :: nu_Spitzer_cond, nu_Spitzer_res
  private :: conductivity_metals, resistivity_metals
  private :: resistivity_tab, resistivity_H, cond_H
  public  :: capacity, resistivity_Spitzer, ion_capacity
  public  :: Omega_ei, resis, cond, initialize_resistivity

  contains

  !=====================================================================================
  !                                    Thermal capacities
  !=====================================================================================

  elemental function capacity(Z, ni, Te)
    ! input : Z  = atomic number in ()
    !         ni = ion density in (/cm^3)
    !         Te = background electron Temperature in (K)
    ! output : capacity = background electron thermal capacity in (erg/K/cm^3)
    !          according to a fit between the perfect electron gaz model 
    !                                 and the Sommerfeld model for solids  
    implicit none
    real(PR), intent(in) :: Z, ni, Te
    real(PR)             :: capacity   
    if (Z == Z_Al) then
      capacity =( ((912._PR*Te)**(-2._PR)) &
      &+ ((1.5_PR * kb * zeff(Z, ni, Te) * ni)**(-2._PR)) )**(-0.5_PR)
    elseif (Z == Z_Cu) then
      capacity =( ((968._PR*Te)**(-2._PR)) &
      &+ ((1.5_PR * kb * zeff(Z, ni, Te) * ni)**(-2._PR)) )**(-0.5_PR)
    elseif (Z == Z_Ta) then
      capacity =( ((5428.8_PR*Te)**(-2._PR)) &
      &+ ((1.5_PR * kb * zeff(Z, ni, Te) * ni)**(-2._PR)) )**(-0.5_PR)
    else
      capacity = 1.5_PR * kb * zeff(Z, ni, Te) * ni
    endif
  end function capacity

  elemental function ion_capacity(Z, ni, Ti)
    ! input : Z  = atomic number in ()
    !         ni = ion density in (/cm^3)
    !         Ti = background ion Temperature in (K)
    ! output : ion_capacity = background ion/phonons thermal capacity in (erg/K/cm^3)
    !          according to the Einstein model
    implicit none
    real(PR), intent(in) :: Z, ni, Ti
    real(PR)             :: x
    real(PR)             :: ion_capacity    
    if (Z == Z_Al) then
      x=Ti/284._PR
      ion_capacity = 1.5_PR * kb * ni &
      &* ((1._PR/x)**2._PR) * exp(1._PR/x) / ((exp(1._PR/x)-1._PR)**2._PR)
    elseif (Z == Z_Cu) then
      x=Ti/278._PR
      ion_capacity = 1.5_PR * kb * ni &
      &* ((1._PR/x)**2._PR) * exp(1._PR/x) / ((exp(1._PR/x)-1._PR)**2._PR)
    elseif (Z == Z_Ta) then
      x=Ti/193._PR
      ion_capacity = 1.5_PR * kb * ni &
      &* ((1._PR/x)**2._PR) * exp(1._PR/x) / ((exp(1._PR/x)-1._PR)**2._PR)
    else
      ion_capacity = 1.5_PR * kb * ni
    end if
  end function ion_capacity

  !=====================================================================================
  !                                Electron collision rates
  !=====================================================================================

  elemental function nu_ei(Z, ni, Te, Ti)
    ! input : Z      = atomic number in ()
    !         ni     = ion density in (/cm^3)
    !         Te/Ti  = background electron/ion Temperature in (K)
    ! output : nu_ei = electron-ion collision frequency in (/s) for plasmas
    !          according to Lee Y. and More R. , Physics of Fluids Vol. 27 p. 1273 (1984)
    implicit none
    real(PR), intent(in) :: Z, ni, Te, Ti
    real(PR)             :: mu, nu_ei
    mu    = Chemical_potential(Z, ni, Te)/(kB*Te)
    nu_ei = ( 2._PR * sqrt(2._PR) * pi * (e**4._PR) ) / ( 3._PR * sqrt(me) )
    nu_ei = nu_ei * (zeff(Z, ni, Te)**2._PR) * ni 
    nu_ei = nu_ei * log_clas_ei(Z, ni, Te, Ti)
    nu_ei = nu_ei / ( (kB*Te)**1.5_PR )
    nu_ei = nu_ei / ( (1._PR + exp(-mu)) * Fermi_integrale(0.5_PR, Z, ni, Te) )
  end function nu_ei

  elemental function nu_eph(Z, ni,Ti)
    ! input : Z  = atomic number in ()
    !         ni = ion density in (/cm^3)
    !         Te = background electron Temperature in (K)
    ! output : nu_eph = electron-phonon collision frequency in (/s) for metals
    !          according to Chimier N. et al., Phys. Rev. B 75, 195124 (2007)
    !          and the electrical resistivity value of the metal at ambiant cond.
    implicit none
    real(PR), intent(in) :: Z, ni, Ti
    real(PR)             :: ks, nu_eph
    if (Z == Z_Al) then
      ks = 1.6714_PR
    elseif (Z == Z_Cu) then
      ks = 0.3764_PR
    elseif (Z == Z_Ta) then
      ks = 4.7442_PR
    else
      ks = 1._PR ! arbitrary value
    end if
    nu_eph = ks * 2._PR * (e**2._PR) * kb * Ti /&
    & ((hbar**2._PR) * sqrt(2._PR* kb * Tfermi(Z, ni) / me))
  end function nu_eph

  elemental function nu_ee(Z, ni, Te)
    ! input : Z  = atomic number in ()
    !         ni = ion density in (/cm^3)
    !         Te = background electron Temperature in (K)
    ! output : nu_ee = electron-electron collision frequency in (/s) for metals
    !          according to a fit of some results found by 
    !          Petrov et al., JETP Letters Vol. 97 p. 20 (2013)
    implicit none
    real(PR), intent(in) :: Z, ni, Te
    real(PR)             :: A_nu, B_nu, nu_ee1, nu_ee2, nu_ee
    if (Z == Z_Al) then
      A_nu = 1.9066_PR ! collisions between s-band electrons and themselves
      B_nu = 0.125_PR  ! collisions between s-band electrons and themselves
    elseif (Z == Z_Cu) then
      !A_nu = 136.81_PR ! collisions between s-band electrons and themselves, only
      !B_nu = 0.63_PR   ! collisions between s-band electrons and themselves, only
      A_nu = 252.47_PR  ! collisions between s-band and s or d-band electrons 
      B_nu = 0.819_PR   ! collisions between s-band and s or d-band electrons
    elseif (Z == Z_Ta) then
      !A_nu = 1.95_PR  ! collisions between s-band electrons and themselves, only
      !B_nu = 0.14_PR  ! collisions between s-band electrons and themselves, only
      A_nu = 13.393_PR ! collisions between s-band and s or d-band electrons 
      B_nu = 0.280_PR  ! collisions between s-band and s or d-band electrons 
    else
      A_nu = 1._PR ! arbitrary value
      B_nu = 0._PR ! arbitrary value
    end if 
    nu_ee1 = A_nu * kb  * (Te**2._PR) / (hbar * Tfermi(Z, ni))
    nu_ee2 = B_nu * kb * sqrt( Te * Tfermi(Z, ni) ) / hbar
    nu_ee = ( (nu_ee1**(-2._PR)) + (nu_ee2**(-2._PR)) )**(-0.5_PR)  
  end function nu_ee

  elemental function nu_c(Z, ni, Te)
    ! input : Z  = atomic number in ()
    !         ni = ion density in (/cm^3)
    !         Te = background electron Temperature in (K)
    ! output : nu_c = collision frequency threshold in (/s) for metals
    !          corresponding to the fact that the electron mean free path cannot exceed
    !          the mean distance between ions
    implicit none
    real(PR), intent(in) :: Z, ni, Te
    real(PR)             :: nu_c
    nu_c = vTh_e(Z, ni, Te) / lambda_i(ni)
  end function nu_c

  !=====================================================================================
  !                           Electron-ion thermal equilibration             
  !=====================================================================================

  elemental function Omega_ei(A, Z, ni, Te, Ti) 
    ! input  : A   = Atomic weight of the material in ()
    !          Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Te  = Electron Temperature in the material in (K)
    !          Ti  = Ion Temperature in the material in (K)
    ! output : Omega_ei = electron-ion or electron-latice coupling factor in
    !          (erg/s/cm^3/K) inspired by
    !          Chen J. K. et al., Journal of laser applications, 17, p. 63 (2005) 
    !          for metals and Brysk H. et al., Plasma phys. 17, 6, p. 2714 (1975) for 
    !          plasmas.
    implicit none
    real(PR), intent(in)                    :: A, Z, ni, Te, Ti
    real(PR)                                :: mi, mmu, G_RT
    real(PR)                                :: Omega_ei_hot, Omega_ei_cold, Omega_ei
    mi = A*mu
    mmu = Chemical_potential(Z, ni, Te)
    Omega_ei_hot = (3._PR*me/mi) * kB * zeff(Z, ni, Te) * ni * nu_ei(Z, ni, Te, Ti)
    ! Metals
    ! (The factor 1/5 has been chosen for patching the solid and plasma states)
    if (Z == Z_Al) then
      G_RT  = 3.e18_PR
      omega_ei_cold = G_RT * (1._PR + ( ( Te + Ti )/(Tfermi(Z, ni)/5._PR) ))
      Omega_ei = ( (omega_ei_cold**(-2._PR)) +  (Omega_ei_hot**(-2._PR)) )**(-0.5_PR)
    elseif (Z == Z_Cu) then
      G_RT  = 1.e18_PR
      omega_ei_cold = G_RT * (1._PR + ( ( Te + Ti )/(Tfermi(Z, ni)/5._PR) ))
      Omega_ei = ( (omega_ei_cold**(-2._PR)) +  (Omega_ei_hot**(-2._PR)) )**(-0.5_PR)
    elseif (Z == Z_Ta) then
      G_RT  = 1.5e18_PR
      omega_ei_cold = G_RT * (1._PR + ( ( Te + Ti )/(Tfermi(Z, ni)/5._PR) ))
      Omega_ei = ( (omega_ei_cold**(-2._PR)) +  (Omega_ei_hot**(-2._PR)) )**(-0.5_PR)
    ! Plasma
    else
      Omega_ei = Omega_ei_hot
    end if
  end function Omega_ei

  !=====================================================================================
  !                            Electrical Resitivity for metals
  !=====================================================================================

  elemental function gama_E(Z, ni, Te)
    ! input  : Z       = Atomic number of the material in ()
    !          ni      = Ion density in the material in (/cm^3)
    ! output : gama_E = Electron-electron collision correction factor for the electrical 
    !                    conductivity in () according to a simple fit of results given by
    !                    Spitzer L. and Härm R., Phys. Rev. 89, 977 (1953)    
    implicit none
    real(PR), intent(in) :: Z, ni, Te
    real(PR)             :: gama_E
    gama_E = (zeff(Z, ni, Te) + 0.9833_PR) / (zeff(Z, ni, Te) + 2.4101_PR)
  end function gama_E

  elemental function nu_Spitzer_res(Z, ni, Te, Ti)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Te  = Electron Temperature in the material in (K)
    !          Ti  = Ion Temperature in the material in (K)
    ! output nu_Spitzer_res = Electron collision frequency in a plasma accounting for the
    !                         electron-electron collision correction factor for the  
    !                         electrical conductivity in (/s)
    implicit none
    real(PR), intent(in) :: Z, ni, Te, Ti
    real(PR)             :: nu_Spitzer_res
    nu_Spitzer_res = nu_ei(Z, ni, Te, Ti) / gama_E(Z, ni, Te)
  end function nu_spitzer_res

  elemental function nu_res(Z, ni, Te, Ti)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Te  = Electron Temperature in the material in (K)
    !          Ti  = Ion Temperature in the material in (K)
    ! output : nu_res = harmonic mean of all electron collision frequencies allowing 
    !                   to relate all regimes for the solid to the plasma state in (/s)  
    implicit none
    real(PR), intent(in) :: Z, ni, Te, Ti
    real(PR)             :: nu_res
    nu_res = (&
    & (( &
    &    nu_eph(Z, ni, Ti) &
    &  + nu_ee(Z, ni, Te)  & 
    &  )**(-2._PR)) &
    & + (nu_c(Z, ni, Te)**(-2._PR)) &
    & + (nu_Spitzer_res(Z, ni, Te, Ti)**(-2._PR))&
    & )**(-0.5_PR)
  end function nu_res

  elemental function A_alpha(Z, ni, Te)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Te  = Electron Temperature in the material in (K)
    ! output : A_alpha = Function without dimension depending on the Fermi integrales given 
    !                    by Lee Y. and More R. , Physics of Fluids Vol. 27 p. 1273 (1984)   
    implicit none
    real(PR), intent(in) :: Z, ni, Te
    real(PR)             :: mu, F_undemi, F_deux, aa, bb, A_alpha
    mu=Chemical_potential(Z, ni, Te)/(kB*Te)
    F_undemi=Fermi_integrale(0.5_PR, Z, ni, Te)
    F_deux = Fermi_integrale(2._PR, Z, ni, Te)
    aa = F_deux / F_undemi
    bb = (1._PR + exp(-mu)) * F_undemi
    A_alpha = (4._PR/3._PR) * aa / bb
  end function A_alpha

  elemental function resistivity_metals(Z, ni, Te, Ti)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Te  = Electron Temperature in the material in (K)
    !          Ti  = Ion Temperature in the material in (K)
    ! output : resistivity_metals = electrical resistivity for metals in (s)
    implicit none
    real(kind=PR), intent(in) :: Z, ni, Te, Ti
    real(kind=PR)             :: resistivity_metals
    resistivity_metals = me * nu_res(Z, ni, Te, Ti) /&
    & (zeff(Z, ni, Te) * ni * (e**2._PR) * A_alpha(Z, ni, Te))
  end function resistivity_metals

  !=====================================================================================
  !                      Thermal electron conductivity for metals
  !=====================================================================================

  elemental function eps_delta_T(Z, ni, Te, Ti)
    ! input  : Z           = Atomic number of the material in ()
    !          ni          = Ion density in the material in (/cm^3)
    ! output : eps_delta_T = Electron-electron collision correction factor for the  
    !                        electron thermal conductivity in () according to 
    !                        Brysk H. et al., Plasma phys. 17, 6, p. 2714 (1975)
    implicit none
    real(PR), intent(in) :: Z, ni, Te, Ti
    real(PR)             :: zn
    real(PR)             :: eps_delta_T
    zn = (2._PR**(-2.5_PR)) * zeff(Z, ni, Te) * log_clas_ei(Z, ni, Te, Ti) &
    &/ log_clas_ee(Z, ni, Te)
    eps_delta_T = (15._PR*pi/256._PR) * ( (45._PR*zn) + (433._PR*(zn**2._PR)) ) &
    & / (9._PR + (151._PR*zn) + (217._PR*(zn**2._PR)) )
  end function eps_delta_T

  elemental function nu_Spitzer_cond(Z, ni, Te, Ti)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Te  = Electron Temperature in the material in (K)
    !          Ti  = Ion Temperature in the material in (K)
    ! output : nu_Spitzer_cond = Electron collision frequency in a plasma accounting 
    !                            for the electron-electron collision correction factor 
    !                            in the expression of the electron thermal conductivity 
    !                            in (/s)
    implicit none
    real(PR), intent(in) :: Z, ni, Te, Ti
    real(PR)             :: nu_Spitzer_cond
    nu_Spitzer_cond = nu_ei(Z, ni, Te, Ti) / eps_delta_T(Z, ni, Te, Ti)
  end function nu_spitzer_cond

  elemental function nu_cond(Z, ni, Te, Ti)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Te  = Electron Temperature in the material in (K)
    !          Ti  = Ion Temperature in the material in (K)
    ! output : nu_cond = harmonic mean of all electron collision frequencies allowing 
    !                    to relate all regimes for the solid to the plasma state in (/s) 
    implicit none
    real(PR), intent(in) :: Z, ni, Te, Ti
    real(PR)             :: nu_cond
    nu_cond = ( ((nu_eph(Z, ni, Ti)+nu_ee(Z, ni, Te) )**(-2._PR)) &
    & + (nu_c(Z, ni, Te)**(-2._PR)) &
    & + (nu_Spitzer_cond(Z, ni, Te, Ti)**(-2._PR)) )**(-0.5_PR)
  end function nu_cond

  elemental function A_beta(Z, ni, Te)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Te  = Electron Temperature in the material in (K)
    ! output : A_beta = Function without dimension depending on the Fermi integrales given 
    !                   by Lee Y. and More R. , Phys. of Fluids Vol. 27 p. 1273 (1984)   
    implicit none
    real(PR), intent(in) :: Z, ni, Te
    real(PR)             :: mu, F_undemi, F_deux, F_trois, F_quatre, aa, bb, cc, dd, A_beta
    mu=Chemical_potential(Z, ni, Te)/(kB*Te)
    F_undemi=Fermi_integrale(0.5_PR, Z, ni, Te)
    F_deux = Fermi_integrale(2._PR, Z, ni, Te)
    F_trois = Fermi_integrale(3._PR, Z, ni, Te)
    F_quatre = Fermi_integrale(4._PR, Z, ni, Te)
    aa = F_quatre / F_undemi
    bb = (1._PR + exp(-mu)) * F_undemi
    cc = (15._PR*F_quatre*F_deux-16._PR*(F_trois**2._PR))
    dd = cc / (15._PR*F_quatre*F_deux)
    A_beta=(20._PR/9._PR) * aa * dd / bb
  end function A_beta

  elemental function gama_Lorenz(Z, ni, Te)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Te  = Electron Temperature in the material in (K)
    ! output : gama_Lorenz = Lorenz factor kappa / (sigma Te) for a Lorenzian plasma
    !                         where kappa is the electron thermal conductivity
    !                               sigma the electrical conductivity and
    !                               Te(K) the electron Temperature 
    implicit none
    real(PR), intent(in) :: Z, ni, Te
    real(PR)             :: gama_Lorenz
    gama_Lorenz = (A_beta(Z, ni, Te)/A_alpha(Z, ni, Te)) * ((kb/e)**2._PR) 
  end function gama_Lorenz

  elemental function conductivity_metals(Z, ni, Te, Ti)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Te  = Electron Temperature in the material in (K)
    !          Ti  = Ion Temperature in the material in (K)
    ! output : conductivity_metals = thermal electron conductivity for metals in 
    !                                (erg/cm/K/s)
    implicit none
    real(kind=PR), intent(in) :: Z, ni, Te, Ti
    real(PR)                  :: AA,BB
    real(kind=PR)             :: conductivity_metals
    AA = gama_Lorenz(Z, ni, Te) * Te * zeff(Z, ni, Te) &
    & * ni * (e**2._PR) * A_alpha(Z, ni, Te)
    BB = (me * nu_cond(Z, ni, Te, Ti) )
    conductivity_metals = AA / BB
  end function conductivity_metals

  !=====================================================================================
  !                      Electrical Resitivity for Hydrogen plasmas
  !=====================================================================================

  elemental function resistivity_Spitzer(Z, ni, Te, Ti)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Te  = Electron Temperature in the material in (K)
    !          Ti  = Ion Temperature in the material in (K)
    ! output : resistivity_Spitzer = electrical resistivity for metals in (s)
    !                                according to 
    !                                Spitzer L. and Härm R., Phys. Rev. 89, 977 (1953) 
    implicit none
    real(PR),intent(in) :: Z, ni, Te, Ti
    real(PR)            :: nu, Zf, resistivity_Spitzer
    Zf = Zeff(Z, ni, Te)
    nu = nu_Spitzer_res(Z, ni, Te, Ti)
    resistivity_Spitzer = me * nu / (Zf * ni * (e**2._PR) * A_alpha(Z, ni, Te))
  end function resistivity_Spitzer

  elemental function cond_Spitzer(Z, ni, temp_e, temp_i)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Te  = Electron Temperature in the material in (K)
    !          Ti  = Ion Temperature in the material in (K)
    ! output : cond_Spitzer = electron thermal resistivity for metals in (erg/cm/K/s)
    !                                according to 
    !                                Spitzer L. and Härm R., Phys. Rev. 89, 977 (1953) 
    implicit none
    real(PR), intent(in) :: Z, ni, temp_e, temp_i
    real(PR)             :: cond_Spitzer
    real(PR)             :: Zf, Lorentz_factor
    Zf = Zeff(Z, ni, temp_e)
    Lorentz_factor = 4._PR * (Zf/(Zf+0.2_PR*log(Zf)+3.44_PR))*&
    ((Zf+2.2_PR)/(Zf+0.7_PR))
    cond_Spitzer = Lorentz_factor * ((kb/e)**2._PR) * temp_e /&
    resistivity_Spitzer(Z, ni, temp_e, temp_i)
  end function cond_Spitzer

  elemental function gama_ii(Z, ni, Ti)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Ti  = Ion Temperature in the material in (K)
    ! output : gama_ii = Dimensionless function given by
    !                     Hubbard W., Astrophys. Jour. 146, 858 (1966)
    implicit none
    real(PR),intent(in) :: Z, ni,Ti
    real(PR)            :: gama_ii
    gama_ii = ((Z * e)**2._PR) / (kb * Ti * lambda_i(ni))
  end function gama_ii

  elemental function G_gama(gama, kappa)
    ! input  : gama = dimensionless function (cf. gama_ii)
    !          kappa = dimensionless wavelength number
    ! output : G_gama = Dimensionless function given by
    !                    Hubbard W., Astrophys. Jour. 146, 858 (1966)
    implicit none
    real(PR), intent(in) :: gama, kappa
    real(PR)             :: G_gama
    ! The expression is valid only 
    ! if kappa>=3/2 and kappa<=2 i.e. only for hydrogen plasmas
    G_gama = (&
    ((log(1._PR + (4._PR*(kappa**2._PR)/(3._PR*gama))))**(-1._PR))**(-2._PR) +&
    (log( ((gama+5._PR)**(5._PR/12._PR)) - 0.5_PR ))**(-2._PR)&
    )**(-0.5_PR)
  end function G_gama

  elemental function S_11(Z, ni, Ti)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Ti  = Ion Temperature in the material in (K)
    ! output : S_11 = diagonal component of transport coefficient tensors
    !                 corresponding to the electrical conductivity in (/s/statcoulomb^2) 
    !                 given by Hubbard W., Astrophys. Jour. 146, 858 (1966)
    implicit none
    real(PR), intent(in) :: Z, ni, Ti
    real(PR)             :: kappa_F, C1
    real(PR)             :: S_11
    C1 = ((2._PR*pi*hbar)**3._PR) * ni /&
    (27._PR * (pi**4._PR) * (me**2._PR)* (e**4._PR) * (Z**2._PR))
    kappa_F = (9._PR * pi * Z / 4._PR)**(1._PR/3._PR)
    S_11 = C1 * (Kappa_F**6._PR) * G_gama(gama_ii(Z, ni, Ti), kappa_F)
  end function S_11

  elemental function resistivity_Hubbard(Z, ni, Ti)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Ti  = Ion Temperature in the material in (K)
    ! output : resistivity_Hubbard = electrical resistivity in (s) 
    !          according to Hubbard W., Astrophys. Jour. 146, 858 (1966)
    implicit none
    real(PR),intent(in) :: Z, ni, Ti
    real(PR)            :: resistivity_Hubbard
    resistivity_hubbard  = (S_11(Z, ni, Ti) * (e**2._PR))**(-1._PR) 
  end function resistivity_Hubbard

  elemental function cond_Hubbard(Z, ni, Te, Ti)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Ti  = Ion Temperature in the material in (K)
    ! output : cond_Hubbard = electron thermal conductivity in (erg/cm/K/s) 
    !          according to Hubbard W., Astrophys. Jour. 146, 858 (1966)
    implicit none
    real(PR), intent(in) :: Z, ni, Te, Ti
    real(PR)             :: cond_Hubbard
    cond_Hubbard = ((pi**2._PR)/3._PR)*(kB**2._PR)*S_11(Z, ni, Ti)*Te
  end function cond_Hubbard

  elemental function resistivity_H(Z, ni, Te, Ti)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Te  = Electron Temperature in the material in (K)
    !          Ti  = Ion Temperature in the material in (K)
    ! output : resistivity_H = electrical resistivity for hydrogen plasma in (s)
    !          assuming an harmonic mean between the Hubbard and Spitzer's models
    implicit none
    real(PR),intent(in) :: Z, ni, Te, Ti
    real(PR)            :: cold, hot, resistivity_H
    cold = resistivity_Hubbard(Z, ni, Ti)**(-2._PR)
    hot  = resistivity_Spitzer(Z, ni, Te, Ti)**(-2._PR)
    resistivity_H = (cold + hot)**(-0.5_PR)
  end function resistivity_H

  elemental function cond_H(Z, ni, Te, Ti)
    ! input  : Z   = Atomic number of the material in ()
    !          ni  = Ion density in the material in (/cm^3)
    !          Te  = Electron Temperature in the material in (K)
    !          Ti  = Ion Temperature in the material in (K)
    ! output : cond_H = electrical resistivity for hydrogen plasma in (erg/cm/s/K)
    !          assuming an harmonic mean between the Hubbard and Spitzer's models
    implicit none
    real(PR),intent(in) :: Z, ni, Te, Ti
    real(PR)            :: cold, hot, cond_H
    cold = cond_Hubbard(Z, ni, Te, Ti)**(2._PR)
    hot  = cond_Spitzer(Z, ni, Te, Ti)**(2._PR)
    cond_H = (cold + hot)**(0.5_PR)
  end function cond_H

  !=====================================================================================
  !                               Tabulated Resitivity
  !=====================================================================================

  subroutine initialize_resistivity(eta_tab,N_eta_tab)
    ! output  : eta_tab = tabulated values the electrical resistivity if chosen by the user
    implicit none
    real(PR), dimension(:,:), allocatable, intent(inout) :: eta_tab
    integer, intent(out)                                 :: N_eta_tab
    integer                                              :: reason, i
    character                                            :: str
    ! Find the number of lines in the file 'sources/user/resistivity_tab.dat' and allocate the
    ! table eta_tab :
    reason = 0
    i = -2
    open (unit=700,file='sources/user/resistivity_tab.dat',form='formatted',status='unknown')      
    do while (reason==0)
      read(700,*,IOSTAT=Reason) str
      i = i + 1
      if (str == '') reason = 1
    end do
    close(700)
    N_eta_tab = i
    allocate(eta_tab(1:N_eta_tab,1:2))
    !  Read the file 'sources/user/resistivity_tab.dat' and define the table eta_tab with the 
    !  corresponding values
    open (unit=40,file='sources/user/resistivity_tab.dat',form='formatted',status='unknown')      
    read(40,*)
    do i =1,N_eta_tab,1
      read(40,*) eta_tab(i,1),eta_tab(i,2)
    end do
    close(40)
  end subroutine initialize_resistivity

  function resistivity_tab(eta_tab, Z, ni, Te, Ti)
    ! input : eta_tab = tabulated values the electrical resistivity if chosen by the user
    !         Z       = atomic number of the material in ()
    !         ni      = ion density of the material in (/cm^3)
    !         Te      = electron temperature in the material in (K)
    !         Ti      = ion temperature in the material in (K)
    implicit none
    real(PR), dimension(1:N_eta_tab,1:2), intent(in) :: eta_tab
    real(PR),intent(in)                              :: Z, ni, Te, Ti
    real(PR)                                         :: resistivity_tab
    integer                                          :: il, ir, im
    logical                                          :: not_found
    real(PR)                                         :: alpha, beta
    if (Te > eta_tab(N_eta_tab,1)) then
      resistivity_tab = resistivity_Spitzer(Z, ni, Te, Ti)
    else if (Te < eta_tab(1,1)) then
      resistivity_tab = eta_tab(1,2)
      ! the corresponding index is found by dichotomy
    else
      not_found = .true.
      il = 1
      ir = N_eta_tab
      im = (il+ir)/2
      do while (not_found)
        if ((eta_tab(im,1)) <= Te) then
          il = im
        else
          ir = im
        end if
        im = (il+ir)/2
        not_found = (ir - il) > 1   
      end do
      if ((im.ne.N_eta_tab).and.(eta_tab(im,2).ne.eta_tab(im+1,2))) then 
        ! linear interpolation
        alpha = eta_tab(im+1,1) - Te
        beta  = Te - eta_tab(im,1)
        resistivity_tab = ((alpha*eta_tab(im,2))+(beta*eta_tab(im+1,2)))&
                        /(eta_tab(im+1,1)-eta_tab(im,1))
      else
        resistivity_tab = eta_tab(im,2)
      end if
    end if
  end function resistivity_tab

  !=====================================================================================
  !                              Electrical Resitivity
  !=====================================================================================

  function resis(eta_tab, Z, ni, Te, Ti)
    ! input  : eta_tab = tabulated values of the electrical resistivity if chosen by the user
    !          Z      = Atomic number of the material in ()
    !          ni     = Ion density in the material in (/cm^3)
    !          Te     = Electron Temperature in the material in (K)
    !          Ti     = Ion Temperature in the material in (K)
    ! output : resis = electrical resistivity of the background in (s)
    implicit none
    real(PR), dimension(1:N_eta_tab,1:2), intent(in) :: eta_tab
    real(PR), intent(in)                             :: Z, ni, Te, Ti
    real(PR)                                         :: resis
    ! Metals
    if ((Z == Z_Al).or.(Z == Z_Cu).or.(Z == Z_Ta)) then
      resis = resistivity_metals(Z, ni, Te, Ti)
    ! Hydrogen plasmas
    else if (Z == 1) then
      resis = resistivity_H(Z, ni, Te, Ti)
    ! Other plasmas
    else 
      if (tabulated_resistivity) then
        resis = resistivity_tab(eta_tab, Z, ni, Te, Ti)
      else
        resis = resistivity_Spitzer(Z, ni, Te, Ti)
      end if
    end if
  end function resis

  !=====================================================================================
  !                            Thermal electron Conductivity
  !=====================================================================================

  elemental function cond(Z, ni, Te, Ti)
    ! input  : eta_tab = tabulated values if chosen by the user
    !          Z       = Atomic number of the material in ()
    !          ni      = Ion density in the material in (/cm^3)
    !          Te      = Electron Temperature in the material in (K)
    !          Ti      = Ion Temperature in the material in (K)
    ! output : cond = thermal electron conductivity of the background in (erg/cm/s/K)
    implicit none
    real(PR), intent(in) :: Z, ni, Te, Ti
    real(PR)             :: cond
    !
    if (tabulated_resistivity) then
      cond = zero
    else
      ! Metals
      if ((Z == Z_Al).or.(Z == Z_Cu).or.(Z == Z_Ta)) then
        cond = conductivity_metals(Z, ni, Te, Ti)
      ! Hydrogen plasmas   
      else if(Z == 1) then
        cond = cond_H(Z, ni, Te, Ti)
      ! Other plasmas
      else
        cond = cond_Spitzer(Z, ni, Te, Ti)
      end if
    end if
  end function cond

end module transport_coef
