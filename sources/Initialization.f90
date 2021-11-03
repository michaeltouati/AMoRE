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
module initialization

  use acuracy
  use constants
  use physics_library
  use fokker_planck_coef
  use input

  implicit none
  public :: calcul_Vb0_Intfeps_nrjmean_Intfx
  public :: initialize_spectrum,function_energy, function_space, function_theta, function_time
  public :: initialize_phi, initialize, initialize_next_step

  contains

  !=========================================================================================
  !  Initialization of angular moments of the injected fast e- beam distribution function
  !                                      at z = 0
  !=========================================================================================

  subroutine initialize_spectrum(eps_tab,N_eps_tab)
    ! output  : eps_tab = tabulated values of the electron beam kinetic energy spectrum if chosen by the user
    implicit none
    real(PR), dimension(:,:), allocatable, intent(out) :: eps_tab
    integer, intent(out)                               :: N_eps_tab
    integer                                            :: reason, i
    character                                          :: str
    ! Find the number of lines in the file 'sources/user/spectrum_tab.dat' and allocate the
    ! table eps_tab :
    reason = 0
    i = -2   
    open (unit=500,file='sources/user/spectrum_tab.dat',form='formatted',status='unknown')  
    do while (reason==0)
      read(500,*,IOSTAT=Reason) str
      i = i + 1
      if (str == '') reason = 1
    end do
    close(500)
    N_eps_tab = i
    allocate(eps_tab(1:N_eps_tab,1:2))
    !  Read the file 'sources/user/spectrum_tab.dat' and define the table eps_tab with the 
    !  corresponding values
    open (unit=400,file='sources/user/spectrum_tab.dat',form='formatted',status='unknown')      
    read(400,*)
    do i =1,N_eps_tab,1
      read(400,*) eps_tab(i,1),eps_tab(i,2)
    end do
    close(400)
    ! Conversion MeV -> keV
    do i =1,N_eps_tab,1
      eps_tab(i,1) = eps_tab(i,1)*1.e3_PR
      eps_tab(i,2) = eps_tab(i,2)/1.e3_PR
    end do
  end subroutine initialize_spectrum

  function spectrum_tab(eps_tab, nrj)
    ! input : eps_tab = tabulated values of the electron beam kinetic energy spectrum if chosen by the user
    !         nrj     = beam electrons kinetic energy in (keV)
    implicit none
    real(PR), dimension(1:N_eps_tab,1:2), intent(in) :: eps_tab
    real(PR)                                         :: nrj
    !
    real(PR)                                         :: spectrum_tab
    integer                                          :: il, ir, im
    logical                                          :: not_found
    real(PR)                                         :: alpha, beta
    if (nrj > eps_tab(N_eps_tab,1)) then
      spectrum_tab = 0.
    else if (nrj < eps_tab(1,1)) then
      spectrum_tab = 0.
    ! the corresponding index is found by dichotomy
    else
      not_found = .true.
      il = 1
      ir = N_eps_tab
      im = (il+ir)/2
      do while (not_found)
        if ((eps_tab(im,1)) <= nrj) then
          il = im
        else
          ir = im
        end if
        im = (il+ir)/2
        not_found = (ir - il) > 1   
      end do
      if ((im.ne.N_eps_tab).and.(eps_tab(im,2).ne.eps_tab(im+1,2))) then 
        ! linear interpolation
        alpha = eps_tab(im+1,1) - nrj
        beta  = nrj - eps_tab(im,1)
        spectrum_tab = ((alpha*eps_tab(im,2))+(beta*eps_tab(im+1,2)))&
                        /(eps_tab(im+1,1)-eps_tab(im,1))
      else
        spectrum_tab = eps_tab(im,2)
      end if
    end if
  end function spectrum_tab

  function function_energy(eps_tab, nrj)
    ! input  : eps_tab = tabulated values of the electron beam kinetic energy spectrum 
    !                    if chosen by the user
    !          nrj = kinetic energy of fast e- in (keV)
    ! output : function_energy = fast e- beam spectrum in (/keV)
    implicit none
    real(PR), dimension(1:N_eps_tab,1:2), intent(in) :: eps_tab
    real(PR),intent(in)                              :: nrj
    real(PR)                                         :: function_energy    
    select case (spectrum)
    case (1)
      function_energy = exp(- (nrj-eps0)**2._PR/(2*(sigma_eps**2._PR)) ) &
      &/ sqrt(2._PR*pi*(sigma_eps**2._PR))
    case default ! case (2)
      function_energy = exp(- nrj / Tb0 ) / Tb0
    case (3)
      function_energy = (exp(- nrj / Tb0 ) / nrj) + ( alpha1 * exp(- nrj / Tb1 ) )
    case (4)
      function_energy = (exp(- nrj / Tb0 ) * (eps1/ nrj) * ((eps0/nrj)**alpha0) ) &
      &               + ( alpha1 * exp(- nrj / Tb1 ) )
    case (5)
      function_energy = spectrum_tab(eps_tab, nrj)
    end select
  end function function_energy

  elemental function function_space(r)
    ! input  : r = transverse position in (microns)
    ! output : function_space = fast e- beam transverse distribution in (/cm)
    implicit none
    real(PR),intent(in)   :: r
    real(PR)              :: function_space
    function_space = exp(- r**2._PR/(2*(sigma_r**2._PR))) / sqrt(2._PR*pi*((sigma_r*microns)**2._PR))
  end function function_space

  elemental function function_theta(r)
    ! input  : r = transverse position in (microns)
    ! output : function_theta = angle between the mean propagation direction of fast e-
    !                           and the longitudinal z-axis in (rad)
    implicit none
    real(PR),intent(in)   :: r
    real(PR)              :: function_theta    
    function_theta = (angle0*pi/180._PR)*tanh(r/sigma_r)    
  end function function_theta

  elemental function function_time(Vb, t)
    ! input  : Vb = injected fast e- beam kinetic energy flux mean velocity in (cm/s)
    !          t  = time in (fs)
    ! output : function_time = fast e- beam longitudinal distribution in (/cm)
    implicit none
    real(PR),intent(in)   :: Vb, t
    real(PR)              :: function_time
    function_time = exp(- (t-t_c)**2._PR/(2*(sigma_t**2._PR)) ) / (Vb*sqrt(2._PR*pi*((sigma_t*fs)**2._PR)))
  end function function_time

  subroutine calcul_Vb0_Intfeps_nrjmean_Intfx(eps_tab, xa, epsa, Vb, Intfeps, nrj_mean, Intfx)
    ! input  : eps_tab  = tabulated values of the electron beam kinetic energy spectrum if chosen by the user
    !          xa       = array containing all transverse position x(i) in (microns)
    !          epsa     = array containing all kinetic energies epsa(l) in (keV)
    ! output : Vb       = injected fast e- beam kinetic energy flux mean velocity in (cm/s)
    !          Intfeps  = Integrale of injected fast e- beam spectrum over all kinetic energies
    !                     in (keV)
    !          nrj_mean = mean kinetic energy of fast electron beam injected in the simulation 
    !                     box at z = 0 in (keV)
    !          Intfx    = Integrale of injected fast e- beam transverse distribution over all 
    !                     transverse position in (keV)
    implicit none
    real(PR), dimension(1:N_eps_tab,1:2), intent(in) :: eps_tab
    real(PR), dimension(1:N_x), intent(in)           :: xa
    real(PR), dimension(1:N_eps), intent(in)         :: epsa
    real(PR), intent(out)                            :: Vb, Intfeps, nrj_mean, Intfx
    real(PR)                                         :: vit_b, a2, omega0, theta
    integer                                          :: l,i
    vit_b    = 0._PR
    Intfeps  = 0._PR
    Intfx    = 0._PR
    nrj_mean = 0._PR
    l = 1
    if (sigma_theta.ne.0._PR) then
      a2 = 1._PR / ((pi*sigma_theta/180._PR)**2._PR) 
      omega0 = (1._PR/tanh(a2)) - (1._PR /a2)
    else
      omega0 = 1._PR
    end if
    do i=1,N_x,1
      Intfx = Intfx + (function_space(xa(i)) * (d_x*microns))
    end do
    do l=1,N_eps,1
      Intfeps = Intfeps + (function_energy(eps_tab,epsa(l)) * d_eps)
      nrj_mean =  nrj_mean + (epsa(l) * function_energy(eps_tab,epsa(l)) * d_eps)
      do i=1,N_x,1
        theta   = function_theta(xa(i))
        vit_b = vit_b + (vit(epsa(l)) * omega0 * cos(theta) * epsa(l) * function_space(xa(i)) &
                      & * function_energy(eps_tab,epsa(l))  * d_eps * (d_x*microns))
      end do
    end do
    nrj_mean = nrj_mean / Intfeps
    nrj_mean = nrj_mean / Intfx
    vit_b    = vit_b    / Intfeps
    vit_b    = vit_b    / Intfx
    Vb       = vit_b    / nrj_mean
  end subroutine calcul_Vb0_Intfeps_nrjmean_Intfx

  subroutine initialize_phi(eps_tab, time, dt, xa, epsa, phin, phinp1, norm)
    ! input  : eps_tab = tabulated values of the electron beam kinetic energy spectrum if chosen by the user
    !          time   = time in (fs)
    !          dt     = time step in (fs)
    !          xa     = array containing all transverse position x(i) in (microns)
    !          epsa   = array containing all kinetic energies epsa(l) in (keV)
    ! output : phin   = table containing the 0th and 1st order angular moment components
    !                   of the fast e- beam injected in () at z = 0 and t=time 
    !          phinp1 = table containing the 0th and 1st order angular moment components
    !                   of the fast e- beam injected in () at z = 0 and t=time+dt
    !          norm   = coefficient of normalization for the distribution function
    !                   in (/cm^3/keV)
    implicit none
    real(PR), dimension(1:N_eps_tab,1:2), intent(in)                                                :: eps_tab
    real(PR), intent(in)                                                                            :: time, dt
    real(PR), dimension(1:N_x), intent(in)                                                          :: xa
    real(PR), dimension(1:N_eps), intent(in)                                                        :: epsa
    real(PR), intent(out)                                                                           :: norm
    real(PR), dimension(forward:backward,psi0:psi1z, -1:N_x+2, -1:N_z+2, -1:N_eps+2), intent(inout) :: phin
    real(PR), dimension(forward:backward,psi0:psi1z, -1:N_x+2, -1:N_z+2, -1:N_eps+2), intent(inout) :: phinp1
    integer   :: l, i
    real(PR)  :: third_dim_norm
    real(PR)  :: N0, Vb, a2
    real(PR)  :: omega0
    real(PR)  :: theta,Int_feps
    real(PR)  :: nrj_mean,Int_fx
    integer   :: shape
    shape = 0 ! shape = 1 for the particular case of L. Volpe hyperbolic secant distribution profile 
    call calcul_Vb0_Intfeps_nrjmean_Intfx(eps_tab,xa,epsa,Vb,Int_feps,nrj_mean,Int_fx)
    ! total number of fast e- beam injected in the simulation box in (/cm) (2D simulation)
    if (shape==1) then
      third_dim_norm = sqrt( pi*Catalan*((sigma_r*microns)**2._PR))
    else  
      third_dim_norm = sqrt(2._PR*pi*((sigma_r*microns)**2._PR))
    end if
    N0 = E_tot*Joules / (nrj_mean*keV*third_dim_norm)
    norm = N0 * function_energy(eps_tab,eps_min)
    norm = norm * function_time(Vb,t_c)
    norm = norm * function_space(0._PR)
    if (time <= source_duration) then       
      if (sigma_theta.ne.0._PR) then
        a2 = 1._PR / ((pi*sigma_theta/180._PR)**2._PR)
        omega0 = (1._PR/tanh(a2)) - (1._PR /a2)
      else
        omega0 = 1._PR
      end if
      !$omp PARALLEL DO PRIVATE(l,i,theta) COLLAPSE(2)
      do l=1,N_eps,1
        do i=1,N_x,1   
          theta   = function_theta(xa(i))
          phin(forward,psi0,i,0,l)     = N0 * function_energy(eps_tab,epsa(l))/Int_feps
          phin(forward,psi0,i,0,l)     = phin(forward,psi0,i,0,l) * function_time(Vb, time)
          phin(forward,psi0,i,0,l)     = phin(forward,psi0,i,0,l) * function_space(xa(i))/Int_fx
          phin(forward,psi0,i,0,l)     = phin(forward,psi0,i,0,l) / norm
          phin(forward,psi1x,i,0,l)    = omega0 * sin(theta) * phin(forward,psi0,i,0,l) 
          phin(forward,psi1z,i,0,l)    = omega0 * cos(theta) * phin(forward,psi0,i,0,l) 
          phinp1(forward,psi0,i,0,l)   = N0 * function_energy(eps_tab,epsa(l))/Int_feps
          phinp1(forward,psi0,i,0,l)   = phinp1(forward,psi0,i,0,l) * function_time(Vb, time + dt)
          phinp1(forward,psi0,i,0,l)   = phinp1(forward,psi0,i,0,l) * function_space(xa(i))/Int_fx
          phinp1(forward,psi0,i,0,l)   = phinp1(forward,psi0,i,0,l) / norm
          phinp1(forward,psi1x,i,0,l)  = omega0 * sin(theta) * phinp1(forward,psi0,i,0,l) 
          phinp1(forward,psi1z,i,0,l)  = omega0 * cos(theta) * phinp1(forward,psi0,i,0,l) 
        end do
      end do
    else
      phin(forward,:,:,0,:)   = 0._PR
      phinp1(forward,:,:,0,:) = 0._PR
    end if
  end subroutine initialize_phi

  !=========================================================================================
  !   Initialization of all physical quantities at t=0 and at the end of each time step
  !=========================================================================================

  subroutine get_A_Z_ni(A1,Z1,ni1)
    implicit none
    real(PR), intent(out) :: A1, Z1, ni1
    select case (Material)
      case (1)
        A1  = A0
        Z1  = Z0
        ni1 = rho/(mu*A0)
      case (2)
        select case (Mat)
          case ('Al')
            ni1 = ni_Al
            A1  = A_Al
            Z1  = Z_Al
          case ('Cu')
            ni1 = ni_Cu
            A1  = A_Cu
            Z1  = Z_Cu
          case ('Ta')
            ni1 = ni_Ta
            A1  = A_Ta
            Z1  = Z_Ta
          case default
            ni1 = rho/(mu*A0)
            A1  = A0
            Z1  = Z0
        end select
      case default
        ni1 = rho/(mu*A0)
        A1  = A0
        Z1  = Z0
    end select
  end subroutine get_A_Z_ni

  subroutine initialize(x_a, z_a, eps_a, ni, A, Z, Ex, Ez, By_n, By_np1,&
                        nb, jbx_n, jbz_n, jrx_n, jrz_n,jbx_np1, jbz_np1, jrx_np1, jrz_np1,&
                        pdepos, plost_col, plost_res, Te, Ti,&
                        nK_Holes, nKalpha, nKbeta,ioniztime,&
                        ud, Udcol, Udres, ue, usf, usb, usu, usd, ub, uel, uma,Sva, Sva_lp1, nua)
    implicit none
    real(PR), dimension(1:N_x), intent(in)                :: x_a
    real(PR), dimension(1:N_z), intent(in)                :: z_a
    real(PR), dimension(1:N_eps), intent(in)              :: eps_a
    real(PR), intent(out)                                 :: Ud, udcol, udres, ue, usf
    real(PR), intent(out)                                 :: usb, usu, usd, ub, uel, uma
    real(PR), intent(out), dimension(1:N_x,1:N_z)         :: Ex, Ez, By_n, By_np1
    real(PR), intent(out), dimension(1:N_x,1:N_z)         :: nb, jbx_n, jbz_n, jrx_n, jrz_n
    real(PR), intent(out), dimension(1:N_x,1:N_z)         :: jbx_np1, jbz_np1, jrx_np1, jrz_np1
    real(PR), intent(out), dimension(1:N_x,1:N_z)         :: pdepos, plost_col, plost_res
    real(PR), intent(out), dimension(1:N_x,1:N_z)         :: nK_Holes, nKalpha, nKbeta, ioniztime
    real(PR), intent(out), dimension(0:N_x+1,0:N_z+1)     :: Te, Ti
    real(PR), intent(out), dimension(0:N_x+1,0:N_z+1)     :: ni, A, Z
    real(PR), dimension(1:N_x,1:N_z,1:N_eps), intent(out) :: Sva, Sva_lp1, nua
    integer                                               :: i, k, l, reason, N_plasma_tab
    integer                                               :: Nx_plasma_tab, Nz_plasma_tab
    real(PR), dimension(:,:), allocatable                 :: plasma_tab
    integer                                               :: il, ir, im, kl, kr, km
    logical                                               :: not_found
    real(PR)                                              :: alpha_x, alpha_z, r, t
    real(PR)                                              :: beta_x, beta_z
    real(PR)                                              :: A1, Z1, ni1
    character                                             :: str
    !
    call get_A_Z_ni(A1, Z1,ni1)
    !
    select case  (Material)
      case (1)
        if (tabulated_plasma) then
          !  Read the file 'sources/user/plasma_tab.dat' and find the number of lines :
          reason = 0
          i = -2
          open (unit=600,file='sources/user/plasma_tab.dat',form='formatted',status='unknown')      
          do while (reason==0)
            read(600,*,IOSTAT=Reason) str
            i = i + 1
            if (str == '') reason = 1
          end do
          close(600)
          N_plasma_tab = i
          allocate(plasma_tab(1:N_plasma_tab,1:4))
          !  Read the file 'sources/user/plasma_tab.dat' and define the table plasma_tab with the 
          !  corresponding values :
          open (unit=700,file='sources/user/plasma_tab.dat',form='formatted',status='unknown')      
          read(700,*)
          do i =1,N_plasma_tab,1
            read(700,*) plasma_tab(i,1),plasma_tab(i,2),plasma_tab(i,3),plasma_tab(i,4)
          end do
          close(700)
          ! find the number of cells on x-axis and deduce the
          !      the number of cells on z-axis in 'sources/user/plasma_tab.dat' :
          i = 1
          reason = 0
          do while (reason==0)
            if ((i.ne.1)) then
              if (plasma_tab(i,1).ne.plasma_tab(i-1,1)) then
                reason =1
              end if
            end if
            i = i + 1
          end do
          Nx_plasma_tab = i - 2
          Nz_plasma_tab = N_plasma_tab / Nx_plasma_tab
          !$omp PARALLEL DO DEFAULT(NONE) &
          !$omp& SHARED(N_z, N_x, Nz_plasma_tab, Nx_plasma_tab) &
          !$omp& SHARED(plasma_tab, x_a, z_a, A, A1, Z, Z1, ni, Te, Ti) &
          !$omp& PRIVATE(k,i,not_found,il,ir,im,kl,kr,km, alpha_x,beta_x,alpha_z,beta_z,r,t) &
          !$omp& COLLAPSE(2)
          do k=1,N_z,1
            do i=1,N_x,1
              ! find by dichotomy the closest cell :
              not_found = .true.
              il = 1
              ir = Nx_plasma_tab
              im = (il+ir)/2
              do while (not_found)
                if ((plasma_tab((im-1)*Nz_plasma_tab+1,1)) <= x_a(i)) then
                  il = im
                else
                  ir = im
                end if
                im = (il+ir)/2
                not_found = (ir - il) > 1   
              end do
              not_found = .true.
              kl = 1
              kr = Nz_plasma_tab
              km = (kl+kr)/2
              do while (not_found)
                if ((plasma_tab((im-1)*Nz_plasma_tab+km,2)) <= z_a(k)) then
                  kl = km
                else
                  kr = km
                end if
                km = (kl+kr)/2
                not_found = (kr - kl) > 1   
              end do
              if ((im.ne.Nx_plasma_tab).and.(km.ne.Nz_plasma_tab)) then
                ! linear interpolation
                alpha_x = (plasma_tab(im*Nz_plasma_tab+km,1) - x_a(i))&
                        / (plasma_tab(im*Nz_plasma_tab+km,1) - plasma_tab((im-1)*Nz_plasma_tab+km,1))
                beta_x  = (x_a(i) - plasma_tab((im-1)*Nz_plasma_tab+km,1))&
                        / (plasma_tab(im*Nz_plasma_tab+km,1) - plasma_tab((im-1)*Nz_plasma_tab+km,1))
                alpha_z = (plasma_tab((im-1)*Nz_plasma_tab+km+1,2) - z_a(k))&
                        / (plasma_tab((im-1)*Nz_plasma_tab+km+1,2) - plasma_tab((im-1)*Nz_plasma_tab+km,2))
                beta_z  = (z_a(k) - plasma_tab((im-1)*Nz_plasma_tab+km,2))&
                        / (plasma_tab((im-1)*Nz_plasma_tab+km+1,2) - plasma_tab((im-1)*Nz_plasma_tab+km,2))
                r = (alpha_x * ((alpha_z * plasma_tab((im-1)*Nz_plasma_tab+km  ,3))&
                               +(beta_z  * plasma_tab((im-1)*Nz_plasma_tab+km+1,3))))&
                  + (beta_x  * ((alpha_z * plasma_tab( im   *Nz_plasma_tab+km  ,3))&
                               +(beta_z  * plasma_tab( im   *Nz_plasma_tab+km+1,3))))
                t = (alpha_x * ((alpha_z * plasma_tab((im-1)*Nz_plasma_tab+km  ,4))&
                               +(beta_z  * plasma_tab((im-1)*Nz_plasma_tab+km+1,4))))&
                  + (beta_x  * ((alpha_z * plasma_tab( im   *Nz_plasma_tab+km  ,4))&
                               +(beta_z  * plasma_tab( im   *Nz_plasma_tab+km+1,4))))
              else if ((im == Nx_plasma_tab).and.(km.ne.Nz_plasma_tab)) then
                alpha_z = (plasma_tab((im-1)*Nz_plasma_tab+km+1,2) - z_a(k))&
                        / (plasma_tab((im-1)*Nz_plasma_tab+km+1,2) - plasma_tab((im-1)*Nz_plasma_tab+km,2))
                beta_z  = (z_a(k) - plasma_tab((im-1)*Nz_plasma_tab+km,2))&
                        / (plasma_tab((im-1)*Nz_plasma_tab+km+1,2) - plasma_tab((im-1)*Nz_plasma_tab+km,2))
                r = (alpha_z * plasma_tab((im-1)*Nz_plasma_tab+km  ,3))&
                  + (beta_z  * plasma_tab((im-1)*Nz_plasma_tab+km+1,3))
                t = (alpha_z * plasma_tab((im-1)*Nz_plasma_tab+km  ,4))&
                  + (beta_z  * plasma_tab((im-1)*Nz_plasma_tab+km+1,4))
              else if ((im.ne.Nx_plasma_tab).and.(km == Nz_plasma_tab)) then
                alpha_x = (plasma_tab(im * Nz_plasma_tab+km,1) - x_a(i))&
                        / (plasma_tab(im * Nz_plasma_tab+km,1) - plasma_tab((im-1)*Nz_plasma_tab+km,1))
                beta_x  = (x_a(i) - plasma_tab((im-1)*Nz_plasma_tab+km,1))&
                        / (plasma_tab(im*Nz_plasma_tab+km,1) - plasma_tab((im-1)*Nz_plasma_tab+km,1))
                r = (alpha_x * plasma_tab((im-1)*Nz_plasma_tab+km,3))&
                  + (beta_x  * plasma_tab( im   *Nz_plasma_tab+km,3))
                t = (alpha_x * plasma_tab((im-1)*Nz_plasma_tab+km,4))&
                  + (beta_x  * plasma_tab( im   *Nz_plasma_tab+km,4))
              else
                r = plasma_tab((im-1)*Nz_plasma_tab+km,3)
                t = plasma_tab((im-1)*Nz_plasma_tab+km,4)
              end if
              A(i,k)  = A1
              Z(i,k)  = Z1
              ni(i,k) = r / (A1 * mu)
              Te(i,k) = t * eV
              Ti(i,k) = t * eV
            end do
          end do
          !$omp END PARALLEL DO
          deallocate(plasma_tab)
          ! guard cells :
          do i = 0,N_x+1,1
            A(i,0)      = A(i,1)
            A(i,N_z+1)  = A(i,N_z)
            Z(i,0)      = Z(i,1)
            Z(i,N_z+1)  = Z(i,N_z)
            ni(i,0)     = ni(i,1)
            ni(i,N_z+1) = ni(i,N_z)
            Te(i,0)     = Te(i,1)
            Te(i,N_z+1) = Te(i,N_z)
            Ti(i,0)     = Ti(i,1)
            Ti(i,N_z+1) = Ti(i,N_z)
          end do
          do k = 1,N_z,1
            A(0,k)      = A(1,k)
            A(N_x+1,k)  = A(N_x,k)
            Z(0,k)      = Z(1,k)
            Z(N_x+1,k)  = Z(N_x,k)
            ni(0,k)     = ni(1,k)
            ni(N_x+1,k) = ni(N_x,k)
            Te(0,k)     = Te(1,k)
            Te(N_x+1,k) = Te(N_x,k)
            Ti(0,k)     = Ti(1,k)
            Ti(N_x+1,k) = Ti(N_x,k)
          end do
        else
          !$omp PARALLEL DO DEFAULT(NONE) &
          !$omp SHARED(N_z, N_x, A, Z, ni, Te, Ti, A1, Z1, ni1, T_ini) & 
          !$omp PRIVATE(i,k) COLLAPSE(2)
          do i=0,N_x+1,1
            do k=0,N_z+1,1
              A(i,k)     = A1
              Z(i,k)     = Z1
              ni(i,k)    = ni1
              Te(i,k)    = T_ini 
              Ti(i,k)    = T_ini
            end do
          end do
          !$omp END PARALLEL DO
        end if
      case (2)
        !$omp PARALLEL DO DEFAULT(NONE) &
        !$omp SHARED(N_x, N_z, A, Z, ni, Te, Ti, A1, Z1, ni1) & 
        !$omp PRIVATE(i,k) COLLAPSE(2)
        do i=0,N_x+1,1
          do k=0,N_z+1,1
            ni(i,k)    = ni1
            A(i,k)     = A1
            Z(i,k)     = Z1
            Te(i,k)    = Tamb 
            Ti(i,k)    = Tamb
          end do
        end do
        !$omp END PARALLEL DO
    end select
    if (Tracer==1) then
      !$omp PARALLEL DO DEFAULT(NONE) &
      !$omp SHARED(N_x, N_z, A, Z, ni, z_a) &
      !$omp SHARED(z_tracer_start, z_tracer_stop, A_tracer, rho_tracer, Z_tracer) & 
      !$omp PRIVATE(i,k) COLLAPSE(2)
      do i=0,N_x+1,1
        do k=1,N_z,1
          if ( (z_a(k) >= z_tracer_start) .and. (z_a(k) < z_tracer_stop) ) then
            A(i,k)     = A_tracer
            ni(i,k)    = rho_tracer / (mu * A_tracer)
            Z(i,k)     = Z_tracer
          end if
        end do
      end do
      !$omp END PARALLEL DO
      ! guard cells
      if (z_a(1) == z_tracer_start) then
        !$omp PARALLEL DO DEFAULT(NONE) &
        !$omp SHARED(N_x, A, Z, ni) & 
        !$omp PRIVATE(i) COLLAPSE(1)
        do i=0,N_x+1,1
          ni(i,0) = ni(i,1)
          A(i,0)  = A(i,1)
          Z(i,0)  = Z(i,1)
        end do
        !$omp END PARALLEL DO
      end if
      if (z_a(N_z) == z_tracer_stop) then
        !$omp PARALLEL DO DEFAULT(NONE) &
        !$omp SHARED(N_x, N_z, A, Z, ni) & 
        !$omp PRIVATE(i) COLLAPSE(1)
        do i=0,N_x+1,1
          ni(i,N_z+1) = ni(i,N_z)
          A(i,N_z+1)  = A(i,N_z)
          Z(i,N_z+1)  = Z(i,N_z)
        end do
        !$omp END PARALLEL DO
      end if
    end if
    !$omp PARALLEL DO DEFAULT(NONE) &
    !$omp SHARED(N_x, N_z, Ex, Ez, By_n, By_np1) &
    !$omp SHARED(nb, jbx_n, jbz_n, jrx_n, jrz_n) &
    !$omp SHARED(jbx_np1, jbz_np1, jrx_np1, jrz_np1) &
    !$omp SHARED(pdepos, plost_col, plost_res) &
    !$omp SHARED(nK_Holes, nKalpha, nKbeta, ioniztime) & 
    !$omp PRIVATE(i,k) COLLAPSE(2)
    do i=1,N_x,1
      do k=1,N_z,1
        Ex(i,k)          = 0._PR
        Ez(i,k)          = 0._PR
        By_n(i,k)        = 0._PR
        By_np1(i,k)      = 0._PR
        nb(i,k)          = 0._PR
        jbx_n(i,k)       = 0._PR
        jbz_n(i,k)       = 0._PR
        jrx_n(i,k)       = 0._PR
        jrz_n(i,k)       = 0._PR           
        jbx_np1(i,k)     = 0._PR
        jbz_np1(i,k)     = 0._PR
        jrx_np1(i,k)     = 0._PR
        jrz_np1(i,k)     = 0._PR          
        pdepos(i,k)      = 0._PR
        plost_col(i,k)   = 0._PR
        plost_res(i,k)   = 0._PR
        nK_Holes(i,k)    = 0._PR
        nKalpha(i,k)     = 0._PR
        nKbeta(i,k)      = 0._PR
        ioniztime(i,k)   = 0._PR
      end do
    end do
    ud    = 0._PR
    udcol = 0._PR
    udres = 0._PR
    ue    = 0._PR
    usf   = 0._PR
    usb   = 0._PR
    usu   = 0._PR
    usd   = 0._PR 
    ub    = 0._PR
    uel   = 0._PR
    uma   = 0._PR
    do l=1,N_eps,1
      do k=1,N_z,1
        do i=1,N_x,1
          Sva(i,k,l)     = S_tot(A(i,k), Z(i,k), ni(i,k), Te(i,k), Ti(i,k), eps_a(l)) &
          &* vit(eps_a(l))
          Sva_lp1(i,k,l) = S_tot(A(i,k), Z(i,k), ni(i,k), Te(i,k), Ti(i,k), eps_a(l)+d_eps) &
          &* vit(eps_a(l)+d_eps)
          nua(i,k,l)     = nu_tot(Z(i,k), ni(i,k), Te(i,k), Ti(i,k), eps_a(l)) 
        end do
      end do
    end do
  end subroutine initialize

  subroutine initialize_next_step(phin, phinp1,&
                                  By_n, By_np1, Ex, Ez,&
                                  jbx_n, jbz_n, jrx_n, jrz_n,&
                                  jbx_np1, jbz_np1, jrx_np1, jrz_np1,&
                                  pdepos, plost_col, plost_res, ioniztime,&
                                  Ub, Uel, Uma)
    implicit none
    real(PR), dimension(forward:backward,psi0:psi1z, -1:N_x+2, -1:N_z+2, -1:N_eps+2), intent(inout) :: phin
    real(PR), dimension(forward:backward,psi0:psi1z, -1:N_x+2, -1:N_z+2, -1:N_eps+2), intent(inout) :: phinp1
    real(PR), intent(inout), dimension(1:N_x,1:N_z)                                :: Ex, Ez
    real(PR), intent(inout), dimension(1:N_x,1:N_z)                                :: By_n, By_np1
    real(PR), intent(inout), dimension(1:N_x,1:N_z)                                :: jrx_n, jrz_n
    real(PR), intent(inout), dimension(1:N_x,1:N_z)                                :: jbx_n, jbz_n
    real(PR), intent(inout), dimension(1:N_x,1:N_z)                                :: jrx_np1, jrz_np1
    real(PR), intent(inout), dimension(1:N_x,1:N_z)                                :: jbx_np1, jbz_np1
    real(PR), intent(inout), dimension(1:N_x,1:N_z)                                :: pdepos, plost_col
    real(PR), intent(inout), dimension(1:N_x,1:N_z)                                :: plost_res
    real(PR), intent(inout), dimension(1:N_x,1:N_z)                                :: ioniztime
    real(PR), intent(inout)                                                        :: ub, uel, uma
    integer                                                                        :: m,i, k, l
    !$omp  PARALLEL DO DEFAULT(NONE) &
    !$omp& SHARED(backward,N_eps, N_x, N_z, phin, phinp1) &
    !$omp& PRIVATE(m,l,k,i) COLLAPSE(4)
    do m=forward,backward,1
      do l=-1,N_eps+2,1
        do k=-1,N_z+2,1
          do i=-1,N_x+2,1
            phin(m,psi0:psi1z,i,k,l)    = phinp1(m,psi0:psi1z,i,k,l)
            phinp1(m,psi0:psi1z,i,k,l)  = 0._PR
          end do
        end do
      end do
    end do
    !$omp PARALLEL DO DEFAULT(NONE) &
    !$omp SHARED(N_x, N_z, Ex, Ez, By_n, By_np1) &
    !$omp SHARED(jbx_n, jbz_n, jrx_n, jrz_n) &
    !$omp SHARED(jbx_np1, jbz_np1, jrx_np1, jrz_np1) &
    !$omp SHARED(pdepos, plost_col, plost_res) &
    !$omp SHARED(ioniztime) & 
    !$omp PRIVATE(i,k) COLLAPSE(2)
    do k=1,N_z,1
      do i=1,N_x,1
        By_n(i,k)      = By_np1(i,k)             
        By_np1(i,k)    = 0._PR
        Ex(i,k)        = 0._PR
        Ez(i,k)        = 0._PR
        jbx_n(i,k)     = jbx_np1(i,k)
        jbz_n(i,k)     = jbz_np1(i,k)
        jbx_np1(i,k)   = 0._PR 
        jbz_np1(i,k)   = 0._PR
        jrx_n(i,k)     = jrx_np1(i,k)
        jrz_n(i,k)     = jrz_np1(i,k)
        jrx_np1(i,k)   = 0._PR 
        jrz_np1(i,k)   = 0._PR
        pdepos(i,k)    = 0._PR
        plost_col(i,k) = 0._PR
        plost_res(i,k) = 0._PR
        ioniztime(i,k) = 0._PR
      end do
    end do
    Ub  = 0._PR
    Uel = 0._PR
    Uma = 0._PR
  end subroutine initialize_next_step

end module initialization