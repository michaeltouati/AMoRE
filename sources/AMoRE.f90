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
program amore

  use acuracy
  use constants
  use input
  use vfp
  use diagnostics
  use initialization
  use em_fields
  use heat_equations

  ! Declaration of variables :
  implicit none
  ! N_t  = iteration number
  ! d_t  = time step computed according to a CFL in (fs); time = time in (fs)
  ! norm = coefficient of normalization for the angular moment components in (/cm^3/keV)
  integer  :: N_t 
  real(PR) :: d_t, time, norm
  ! z/x_a = array containing the spatial grid on z/x-axis in (microns) 
  ! eps_a = array containing the kinetic energy grid in (keV)
  real(PR), dimension(:), allocatable :: z_a, x_a, eps_a
  ! Sv_a/Sv_a_lp1  = table containing the total stopping power
  !                  of fast e- with kinetic energy eps(l)/eps(l+1) multiplied by their 
  !                  corresponding velocity in (erg/s)
  ! nua            = table containing the total angular collision rate of fast e-
  !                  with kinetic energy eps(l) in (/s)
  real(PR), dimension(:,:,:), allocatable :: Sv_a, Sv_a_lp1, nua
  ! phi_n/np1 = table containing the two first angular moment components of 
  !             the distribution function (Psi_0 Psi_1x Psi_1z) at iteration N_t/N_t+1 in ()
  real(PR), dimension(:,:,:,:,:), allocatable :: phi_n, phi_np1
  ! dF_dx/z/eps = table containing the derivatives along x/z/eps of angular moment equations
  real(PR), dimension(:,:,:,:,:), allocatable :: dF_dx, dF_dz, dF_deps
  ! gama_e/b/nu = table containing the source term of angular moment equations due to the
  !               electric field/magnetic field/angular collisions
  real(PR), dimension(:,:,:,:,:), allocatable :: gama_e, gama_b, gama_nu
  ! nb           = table containing the fast e- beam density in (/cm^3) computed according
  !                 to the 0th order angular moment of the distribution function
  ! jbx/z_n/np1 = table containing the fast e- beam current density component on x/z-axis
  !                 at iteration N_t/N_t +1 in (statAmpere/cm^2) computed according to 
  !                 the 1st order angular moment of the distribution function
  ! jrx/z_n/np1 = table containing the return current density component on x/z-axis
  !                 at iteration N_t/N_t +1 in (statAmpere/cm^2) computed according to 
  !                 the Ampere's equation
  real(PR), dimension(:,:), allocatable :: nb,  jbx_n, jbz_n, jbx_np1, jbz_np1
  real(PR), dimension(:,:), allocatable :: jrx_n, jrz_n, jrx_np1, jrz_np1
  ! Ex/z = table containing the self-generated electric field component on x/z-axis 
  !        at iteration N_t in (statVolt/cm) computed according to the Ohm's law
  ! By_n/By_np1 = table containing the self-generated magnetic field component on y-axis 
  !               at iteration N_t/N_t+1 in (Gauss) computed according to the Faraday,
  !               Ohm and Ampere's equations without displacement current
  real(PR), dimension(:,:), allocatable :: Ex, Ez, By_n, By_np1
  ! Kalpha_tab    = table containing data concerning Kalpha photon emission of radiations
  !                 (see file 'Kalpha_tab.dat')
  ! nK_holes     = table containing the K-shell hole density in (/cm^3) induced by impact
  !                 ionization with fast e-
  ! nKalpha/beta = table containing time integrated density of Kalpha/beta photons in 
  !                 (/cm^3/sr) emitted due to K-shell e- vacancies
  ! ioniz_rate    = table containing K-shell e- ionization rate in (/s) due to collisions 
  !                 with fast e- 
  real(PR), dimension(1:79,1:7)         :: Kalpha_tab
  real(PR), dimension(:,:), allocatable :: nK_holes, nKalpha, nKbeta, ioniz_rate
  ! A  = table containing the atomic weight of the material in ()
  ! Z  = table containing the atomic number of the material in ()
  ! ni = table containing the ion density in the material in (/cm^3)
  real(PR), dimension(:,:), allocatable :: A, Z, ni 
  ! p_lost_col = density of power lost by the fast e- beam in (erg/cm^3/s) due to collisions
  ! p_lost_res = density of power lost by the fast e- beam in (erg/cm^3/s) due to the 
  !              resistive electric field
  ! p_depos    = density of power lost by the fast e- beam and deposited on the background
  !              e- in the material in (erg/s/cm^3) 
  real(PR), dimension(:,:), allocatable :: p_depos, p_lost_col, p_lost_res, T_e, T_i
  ! Ue    = time integrated fast e- beam kinetic energy injected in the simulation box 
  !         at iteration N_t in (J)
  ! Ub    = instantaneous beam energy in the simulation box at iteration N_t in (J)
  ! Udcol = time integrated energy lost by fast e- due to collisions at iteration N_t in (J) 
  ! Udres = time integrated energy lost by fast e- due to electric fields at iteration N_t 
  !         in (J)
  ! Uel   = instantaneous electric energy in the simulation box at iteration N_t in (J)
  ! Uma   = instantaneous magnetic energy in the simulation box at iteration N_t in (J)
  ! Usf   = time integrated fast e- beam energy escaping from the simulation box at z = L_z
  !         and iteration N_t in (J)
  ! Usb   = time integrated fast e- beam energy escaping from the simulation box at z = 0
  !         and iteration N_t in (J)
  ! Usu   = time integrated fast e- beam energy escaping from the simulation box at 
  !         x = +L_x/2 and iteration N_t in (J)
  ! Usd   = time integrated fast e- beam energy escaping from the simulation box at 
  !         x = -L_x/2 and iteration N_t in (J)
  real(PR) :: Ue, Ub, Ud, Ud_col, Ud_res
  real(PR) :: Uel, Uma, Usf, Usb, Usu, Usd 
  ! eta_tab = table containing the electrical resistivity if the user want to use the 
  !           tabulated electrical resistivity 'resistivity.dat' 
  real(PR), dimension(:,:), allocatable :: eta_tab, eps_tab
  ! diag_condition = logical defined to ensure to write computation results in .txt 
  !                  files every Delta_t_diag as defined by the user
  logical :: diag_condition

  ! Read the simulation parameter
  call read_init_parameters  
  print *, 'Initialization of the simulation'
  write(*,*)' '
  ! Initialization of diagnostics
  call initialize_spectrum(eps_tab, N_eps_tab)
  call initialize_resistivity(eta_tab, N_eta_tab)
  call plasma_diagnostics(eta_tab)
  call initialize_diagnostics(Kalpha_tab)

  ! Initialization
  N_t  = 1
  time = 0._PR

  allocate(z_a(1:N_z), x_a(1:N_x), eps_a(1:N_eps),&
           Sv_a(1:N_x,1:N_z,1:N_eps), Sv_a_lp1(1:N_x,1:N_z,1:N_eps), nUa(1:N_x,1:N_z,1:N_eps),&
           A(0:N_x+1,0:N_z+1), Z(0:N_x+1,0:N_z+1), ni(0:N_x+1,0:N_z+1),&
           Ex(1:N_x,1:N_z), Ez(1:N_x,1:N_z),By_n(1:N_x,1:N_z), By_np1(1:N_x,1:N_z),&
           nb(1:N_x,1:N_z), jbx_n(1:N_x,1:N_z), jbz_n(1:N_x,1:N_z),&
           jbx_np1(1:N_x,1:N_z), jbz_np1(1:N_x,1:N_z),&
           jrx_n(1:N_x,1:N_z), jrz_n(1:N_x,1:N_z),jrx_np1(1:N_x,1:N_z), jrz_np1(1:N_x,1:N_z),&
           p_depos(1:N_x,1:N_z), p_lost_col(1:N_x,1:N_z), p_lost_res(1:N_x,1:N_z),&
           T_e(0:N_x+1,0:N_z+1), T_i(0:N_x+1,0:N_z+1),&
           nK_holes(1:N_x,1:N_z), nKalpha(1:N_x,1:N_z), nKbeta(1:N_x,1:N_z),&
           ioniz_rate(1:N_x,1:N_z),&
           phi_n(forward:backward  ,psi0:psi1z, -1:N_x+2, -1:N_z+2, -1:N_eps+2),&
           phi_np1(forward:backward,psi0:psi1z, -1:N_x+2, -1:N_z+2, -1:N_eps+2))

  ! Definition of space and kinetic energy grids :
  call grid(z_a,x_a,eps_a)
  ! Diagnostics of the initialized electron beam
  call electron_diagnostics(eps_tab,x_a,eps_a)
  ! Initialize arrays
  call initialize(x_a, z_a, eps_a, ni, A, Z, Ex, Ez, By_n, By_np1,&
                  nb, jbx_n, jbz_n, jrx_n, jrz_n,&
                  jbx_np1, jbz_np1, jrx_np1, jrz_np1,&
                  p_depos, p_lost_col, p_lost_res, T_e, T_i,&
                  nK_holes, nKalpha, nKbeta, ioniz_rate,&
                  Ud, Ud_col, Ud_res, Ue, Usf, Usb, Usu, Usd, Ub, Uel, Uma,&
                  Sv_a, Sv_a_lp1, nUa)

  ! Time loop :
  do while (time <= L_t) 

  ! Update of self-generated electromagnetic fields :
  call magnetic_field(eta_tab, Z, ni, T_e, T_i,&
                      jbx_n, jbz_n, d_t, By_n, By_np1)       
  call electric_field(eta_tab, Z, ni, T_e, T_i,&
                      jrx_n, jrz_n, Ex, Ez)

  ! Computation of the CFL :
  call cfl_scheme(eps_a, Sv_a, nUa, Ex, Ez, By_n, d_t)  

  ! Injection of the fast e- beam :
  ! (initialize phi_n and phi_np1 at z = 0 according to the simulation parameters)
  call initialize_phi(eps_tab, time, d_t, x_a, eps_a,&
                      phi_n, phi_np1, norm)

  ! Update phi_n depending on the boundary conditions :
  call boundary_cond(norm, time, d_t, eps_a, Ue, Usf, Usb, Usu, Usd, phi_n)

  ! Update hydrodynamic quantities   
  call hydrodynamic_quantities(d_t, norm, Sv_a, phi_n, &
                               A, Z, ni, T_e, T_i, &
                               Ex, Ez, jrx_n, jrz_n, &
                               nb, jbx_n, jbz_n, &
                               p_depos, p_lost_col, p_lost_res)

  ! Allocation of temporary arrays for computing the angular moment equations :
  allocate(dF_dx(forward:backward,psi0:psi1z, 1:N_x, 1:N_z, 1:N_eps),&
  dF_dz(forward:backward,psi0:psi1z, 1:N_x, 1:N_z, 1:N_eps),&
  dF_deps(forward:backward,psi0:psi1z, 1:N_x, 1:N_z, 1:N_eps),&
  gama_e(forward:backward,psi0:psi1z, 1:N_x, 1:N_z, 1:N_eps),&
  gama_b(forward:backward,psi0:psi1z, 1:N_x, 1:N_z, 1:N_eps),&
  gama_nu(forward:backward,psi0:psi1z, 1:N_x, 1:N_z, 1:N_eps))

  ! Compute the spatial derivatives in the z and x directions : 
  call spatial_derivatives(eps_a,phi_n,dF_dz,dF_dx)    

  ! Compute the kinetic energy derivative :
  call energy_derivative(Ex,Ez,eps_a,Sv_a,Sv_a_lp1,phi_n,dF_deps)

  ! Compute the source terms :
  call source_terms(Ex,Ez,By_n,eps_a,nUa,phi_n,gama_E,gama_B,gama_nu)

  ! Update the new angular moments at iteration N_t+1 :
  call update(phi_n, d_t, dF_dx, dF_dz, dF_deps, gama_e, gama_b, gama_nu,&
  Sv_a, Sv_a_lp1, nUa, phi_np1)

  ! Deallocation of temporary arrays for computing the angular moment equations :
  deallocate(dF_dx,dF_dz,dF_deps,gama_e,gama_b,gama_nu)

  ! Update the fast e- beam and return current densities :
  call currents(norm,eps_a,phi_np1,By_np1,jbx_np1,jbz_np1,jrx_np1,jrz_np1)

  ! Diagnostics :
  if (Kalpha) then
    ! -Update the time-integrated densities Kalpha and Kbeta photons emitted
    call Kalpha_emission(Kalpha_tab, norm, d_t, Z, ni, eps_a, phi_n,&
    nK_holes, nKalpha, nKbeta, ioniz_rate) 
  end if
  diag_condition = abs(mod(time,Delta_t_diag)) < d_t
  diag_condition = diag_condition.or.(time.eq.0._PR) 
  diag_condition = diag_condition.or.((time <= L_t).and.((time+d_t) > L_t))
  if (diag_condition) then
    ! -write the results in text files
    call write_results(eta_tab, A, Z, ni, nb, jbx_n, jbz_n, jrx_n, jrz_n,&
                       Ex, Ez, By_n, &
                       p_depos, T_e, T_i,&
                       nKalpha, nKbeta, ioniz_rate,&
                       phi_n, norm, time, x_a, z_a, eps_a)
  end if
  ! -Energy balance 
  !  (write the results in text files and print the global balance in the console)
  call energy_balance(diag_condition, N_t, time, d_t, norm, eps_a, phi_n,&
                      p_depos, p_lost_col, p_lost_res, Ex, Ez, By_n,&
                      Ue, Ub, Ud, Ud_col, Ud_res, Uel, Uma,&
                      Usf, Usb, Usu, Usd)

  ! Initialization for the next iteration :
  N_t = N_t + 1
  time = time + d_t 
  call initialize_next_step(phi_n, phi_np1,&
                            By_n, By_np1,&
                            Ex, Ez,&
                            jbx_n, jbz_n,jrx_n, jrz_n,&
                            jbx_np1, jbz_np1,jrx_np1, jrz_np1,&
                            p_depos, p_lost_col, p_lost_res, ioniz_rate,&
                            Ub, Uel, Uma)     

  ! End of time loop :
  end do  

  deallocate(z_a,x_a,eps_a,Sv_a,Sv_a_lp1,nua,eta_tab,&
             Ex, Ez, nb, jbx_n, jbz_n, jbx_np1, jbz_np1, By_n, By_np1,&
             jrx_n, jrz_n, jrx_np1, jrz_np1,p_depos, p_lost_col, p_lost_res, T_e, T_i,&
             nK_holes, nKalpha, nKbeta, ioniz_rate, phi_n, phi_np1)
  !
end program amore