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
module input

use acuracy
use constants
use physics_library
use omp_lib

implicit none
public :: read_init_parameters

!=====================================================================================
!                      Public variables known in the whole code 
!                 (see the file init_parameters for deenditions)
!=====================================================================================

character(len=60), public :: simu
real(PR)         , public :: sigma_t, L_t
real(PR)         , public :: sigma_r
real(PR)         , public :: d_x, L_x
real(PR)         , public :: d_z, L_z
integer          , public :: spectrum
real(PR)         , public :: E_tot, eps0, eps1, sigma_eps, Tb0, alpha0, alpha1, Tb1
real(PR)         , public :: eps_min, d_eps, L_eps
real(PR)         , public :: sigma_theta, angle0
integer          , public :: Material 
Character(len=2) , public :: Mat
real(PR)         , public :: Z0, A0, rho, T_ini
integer          , public :: Tracer
real(PR)         , public :: z_tracer_start, z_tracer_stop, Z_tracer, A_tracer, rho_tracer
real(PR)         , public :: Delta_t_diag
real(PR)         , public :: cfl
integer          , public :: hll_order, N_th
logical          , public :: coll_implicit_scheme, thermal_conduction
logical          , public :: bitemperature, magnetic_diffusion
logical          , public :: Kalpha
logical          , public :: tabulated_resistivity
logical          , public :: tabulated_plasma
real(PR)         , public :: sim_box_thickness, t_c, source_duration
integer          , public :: N_x, N_z, N_eps, N_threads, N_eta_tab, N_eps_tab
integer          , public :: backward
logical          , public :: irradiated_side_refluxing

contains

subroutine read_init_parameters
  implicit none
  character(len=80)         :: str
  integer                   :: i, istr
  integer                   :: N_th_max
  real(PR)                  :: ni0
  logical                   :: Tmin_alert, Tmin_alert2
  ! read the file 'init_parameters'
  open (unit=1, file='input-deck',&
  & position='rewind', access='sequential',&
  & form='formatted', action='read',status='old')
  read: do
    call get_str(str)
    if (str == 'end') exit read
    if (str(1:1) /= '#') cycle read
    istr = 80
    do i = 1, len(str)
      if (str(i:i) == ' ') then
        istr = i
        exit
      end if
    end do
    select case (str(1:istr))
      case ('#simu')
        simu=get_char(str(istr+1:))
      case ('#hll_order')
        hll_order = get_integer(str(istr+1:))
      case ('#implicit_coll')  
        coll_implicit_scheme = get_logical(str(istr+1:))
      case ('#cfl')
        cfl = get_real(str(istr+1:))  
      case ('#bi_temp')  
        bitemperature = get_logical(str(istr+1:))
      case ('#magnetic_diff')
        magnetic_diffusion = get_logical(str(istr+1:)) 
        thermal_conduction = get_logical(str(istr+1:))
      case ('#N_threads')
        N_th = get_integer(str(istr+1:))
      case ('#backward')
        backward = get_integer(str(istr+1:))
      case ('#irradiated_side_refluxing')
        irradiated_side_refluxing = get_logical(str(istr+1:))
      case ('#L_t')
        L_t = get_real(str(istr+1:))
      case ('#d_z')
        d_z = get_real(str(istr+1:))
      case ('#d_x')
        d_x = get_real(str(istr+1:))
      case ('#L_z')
        L_z = get_real(str(istr+1:))
      case ('#L_x')
        L_x = get_real(str(istr+1:))
      case ('#d_eps')
        d_eps = get_real(str(istr+1:))
      case ('#eps_min')
        eps_min = get_real(str(istr+1:))
      case ('#L_eps')
        L_eps = get_real(str(istr+1:))
      case ('#E_tot')
        E_tot = get_real(str(istr+1:))
      case ('#Delta_t')
        sigma_t = get_real(str(istr+1:))
        sigma_t = sigma_t / sqrt(8._PR*log(2._PR))
      case ('#Delta_r')
        sigma_r = get_real(str(istr+1:))
        sigma_r = sigma_r / sqrt(8._PR*log(2._PR))
      case ('#spectrum')
        spectrum = get_integer(str(istr+1:))
      case ('#eps0')
        eps0 = get_real(str(istr+1:))
      case ('#eps1')
        eps1 = get_real(str(istr+1:))
      case ('#sigma_eps')
        sigma_eps = get_real(str(istr+1:))
      case ('#Tb0')
        Tb0 = get_real(str(istr+1:))
      case ('#alpha0')
        alpha0 = get_real(str(istr+1:))
      case ('#alpha1')
        alpha1 = get_real(str(istr+1:))
      case ('#Tb1')
        Tb1 = get_real(str(istr+1:))
      case ('#angle0')
        angle0 = get_real(str(istr+1:))
      case ('#Delta_theta')
        sigma_theta = get_real(str(istr+1:))
        sigma_theta = sigma_theta / sqrt(8._PR*log(2._PR))
      case('#tabulated_resistivity')
        tabulated_resistivity = get_logical(str(istr+1:))
      case ('#Material')
        Material = get_integer(str(istr+1:))
      case ('#Mat')
        Mat = trim(get_char(str(istr+1:)))
      case ('#tabulated_plasma')
        tabulated_plasma = get_logical(str(istr+1:)) 
      case ('#Z0')
        Z0 = get_real(str(istr+1:))
      case ('#A0')
        A0 = get_real(str(istr+1:))
      case ('#rho')
        rho = get_real(str(istr+1:))
      case ('#T_ini')
        T_ini = get_real(str(istr+1:)) 
        T_ini = T_ini*eV
      case ('#Tracer')
        Tracer = get_integer(str(istr+1:))
      case ('#z_tracer_start')
        z_tracer_start = get_real(str(istr+1:))
      case ('#z_tracer_stop')
        z_tracer_stop  = get_real(str(istr+1:))
      case ('#Z_tracer')
        Z_tracer       = get_real(str(istr+1:))
      case ('#A_tracer')
        A_tracer       = get_real(str(istr+1:))
      case ('#rho_tracer')
        rho_tracer     = get_real(str(istr+1:))
      case ('#Kalpha')
        Kalpha = get_logical(str(istr+1:)) 
      case ('#Delta_t_diag')
        Delta_t_diag = get_real(str(istr+1:))
    end select
  end do read
  close(1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Deduced simu. parameters !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Time in (fs) at which the fast e- number injected in the simulation box is maximum :
  t_c = 3.5_PR*sigma_t

  ! Time duration of fast e- in the simulation box from t = 0 to source_duration in (fs) :
  source_duration = 7._PR*sigma_t

  ! Simulation box thickness in (cm) :
  sim_box_thickness = sigma_r*microns*sqrt(2._PR*pi)

  ! Number of grid cells on x-axis :
  N_x   = floor(L_x/d_x)

  ! Number of grid cells on z-axis :
  N_z   = floor(L_z/d_z)

  ! Number of grid cell on eps-axis (fast e- Kinetic energy) :
  N_eps = floor(L_eps/d_eps)
  
  ! Maximum number of available CPUs on the machine :
  N_th_max = omp_get_max_threads()
  N_threads = N_th_max
  if (N_th.gt.0) N_threads = N_th
  call omp_set_num_threads(N_threads)

  ! Target ion density
  ni0 = rho / (mu * A0)

  ! Alert initial temperature
  Tmin_alert = (Material == 1).and.(T_ini <= Tfermi(Z0, ni0)).and.(Z0.ne.1._PR)
  Tmin_alert = Tmin_alert.and.((Z0.ne.13._PR).and.(Z0.ne.29._PR).and.(Z0.ne.73._PR))
  Tmin_alert = Tmin_alert.and.(tabulated_resistivity.eqv..false.)
  Tmin_alert2 = (((Tracer == 1).and.(T_ini <= Tfermi(Z_tracer, rho_tracer/(A_tracer*mu))).and.(Z_tracer.ne.1)))
  Tmin_alert2 = Tmin_alert2.and.((Z_tracer.ne.13._PR).and.(Z_tracer.ne.29._PR).and.(Z_tracer.ne.73._PR))
  if (Tmin_alert.or.Tmin_alert2) then
    write(*,*)'The plasma temperature T_ini must be greater than the Fermi temperature :'
    if (Tmin_alert) write(*,*)'because the electrical resistivity is not yet implemented for this plasma'
    if (Tmin_alert2) write(*,*)'because the electrical resistivity is not yet implemented for this tracer'
    write(*,*) 'T_ini =',T_ini,'eV while T_F = ',Tfermi(Z0, ni0)/eV,' eV'
    stop
  end if

  ! Alert bitemperature MHD equations
  if (((Material == 1).or.tabulated_resistivity).and.(bitemperature.eqv..true.)) then 
    write(*,*)'For plasmas and/or a tabulated resistivity, Te different from Ti is not implemented'
    write(*,*)'bi_temp must be set .false.'
    stop
  end if

  ! Alert magnetic field diffusion
  ! if (tabulated_resistivity.and.magnetic_diffusion) then
  !   write(*,*)'For a tabulated resistivity, the magnetic field diffusion is not yet implemented'
  !   write(*,*)'magnetic_diff must be set .false.'
  !   stop
  ! end if

  ! Alerts refluxing 
  if ((backward .ne. 1) .and. (backward .ne. 2)) then
    write(*,*)'Target rear side fast electron refluxing parameter error'
    write(*,*)'backward should be 1 (no specular reflection) or 2 (specular reflection)!'
    stop
  end if
  if ((irradiated_side_refluxing .eqv. .true.) .and. (backward .eq. 1)) then
    irradiated_side_refluxing = .false.
    write(*,*)'Target irradiated side fast electron refluxing parameter error'
    write(*,*)'It can be taken into account only if there is electron refluxing'
    write(*,*)'at the target rear side .i.e. only if backward = 2'
    stop
  end if

  write(*,*)'-----------------------------------------'
  write(*,*)'Recapitulation of simulation parameters :'
  write(*,*)'-----------------------------------------'
  write(*,*)'* Simulation :'
  write(*,*)'  simu = ',trim(simu)
  write(*,*)'-----------------------------------------'
  write(*,*)'* Number of OpenMP threads :'
  write(*,'(A,1I4)')'   N_threads = ',N_threads
  write(*,*)'-----------------------------------------'
  write(*,*)'* Simulation properties : '
  write(*,'(A,1I2)')'   hll_order                 = ',hll_order
  write(*,'(A,1L2)')'   coll_implicit_scheme      = ',coll_implicit_scheme
  write(*,'(A,1L2)')'   bitemperature             = ',bitemperature
  write(*,'(A,1L2)')'   magnetic_diffusion        = ',magnetic_diffusion
  write(*,'(A,1L2)')'   thermal_conduction        = ',thermal_conduction
  write(*,'(A,1I2)')'   backward                  = ',backward
  write(*,'(A,1L2)')'   irradiated_side_refluxing = ',irradiated_side_refluxing
  write(*,'(A,1L2)')'   Kalpha                    = ',Kalpha
  write(*,*)'-----------------------------------------'
  write(*,*)'* Simulation time properties : '
  write(*,'(A,1E21.14)')'   cfl          = ',cfl
  write(*,'(A,1E21.14)')'   Delta_t_diag = ',Delta_t_diag
  write(*,'(A,1E21.14)')'   L_t          = ',L_t
  write(*,*)'-----------------------------------------'
  write(*,*)'* Relativistic electron beam properties :'
  write(*,'(A,1E21.14)')'   E_tot       = ',E_tot
  write(*,'(A,1E21.14)')'   sigma_t     = ',sigma_t
  write(*,'(A,1E21.14)')'   sigma_r     = ',sigma_r
  write(*,'(A,1E21.14)')'   angle0      = ',angle0
  write(*,'(A,1E21.14)')'   sigma_theta = ',sigma_theta
  write(*,'(A,1I2)')'   spectrum    = ',spectrum
  select case (spectrum)
    case (1)
      write(*,'(A,1E21.14)')'   eps0        = ',eps0
      write(*,'(A,1E21.14)')'   sigma_eps   = ',sigma_eps
    case (2)
      write(*,'(A,1E21.14)')'   Tb0         = ',Tb0
    case (3)
      write(*,'(A,1E21.14)')'   Tb0         = ',Tb0
      write(*,'(A,1E21.14)')'   alpha1      = ',alpha1
      write(*,'(A,1E21.14)')'   Tb1         = ',Tb1
    case (4)
      write(*,'(A,1E21.14)')'   alpha0      = ',alpha0
      write(*,'(A,1E21.14)')'   eps0        = ',eps0
      write(*,'(A,1E21.14)')'   Tb0         = ',Tb0
      write(*,'(A,1E21.14)')'   alpha1      = ',alpha1
      write(*,'(A,1E21.14)')'   eps1        = ',eps1
      write(*,'(A,1E21.14)')'   Tb1         = ',Tb1
    case (5)
      write(*,*)'   tabulated'
  end select
  write(*,*)'-----------------------------------------'
  write(*,*)'* Simulation box properties : '
  write(*,'(A,1E21.14)')'   d_z     = ',d_z
  write(*,'(A,1E21.14)')'   L_z     = ',L_z
  write(*,'(A,1E21.14)')'   d_x     = ',d_x
  write(*,'(A,1E21.14)')'   L_x     = ',L_x
  write(*,'(A,1E21.14)')'   eps_min = ',eps_min
  write(*,'(A,1E21.14)')'   d_eps   = ',d_eps
  write(*,'(A,1E21.14)')'   L_eps   = ',L_eps
  write(*,*)'-----------------------------------------'
  write(*,*)'* Target properties : '
  write(*,'(A,1I2)')'   Material              = ',Material
  select case (Material)
    case (1)
      write(*,'(A,1E21.14)')'   Z0                    = ',Z0
      write(*,'(A,1E21.14)')'   A0                    = ',A0
      write(*,'(A,1E21.14)')'   rho                   = ',rho
      write(*,'(A,1E21.14)')'   T_ini                 = ',T_ini
      write(*,'(A,1L2)')'   tabulated_plasma      = ',tabulated_plasma
      if (tabulated_plasma) then
        write(*,*)'   material properties diagnostics have been done with :'
        write(*,*)'   rho = ', rho, 'g/cm^3'
      end if
    case (2)
      write(*,'(A,A)')'   Mat                   = ',Mat
  end select
  write(*,'(A,1I2)')'   Tracer                = ',Tracer
  if (Tracer == 1) then
    write(*,'(A,1E21.14)')'   z_tracer_start        = ',z_tracer_start
    write(*,'(A,1E21.14)')'   z_tracer_stop         = ',z_tracer_stop
    write(*,'(A,1E21.14)')'   Z_tracer              = ',Z_tracer
    write(*,'(A,1E21.14)')'   A_tracer              = ',A_tracer
    write(*,'(A,1E21.14)')'   rho_tracer            = ',rho_tracer
  end if
  write(*,'(A,1L2)')'   tabulated_resistivity = ',tabulated_resistivity
  write(*,*)'-----------------------------------------'
  write(*,*)'* Deduced simulation parameters :'
  write(*,'(A,1I6)')'   N_z   = ',N_z
  write(*,'(A,1I6)')'   N_x   = ',N_x
  write(*,'(A,1I6)')'   N_eps = ',N_eps
  write(*,*)'-----------------------------------------'
  write(*,*)' '
end subroutine read_init_parameters

subroutine get_str(str)
  character(len=*), intent(inout) :: str
  read (1,'(A)',end=10) str
  str = adjustl(str)
  return
  continue ; 10 str = 'end'
end subroutine get_str

function get_logical(str)
  character(len=*), intent(in) :: str
  logical :: get_logical
  read (str, *) get_logical
  return
end function get_logical

function get_integer(str)
  character(len=*), intent(in) :: str
  integer :: get_integer
  read (str, *) get_integer
  return
end function get_integer

function get_real(str)
  character(len=*), intent(in) :: str
  real(PR) :: get_real
  read (str, *) get_real
  return
end function get_real

function get_char(str)
  character(len=*), intent(in) :: str
  Character(len=60)            :: get_char
  read (str, *) get_char
  get_char = adjustl(get_char)
  return
end function get_char

end module input