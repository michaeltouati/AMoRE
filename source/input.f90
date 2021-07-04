module input

use constants
use physics_library
use omp_lib

implicit none
public :: read_init_parameters

!=====================================================================================
!                      Public variables known in the whole code 
!                 (see the file init_parameters for deenditions)
!=====================================================================================

integer, public  :: shape
real(PR), public :: sigma_t, L_t
real(PR), public :: sigma_r
real(PR), public :: d_x, L_x
real(PR), public :: d_z, L_z
integer, public  :: spectrum
real(PR),public  :: E_tot, eps0, eps1, sigma_eps, Tb0, alpha0, alpha1, Tb1
real(PR),public  :: eps_min, d_eps, L_eps
real(PR), public :: sigma_theta, angle0
integer, public  :: Material 
Character(len=2) :: Mat
real(PR), public :: Z0, A0, rho, T_ini
integer, public  :: Tracer
real(PR), public :: z_tracer_start, z_tracer_stop, Z_tracer, A_tracer, rho_tracer
real(PR), public :: Delta_t_diag
real(PR), public :: cfl
integer, public  :: hll_order, N_th
logical, public  :: coll_implicit_scheme, thermal_conduction
logical, public  :: bitemperature, magnetic_diffusion
logical, public  :: Kalpha
logical, public  :: tabulated_resistivity
logical, public  :: tabulated_plasma
real(PR), public :: sim_box_thickness, t_c, source_duration
integer, public  :: N_x, N_z, N_eps, N_threads, N_eta_tab, N_eps_tab
integer, public  :: backward
logical, public  :: irradiated_side_refluxing

contains

subroutine read_init_parameters
implicit none
Character(len=80)         :: str
Integer                   :: i, istr
integer                   :: N_th_max
real(PR)                  :: ni0
logical                   :: Tmin_alert
! read the file 'init_parameters'
open (unit=1, file='init_parameters',&
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
case ('#Shape')
shape = get_integer(str(istr+1:))
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
Mat = get_char(str(istr+1:))
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
ni0 = rho / (mu * A0)
Tmin_alert = (Material == 1).and.(T_ini <= Tfermi(Z0, ni0)).and.(Z0.ne.1)
Tmin_alert = Tmin_alert.and.((Z0.ne.13).and.(Z0.ne.29).and.(Z0.ne.73))
Tmin_alert = Tmin_alert.and.(tabulated_resistivity.eqv..false.)
if (Tmin_alert) then
print*,'The plasma temperature T_ini must be greater than the Fermi temperature :'
print*, 'T_F = ',Tfermi(Z0, ni0)/eV,'eV'
stop
end if
if ((Material == 1).or.(tabulated_resistivity)) then
bitemperature = .false.   
print*,'bi_temp has been set .false. (=> Te = Ti)'
end if
if (tabulated_resistivity) then
magnetic_diffusion = .false.   
print*,'magnetic_diff has been set .false.'
end if
if (tabulated_plasma) then
print*,'material properties diagnostics have been done with :'
print*,'rho = ', rho, 'g/cm^3'
end if
if ((backward .ne. 1) .and. (backward .ne. 2)) then
backward = 1
print*,'fast electron refluxing parameter error (input deck variable backward distinct from 1 or 2)'
print*,'=> No specular reflection : backward = 1 has been imposed!'
end if
if ((irradiated_side_refluxing .eqv. .true.) .and. (backward .eq. 1)) then
irradiated_side_refluxing = .false.
print*,'Fast electron refluxing at the target irradiated side can be taken into account'
print*,'only if their refluxing at the target rear side is taken into account .i.e. if backward = 2'
end if
end subroutine read_init_parameters

subroutine get_str(str)
character(len=*), intent(inout) :: str
read (1,'(A)',end=10) str
str = adjustl(str)
return
continue ; 10 str = 'end'
end subroutine get_str

elemental function get_logical(str)
character(len=*), intent(in) :: str
logical :: get_logical
read (str, *) get_logical
return
end function get_logical

elemental function get_integer(str)
character(len=*), intent(in) :: str
integer :: get_integer
read (str, *) get_integer
return
end function get_integer

elemental function get_real(str)
character(len=*), intent(in) :: str
real(PR) :: get_real
read (str, *) get_real
return
end function get_real

elemental function get_char(str)
character(len=*), intent(in) :: str
Character(len=2) :: get_char
read (str, *) get_char
return
end function get_char

end module input