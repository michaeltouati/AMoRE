module constants

implicit none
! Number of bytes occupied in the memory by one real number :
integer, parameter        :: PR = 8 
! Rounding error for diagnostics :                                                
real(PR),parameter,public :: zero = 1.E-15_PR

!=====================================================================================
!                              Mathematical constants
!=====================================================================================

real(PR),parameter,public :: pi      = 3.14159
real(PR),parameter,public :: Catalan = 0.9159656

!=====================================================================================
!                          Physical_constants constants
!=====================================================================================

! Atomic mass unit [g] :
real(PR),parameter,public :: mu   = 1.6605E-24_PR      
! Proton mass unit [g] :                     
real(PR),parameter,public :: mp   = 1.6726E-24_PR           
! Electron mass [g] :                
real(PR),parameter,public :: me   = 9.1094E-28_PR
! Elementary charge [statcoulomb] :                  
real(PR),parameter,public :: e    = 4.8032E-10_PR
! Speed of light in vacuum [cm/s] :                           
real(PR),parameter,public :: c    = 2.9979E10_PR
! me c^2 en keV :                             
real(PR),parameter,public :: mec2 = me*(c**2._PR)*1.E-10_PR / 1.6022E-19_PR
! Boltzmann constant [erg/K] :  
real(PR),parameter,public :: kb   = 1.3807E-16_PR
! Planck constant [erg.s] :                           
real(PR),parameter,public :: h    = 6.6261E-27_PR                           
real(PR),parameter,public :: hbar = h / (2*pi)
! Fine Structure constant :
real(PR),parameter,public :: alpha = 1._PR / 137.035999_PR                   

!=====================================================================================
!                              Conversion factors
!=====================================================================================

! 1 eV = 11604 K :
real(PR),parameter,public   :: eV = 1.6022E-19_PR / 1.3807E-23_PR
! 1 keV in erg :             
real(PR),parameter,public   :: keV = 1.E10_PR * 1.60022E-19_PR
! 1 fs in s :                
real(PR), parameter, public :: fs = 1.E-15_PR
! 1 micron in cm :
real(PR), parameter, public :: microns = 1.E-4_PR
! Ambiant temperature [K] :
real(PR), parameter, public :: Tamb=300
! 1 J in erg :
real(PR), parameter, public :: Joules = 1.E7_PR                              

!=====================================================================================
!                 indices of the zero and first order angular moments
!=====================================================================================

integer, parameter,public :: psi0   = 1
integer, parameter,public :: psi1x  = 2
integer, parameter,public :: psi1z  = 3

!=====================================================================================
!          indices of the laser-generated electron beam and the refluxing one 
!=====================================================================================

integer, parameter,public :: forward  = 1

!=====================================================================================
!                              Solid properties
!=====================================================================================
! Atomic weight :
real(PR), parameter, public    :: A_Al = 26.9815_PR
real(PR), parameter, public    :: A_Cu = 63.546_PR
real(PR), parameter, public    :: A_Ta = 180.9479_PR
real(PR), parameter, public    :: A_H  = 1.0079_PR
real(PR), parameter, public    :: A_C  = 12.011_PR
real(PR), parameter, public    :: A_CH = 6.50945_PR
! Atomic number :
real(PR), parameter, public    :: Z_Al = 13._PR
real(PR), parameter, public    :: Z_Cu = 29._PR
real(PR), parameter, public    :: Z_Ta = 73._PR
real(PR), parameter, public    :: Z_H  = 1._PR
real(PR), parameter, public    :: Z_C  = 6._PR
real(PR), parameter, public    :: Z_CH = 3.5_PR
! Ion density :
real(PR), parameter, public    :: ni_Al = 2.6989_PR / (mu * A_Al)
real(PR), parameter, public    :: ni_Cu = 8.96_PR   / (mu * A_Cu)
real(PR), parameter, public    :: ni_Ta = 16.69_PR  / (mu * A_Ta)
real(PR), parameter, public    :: ni_H = 10._PR     / (mu * A_H)
real(PR), parameter, public    :: ni_C = 2.1_PR     / (mu * A_C)
real(PR), parameter, public    :: ni_CH = 1.05_PR   / (mu * A_CH)

end module constants