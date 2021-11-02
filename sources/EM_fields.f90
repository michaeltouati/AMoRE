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
module em_fields

  use acuracy
  use constants
  use transport_coef
  use input

  implicit none
  private :: matrix_vector_product, scalar_product
  private :: create_diffusion_matrix, Jacobi_preconditioner, conj_grad 
  public  :: magnetic_field, electric_field

  contains

  subroutine magnetic_field(eta_tab, Zn, ni, Te, Ti, jbx, jbz, dt, By_n, By_np1)    
    ! input  : eta_tab = tabulated resistivity if provided by the user 
    !                    (tabulated_resistivity = .true.)
    !          Zn      = table containing the atomic number in ()
    !          ni      = table containing the ion density in (/cm^3)
    !          Te(i)   = table containing the background electron (ion) temperature at t(n) 
    !                    in (K)
    !          jbx(z)  = table containing the beam current density component on x(z)-axis
    !                    at t(n) in (statAmpere/cm)
    !          dt      = time step in (fs)
    !          By_n    = table containing the magnetic field at time t(n) in (gauss)
    ! output : By_np1  = table containing the magnetic field at time t(n+1) in (gauss)
    implicit none
    real(PR), intent(in)                                    :: dt
    real(kind=PR), dimension(1:N_eta_tab,1:2), intent(in)   :: eta_tab
    real(PR), dimension(0:N_x+1,0:N_z+1), intent(in)        :: Zn, ni, Te, Ti
    real(PR), dimension(1:N_x,1:N_z), intent(in)            :: jbx, jbz, By_n
    real(PR), dimension(1:N_x,1:N_z), intent(out)           :: By_np1
    real(PR),dimension(1:N_x,1:N_z)                         :: eta
    real(PR),dimension(1:N_x*N_z)                           :: diag1,diag2,diag3,diag4,diag5
    real(PR),dimension(1:N_x*N_z)                           :: x, sol
    integer                                                 :: i, k, l
    real(PR)                                                :: Bcour, Bres, Bcross
    real(PR)                                                :: ne_i_k, ne_i_km1, ne_i_kp1
    real(PR)                                                :: ne_im1_k, ne_ip1_k
    real(PR)                                                :: eta_i_k, eta_im1_k, eta_ip1_k
    real(PR)                                                :: eta_i_km1, eta_i_kp1
    !$omp  PARALLEL DO DEFAULT(SHARED) &
    !$omp& PRIVATE(i,k,Bcour,Bres,Bcross) &
    !$omp& PRIVATE(eta_i_k, eta_im1_k, eta_ip1_k, eta_i_km1, eta_i_kp1) &
    !$omp& PRIVATE(ne_i_k, ne_i_km1, ne_i_kp1, ne_im1_k, ne_ip1_k) &
    !$omp& COLLAPSE(2)
    do k = 2,N_z-1,1
      do i=2,N_x-1,1   
        ! Electrical resistivity values
        eta_i_k   = resis(eta_tab, Zn(i,k),   ni(i,k),   Te(i,k),   Ti(i,k))
        eta_i_kp1 = resis(eta_tab, Zn(i,k+1), ni(i,k+1), Te(i,k+1), Ti(i,k+1))
        eta_i_km1 = resis(eta_tab, Zn(i,k-1), ni(i,k-1), Te(i,k-1), Ti(i,k-1))
        eta_ip1_k = resis(eta_tab, Zn(i+1,k), ni(i+1,k), Te(i+1,k), Ti(i+1,k))
        eta_im1_k = resis(eta_tab, Zn(i-1,k), ni(i-1,k), Te(i-1,k), Ti(i-1,k))
        ! Background electron density values
        ne_i_k   = zeff(Zn(i,k),   ni(i,k),   Te(i,k))   * ni(i,  k)
        ne_ip1_k = zeff(Zn(i+1,k), ni(i+1,k), Te(i+1,k)) * ni(i+1,k)
        ne_im1_k = zeff(Zn(i-1,k), ni(i-1,k), Te(i-1,k)) * ni(i-1,k)
        ne_i_kp1 = zeff(Zn(i,k+1), ni(i,k+1), Te(i,k+1)) * ni(i,k+1)
        ne_i_km1 = zeff(Zn(i,k-1), ni(i,k-1), Te(i,k-1)) * ni(i,k-1)
        ! Contribution due to the curl of the electron beam current density :
        Bcour = dt*fs*c*eta_i_k*(&
        &((jbx(i,k+1)-jbx(i,k-1))/(2._PR*d_z*microns))&
        &-&
        &((jbz(i+1,k)-jbz(i-1,k))/(2._PR*d_x*microns))&
        &)
        ! Contribution due to the resistivity gradients :
        Bres = dt*fs*c*(&
        &(jbx(i,k)*(eta_i_kp1-eta_i_km1)/(2._PR*d_z*microns))&
        &-&
        &(jbz(i,k)*(eta_ip1_k-eta_im1_k)/(2._PR*d_x*microns))&
        &)          
        ! Contribution due to the temperature-density crossed gradients :
        Bcross = dt*fs*(kb*c/ne_i_k)*(&
        &((ne_ip1_k-ne_im1_k)*(Te(i,k+1)-Te(i,k-1))/(4._PR*d_x*d_z*(microns**2._PR)))&
        &-&
        &((ne_i_kp1-ne_i_km1)*(Te(i+1,k)-Te(i-1,k))/(4._PR*d_x*d_z*(microns**2._PR)))&
        &)
        ! Magnetic field source terms :
        By_np1(i,k) = By_n(i,k) + Bcour + Bres + Bcross         
      end do
    end do
    !$omp END PARALLEL DO
    ! Boundary conditions :
    ! - the value of the magnetic field in the longitudinal guard cells are taken equal to the
    !   value in the closest cell inside the simulation box
    ! - the value of the magnetic field in the transverse guard cells are taken equal to 0
    !   assuming the transverse dimension of the simulation box have been chosen sufficiently
    !   large
    ! x = +/- L_x/2 :
    do k = 1,N_z,1
      By_np1(1,k)       =  By_np1(2,k)
      By_np1(N_x,k)     =  By_np1(N_x-1,k)
    end do
    !z = 0 and z = L_z :
    do i = 1,N_x,1
      By_np1(i,1)       =  By_np1(i,2)
      By_np1(i,N_z)     =  By_np1(i,N_z-1)
    end do
    ! Computation of the magnetic field diffusion according to the second order implicit
    ! scheme inverted according to the conjugated gradients algorithm :
    if (magnetic_diffusion) then    
      ! 1) Creation of the diffusion matrix and the second member vector
      do l=1,N_x*N_z,1
        k = int((l-1)/N_x)+1
        i = l-(k-1)*N_x
        x(l)       = By_np1(i,k)
        eta(i,k)   = resis(eta_tab, Zn(i,k), ni(i,k), Te(i,k), Ti(i,k))
      end do
      call create_diffusion_matrix(dt,eta,diag1,diag2,diag3,diag4,diag5)    
      ! Jacobi's Preconditioner    
      call Jacobi_preconditioner(diag1,diag2,diag3,diag4,diag5,x)    
      ! Conjugated gradients algorithm    
      call conj_grad(diag1,diag2,diag3,diag4,diag5,x,sol,zero)
      do l=1,N_x*N_z,1
        k = int((l-1)/N_x)+1
        i = l-(k-1)*N_x
        By_np1(i,k) = sol(l)
      end do
    end if
  end subroutine magnetic_field

  subroutine electric_field(eta_tab, Zn, ni, Te, Ti, jrx, jrz, Ex, Ez)    
    ! input  : eta_tab = tabulated resistivity if provided by the user 
    !                    (tabulated_resistivity = .true.)
    !          Zn      = table containing the atomic number in ()
    !          ni      = table containing the ion density in (/cm^3)
    !          Te(i)   = table containing the background electron (ion) temperature at t(n) 
    !                    in (K)
    !          jrx(z)  = table containing the return current density component on x(z)-axis
    !                    at t(n) in (statAmpere/cm)
    ! output : Ex(z)  = table containing the electric field component on x(z)-axis at time t(n)  
    !                   in (statVolt/cm)
    implicit none
    real(kind=PR), dimension(1:N_eta_tab,1:2), intent(in)  :: eta_tab
    real(PR), dimension(0:N_x+1,0:N_z+1), intent(in)       :: Zn, ni, Te, Ti
    real(PR), dimension(1:N_x,1:N_z), intent(in)           :: jrx, jrz
    real(PR), dimension(1:N_x,1:N_z), intent(out)          :: Ex, Ez
    real(PR)                                               :: eta_temp, ne, Pe_mup1, Pe_mum1
    integer                                                :: i, k
    real(PR)                                               :: Excour, Ezcour, Expres, Ezpres
    !$omp  PARALLEL DO DEFAULT(SHARED) &
    !$omp& PRIVATE(i, k, eta_temp, ne, Pe_mup1, Pe_mum1, Excour, Ezcour, Expres, Ezpres) &
    !$omp& COLLAPSE(2)
    do k= 2,N_z-1,1
      do i=2,N_x-1,1
        ! Temporary values
        eta_temp = resis(eta_tab, Zn(i,k), ni(i,k), Te(i,k), Ti(i,k))
        ! Contribution due to the Ohm effect
        Excour = eta_temp*jrx(i,k)
        Ezcour = eta_temp*jrz(i,k)        
        ! temporary values for mu = i
        ne      = zeff(Zn(i,k), ni(i,k), Te(i,k)) * ni(i,k)
        Pe_mup1 = zeff(Zn(i+1,k), ni(i+1,k), Te(i+1,k)) * ni(i+1,k) * kb * Te(i+1,k)
        Pe_mum1 = zeff(Zn(i-1,k), ni(i-1,k), Te(i-1,k)) * ni(i-1,k) * kb * Te(i-1,k)
        ! Contribution due to the pressure gradients in the x-direction
        Expres = - (Pe_mup1 - Pe_mum1) / (2._PR * d_x*microns * e * ne)          
        ! temporary values mu = k
        Pe_mup1 = zeff(Zn(i,k+1), ni(i,k+1), Te(i,k+1)) * ni(i,k+1) * kb * Te(i,k+1)
        Pe_mum1 = zeff(Zn(i,k-1), ni(i,k-1), Te(i,k-1)) * ni(i,k-1) * kb * Te(i,k-1)
        ! Contribution due to the pressure gradients in the z-direction
        Ezpres = - (Pe_mup1-Pe_mum1)/ (2._PR * d_z*microns * e * ne)        
        ! Electric field
        Ex(i,k) =  Excour + Expres
        Ez(i,k) =  Ezcour + Ezpres         
      end do
    end do
    !$omp END PARALLEL DO
    ! Boundary conditions   
    ! x = -Lx/2  and x = Lx/2
    do k = 1,N_z,1 
      Ex(1  ,k)  = Ex(2    ,k)
      Ez(1  ,k)  = Ez(2    ,k)
      Ex(N_x,k)  = Ex(N_x-1,k)
      Ez(N_x,k)  = Ez(N_x-1,k)   
    end do
    ! z = 0 and z = Lz
    do i = 1,N_x,1 
      Ex(i,1)    = Ex(i,2)
      Ez(i,1)    = Ez(i,2)
      Ex(i,N_z)  = Ex(i,N_z-1)
      Ez(i,N_z)  = Ez(i,N_z-1)
    end do
  end subroutine electric_field

  !========================================================================================
  !           Create the second order implicit scheme matrix to be inverted
  !========================================================================================

  pure subroutine create_diffusion_matrix(dt,eta,diag1,diag2,diag3,diag4,diag5)
    ! input  : dt  : time step in (fs)
    !          eta : table containing the electrical resistivity values in (s)
    ! output : block-tridiagonal diffusion matrix : | d1  d4  ..  (0) ..  d5  ..  (0) | 
    !                                               | d2                              |
    !                                               | :                               |
    !                                           A = | (0)                             |
    !                                               | :                               |
    !                                               | d3                              |
    !                                               | :                               |
    !                                               | (0)                             |       
    implicit none
    real(PR), intent(in)                           :: dt
    real(PR),dimension(1:N_x,1:N_z), intent(in)    :: eta
    real(PR),dimension(1:N_x*N_z), intent(out)     :: diag1, diag2, diag3, diag4, diag5
    integer                                        :: l,i,k,N
    real(PR)                                       :: eta_iph_k, eta_imh_k
    real(PR)                                       :: eta_i_kph, eta_i_kmh 
    N=N_x*N_z
    do l=1,N,1
      k = int((l-1)/N_x)+1
      i = l-(k-1)*N_x
      ! Diffusion coefficient values at the different interfaces :
      if (i.ne.N_x) then
        eta_iph_k = 2._PR * eta(i+1,k) * eta(i,k)   / (eta(i+1,k) + eta(i,k)  )
      else
        eta_iph_k = eta(i,k)
      end if
      if (i.ne.1) then
        eta_imh_k = 2._PR * eta(i,k)   * eta(i-1,k) / (eta(i,k)   + eta(i-1,k))
      else
        eta_imh_k = eta(i,k)
      end if
      if (k.ne.N_z) then
        eta_i_kph = 2._PR * eta(i,k+1) * eta(i,k)   / (eta(i,k+1) + eta(i,k)  )
      else
        eta_i_kph = eta(i,k)
      end if
      if (k.ne.1) then
        eta_i_kmh = 2._PR * eta(i,k)   * eta(i,k-1) / (eta(i,k)   + eta(i,k-1))
      else
        eta_i_kmh = eta(i,k)
      end if
      !d1 :
      diag1(l) = 1._PR + ((eta_iph_k+eta_imh_k)*(dt*fs)*(c**2._PR)/(4._PR*pi*((d_x*microns)**2._PR)))&
                      &+ ((eta_i_kph+eta_i_kmh)*(dt*fs)*(c**2._PR)/(4._PR*pi*((d_z*microns)**2._PR)))
      !d2 :
      if (i.ne.N_x) then
        diag2(l) = - eta_iph_k*(dt*fs)*(c**2._PR)/(4._PR*pi*((d_x*microns)**2._PR))
      else
        diag2(l) = 0._PR
      end if
      !d3 :   
      diag3(l) = -eta_i_kph*(dt*fs)*(c**2._PR)/(4._PR*pi*((d_z*microns)**2._PR))
      !d4 :
      if (i.ne.1) then
        diag4(l) = -eta_imh_k*(dt*fs)*(c**2._PR)/(4._PR*pi*((d_x*microns)**2._PR))
      else
        diag4(l)=0._PR
      end if
      !d5 :
      diag5(l) = - eta_imh_k*(dt*fs)*(c**2._PR)/(4._PR*pi*((d_z*microns)**2._PR))      
    end do  
  end subroutine create_diffusion_matrix

  pure subroutine Jacobi_preconditioner(diag1,diag2,diag3,diag4,diag5,vector)
    ! in(out)put : diag1... = block-tridiagonal diffusion matrix A
    !                         (see the subroutine create_diffusion_matrix)
    !              vector   = second member of the linear equation to be inverted
    implicit none
    real(PR),dimension(1:N_x*N_z),intent(inout)   :: diag1, diag2, diag3, diag4, diag5
    real(PR),dimension(1:N_x*N_z),intent(inout)   :: vector
    integer                                       :: N, i, k, l
    N = N_x * N_z
    do l=1,N,1
      k = int((l-1)/N_x)+1
      i = l-(k-1)*N_x
      if (i == N_x) then
        diag2(l) = 0._PR
      else
        diag2(l) = diag2(l)/diag1(l+1) 
      end if      
      if (k <= (N_z-1)) then
        diag3(l) = diag3(l)/diag1(l+N_x)
      end if
      if (i == 1) then
        diag4(l) = 0._PR
      else
        diag4(l) = diag4(l)/diag1(l-1) 
      end if      
      if (k >= 2) then                                                                         
        diag5(l) = diag5(l)/diag1(l-N_x) 
      end if       
      vector(l) = vector(l)/diag1(l)
      diag1(l)  = 1._PR
    end do
    ! boundary conditions at z = 0 and z = L_z :
    do l=1,N,1
      k = int((l-1)/N_x)+1       
      if (k == 1)   diag5(l) = diag5(l+N_x)                         
      if (k == N_z) diag3(l) = diag3(l-N_x)
    end do
  end subroutine Jacobi_preconditioner

  !========================================================================================
  !                             Matrix-vector and scalar products
  !========================================================================================

  pure subroutine matrix_vector_product(diag1,diag2,diag3,diag4,diag5,vector_in,vector_out)
    ! input  : diag1...   = block-tridiagonal diffusion matrix A
    !                       (see the subroutine create_diffusion_matrix)
    !          vector_in  = vector x
    ! output : vector out = A.x
    implicit none
    real(PR),dimension(1:N_x*N_z),intent(in)      :: diag1, diag2, diag3, diag4, diag5
    real(PR),dimension(1:N_x*N_z),intent(in)      :: vector_in
    real(PR),dimension(1:N_x*N_z),intent(out)     :: vector_out
    integer                                       :: N, i, k, l
    N=N_x*N_z
    do l=1,N,1
      k = int((l-1)/N_x)+1
      i = l-(k-1)*N_x
      if ((k >= 2).and.(k <= (N_z-1)).and.(i >= 2).and.(i <= (N_x-1))) then
        vector_out(l) = diag3(l) * vector_in(l+N_x)&
        + diag2(l) * vector_in(l+1)&
        + diag1(l) * vector_in(l)&
        + diag4(l) * vector_in(l-1)&
        + diag5(l) * vector_in(l-N_x)
      else if ((k >= 2).and.(k <= (N_z-1)).and.(i == 1)) then
        vector_out(l) = diag3(l) * vector_in(l+N_x)&
        + diag2(l) * vector_in(l+1)&
        + diag1(l) * vector_in(l)&
        + diag5(l) * vector_in(l-N_x)
      else if ((k >= 2).and.(k <= (N_z-1)).and.(i == N_x)) then
        vector_out(l) = diag3(l) * vector_in(l+N_x)&
        + diag1(l) * vector_in(l)&
        + diag4(l) * vector_in(l-1)&
        + diag5(l) * vector_in(l-N_x)
      else if ((k == 1).and.(i >= 2).and.(i <= (N_x-1))) then
        vector_out(l) = diag3(l) * vector_in(l+N_x)&
        + diag2(l) * vector_in(l+1)&
        + diag1(l) * vector_in(l)&
        + diag4(l) * vector_in(l-1)&
        + diag5(l) * vector_in(l)             ! boundary condition
      else if ((k == 1).and.(i == 1)) then
        vector_out(l) = diag3(l) * vector_in(l+N_x)&
        + diag2(l) * vector_in(l+1)&
        + diag1(l) * vector_in(l)&
        + diag5(l) * vector_in(l)             ! boundary condition
      else if ((k == 1).and.(i == N_x)) then
        vector_out(l) = diag3(l) * vector_in(l+N_x)&
        + diag1(l) * vector_in(l)&
        + diag4(l) * vector_in(l-1)&
        + diag5(l) * vector_in(l)             ! boundary condition
      else if ((k == N_z).and.(i >= 2).and.(i <= (N_x-1))) then
        vector_out(l) = diag2(l) * vector_in(l+1)&
        + diag1(l) * vector_in(l)&
        + diag4(l) * vector_in(l-1)&
        + diag5(l) * vector_in(l-N_x)&
        + diag3(l) * vector_in(l)             ! boundary condition
      else if ((k == N_z).and.(i == 1)) then
        vector_out(l) = diag2(l) * vector_in(l+1)&
        + diag1(l) * vector_in(l)&
        + diag5(l) * vector_in(l-N_x)&
        + diag3(l) * vector_in(l)             ! boundary condition
      else if ((k == N_z).and.(i == N_x)) then
        vector_out(l) = diag1(l) * vector_in(l)&
        + diag4(l) * vector_in(l-1)&
        + diag5(l) * vector_in(l-N_x)&
        + diag3(l) * vector_in(l)             ! boundary condition
      end if
    end do
  end subroutine matrix_vector_product

  pure subroutine scalar_product(Vector_1,Vector_2,res)   
    ! input  : vector_1 and vector_2
    ! output : res = scalar product (vector_1,vector_2)
    implicit none
    real(PR),dimension(1:N_x*N_z),intent(in) :: vector_1, vector_2
    real(PR),intent(out)                       :: res
    integer                                    :: N    
    N=N_x*N_z
    res = sum(vector_1(1:N)*vector_2(1:N))  
  end subroutine scalar_product

  !========================================================================================
  !                            Conjugated gradients algorithm
  !========================================================================================

  subroutine conj_grad(diag1,diag2,diag3,diag4,diag5,sm,sol,error)    
    ! input  : diag1...  = block-tridiagonal diffusion matrix A
    !                      (see the subroutine create_diffusion_matrix)
    !          sm        = second member y of the linear equation to be inverted
    !                      A.x = y
    !          error     = residue wanted
    ! output : sol       = solution x
    implicit none
    real(PR),dimension(1:N_x*N_z),intent(in)      :: diag1, diag2, diag3, diag4, diag5
    real(PR),dimension(1:N_x*N_z),intent(in)      :: sm
    real(PR),intent(in)                           :: error
    real(PR),dimension(1:N_x*N_z),intent(out)     :: sol
    real(PR),dimension(1:N_x*N_z)                 :: y, z, rest, dir
    real(PR)                                      :: coef, sp_num, sp_den, rest_norm, error_norm
    integer                                       :: m, ex_wh
    !initialization for the resolution of Ax=b where A=(diag1,...), x=B(n+1) et b=B(n)
    ! sol = B(n+1) :
    sol  = sm
    call  matrix_vector_product(diag1,diag2,diag3,diag4,diag5,sol,y)
    ! rest = b - A.x :
    rest = sm - y
    dir = rest
    m = 0  
    ! compute the wanted residue error_norm :
    call scalar_product(rest,rest,rest_norm)
    rest_norm = sqrt(rest_norm)
    call scalar_product(sm,sm,error_norm)
    error_norm = sqrt(error_norm)
    error_norm = error*error_norm
    ! conjugated gradients algorithm :
    ex_wh=0
    do while ((rest_norm >= error_norm).and.(ex_wh == 0))
      ! 1)   Compute <rest_i,rest_i>/<A,dir_i,dir_i> :
      ! 1.1) Compute <rest_i,rest_i> :
      call scalar_product(rest,rest,sp_num)
      ! 1.2) Compute <A.dir_i,dir_i> :
      call matrix_vector_product(diag1,diag2,diag3,diag4,diag5,dir,z)
      call scalar_product(dir,z,sp_den)
      if (abs(sp_den) >= error) then ! ensure not dividing by 0 if <A.dir_i,dir_i> = 0!
        ! 1.3) Compute <rest_i,rest_i>/<A.dir_i,dir_i> :
        coef = sp_num / sp_den
        ! 1.4) Conpute sol_i+1 = sol_i + (<rest_i,rest_i>/<A.dir_i,dir_i>) dir_i :
        sol = sol + coef * dir  
        ! 1.4) Conpute rest_i+1 = rest_i - (<rest_i,rest_i>/<A.dir_i,dir_i>) A.dir_i :
        rest = rest - coef * z  
        ! 2)   Compute the coefficient for dir_i+1 : <rest_i+1,rest_i+1>/ <rest_i,rest_i> :
        ! 2.1) Compute <rest_i,rest_i> :
        sp_den = sp_num           
        ! 2.2) Compute  <rest_i+1,rest_i+1> :
        call scalar_product(rest,rest,sp_num)          
        coef = sp_num/sp_den
        ! 2.3) Compute  dir_i+1 = rest_i+1 + (<rest_i+1,rest_i+1>/ <rest_i,rest_i>) dir_i :
        dir = rest + coef * dir
        ! 3) Compute the residue foor the next step r_i+1 :
        call scalar_product(rest,rest,rest_norm)
        rest_norm = sqrt(rest_norm)
        m=m+1
      else
        ex_wh=1
      end if
    end do
  end subroutine conj_grad

end module em_fields