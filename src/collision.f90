module collision
  use Types, only: wp=>dp
  use soot_params
  use constants
  implicit none

  public :: getKernel

  real(wp) :: P_in, T_in, rho_in, m_mono_in, v_mono_in

  abstract interface
     function cvg_kernel(i, j) result(beta)
       use Types, only: wp=>dp
       integer, intent(in) :: i,j
       real(wp) :: beta
     end function cvg_kernel
  end interface

contains

  function getKernel(T, P, rho, m_mono) result(kernel)
    real(wp), intent(in) :: T, P, rho, m_mono
    procedure(cvg_kernel), pointer :: kernel
    T_in = T
    P_in = P
    rho_in = rho
    m_mono_in = m_mono
    v_mono_in = m_mono_in/rho_in
    kernel => kernel_fuchs 
  end function getKernel

  function kernel_fuchs (i, j) result(beta)
    integer, intent(in) :: i, j
    real(wp) :: beta

    real(wp) :: mu, mfp, rpi, rpj, kni, knj, a, b, c
    real(wp) :: Cc_i, Cc_j, Bi, Bj, Di, Dj, mi, mj
    real(wp) :: ci, cj, tau_i, tau_j, lbi, lbj
    real(wp) :: deltai, deltaj, r_mean, beta_corr

    mu = MUREF*(TREF+110.4_wp)/(T_in+110.4_wp)*(T_in/TREF)**(1.5_wp)  !corrected viscosity
    mfp = MFPREF*(T_in/TREF)*(PREF/P_in)*(1_wp+110.4_wp/TREF)/(1.0_wp+110.4_wp/T_in) !corrected mean free path
    
    rpi = ((i*m_mono_in)/rho_in*6.0_wp/PI)**(1.0_wp/3.0_wp)/2.0_wp
    rpj = ((j*m_mono_in)/rho_in*6_wp/PI)**(1.0_wp/3.0_wp)/2.0_wp !size in m
    
    kni = mfp/rpi;
    knj = mfp/rpj;
    
    a = 1.257_wp !Cunningham correction coefficients
    b = 0.4_wp
    c = 1.1_wp
    
    Cc_i = 1.0_wp + kni*(a + b*exp(-c/kni));
    Cc_j = 1.0_wp + knj*(a + b*exp(-c/knj));
    
    Bi = Cc_i/(6.0_wp*PI*mu*rpi) !s/kg
    Bj = Cc_j/(6.0_wp*PI*mu*rpj)
    Di = Bi*BOLTZ*T_in; !m^2/s, close to Table 2.1 in Friedlander's book
    Dj = Bj*BOLTZ*T_in;
    
    mi = m_mono_in*i !mass
    mj = m_mono_in*j !mass
    
    ci = (8.0_wp*BOLTZ*T_in/(PI*mi))**0.5_wp ! m/s, mean thermal velocity of the particles
    cj = (8.0_wp*BOLTZ*T_in/(PI*mj))**0.5_wp
    
    !mean traveling time
    tau_i = (i*m_mono_in*Cc_i)/6.0_wp/PI/mu/rpi;
    tau_j = (j*m_mono_in*Cc_j)/6.0_wp/PI/mu/rpj;
    
    !'mean free path' of the particles
    lbi = tau_i*ci;
    lbj = tau_j*cj;
    
    !detla values
    deltai = 1.0_wp/6.0_wp/rpi/lbi*((2.0_wp*rpi + lbi)**(3.0_wp) &
         &- (4.0_wp*rpi**(2.0_wp)+lbi**(2.0_wp))**(1.5_wp))-2.0_wp*rpi;
    deltaj = 1.0_wp/6.0_wp/rpj/lbj*((2.0_wp*rpj + lbj)**(3.0_wp) &
         &- (4.0_wp*rpj**(2.0_wp)+lbj**(2.0_wp))**(1.5_wp))-2.0_wp*rpj;
    
    r_mean = (rpi+rpj)/2.0_wp;
    beta_corr = r_mean/(r_mean + sqrt(deltai**2.0_wp + deltaj**2.0_wp)/2.0_wp) + 4.0_wp*(Di+Dj)/2.0_wp/sqrt(ci**(2.0_wp)+cj**(2.0_wp))/r_mean;
    beta = 2.0_wp*8.0_wp*PI*(rpi+rpj)/2.0_wp*(Di+Dj)/2.0_wp/beta_corr*1.0e6_wp;
  end function kernel_fuchs




  

end module collision
