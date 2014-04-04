!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine slater (rs, ex, vx)
  !-----------------------------------------------------------------------
  !        Slater exchange with alpha=2/3
  !
  implicit none
  real(8) :: rs, ex, vx
  real(8), parameter  :: f= -0.687247939924714d0, alpha = 2.0d0/3.0d0
  ! f = -9/8*(3/2pi)^(2/3)
  !
  ex = f * alpha / rs
  vx = 4.d0 / 3.d0 * f * alpha / rs
  !
  return
end subroutine slater
!
!-----------------------------------------------------------------------
subroutine slater1(rs, ex, vx)
  !-----------------------------------------------------------------------
  !        Slater exchange with alpha=1, corresponding to -1.374/r_s Ry
  !        used to recover old results
  !
  implicit none
  real(8) :: rs, ex, vx
  real(8), parameter  :: f= -0.687247939924714d0, alpha = 1.0d0
  !
  ex = f * alpha / rs
  vx = 4.d0 / 3.d0 * f * alpha / rs
  !
  return
end subroutine slater1
!
!-----------------------------------------------------------------------
subroutine slater_rxc (rs, ex, vx)
  !-----------------------------------------------------------------------
  !        Slater exchange with alpha=2/3 and Relativistic exchange
  !
  IMPLICIT none
  REAL(8), PARAMETER :: pi     = 3.14159265358979323846_8 
  REAL(8), PARAMETER :: tpi    = 2.0_8 * pi
  REAL(8), PARAMETER :: HARTREE_SI       = 4.35974394E-18_8   ! J
  REAL(8), PARAMETER :: C_SI             = 2.99792458E+8_8    ! m sec^-1
  REAL(8), PARAMETER :: BOHR_RADIUS_SI   = 0.52917720859E-10_8 ! m
  REAL(8), PARAMETER :: H_PLANCK_SI      = 6.62606896E-34_8   ! J s
  REAL(8), PARAMETER :: AU_SEC           = H_PLANCK_SI/tpi/HARTREE_SI
  REAL(8), PARAMETER :: C_AU             = C_SI / BOHR_RADIUS_SI * AU_SEC
  real (8):: rs, ex, vx
  !
  real(8), PARAMETER :: ZERO=0.D0, ONE=1.D0, PFIVE=.5D0, &
       OPF=1.5D0 !, C014=0.014D0
  real (8):: trd, ftrd, tftm, a0, alp, z, fz, fzp, vxp, xp, &
       beta, sb, alb, c014
  !
  TRD = ONE/3.d0
  FTRD = 4.d0*TRD
  TFTM = 2**FTRD-2.d0
  A0 = (4.d0/(9.d0*PI))**TRD
  C014= 1.0_8/a0/c_au
  
  !      X-alpha parameter:
  ALP = 2.d0 * TRD
  
  Z = ZERO
  FZ = ZERO
  FZP = ZERO
  
  VXP = -3.d0*ALP/(2.d0*PI*A0*RS)
  XP = 3.d0*VXP/4.d0
  BETA = C014/RS
  SB = SQRT(1.d0+BETA*BETA)
  ALB = LOG(BETA+SB)
  VXP = VXP * (-PFIVE + OPF * ALB / (BETA*SB))
  XP = XP * (ONE-OPF*((BETA*SB-ALB)/BETA**2)**2)
  !  VXF = 2**TRD*VXP
  !  EXF = 2**TRD*XP
  VX = VXP
  EX = XP
END SUBROUTINE slater_rxc

!
!-----------------------------------------------------------------------
  subroutine slaterKZK (rs, ex, vx, vol)
  !-----------------------------------------------------------------------
  !        Slater exchange with alpha=2/3, Kwee, Zhang and Krakauer KE
  !        correction
  !
  implicit none
  real(8) :: rs, ex, vx, dL, vol, ga, pi, a0
  real(8), parameter  ::  a1 = -2.2037d0, &
              a2 = 0.4710d0, a3 = -0.015d0, ry2h = 0.5d0
  real(8), parameter  :: f= -0.687247939924714d0, alpha = 2.0d0/3.0d0
  ! f = -9/8*(3/2pi)^(2/3)
  !
  pi = 4.d0 * atan(1.d0)
  a0 = f * alpha * 2.d0

  dL = vol**(1.d0/3.d0)
  ga = 0.5d0 * dL *(3.d0 /pi)**(1.d0/3.d0)
  !
  if ( rs .le. ga) then
   ex = a0 / rs + a1 * rs / dL**2.d0 + a2 * rs**2.d0 / dL**3.d0
   vx = (4.d0 * a0 / rs + 2.d0 * a1 * rs / dL**2.d0 + &
        a2 * rs**2.d0 / dL**3.d0 ) / 3.d0
  else
    ex = a0 / ga + a1 * ga / dL**2.d0 + a2 * ga**2.d0 / dL**3.d0 ! solids
    vx = ex
!   ex = a3 * dL**5.d0 / rs**6.d0     ! molecules
!   vx = 3.d0 * ex  
  endif

  ex = ry2h * ex    ! Ry to Hartree
  vx = ry2h * vx
  !
  return
end subroutine slaterKZK
!
!-----------------------------------------------------------------------
subroutine pz (rs, iflag, ec, vc)
  !-----------------------------------------------------------------------
  !     LDA parameterization from Monte Carlo data
  !     iflag=1: J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
  !     iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
  !
  implicit none
  real(8) :: rs, ec, vc
  integer :: iflag
  !
  real(8) :: a (2), b (2), c (2), d (2), gc (2), b1 (2), b2 (2)
  real(8) :: lnrs, rs12, ox, dox
  !
  data a / 0.0311d0, 0.031091d0 /, b / -0.048d0, -0.046644d0 /, &
       c / 0.0020d0, 0.00419d0 /, d / -0.0116d0, -0.00983d0 /
  data gc / -0.1423d0, -0.103756d0 /, b1 / 1.0529d0, 0.56371d0 /, &
       b2 / 0.3334d0, 0.27358d0 /
  !
  if (rs.lt.1.0d0) then
     ! high density formula
     lnrs = log (rs)
     ec = a (iflag) * lnrs + b (iflag) + c (iflag) * rs * lnrs + d ( &
          iflag) * rs
     vc = a (iflag) * lnrs + (b (iflag) - a (iflag) / 3.d0) + 2.d0 / &
          3.d0 * c (iflag) * rs * lnrs + (2.d0 * d (iflag) - c (iflag) ) &
          / 3.d0 * rs
  else
     ! interpolation formula
     rs12 = sqrt (rs)
     ox = 1.d0 + b1 (iflag) * rs12 + b2 (iflag) * rs
     dox = 1.d0 + 7.d0 / 6.d0 * b1 (iflag) * rs12 + 4.d0 / 3.d0 * &
          b2 (iflag) * rs
     ec = gc (iflag) / ox
     vc = ec * dox / ox
  endif
  !
  return
end subroutine pz
!
!-----------------------------------------------------------------------
subroutine pzKZK (rs, ec, vc, vol)
  !-----------------------------------------------------------------------
  !     LDA parameterization from Monte Carlo data
  !     iflag=1: J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
  !     iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
  !
  implicit none
  real(8) :: rs, ec, vc, ec0 (2), vc0(2), ec0p 
  integer :: iflag, kr
  !
  real(8) :: a (2), b (2), c (2), d (2), gc (2), b1 (2), b2 (2)
  real(8) :: lnrs, rs12, ox, dox, lnrsk, rsk
  real(8) :: a1, grs, g1, g2, g3, g4, dL, vol, gh, gl, grsp
  real(8) :: f3, f2, f1, f0, pi
  real(8) :: D1, D2, D3, P1, P2, ry2h
  !
  data a / 0.0311d0, 0.031091d0 /, b / -0.048d0, -0.046644d0 /, &
       c / 0.0020d0, 0.00419d0 /, d / -0.0116d0, -0.00983d0 /
  data gc / -0.1423d0, -0.103756d0 /, b1 / 1.0529d0, 0.56371d0 /, &
       b2 / 0.3334d0, 0.27358d0 /
  data a1 / -2.2037 /, g1 / 0.1182 /, g2 / 1.1656 /, g3 / -5.2884 /, &
       g4 / -1.1233 /
  data ry2h / 0.5d0 /
  !
  iflag = 1
  pi = 4.d0 * atan(1.d0)
  dL = vol**(1.d0/3.d0)
  gh = 0.5d0 * dL / (2.d0 * pi)**(1.d0/3.d0)
  gl = dL * (3.d0 / 2.d0 / pi)**(1.d0/3.d0)

  rsk = gh
  do kr = 1, 2
  lnrsk = log (rsk)
  if (rsk.lt.1.0d0) then
     ! high density formula
    ec0(kr) = a(iflag) *lnrsk + b(iflag) + c(iflag) * rsk * lnrsk + d( &
          iflag) * rsk
    vc0(kr) = a(iflag) * lnrsk + (b(iflag) - a(iflag) / 3.d0) + 2.d0 / &
      3.d0 * c (iflag) * rsk * lnrsk + (2.d0 * d (iflag) - c (iflag) ) &
          / 3.d0 * rsk
  else
     ! interpolation formula
     rs12 = sqrt (rsk)
     ox = 1.d0 + b1 (iflag) * rs12 + b2 (iflag) * rsk
     dox = 1.d0 + 7.d0 / 6.d0 * b1 (iflag) * rs12 + 4.d0 / 3.d0 * &
          b2 (iflag) * rsk
     ec0(kr) = gc (iflag) / ox
     vc0(kr) = ec0(kr) * dox / ox
  endif
  !
   grs  = g1 * rsk * lnrsk + g2 * rsk + g3 * rsk**1.5d0 + g4 * rsk**2.d0
   grsp = g1 * lnrsk + g1 + g2 + 1.5d0 * g3 * rsk**0.5d0 + & 
        2.d0 * g4 * rsk
   ec0(kr)  = ec0(kr) + (-a1 * rsk / dL**2.d0 + grs / dL**3.d0) * ry2h
   vc0(kr)  = vc0(kr) + (-2.d0 * a1 * rsk / dL**2.d0 / 3.d0 + &
           grs / dL**3.d0 -  grsp * rsk / 3.d0 / dL**3.d0) * ry2h
  !
  rsk = rs
  enddo

  lnrs = log (rs)
  if (rs .le. gh) then
   ec = ec0(2)
   vc = vc0(2)
  else
     if ( rs .le. gl) then
        ec0p = 3.d0 * (ec0(1) - vc0(1)) / gh
        P1 = 3.d0 *  ec0(1) - gh * ec0p
        P2 = ec0p
        D1 = gl - gh
        D2 = gl**2.d0 - gh**2.d0
        D3 = gl**3.d0 - gh**3.d0
        f2 = 2.d0 * gl**2.d0 * P2 * D1 + D2 * P1
        f2 = f2 / (-(2.d0*gl*D1)**2.d0 + 4.d0*gl*D1*D2 - D2**2.d0 )
        f3 = - (P2 + 2.d0*D1*f2) / (3.d0 * D2)
        f1 = - (P1 + D2 * f2) / (2.d0 * D1)
        f0 = - gl * (gl * f2 + 2.d0 * f1) / 3.d0
        !
        ec = f3 * rs**3.d0 + f2 * rs**2.d0 + f1 * rs + f0
        vc = f2 * rs**2.d0 / 3.d0 + f1 * 2.d0 * rs / 3.d0 + f0
     else
        ec = 0.d0
        vc = 0.d0
     endif
    endif
  !
  return
end subroutine pzKZK
!
!-----------------------------------------------------------------------
subroutine vwn (rs, ec, vc)
  !-----------------------------------------------------------------------
  !     S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980)
  !
  implicit none
  real(8) :: rs, ec, vc
  real(8) :: a, b, c, x0
  parameter (a = 0.0310907d0, b = 3.72744d0, c = 12.9352d0, x0 = -0.10498d0)
  real(8) :: q, f1, f2, f3, rs12, fx, qx, tx, tt
  !
  q = sqrt (4.d0 * c - b * b)
  f1 = 2.d0 * b / q
  f2 = b * x0 / (x0 * x0 + b * x0 + c)
  f3 = 2.d0 * (2.d0 * x0 + b) / q
  rs12 = sqrt (rs)
  fx = rs + b * rs12 + c
  qx = atan (q / (2.d0 * rs12 + b) )
  ec = a * (log (rs / fx) + f1 * qx - f2 * (log ( (rs12 - x0) **2 / &
       fx) + f3 * qx) )
  tx = 2.d0 * rs12 + b
  tt = tx * tx + q * q
  vc = ec - rs12 * a / 6.d0 * (2.d0 / rs12 - tx / fx - 4.d0 * b / &
       tt - f2 * (2.d0 / (rs12 - x0) - tx / fx - 4.d0 * (2.d0 * x0 + b) &
       / tt) )
  !
  return
end subroutine vwn
!-----------------------------------------------------------------------
subroutine lyp (rs, ec, vc)
  !-----------------------------------------------------------------------
  !     C. Lee, W. Yang, and R.G. Parr, PRB 37, 785 (1988)
  !     LDA part only
  !
  implicit none
  real(8) :: rs, ec, vc
  real(8) :: a, b, c, d, pi43
  parameter (a = 0.04918d0, b = 0.132d0 * 2.87123400018819108d0)
  ! pi43 = (4pi/3)^(1/3)
  parameter (pi43 = 1.61199195401647d0, c = 0.2533d0 * pi43, d = &
       0.349d0 * pi43)
  real(8) :: ecrs, ox
  !
  ecrs = b * exp ( - c * rs)
  ox = 1.d0 / (1.d0 + d * rs)
  ec = - a * ox * (1.d0 + ecrs)
  vc = ec - rs / 3.d0 * a * ox * (d * ox + ecrs * (d * ox + c) )
  !
  return
end subroutine lyp
!
!-----------------------------------------------------------------------
subroutine pw (rs, iflag, ec, vc)
  !-----------------------------------------------------------------------
  !     iflag=1: J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
  !     iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
  !
  implicit none
  real(8) :: rs, ec, vc
  integer :: iflag
  !
  real(8) :: a, b1, b2, c0, c1, c2, c3, d0, d1
  parameter (a = 0.031091d0, b1 = 7.5957d0, b2 = 3.5876d0, c0 = a, &
       c1 = 0.046644d0, c2 = 0.00664d0, c3 = 0.01043d0, d0 = 0.4335d0, &
       d1 = 1.4408d0)
  real(8) :: lnrs, rs12, rs32, rs2, om, dom, olog
  real(8) :: a1 (2), b3 (2), b4 (2)
  data a1 / 0.21370d0, 0.026481d0 /, b3 / 1.6382d0, -0.46647d0 /, &
       b4 / 0.49294d0, 0.13354d0 /
  !
  ! high- and low-density formulae implemented but not used in PW case
  ! (reason: inconsistencies in PBE/PW91 functionals)
  !
  if (rs.lt.1d0.and.iflag.eq.2) then
     ! high density formula
     lnrs = log (rs)
     ec = c0 * lnrs - c1 + c2 * rs * lnrs - c3 * rs
     vc = c0 * lnrs - (c1 + c0 / 3.d0) + 2.d0 / 3.d0 * c2 * rs * &
          lnrs - (2.d0 * c3 + c2) / 3.d0 * rs
  elseif (rs.gt.100.d0.and.iflag.eq.2) then
     ! low density formula
     ec = - d0 / rs + d1 / rs**1.5d0
     vc = - 4.d0 / 3.d0 * d0 / rs + 1.5d0 * d1 / rs**1.5d0
  else
     ! interpolation formula
     rs12 = sqrt (rs)
     rs32 = rs * rs12
     rs2 = rs**2
     om = 2.d0 * a * (b1 * rs12 + b2 * rs + b3 (iflag) * rs32 + b4 ( &
          iflag) * rs2)
     dom = 2.d0 * a * (0.5d0 * b1 * rs12 + b2 * rs + 1.5d0 * b3 ( &
          iflag) * rs32 + 2.d0 * b4 (iflag) * rs2)
     olog = log (1.d0 + 1.0d0 / om)
     ec = - 2.d0 * a * (1.d0 + a1 (iflag) * rs) * olog
     vc = - 2.d0 * a * (1.d0 + 2.d0 / 3.d0 * a1 (iflag) * rs) &
          * olog - 2.d0 / 3.d0 * a * (1.d0 + a1 (iflag) * rs) * dom / &
          (om * (om + 1.d0) )
  endif
  !
  return
end subroutine pw
!
!-----------------------------------------------------------------------
subroutine wigner (rs, ec, vc)
  !-----------------------------------------------------------------------
  !        Wigner correlation
  !
  implicit none
  real(8) :: rs, ec, vc
  real(8) :: pi34, rho13
  parameter (pi34 = 0.6203504908994d0)
  ! pi34=(3/4pi)^(1/3), rho13=rho^(1/3)
  !
  rho13 = pi34 / rs
  vc = - rho13 * ( (0.943656d0 + 8.8963d0 * rho13) / (1.d0 + &
       12.57d0 * rho13) **2)
  ec = - 0.738d0 * rho13 * (0.959d0 / (1.d0 + 12.57d0 * rho13) )
  !
  return
end subroutine wigner
!
!-----------------------------------------------------------------------
subroutine hl (rs, ec, vc)
  !-----------------------------------------------------------------------
  !     L. Hedin and  B.I. Lundqvist,  J. Phys. C 4, 2064 (1971)
  !
  implicit none
  real(8) :: rs, ec, vc
  real(8) :: a, x
  !
  a = log (1.0d0 + 21.d0 / rs)
  x = rs / 21.0d0
  ec = a + (x**3 * a - x * x) + x / 2.d0 - 1.0d0 / 3.0d0
  ec = - 0.0225d0 * ec
  vc = - 0.0225d0 * a
  !
  return
end subroutine hl
!
!-----------------------------------------------------------------------
subroutine gl (rs, ec, vc)
  !-----------------------------------------------------------------------
  !  O. Gunnarsson and B. I. Lundqvist, PRB 13, 4274 (1976)
  !
  implicit none
  real(8) :: rs, vc, ec
  real(8) :: c, r, x
  parameter (c = 0.0333d0, r = 11.4d0)
  ! c=0.0203, r=15.9 for the paramagnetic case
  !
  x = rs / r
  vc = - c * log (1.d0 + 1.d0 / x)
  ec = - c * ( (1.d0 + x**3) * log (1.d0 + 1.d0 / x) - 1.0d0 / &
       3.0d0 + x * (0.5d0 - x) )
  !
  return
end subroutine gl
!
!-----------------------------------------------------------------------
subroutine becke88 (rho, grho, sx, v1x, v2x)
  !-----------------------------------------------------------------------
  ! Becke exchange: A.D. Becke, PRA 38, 3098 (1988)
  ! only gradient-corrected part, no Slater term included
  !
  implicit none
  real(8) :: rho, grho, sx, v1x, v2x
  real(8) :: beta, third, two13
  parameter (beta = 0.0042d0)
  parameter (third = 1.d0 / 3.d0, two13 = 1.259921049894873d0)
  ! two13 = 2^(1/3)
  real(8) :: rho13, rho43, xs, xs2, sa2b8, shm1, dd, dd2, ee
  !
  rho13 = rho**third
  rho43 = rho13**4
  xs = two13 * sqrt (grho) / rho43
  xs2 = xs * xs
  sa2b8 = sqrt (1.0d0 + xs2)
  shm1 = log (xs + sa2b8)
  dd = 1.0d0 + 6.0d0 * beta * xs * shm1
  dd2 = dd * dd
  ee = 6.0d0 * beta * xs2 / sa2b8 - 1.d0
  sx = two13 * grho / rho43 * ( - beta / dd)
  v1x = - (4.d0 / 3.d0) / two13 * xs2 * beta * rho13 * ee / dd2
  v2x = two13 * beta * (ee-dd) / (rho43 * dd2)
  !
  return
end subroutine becke88
!
!-----------------------------------------------------------------------
subroutine ggax (rho, grho, sx, v1x, v2x)
  !-----------------------------------------------------------------------
  ! Perdew-Wang GGA (PW91), exchange part:
  ! J.P. Perdew et al.,PRB 46, 6671 (1992)
  !
  implicit none
  real(8) :: rho, grho, sx, v1x, v2x
  real(8) :: f1, f2, f3, f4, f5
  parameter (f1 = 0.19645d0, f2 = 7.7956d0, f3 = 0.2743d0, f4 = &
       0.1508d0, f5 = 0.004d0)
  real(8) :: fp1, fp2
  parameter (fp1 = -0.019292021296426d0, fp2 = 0.161620459673995d0)
  ! fp1 = -3/(16 pi)*(3 pi^2)^(-1/3)
  ! fp2 = (1/2)(3 pi^2)**(-1/3)
  real(8) :: rhom43, s, s2, s3, s4, exps, as, sa2b8, shm1, bs, das, &
       dbs, dls
  !
  rhom43 = rho** ( - 4.d0 / 3.d0)
  s = fp2 * sqrt (grho) * rhom43
  s2 = s * s
  s3 = s2 * s
  s4 = s2 * s2
  exps = f4 * exp ( - 100.d0 * s2)
  as = f3 - exps - f5 * s2
  sa2b8 = sqrt (1.0d0 + f2 * f2 * s2)
  shm1 = log (f2 * s + sa2b8)
  bs = 1.d0 + f1 * s * shm1 + f5 * s4
  das = (200.d0 * exps - 2.d0 * f5) * s
  dbs = f1 * (shm1 + f2 * s / sa2b8) + 4.d0 * f5 * s3
  dls = (das / as - dbs / bs)
  sx = fp1 * grho * rhom43 * as / bs
  v1x = - 4.d0 / 3.d0 * sx / rho * (1.d0 + s * dls)
  v2x = fp1 * rhom43 * as / bs * (2.d0 + s * dls)
  !
  return
end subroutine ggax

!
!-----------------------------------------------------------------------
subroutine rPW86 (rho, grho, sx, v1x, v2x)
  !-----------------------------------------------------------------------
  ! PRB 33, 8800 (1986) and J. Chem. Theory comp. 5, 2754 (2009)
  !
  implicit none

  real(8), intent(in) :: rho, grho
  real(8), intent(out) :: sx, v1x, v2x
  real(8) :: s, s_2, s_3, s_4, s_5, s_6, fs, grad_rho, df_ds
  real(8) :: a, b, c, s_prefactor, Ax, four_thirds
  parameter( a = 1.851d0, b = 17.33d0, c = 0.163d0, s_prefactor = 6.18733545256027d0, &
       Ax = -0.738558766382022d0, four_thirds = 4.d0/3.d0)

  grad_rho = sqrt(grho)

  s = grad_rho/(s_prefactor*rho**(four_thirds))

  s_2 = s**2
  s_3 = s_2 * s
  s_4 = s_2**2
  s_5 = s_3 * s_2
  s_6 = s_2 * s_4

  !! Calculation of energy
  fs = (1 + a*s_2 + b*s_4 + c*s_6)**(1.d0/15.d0)
  sx = Ax * rho**(four_thirds) * fs


  !! Calculation of the potential
  df_ds = (1.d0/(15.d0*fs**(14)))*(2*a*s + 4*b*s_3 + 6*c*s_5)


  v1x = Ax*(four_thirds)*(rho**(1.d0/3.d0)*fs &
       -grad_rho/(s_prefactor * rho)*df_ds)

  v2x = Ax * df_ds/(s_prefactor*grad_rho)

end subroutine rPW86

!
!---------------------------------------------------------------
subroutine c09x (rho, grho, sx, v1x, v2x)
  !---------------------------------------------------------------
  ! Cooper '09 exchange for vdW-DF (without Slater exchange):
  ! V. R. Cooper, Phys. Rev. B 81, 161104(R) (2010)
  !
  ! Developed thanks to the contribution of 
  ! Ikutaro Hamada - ikutaro@wpi-aimr.tohoku.ac.jp
  ! WPI-Advanced Institute of Materials Research, Tohoku University  
  !
  implicit none
  REAL(8), PARAMETER :: pi     = 3.14159265358979323846_8 
  REAL(8), PARAMETER :: tpi    = 2.0_8 * pi
  REAL(8), PARAMETER :: HARTREE_SI       = 4.35974394E-18_8   ! J
  REAL(8), PARAMETER :: C_SI             = 2.99792458E+8_8    ! m sec^-1
  REAL(8), PARAMETER :: BOHR_RADIUS_SI   = 0.52917720859E-10_8 ! m
  REAL(8), PARAMETER :: H_PLANCK_SI      = 6.62606896E-34_8   ! J s
  REAL(8), PARAMETER :: AU_SEC           = H_PLANCK_SI/tpi/HARTREE_SI
  REAL(8), PARAMETER :: C_AU             = C_SI / BOHR_RADIUS_SI * AU_SEC
  real(8) :: rho, grho, sx, v1x, v2x
  ! input: charge and squared gradient
  ! output: energy
  ! output: potential
  ! local variables
  real(8) :: kf, agrho, s1, s2, ds, dsg, exunif, fx
  ! (3*pi2*|rho|)^(1/3)
  ! |grho|
  ! |grho|/(2*kf*|rho|)
  ! s^2
  ! n*ds/dn
  ! n*ds/d(gn)
  ! exchange energy LDA part
  ! exchange energy gradient part
  real(8) :: dxunif, dfx, f1, f2, f3, dfx1, dfx2
  ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
  real(8) :: third, c1, c2, c5
  parameter (third = 1.d0 / 3.d0, c1 = 0.75d0 / pi , &
       c2 = 3.093667726280136d0, c5 = 4.d0 * third)
  ! parameters of the functional
  real(8) :: kappa, mu, alpha
  data kappa / 1.245d0  /, &
       mu  / 0.0617d0 /, &
       alpha / 0.0483d0 /
  !
  agrho = sqrt (grho)
  kf = c2 * rho**third
  dsg = 0.5d0 / kf
  s1 = agrho * dsg / rho
  s2 = s1 * s1
  ds = - c5 * s1
  !
  !   Energy
  !
  f1 = exp( - alpha * s2 )
  f2 = exp( - alpha * s2 / 2.0d0 )
  f3 = mu * s2 * f1
  fx = f3 + kappa * ( 1.0d0 - f2 )
  exunif = - c1 * kf
  sx = exunif * fx
  !
  !   Potential
  !
  dxunif = exunif * third
  dfx1 = 2.0d0 * mu * s1 * ( 1.0d0 - alpha * s2 ) * f1
  dfx2 = kappa * alpha * s1 * f2 
  dfx = dfx1 + dfx2 
  v1x = sx + dxunif * fx + exunif * dfx * ds
  v2x = exunif * dfx * dsg / agrho

  sx = sx * rho
  return
end subroutine c09x

!
!-----------------------------------------------------------------------
subroutine perdew86 (rho, grho, sc, v1c, v2c)
  !-----------------------------------------------------------------------
  ! Perdew gradient correction on correlation: PRB 33, 8822 (1986)
  !
  implicit none
  real(8) :: rho, grho, sc, v1c, v2c
  real(8) :: p1, p2, p3, p4, pc1, pc2, pci
  parameter (p1 = 0.023266d0, p2 = 7.389d-6, p3 = 8.723d0, p4 = &
       0.472d0)
  parameter (pc1 = 0.001667d0, pc2 = 0.002568d0, pci = pc1 + pc2)
  real(8) :: third, pi34
  parameter (third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0)
  ! pi34=(3/4pi)^(1/3)
  real(8) :: rho13, rho43, rs, rs2, rs3, cna, cnb, cn, drs
  real(8) :: dcna, dcnb, dcn, phi, ephi
  !
  rho13 = rho**third
  rho43 = rho13**4
  rs = pi34 / rho13
  rs2 = rs * rs
  rs3 = rs * rs2
  cna = pc2 + p1 * rs + p2 * rs2
  cnb = 1.d0 + p3 * rs + p4 * rs2 + 1.d4 * p2 * rs3
  cn = pc1 + cna / cnb
  drs = - third * pi34 / rho43
  dcna = (p1 + 2.d0 * p2 * rs) * drs
  dcnb = (p3 + 2.d0 * p4 * rs + 3.d4 * p2 * rs2) * drs
  dcn = dcna / cnb - cna / (cnb * cnb) * dcnb
  phi = 0.192d0 * pci / cn * sqrt (grho) * rho** ( - 7.d0 / 6.d0)
  ! SdG: in the original paper 1.745*0.11=0.19195 is used
  ephi = exp ( - phi)
  sc = grho / rho43 * cn * ephi
  v1c = sc * ( (1.d0 + phi) * dcn / cn - ( (4.d0 / 3.d0) - (7.d0 / &
       6.d0) * phi) / rho)
  v2c = cn * ephi / rho43 * (2.d0 - phi)
  !
  return
end subroutine perdew86
!
!-----------------------------------------------------------------------
subroutine glyp (rho, grho, sc, v1c, v2c)
  !-----------------------------------------------------------------------
  ! Lee Yang Parr: gradient correction part
  !
  implicit none
  real(8) :: rho, grho, sc, v1c, v2c
  real(8) :: a, b, c, d
  parameter (a = 0.04918d0, b = 0.132d0, c = 0.2533d0, d = 0.349d0)
  real(8) :: rhom13, rhom43, rhom53, om, xl, ff, dom, dxl
  !
  rhom13 = rho** ( - 1.d0 / 3.d0)
  om = exp ( - c * rhom13) / (1.d0 + d * rhom13)
  xl = 1.d0 + (7.d0 / 3.d0) * (c * rhom13 + d * rhom13 / (1.d0 + d * &
       rhom13) )
  ff = a * b * grho / 24.d0
  rhom53 = rhom13**5
  sc = ff * rhom53 * om * xl
  dom = - om * (c + d+c * d * rhom13) / (1.d0 + d * rhom13)
  dxl = (7.d0 / 3.d0) * (c + d+2.d0 * c * d * rhom13 + c * d * d * &
       rhom13**2) / (1.d0 + d * rhom13) **2
  rhom43 = rhom13**4
  v1c = - ff * rhom43 / 3.d0 * (5.d0 * rhom43 * om * xl + rhom53 * &
       dom * xl + rhom53 * om * dxl)
  v2c = 2.d0 * sc / grho
  !
  return
end subroutine glyp
!
!-----------------------------------------------------------------------
subroutine ggac (rho, grho, sc, v1c, v2c)
  !-----------------------------------------------------------------------
  ! Perdew-Wang GGA (PW91) correlation part
  !
  implicit none
  real(8) :: rho, grho, sc, v1c, v2c
  real(8) :: al, pa, pb, pc, pd, cx, cxc0, cc0
  parameter (al = 0.09d0, pa = 0.023266d0, pb = 7.389d-6, pc = &
       8.723d0, pd = 0.472d0)
  parameter (cx = -0.001667d0, cxc0 = 0.002568d0, cc0 = - cx + cxc0)
  real(8) :: third, pi34, nu, be, xkf, xks
  parameter (third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0)
  parameter (nu = 15.755920349483144d0, be = nu * cc0)
  parameter (xkf = 1.919158292677513d0, xks = 1.128379167095513d0)
  ! pi34=(3/4pi)^(1/3),  nu=(16/pi)*(3 pi^2)^(1/3)
  ! xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
  real(8) :: kf, ks, rs, rs2, rs3, ec, vc, t, expe, af, bf, y, xy, &
       qy, s1
  real(8) :: h0, dh0, ddh0, ee, cn, dcn, cna, dcna, cnb, dcnb, h1, &
       dh1, ddh1
  !
  rs = pi34 / rho**third
  rs2 = rs * rs
  rs3 = rs * rs2
  call pw (rs, 1, ec, vc)
  kf = xkf / rs
  ks = xks * sqrt (kf)
  t = sqrt (grho) / (2.d0 * ks * rho)
  expe = exp ( - 2.d0 * al * ec / (be * be) )
  af = 2.d0 * al / be * (1.d0 / (expe-1.d0) )
  bf = expe * (vc - ec)
  y = af * t * t
  xy = (1.d0 + y) / (1.d0 + y + y * y)
  qy = y * y * (2.d0 + y) / (1.d0 + y + y * y) **2
  s1 = 1.d0 + 2.d0 * al / be * t * t * xy
  h0 = be * be / (2.d0 * al) * log (s1)
  dh0 = be * t * t / s1 * ( - 7.d0 / 3.d0 * xy - qy * (af * bf / &
       be-7.d0 / 3.d0) )
  ddh0 = be / (2.d0 * ks * ks * rho) * (xy - qy) / s1
  ee = - 100.d0 * (ks / kf * t) **2
  cna = cxc0 + pa * rs + pb * rs2
  dcna = pa * rs + 2.d0 * pb * rs2
  cnb = 1.d0 + pc * rs + pd * rs2 + 1.d4 * pb * rs3
  dcnb = pc * rs + 2.d0 * pd * rs2 + 3.d4 * pb * rs3
  cn = cna / cnb - cx
  dcn = dcna / cnb - cna * dcnb / (cnb * cnb)
  h1 = nu * (cn - cc0 - 3.d0 / 7.d0 * cx) * t * t * exp (ee)
  dh1 = - third * (h1 * (7.d0 + 8.d0 * ee) + nu * t * t * exp (ee) &
       * dcn)
  ddh1 = 2.d0 * h1 * (1.d0 + ee) * rho / grho
  sc = rho * (h0 + h1)
  v1c = h0 + h1 + dh0 + dh1
  v2c = ddh0 + ddh1
  !
  return
end subroutine ggac
!
!---------------------------------------------------------------
subroutine pbex (rho, grho, iflag, sx, v1x, v2x)
  !---------------------------------------------------------------
  !
  ! PBE exchange (without Slater exchange):
  ! iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
  ! iflag=2  "revised' PBE: Y. Zhang et al., PRL 80, 890 (1998)
  ! iflag=3  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
  !
  implicit none
  REAL(8), PARAMETER :: pi     = 3.14159265358979323846_8 
  REAL(8), PARAMETER :: tpi    = 2.0_8 * pi
  REAL(8), PARAMETER :: HARTREE_SI       = 4.35974394E-18_8   ! J
  REAL(8), PARAMETER :: C_SI             = 2.99792458E+8_8    ! m sec^-1
  REAL(8), PARAMETER :: BOHR_RADIUS_SI   = 0.52917720859E-10_8 ! m
  REAL(8), PARAMETER :: H_PLANCK_SI      = 6.62606896E-34_8   ! J s
  REAL(8), PARAMETER :: AU_SEC           = H_PLANCK_SI/tpi/HARTREE_SI
  REAL(8), PARAMETER :: C_AU             = C_SI / BOHR_RADIUS_SI * AU_SEC
  real(8) :: rho, grho, sx, v1x, v2x
  ! input: charge and squared gradient
  ! output: energy
  ! output: potential
  integer :: iflag
  ! local variables
  real(8) :: kf, agrho, s1, s2, ds, dsg, exunif, fx
  ! (3*pi2*|rho|)^(1/3)
  ! |grho|
  ! |grho|/(2*kf*|rho|)
  ! s^2
  ! n*ds/dn
  ! n*ds/d(gn)
  ! exchange energy LDA part
  ! exchange energy gradient part
  real(8) :: dxunif, dfx, f1, f2, f3, dfx1
  ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
  real(8) :: third, c1, c2, c5
  parameter (third = 1.d0 / 3.d0, c1 = 0.75d0 / pi , &
       c2 = 3.093667726280136d0, c5 = 4.d0 * third)
  ! parameters of the functional
  real(8) :: k (3), mu(3)
  data k / 0.804d0, 1.2450D0, 0.804d0 /, &
       mu/ 0.21951d0, 0.21951d0, 0.12345679012345679012d0  /
  !
  agrho = sqrt (grho)
  kf = c2 * rho**third
  dsg = 0.5d0 / kf
  s1 = agrho * dsg / rho
  s2 = s1 * s1
  ds = - c5 * s1
  !
  !   Energy
  !
  f1 = s2 * mu(iflag) / k (iflag)
  f2 = 1.d0 + f1
  f3 = k (iflag) / f2
  fx = k (iflag) - f3
  exunif = - c1 * kf
  sx = exunif * fx
  !
  !   Potential
  !
  dxunif = exunif * third
  dfx1 = f2 * f2
  dfx = 2.d0 * mu(iflag) * s1 / dfx1
  v1x = sx + dxunif * fx + exunif * dfx * ds
  v2x = exunif * dfx * dsg / agrho

  sx = sx * rho
  return
end subroutine pbex
!
!---------------------------------------------------------------
subroutine pbex_vec (rho, grho, iflag, sx, v1x, v2x, length, small)
  !---------------------------------------------------------------
  !
  ! PBE exchange (without Slater exchange):
  ! iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
  ! iflag=2  "revised' PBE: Y. Zhang et al., PRL 80, 890 (1998)
  ! iflag=3  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
  !
  implicit none
  REAL(8), PARAMETER :: pi     = 3.14159265358979323846_8 
  REAL(8), PARAMETER :: tpi    = 2.0_8 * pi
  REAL(8), PARAMETER :: HARTREE_SI       = 4.35974394E-18_8   ! J
  REAL(8), PARAMETER :: C_SI             = 2.99792458E+8_8    ! m sec^-1
  REAL(8), PARAMETER :: BOHR_RADIUS_SI   = 0.52917720859E-10_8 ! m
  REAL(8), PARAMETER :: H_PLANCK_SI      = 6.62606896E-34_8   ! J s
  REAL(8), PARAMETER :: AU_SEC           = H_PLANCK_SI/tpi/HARTREE_SI
  REAL(8), PARAMETER :: C_AU             = C_SI / BOHR_RADIUS_SI * AU_SEC
  integer, intent(in)   :: length
  integer, intent(in)   :: iflag
  real(8), intent(in)  :: small
  real(8), intent(in)  :: rho(length), grho(length)
  real(8), intent(out) :: sx(length), v1x(length), v2x(length)
  ! local variables
  integer :: i
  real(8) :: kf, agrho, s1, dsg, exunif, fx
  ! (3*pi2*|rho|)^(1/3)
  ! |grho|
  ! |grho|/(2*kf*|rho|)
  ! n*ds/d(gn)
  ! exchange energy LDA part
  ! exchange energy gradient part
  real(8) :: dfx, f1, f2
  ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
  real(8) :: third, c1, c2, c5
  parameter (third = 1.0_8 / 3.0_8, c1 = 0.75_8 / pi , &
       c2 = 3.093667726280136_8, c5 = 4.0_8 * third)
  ! parameters of the functional
  real(8) :: k (3), mu(3)
  data k / 0.804_8, 1.245_8, 0.804_8 /, &
       mu/ 0.21951_8, 0.21951_8, 0.12345679012345679012_8  /
  !
  do i=1,length
     if ((rho(i).gt.small).and.(grho(i).gt.small**2)) then
        agrho = sqrt(grho(i))
        kf = c2 * rho(i)**third
        dsg = 0.5_8 / kf
        s1 = agrho * dsg / rho(i)
        !
        !   Energy
        f1 = s1*s1 * mu(iflag) / k(iflag)
        f2 = 1.0_8 / (1.0_8 + f1)
        fx = k(iflag) * (1.0_8 - f2)
        exunif = - c1 * kf
        sx(i) = exunif * fx
        !
        !   Potential
        dfx = 2.0_8 * mu(iflag) * s1 *f2*f2

        v1x(i) = sx(i) + exunif * (third * fx - c5 * dfx * s1)
        v2x(i) = exunif * dfx * dsg / agrho
        sx(i) = sx(i) * rho(i)

     else
        v1x(i) = 0.0_8
        v2x(i) = 0.0_8
        sx(i) = 0.0_8
     end if
  end do

end subroutine pbex_vec
!
!---------------------------------------------------------------
subroutine pbec (rho, grho, iflag, sc, v1c, v2c)
  !---------------------------------------------------------------
  !
  ! PBE correlation (without LDA part)
  ! iflag=1: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
  ! iflag=2: J.P.Perdew et al., PRL 100, 136406 (2008).
  !
  implicit none
  integer, intent(in) :: iflag
  real(8) :: rho, grho, sc, v1c, v2c
  real(8) :: ga, be (2)
  parameter (ga = 0.031091d0)
  data be / 0.066725d0, 0.046d0 /
  real(8) :: third, pi34, xkf, xks
  parameter (third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0)
  parameter (xkf = 1.919158292677513d0, xks = 1.128379167095513d0)
  ! pi34=(3/4pi)^(1/3), xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
  real(8) :: kf, ks, rs, ec, vc, t, expe, af, bf, y, xy, qy
  real(8) :: s1, h0, dh0, ddh0
  !
  rs = pi34 / rho**third
  call pw (rs, 1, ec, vc)
  kf = xkf / rs
  ks = xks * sqrt (kf)
  t = sqrt (grho) / (2.d0 * ks * rho)
  expe = exp ( - ec / ga)
  af = be(iflag) / ga * (1.d0 / (expe-1.d0) )
  bf = expe * (vc - ec)
  y = af * t * t
  xy = (1.d0 + y) / (1.d0 + y + y * y)
  qy = y * y * (2.d0 + y) / (1.d0 + y + y * y) **2
  s1 = 1.d0 + be(iflag) / ga * t * t * xy
  h0 = ga * log (s1)
  dh0 = be(iflag) * t * t / s1 * ( - 7.d0 / 3.d0 * xy - qy * (af * bf / &
       be(iflag)-7.d0 / 3.d0) )
  ddh0 = be(iflag) / (2.d0 * ks * ks * rho) * (xy - qy) / s1
  sc = rho * h0
  v1c = h0 + dh0
  v2c = ddh0
  !
  return
end subroutine pbec

!     ==================================================================
subroutine hcth(rho,grho,sx,v1x,v2x)
  !     ==================================================================
  !     HCTH/120, JCP 109, p. 6264 (1998)
  !     Parameters set-up after N.L. Doltsisnis & M. Sprik (1999)
  !     Present release: Mauro Boero, Tsukuba, 11/05/2004
  !--------------------------------------------------------------------------
  !     rhoa = rhob = 0.5 * rho
  !     grho is the SQUARE of the gradient of rho! --> gr=sqrt(grho)
  !     sx  : total exchange correlation energy at point r
  !     v1x : d(sx)/drho  (eq. dfdra = dfdrb in original)
  !     v2x : 1/gr*d(sx)/d(gr) (eq. 0.5 * dfdza = 0.5 * dfdzb in original)
  !--------------------------------------------------------------------------
  implicit none
  REAL(8), PARAMETER :: pi     = 3.14159265358979323846_8 
  REAL(8), PARAMETER :: tpi    = 2.0_8 * pi
  REAL(8), PARAMETER :: HARTREE_SI       = 4.35974394E-18_8   ! J
  REAL(8), PARAMETER :: C_SI             = 2.99792458E+8_8    ! m sec^-1
  REAL(8), PARAMETER :: BOHR_RADIUS_SI   = 0.52917720859E-10_8 ! m
  REAL(8), PARAMETER :: H_PLANCK_SI      = 6.62606896E-34_8   ! J s
  REAL(8), PARAMETER :: AU_SEC           = H_PLANCK_SI/tpi/HARTREE_SI
  REAL(8), PARAMETER :: C_AU             = C_SI / BOHR_RADIUS_SI * AU_SEC
  real(8) :: rho, grho, sx, v1x, v2x

  real(8), parameter :: o3=1.0d0/3.0d0, o34=4.0d0/3.0d0, fr83=8.d0/3.d0
  real(8) :: cg0(6), cg1(6), caa(6), cab(6), cx(6)
  real(8) :: r3q2, r3pi, gr, rho_o3, rho_o34, xa, xa2, ra, rab, &
       dra_drho, drab_drho, g, dg, era1, dera1_dra, erab0, derab0_drab, &
       ex, dex_drho, uaa, uab, ux, ffaa, ffab,  dffaa_drho, dffab_drho,&
       denaa, denab, denx, f83rho, bygr, gaa, gab, gx, taa, tab, txx, &
       dgaa_drho, dgab_drho, dgx_drho, dgaa_dgr, dgab_dgr, dgx_dgr
  !
  r3q2=2.d0**(-o3)
  r3pi=(3.d0/pi)**o3
  !.....coefficients for pw correlation......................................
  cg0(1)= 0.031091d0
  cg0(2)= 0.213700d0
  cg0(3)= 7.595700d0
  cg0(4)= 3.587600d0
  cg0(5)= 1.638200d0
  cg0(6)= 0.492940d0
  cg1(1)= 0.015545d0
  cg1(2)= 0.205480d0
  cg1(3)=14.118900d0
  cg1(4)= 6.197700d0
  cg1(5)= 3.366200d0
  cg1(6)= 0.625170d0
  !......hcth-19-4.....................................
  caa(1)=  0.489508d+00
  caa(2)= -0.260699d+00
  caa(3)=  0.432917d+00
  caa(4)= -0.199247d+01
  caa(5)=  0.248531d+01
  caa(6)=  0.200000d+00
  cab(1)=  0.514730d+00
  cab(2)=  0.692982d+01
  cab(3)= -0.247073d+02
  cab(4)=  0.231098d+02
  cab(5)= -0.113234d+02
  cab(6)=  0.006000d+00
  cx(1) =  0.109163d+01
  cx(2) = -0.747215d+00
  cx(3) =  0.507833d+01
  cx(4) = -0.410746d+01
  cx(5) =  0.117173d+01
  cx(6)=   0.004000d+00
  !...........................................................................
  gr=DSQRT(grho)
  rho_o3=rho**(o3)
  rho_o34=rho**(o34)
  xa=1.25992105d0*gr/rho_o34
  xa2=xa*xa
  ra=0.781592642d0/rho_o3
  rab=r3q2*ra
  dra_drho=-0.260530881d0/rho_o34
  drab_drho=r3q2*dra_drho
  call pwcorr(ra,cg1,g,dg)
  era1=g
  dera1_dra=dg
  call pwcorr(rab,cg0,g,dg)
  erab0=g
  derab0_drab=dg
  ex=-0.75d0*r3pi*rho_o34
  dex_drho=-r3pi*rho_o3
  uaa=caa(6)*xa2
  uaa=uaa/(1.0d0+uaa)
  uab=cab(6)*xa2
  uab=uab/(1.0d0+uab)
  ux=cx(6)*xa2
  ux=ux/(1.0d0+ux)
  ffaa=rho*era1
  ffab=rho*erab0-ffaa
  dffaa_drho=era1+rho*dera1_dra*dra_drho
  dffab_drho=erab0+rho*derab0_drab*drab_drho-dffaa_drho
  ! mb-> i-loop removed
  denaa=1.d0/(1.0d0+caa(6)*xa2)
  denab=1.d0/(1.0d0+cab(6)*xa2)
  denx =1.d0/(1.0d0+cx(6)*xa2)
  f83rho=fr83/rho
  bygr=2.0d0/gr
  gaa=caa(1)+uaa*(caa(2)+uaa*(caa(3)+uaa*(caa(4)+uaa*caa(5))))
  gab=cab(1)+uab*(cab(2)+uab*(cab(3)+uab*(cab(4)+uab*cab(5))))
  gx=cx(1)+ux*(cx(2)+ux*(cx(3)+ux*(cx(4)+ux*cx(5))))
  taa=denaa*uaa*(caa(2)+uaa*(2.d0*caa(3)+uaa &
      *(3.d0*caa(4)+uaa*4.d0*caa(5))))
  tab=denab*uab*(cab(2)+uab*(2.d0*cab(3)+uab &
       *(3.d0*cab(4)+uab*4.d0*cab(5))))
  txx=denx*ux*(cx(2)+ux*(2.d0*cx(3)+ux &
       *(3.d0*cx(4)+ux*4.d0*cx(5))))
  dgaa_drho=-f83rho*taa
  dgab_drho=-f83rho*tab
  dgx_drho=-f83rho*txx
  dgaa_dgr=bygr*taa
  dgab_dgr=bygr*tab
  dgx_dgr=bygr*txx
  ! mb
  sx=ex*gx+ffaa*gaa+ffab*gab
  v1x=dex_drho*gx+ex*dgx_drho &
       +dffaa_drho*gaa+ffaa*dgaa_drho &
       +dffab_drho*gab+ffab*dgab_drho
  v2x=(ex*dgx_dgr+ffaa*dgaa_dgr+ffab*dgab_dgr)/gr
  return
end subroutine hcth
!-------------------------------------------------------------------=
subroutine pwcorr(r,c,g,dg)
  implicit none
  real(8) :: r, g, dg, c(6)
  real(8) :: r12, r32, r2, rb, drb, sb

  r12=dsqrt(r)
  r32=r*r12
  r2=r*r
  rb=c(3)*r12+c(4)*r+c(5)*r32+c(6)*r2
  sb=1.0d0+1.0d0/(2.0d0*c(1)*rb)
  g=-2.0d0*c(1)*(1.0d0+c(2)*r)*dlog(sb)
  drb=c(3)/(2.0d0*r12)+c(4)+1.5d0*c(5)*r12+2.0d0*c(6)*r
  dg=(1.0d0+c(2)*r)*drb/(rb*rb*sb)-2.0d0*c(1)*c(2)*dlog(sb)

  return
end subroutine pwcorr
!-----------------------------------------------------------------------------
!     ==================================================================
subroutine optx(rho,grho,sx,v1x,v2x)
!     OPTX, Handy et al. JCP 116, p. 5411 (2002) and refs. therein
!     Present release: Mauro Boero, Tsukuba, 10/9/2002
!--------------------------------------------------------------------------
!     rhoa = rhob = 0.5 * rho in LDA implementation
!     grho is the SQUARE of the gradient of rho! --> gr=sqrt(grho)
!     sx  : total exchange correlation energy at point r
!     v1x : d(sx)/drho
!     v2x : 1/gr*d(sx)/d(gr)
!--------------------------------------------------------------------------
  implicit none
  real(8) :: rho, grho, sx, v1x, v2x

  real(8), parameter :: small=1.D-30, smal2=1.D-10
!.......coefficients and exponents....................
  real(8), parameter :: o43=4.0d0/3.0d0, two13=1.259921049894873D0, &
       two53=3.174802103936399D0, gam=0.006D0, a1cx=0.9784571170284421D0,&
       a2=1.43169D0 
  real(8) :: gr, rho43, xa, gamx2, uden, uu
  !.......OPTX in compact form..........................
  if(rho <= small) then
     sx=0.0D0
     v1x=0.0D0
     v2x=0.0D0
  else
     gr = max(grho,SMAL2)
     rho43=rho**o43
     xa=two13*DSQRT(gr)/rho43
     gamx2=gam*xa*xa
     uden=1.d+00/(1.d+00+gamx2)
     uu=a2*gamx2*gamx2*uden*uden
     uden=rho43*uu*uden
     sx=-rho43*(a1cx+uu)/two13
     v1x=o43*(sx+two53*uden)/rho
     v2x=-two53*uden/gr
  endif
  return
end subroutine optx
!
!---------------------------------------------------------------
subroutine wcx (rho, grho, sx, v1x, v2x)
  !---------------------------------------------------------------
  !
  ! Wu-Cohen exchange (without Slater exchange):
  ! Z. Wu and R. E. Cohen, PRB 73, 235116 (2006)
  !
  implicit none
  REAL(8), PARAMETER :: pi     = 3.14159265358979323846_8 
  REAL(8), PARAMETER :: tpi    = 2.0_8 * pi
  REAL(8), PARAMETER :: HARTREE_SI       = 4.35974394E-18_8   ! J
  REAL(8), PARAMETER :: C_SI             = 2.99792458E+8_8    ! m sec^-1
  REAL(8), PARAMETER :: BOHR_RADIUS_SI   = 0.52917720859E-10_8 ! m
  REAL(8), PARAMETER :: H_PLANCK_SI      = 6.62606896E-34_8   ! J s
  REAL(8), PARAMETER :: AU_SEC           = H_PLANCK_SI/tpi/HARTREE_SI
  REAL(8), PARAMETER :: C_AU             = C_SI / BOHR_RADIUS_SI * AU_SEC
  real(8) :: rho, grho, sx, v1x, v2x
  ! input: charge and squared gradient
  ! output: energy
  ! output: potential
  ! local variables
  real(8) :: kf, agrho, s1, s2, es2, ds, dsg, exunif, fx
  ! (3*pi2*|rho|)^(1/3)
  ! |grho|
  ! |grho|/(2*kf*|rho|)
  ! s^2
  ! n*ds/dn
  ! n*ds/d(gn)
  ! exchange energy LDA part
  ! exchange energy gradient part
  real(8) :: dxunif, dfx, f1, f2, f3, dfx1, x1, x2, x3, &
              dxds1, dxds2, dxds3
  ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
  real(8) :: third, c1, c2, c5, c6, teneightyone
  parameter (third = 1.d0 / 3.d0, c1 = 0.75d0 / pi , &
       c2 = 3.093667726280136d0, c5 = 4.d0 * third, &
       teneightyone = 0.123456790123d0)
  ! parameters of the functional
  real(8) :: k, mu, cwc
  parameter (k = 0.804d0, mu = 0.2195149727645171d0, cwc = 0.00793746933516d0) 
  !
  agrho = sqrt (grho)
  kf = c2 * rho**third
  dsg = 0.5d0 / kf
  s1 = agrho * dsg / rho
  s2 = s1 * s1
  es2 = exp(-s2)
  ds = - c5 * s1
  !
  !   Energy
  !
  ! x = 10/81 s^2 + (mu - 10/81) s^2 e^-s^2 + ln (1 + c s^4)
  x1 = teneightyone * s2 
  x2 = (mu - teneightyone) * s2 * es2
  x3 = log(1.d0 + cwc * s2 * s2)
  f1 = (x1 + x2 + x3) / k
  f2 = 1.d0 + f1
  f3 = k / f2
  fx = k - f3
  exunif = - c1 * kf
  sx = exunif * fx
  !
  !   Potential
  !
  dxunif = exunif * third
  dfx1 = f2 * f2
  dxds1 = teneightyone
  dxds2 = (mu - teneightyone) * es2 * (1.d0 - s2)
  dxds3 = 2.d0 * cwc * s2 / (1.d0 + cwc * s2 *s2)
  dfx = 2.d0 * s1 * (dxds1 + dxds2 + dxds3) / dfx1
  v1x = sx + dxunif * fx + exunif * dfx * ds
  v2x = exunif * dfx * dsg / agrho

  sx = sx * rho
  return
end subroutine wcx
!
!-----------------------------------------------------------------------
function dpz (rs, iflg)
  !-----------------------------------------------------------------------
  !  derivative of the correlation potential with respect to local density
  !  Perdew and Zunger parameterization of the Ceperley-Alder functional
  !
  !
  implicit none
  REAL(8), PARAMETER :: pi     = 3.14159265358979323846_8 
  REAL(8), PARAMETER :: tpi    = 2.0_8 * pi
  REAL(8), PARAMETER :: fpi    = 4.0_8 * pi
  REAL(8), PARAMETER :: HARTREE_SI       = 4.35974394E-18_8   ! J
  REAL(8), PARAMETER :: C_SI             = 2.99792458E+8_8    ! m sec^-1
  REAL(8), PARAMETER :: BOHR_RADIUS_SI   = 0.52917720859E-10_8 ! m
  REAL(8), PARAMETER :: H_PLANCK_SI      = 6.62606896E-34_8   ! J s
  REAL(8), PARAMETER :: AU_SEC           = H_PLANCK_SI/tpi/HARTREE_SI
  REAL(8), PARAMETER :: C_AU             = C_SI / BOHR_RADIUS_SI * AU_SEC
  !
  real(8), intent (in) :: rs
  integer, intent(in) :: iflg
  real(8) :: dpz
  !
  !  local variables
  !  a,b,c,d,gc,b1,b2 are the parameters defining the functional
  !
  real(8), parameter :: a = 0.0311d0, b = -0.048d0, c = 0.0020d0, &
       d = -0.0116d0, gc = -0.1423d0, b1 = 1.0529d0, b2 = 0.3334d0,&
       a1 = 7.0d0 * b1 / 6.d0, a2 = 4.d0 * b2 / 3.d0
  real(8) :: x, den, dmx, dmrs
  !
  !
  if (iflg == 1) then
     dmrs = a / rs + 2.d0 / 3.d0 * c * (log (rs) + 1.d0) + &
          (2.d0 * d-c) / 3.d0
  else
     x = sqrt (rs)
     den = 1.d0 + x * (b1 + x * b2)
     dmx = gc * ( (a1 + 2.d0 * a2 * x) * den - 2.d0 * (b1 + 2.d0 * &
          b2 * x) * (1.d0 + x * (a1 + x * a2) ) ) / den**3
     dmrs = 0.5d0 * dmx / x
  endif
  !
  dpz = - fpi * rs**4.d0 / 9.d0 * dmrs
  return
  !
end function dpz
