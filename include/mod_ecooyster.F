!
! Fortran 95 implementation of `EcoOyster', based on the set of equations shared
!   by Mark J. Brush (brush@vims.edu, pers.comm., 20210915, 20220624, 20230710).
!   The implementation is by pst-laurent@vims.edu and was last modified on
!   20230824.
!   The code is implemented in such a way that it can be used as a stand-alone
!   Fortran application, or compiled and executed within ROMS, with relatively
!   few modifications. Blocks that are specific to execution within ROMS are
!   identified with #ifdef ECOOYSTER, while blocks specific to the stand-alone
!   execution are identified with #ifndef ECOOYSTER ("if-not-defined").
!
#include "cppdefs.h" ! CAREFUL, un-comment this line if running inside ROMS.
      module mod_ecooyster
#ifdef ECOOYSTER
      use mod_kinds
      use mod_param, only: NDBio2d
#endif
      implicit none
      private
!
#ifndef ECOOYSTER
!     Define parameters and variables that in ROMS are defined: (1) from reading
!       a .in file, (2) from a module called with a `USE' statement, or (3) from
!       a header file called with a #include statement. In the ROMS code, these
!       modules/headers are called in nearly every procedures. We mimic this
!       ubiquity (and avoid repetition) by making them globally accessible.
!
      integer,    parameter ::                                          &
!       For historical reasons, ROMS uses `r8' to represent the working precision
!         of floating point numbers, regardless of what this precision is.
!         Here, `r8' corresponds to a 4-bytes number and would match a ROMS
!         calculation that has SINGLE_PRECISION activated.
!       `dp' corresponds to a true 8-bytes number (double precision).
     &  dp         = selected_real_kind( P = 12 ),                      &
     &  r8         = selected_real_kind( P =  6 ),                      &
     &  Ngrids     = 1,                                                 &
     &  Lm(Ngrids) = 1,                                                 &
     &  Mm(Ngrids) = 1,                                                 &
     &  N( Ngrids) = 1,                                                 &
     &  NAT        = 2,                                                 &
     &  NPT        = 0,                                                 &
     &  NCS        = 0,                                                 &
     &  NNS        = 4,                                                 &
     &  NBT        = 17,                                                &
     &  NDBio2d    = 19,                                                &
     &  NT(Ngrids) = NAT + NPT + NCS + NNS + NBT,                       &
     &  Istr       = 1,                                                 &
     &  Iend       = Lm(1),                                             &
     &  Jstr       = 1,                                                 &
     &  Jend       = Mm(1),                                             &
     &  IminS      = Istr - 3
      integer,    parameter ::                                          &
     &  ImaxS      = Iend + 3,                                          &
     &  JminS      = Jstr - 3,                                          &
     &  JmaxS      = Jend + 3,                                          &
     &  I_padd     = floor( real(Lm(1) + 2) / 2. )                      &
     &             - floor( real(Lm(1) + 1) / 2. ),                     &
     &  J_padd     = floor( real(Mm(1) + 2) / 2. )                      &
     &             - floor( real(Mm(1) + 1) / 2. ),                     &
     &  Im(Ngrids) = Lm(1) + I_padd,                                    &
     &  Jm(Ngrids) = Mm(1) + J_padd,                                    &
     &  Imin       = 0,                                                 &
     &  Imax       = Im(1) + 1,                                         &
     &  Jmin       = 0,                                                 &
     &  Jmax       = Jm(1) + 1
      integer,    parameter ::                                          &
     &  LBi        = Imin,                                              &
     &  UBi        = Imax,                                              &
     &  LBj        = Jmin,                                              &
     &  UBj        = Jmax,                                              &
     &  itemp      = 1,                                                 &
     &  isalt      = 2,                                                 &
     &  iPhyt      = NAT + NPT + NCS + NNS +  4,                        &
     &  iZoop      = NAT + NPT + NCS + NNS +  5,                        &
     &  iLDeC      = NAT + NPT + NCS + NNS +  8,                        &
     &  iSDeC      = NAT + NPT + NCS + NNS +  9,                        &
     &  iTIC_      = NAT + NPT + NCS + NNS + 10,                        &
     &  iTAlk      = NAT + NPT + NCS + NNS + 11,                        &
     &  iOxyg      = NAT + NPT + NCS + NNS + 12
      integer,    parameter ::                                          &
     &  iwshd      = 11,                                                &
     &  iwtid      = 12,                                                &
     &  iwgon      = 13,                                                &
     &  iwsht      = 14,                                                &
     &  iwtit      = 15,                                                &
     &  ihshd      = 16,                                                &
     &  iherd      = 17,                                                &
     &  ihsht      = 18,                                                &
     &  ihert      = 19
!
      real( r8 ), parameter ::                                          &
     &  sec2day    = 1._r8 / (60._r8 * 60._r8 * 24._r8)
!
!     These are initialized inside the procedure `initialize_variables':
      integer    :: ng, tile, ntimes, nstp( Ngrids ), nnew( Ngrids )
!
      real( r8 ) ::                                                     &
     &  t(         LBi : UBi, LBj : UBj, N(1), 3, NT(1) ),              &
     &  rmask(     LBi : UBi, LBj : UBj                 ),              &
     &  rmask_wet( LBi : UBi, LBj : UBj                 ),              &
     &  pH(        LBi : UBi, LBj : UBj                 ),              &
     &  h(         LBi : UBi, LBj : UBj                 ),              &
     &  ZooCN( Ngrids ), PhyCN( Ngrids ), dt( Ngrids )
      real( dp ) :: dstart, time_ref
!
      public run_stand_alone
#else
      public ecooyster
#endif
!
!     Parameters and variables that are specific to ecooyster:
!
      integer,    parameter ::                                          &
!       Number of oyster classes (>= 1). The original set of equations had two
!         classes corresponding to small/large oysters. The latest set has two
!         classes corresponding to diploids (j=1) and triploids (j=2).
     &  ncla        = 2,                                                &
!       Number of iterations in the Euler backward temporal discretization:
     &  OysIter     = 1
!
      real( r8 ), parameter ::                                          &
     &  epsi        = 1.e-9,                                            &
!       p_j,a_j,b_j,c_j are parameters specific to the oyster classes.
!    &  p_j( ncla ) = (/  0.25, 0.25 /),                                &!Blachman20230713
!    &  p_j( ncla ) = (/  0.20, 0.20 /),                                &!Brush20230710
     &  p_j( ncla ) = (/  0.15, 0.15 /),                                &!Czajka20230903
     &  a_j( ncla ) = (/ 158.2, 158.2/),                                &
     &  b_j( ncla ) = (/ 1.052, 1.052/),                                &
     &  c_j( ncla ) = (/ 0.932, 0.932/),                                &
!       Assimilation efficiency (0<=asef<=1); balance assumed lost to feces.
     &  asef        = 0.75,                                             &
!       Activity respiration, i.e. the fraction `rfra' (0<=rfra<=1) of
!         assimilation that is lost to respiration to support filtering activity
!         of the organism.
     &  rfra        = 0.10,                                             &
!       Carbon to dry weight ratio of oyster tissue (g/g):
     &  c_dw        = 0.45
      real( dp ), parameter ::                                          &
     &  e_dp = tiny( 0._dp ),                                           &
!     Half-saturation coefficient for Michaelis-Menten function inside
!       resorption equation (relevant for gonadal tissue of diploids).
     &  kgon = 2.e-3_dp !psl20230714
!    &  kgon = 1.e-2_dp
!
#ifndef ECOOYSTER
!     When running in stand-alone mode, DiaBio2d is an automatic array of fixed
!       dimensions (recall that LBi and friends are parameters in this case.)
!       Because the dimensions are fixed, the `save' attribute can be used.
      real( r8 ), save ::                                               &
     &  DiaBio2d( LBi : UBi, LBj : UBj, NDbio2d )
#endif
!
      contains
!
      subroutine ecooyster( ng, tile, myit )
#ifdef ECOOYSTER
      use mod_param
      use mod_diags
      use mod_forces
      use mod_grid
      use mod_ncparam
      use mod_ocean
      use mod_stepping
#endif
      implicit none
      integer, intent(in) :: ng, tile, myit
#ifdef ECOOYSTER
#include "tile.h"
#endif
!
      call ecooyster_tile( ng, tile, myit,                              &
     &                     LBi, UBi, LBj, UBj, N(ng), NT(ng),           &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp(ng), nnew(ng),                          &
#ifdef ECOOYSTER
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % rmask_wet,                        &
     &                     GRID(ng) % h,                                &
     &                     DIAGS(ng) % DiaBio2d,                        &
     &                     OCEAN(ng) % pH,                              &
     &                     OCEAN(ng) % t )
#else
     &                     rmask,                                       &
     &                     rmask_wet,                                   &
     &                     h,                                           &
     &                     DiaBio2d,                                    &
     &                     t )
#endif
      end subroutine ecooyster
!
      subroutine ecooyster_tile( ng, tile, myit,                        &
     &                     LBi, UBi, LBj, UBj, UBk, UBt,                &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     nstp, nnew,                                  &
     &                     rmask,                                       &
     &                     rmask_wet,                                   &
     &                     h,                                           &
     &                     DiaBio2d,                                    &
#ifdef ECOOYSTER
     &                     pH,                                          &
#endif
     &                     t )
#ifdef ECOOYSTER
      use biology_mod,   only: pCO2_water_RZ ! From ecb.h.
      use mod_param
      use mod_biology
      use mod_ncparam
      use mod_scalars
      use dateclock_mod, only: caldate
#endif
      implicit none
      integer,  intent(in   ) :: ng, tile, myit
      integer,  intent(in   ) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer,  intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      integer,  intent(in   ) :: nstp, nnew
      real(r8), intent(in   ) :: rmask(     LBi : UBi, LBj : UBj )
      real(r8), intent(in   ) :: rmask_wet( LBi : UBi, LBj : UBj )
      real(r8), intent(in   ) :: h(         LBi : UBi, LBj : UBj )
      real(r8), intent(inout) :: DiaBio2d(LBi : UBi, LBj : UBj, NDbio2d)
#ifdef ECOOYSTER
      real(r8), intent(inout) :: pH(      LBi : UBi, LBj : UBj         )
#endif
      real(r8), intent(inout) :: t( LBi : UBi, LBj : UBj, UBk, 3, UBt )
      integer,  parameter     :: DoNewton = 0
!
      integer    :: i, j, k, Iter, itrc, jcla, mm_i
      real( r8 ) :: dtdays, cff1, cff2, cff3, poc, tss, k_d, cchl, phyc,&
     &              f_te, f_sa, f_ts, f_o2, f_ch, reff, repr, reso,     &
     &  pCO2(    IminS : ImaxS                ),                        &
     &  cass(    IminS : ImaxS                ),                        &
     &  Bio(     IminS : ImaxS, N(ng), NT(ng) )
      real( dp ) ::                                                     &
     &  wtol( ncla ), wtne( ncla ), spaw( ncla ), eros, ctim,           &
     &  exte( LBi : UBi, ncla ),                                        &
!     These 5 are automatic arrays and thus their values are not guaranteed to
!       be preserved from one time-step to the next. (Automatic arrays cannot
!       have the `save' attribute.) We rely on DIAGS(ng)%DiaBio2d to store their
!       values across time-steps and to make the EcoOyster state variables
!       available outside of this module (e.g., for writing outputs.) One draw-
!       back is that DiaBio2d is type r8 and thus we lose information when
!       running in SINGLE_PRECISION.
     &  wshe( LBi : UBi, LBj : UBj, ncla ),                             &
     &  wtis( LBi : UBi, LBj : UBj, ncla ),                             &
     &  wgon( LBi : UBi, LBj : UBj, ncla ),                             &
     &  hshe( LBi : UBi, LBj : UBj, ncla ),                             &
     &  hero( LBi : UBi, LBj : UBj, ncla )
!
#ifdef ECOOYSTER
#include "set_bounds.h"
#endif
!
!     Transfer oyster state variables into the scratch automatic arrays:
      jcla = 1 ! Diploids.
      do j = Jstr, Jend
        do i = Istr, Iend
          wshe( i, j, jcla ) = real( DiaBio2d( i, j, iwshd ), dp )
          wtis( i, j, jcla ) = real( DiaBio2d( i, j, iwtid ), dp )
          wgon( i, j, jcla ) = real( DiaBio2d( i, j, iwgon ), dp )
          hshe( i, j, jcla ) = real( DiaBio2d( i, j, ihshd ), dp )
          hero( i, j, jcla ) = real( DiaBio2d( i, j, iherd ), dp )
        end do
      end do
      jcla = 2 ! Triploids.
      do j = Jstr, Jend
        do i = Istr, Iend
          wshe( i, j, jcla ) = real( DiaBio2d( i, j, iwsht ), dp )
          wtis( i, j, jcla ) = real( DiaBio2d( i, j, iwtit ), dp )
          wgon( i, j, jcla ) = 0._dp
          hshe( i, j, jcla ) = real( DiaBio2d( i, j, ihsht ), dp )
          hero( i, j, jcla ) = real( DiaBio2d( i, j, ihert ), dp )
        end do
      end do
!
      dtdays = dt(ng) * sec2day / real( OysIter, r8 )
!
!     Obtain the calendar month corresponding to the current time.
!     `CurrentTime' should be time in days since time_ref (precision dp).
      ctim = dstart + real( myit, dp ) * real( dt(ng), dp )             &
     &     / 60._dp / 60._dp / 24._dp
      call caldate( CurrentTime = ctim, mm_i = mm_i )
!
      J_LOOP : do j = Jstr, Jend
        do itrc = 1, NT(ng)
          do k = 1, N(ng)
            do i = Istr, Iend
              Bio( i, k, itrc ) = max( t(i, j, k, nstp, itrc), 0._r8 )
            end do
          end do
        end do
!
        ITER_LOOP: do Iter = 1, OysIter
          do i = Istr, Iend
!
!           Compute functions that are independent of oysters' class.
!
            f_te = exp( - 0.006_r8 * (Bio(i, 1, itemp) - 27._r8)**2 )
!
            f_sa = 0._r8
            if ( Bio(i, 1, isalt) >  5._r8 ) then
              f_sa = 0.0926_r8 * Bio(i, 1, isalt) - 0.139_r8
            end if
            if ( Bio(i, 1, isalt) > 12._r8 ) f_sa = 1._r8
!
!           Compute POC in g-C/m3 <=> mg-C/L. In her email of Sep.2,2022, Sara
!             Blachman indicated that `POC' in their EcoOyster implementation is
!             always derived from chlorophyll-a, even when actual POC
!             measurements are available. For consistency, we define `phyc' as
!             phytoplankton-carbon, and `poc' as the full POC (used in the
!             calculation of TSS):
            phyc =   Bio(i, 1, iPhyt) * PhyCN(ng) * 12._r8 * 1.e-3_r8
            poc  = ( Bio(i, 1, iZoop) * ZooCN(ng)                       &
     &             + Bio(i, 1, iPhyt) * PhyCN(ng)                       &
     &             + Bio(i, 1, iSDeC)                                   &
     &             + Bio(i, 1, iLDeC) ) * 12._r8 * 1.e-3_r8
!
!           Turner et al.: TSS=VSS+FSS, VSS=POC*2.9, FSS=0.35*VSS+ISS,
!           Convert POC to VSS (g-C/m3 to g/m3) by using the factor 2.9 of
!           Cerco & Noel (WQSTM, 2017) and then add 0.35*VSS to represent
!             the organic portion of FSS.
            tss = poc * 2.9_r8 * 1.35_r8 ! This is VSS + 0.35 * VSS.
!
!           Convert sand_xx from kg/m3 to g/m3<=>mg/L:
            do itrc = NAT + NPT + NCS + 1, NAT + NPT + NCS + NNS
              tss = tss + max( Bio(i,1,itrc), 0._r8 ) * 1.e3_r8
            end do
!
!           Compute coefficient of diffuse light attenuation (m**(-1)) based on
!             Turner et al. 2021:
            k_d  = max( 0.6_r8,                                         &
     &            0.92_r8 + 0.079_r8 * tss - 0.037_r8 * Bio(i,1,isalt) )
!
!           Carbon to chlorophyll ratio (g-C per g of chl), Cerco & Noel 2004.
            if ( mm_i >= 2 .and. mm_i <= 5 ) then ! Feb-May, spring bloom.
              cchl = 30._r8 + 150._r8 * exp( - 1.18_r8 * k_d )
            else                                  ! Outside of spring.
              cchl = 30._r8 +  90._r8 * exp( - 1.19_r8 * k_d )
            end if
!CAREFUL, YOU ARE TEMPORARILY HARDCODING cchl TO MATCH THE INPUT FILE FROM SARA
!  BLACHMAN (20230713)
!           cchl = 57.5983_r8
!
            f_ts = 1._r8
            if     ( tss > 25._r8 ) then
              f_ts = 10.364_r8 * ( log( tss ) )**(- 2.0477_r8)
!           psl20221201: The piecewise function of EcoOyster is generating
!             unphysical discontinuities wherever TSS drops below 4mg/L. After
!             an email exchange between Catherine and Mark, we decided to remove
!             the sharp drop from the function `f_ts'.
!           elseif ( tss <  4._r8 ) then
!             f_ts = 0.1_r8
            end if
!
!           Convert O2 from mmol-O2/m3 to mg-O2/L:
            f_o2 = Bio(i, 1, iOxyg) * 32.e-3_r8
            f_o2 = (1.75_r8 - f_o2) / (1.75_r8 - 1.5_r8)
            f_o2 = 1._r8 / ( 1._r8 + exp( 1.1_r8 * f_o2 ) )
!
!           Compute chlorophyll-a (in mg-chla/m3) from POC (in g-C/m3) assuming
!             a mass ratio CChl. The original value suggested by Mark on
!             20210915 was 64:1, but Sara used 57.5:1 in her input file for
!             station LE4.3.
!           The original function produces non-physical negative values if
!             Chla<2.2, and thus a limiter (f_ch > 0) is imposed.
!
            f_ch = 1.e3_r8 * phyc / cchl
!           f_ch = 1.e3_r8 * poc  / cchl
!           f_ch = 1.e3_r8 * poc * 0.017391_r8 ! 1/57.5=0.017391 (Sara)
!           f_ch = 1.e3_r8 * poc * 0.015625_r8 ! 1/64  =0.015625 (Mark)
!
            f_ch = 1._r8 - exp( - 0.26_r8 * ( f_ch - 2.2_r8 ) )
            f_ch = max( f_ch, 0._r8 )
!
!           Reproduction efficiency `reff':
            if ( mm_i >= 1 .and. mm_i <= 6 ) then ! Jan-Jun.
              reff = min( 1._r8,                                        &
     &               max(0._r8, 0.054_r8 * Bio(i, 1, itemp) - 0.729_r8))
            else                                  ! Jul-Dec.
              reff = min( 1._r8,                                        &
     &               max(0._r8, 0.047_r8 * Bio(i, 1, itemp) - 0.809_r8))
            end if
!
!           Tissue assimilation and respiration, spawning:
            do jcla = 1, ncla
!             Compute the total tissue weight (somatic+gonadal) at the `old'
!               time-step:
              wtol( jcla ) = max( wtis(i, j, jcla)                      &
     &                          + wgon(i, j, jcla), real(epsi, dp) )
!
              cff1 = 0.17_r8 * f_te * f_sa * f_ts * f_o2 * f_ch         &
     &             * p_j( jcla ) * poc   / c_dw * asef * (1._r8 - rfra) & !psl20230524
     &             - 0.0095_r8                                          &
     &             * exp( 0.069_r8 * (Bio(i, 1, itemp) - 20._r8) )
              cff1 = real( real( cff1,   dp ) / wtol( jcla )**0.25_dp   &
     &                   * real( dtdays, dp ), r8 )
              spaw( jcla ) = wgon(i, j, jcla) / wtol( jcla )
              spaw( jcla ) = max( spaw( jcla ) - 0.2_dp, 0._dp )        &
     &                     * wgon(i, j, jcla)
!
!             Hayes & Menzel 1981, spawning only occurs when T>=25C:
!             spaw( jcla ) = spaw( jcla )                               &
!    &                     * real(floor(Bio(i, 1, itemp) / 25._r8), dp)
!
              if ( jcla == 2 ) then
                spaw( jcla ) = 0._dp
              end if
              wtne( jcla ) = (wtol(jcla) - spaw(jcla) *real(dtdays, dp))&
     &                     / (1._dp - real(cff1, dp))
              wtne( jcla ) = max( wtne( jcla ), real(epsi, dp) )
            end do
!
!           Update gonadal tissue:
            do jcla = 1, ncla
!             Temporarily store Assim-Resp inside `cff2':
              cff2 = 0.17_r8 * f_te * f_sa * f_ts * f_o2 * f_ch         &
     &             * p_j( jcla ) * poc / c_dw * asef * (1._r8 - rfra)   &
     &             - 0.0095_r8                                          &
     &             * exp( 0.069_r8 * (Bio(i, 1, itemp) - 20._r8) )
              cff2 = real( real(cff2, dp) * wtne( jcla )**0.75_dp, r8 )
!
!             Compute shell extension (mm/days):
              exte( i, jcla ) = (real( cff2,   dp ) - spaw( jcla ))     &
     &                        *  real( dtdays, dp ) ! (d wtotal/dt) x dt.
              if     ( jcla == 1 ) then
                exte( i, jcla ) = max( 0._dp, 73.85_dp  * ( wtne( jcla )&
     &                          + exte( i, jcla ) )**0.268              &
     &                          - hshe( i, j, jcla ) - hero( i, j, jcla)&
     &                          ) / real( dtdays, dp )
              elseif ( jcla == 2 ) then
                exte( i, jcla ) = max( 0._dp, 63.243_dp * ( wtne( jcla )&
     &                          + exte( i, jcla ) )**0.4187             &
     &                          - hshe( i, j, jcla ) - hero( i, j, jcla)&
     &                          ) / real( dtdays, dp )
              end if
!
              repr = max( 0._r8, cff2 ) * reff
              reso = real( wgon(i, j, jcla)                             &
     &                  / (wgon(i, j, jcla) + kgon), r8 )               &
     &             * max( 0._r8, - cff2 )
!
              cff2 = real( wgon(i, j, jcla) / wtol(jcla), r8 )
              cff2 = dtdays * max( 0._r8 , cff2 - 0.2_r8 )
              if ( jcla == 2 ) then
                repr = 0._r8
                reso = 0._r8
                cff2 = 0._r8
              end if
              wgon(i, j, jcla) = ( wgon(i, j, jcla)                     &
     &                           + real( dtdays * repr, dp )            &
     &                           - real( dtdays * reso, dp ) )          &
     &                         / (1._dp + real(cff2, dp))
            end do
!
!           Update somatic tissue weight:
            do jcla = 1, ncla
              wtis(i, j, jcla) = wtne( jcla ) - wgon(i, j, jcla)
            end do
!
          end do ! i loop.
!
          k = 1 ! Select the rho vertical level that's closest to the sea floor.
!         Compute pH from salt,temp,TIC,TAlk, for the current `j' and the whole
!           `i' dimension.
          call pCO2_water_RZ( Istr, Iend, LBi, UBi, LBj, UBj,           &
     &                        IminS, ImaxS, j, DoNewton,                &
     &                        rmask,                                    &
     &                        Bio(IminS:,k,itemp), Bio(IminS:,k,isalt), &
     &                        Bio(IminS:,k,iTIC_), Bio(IminS:,k,iTAlk), &
     &                        pH, pCO2 )
!
!         Compute calcite saturation state from pH,temp,salt,TIC,pressure:
          call calcite_saturation_state(                                &
     &           Istr, Iend, LBi, UBi, LBj, UBj,                        &
     &           IminS, ImaxS, j,                                       &
     &           rmask,                                                 &
     &           h,                                                     &
     &           Bio(IminS:,k,itemp), Bio(IminS:,k,isalt),              &
     &           Bio(IminS:,k,iTIC_), pH, cass )
!
          do i = Istr, Iend
!           Net calcification:
            do jcla = 1, ncla
              cff3 = 1._r8 - exp( - b_j(jcla) * (cass(i) - c_j(jcla)) )
              cff3 = a_j(jcla) * cff3 * 1.e-3_r8                        &
     &             * exp( 0.0271_r8 * (Bio(i, 1, itemp) - 25._r8) )
              cff3 = cff3 * dtdays
!
              wshe(i, j, jcla) =  wshe(i, j, jcla)                      &
     &                         + (wtis(i, j, jcla) + wgon(i, j, jcla))  &
     &                         * real( cff3, dp )
              wshe(i, j, jcla) = max( wshe(i, j, jcla), 0._dp )
!
!             Compute erosion of shell height due to dissolution (mm/day):
              eros = real( cff3, dp )                                   &
     &             * (wtis(i, j, jcla) + wgon(i, j, jcla)) ! (d wshe/dt) x dt.
!             Bound `eros' for cases where `wshe' -> 0:
              eros = max( eros, - wshe(i, j, jcla) + epsi )
              eros = max( 0._dp,                                        &
     &                    26.24_dp *  wshe(i, j, jcla        )**0.2594  &
     &                  - 26.24_dp * (wshe(i, j, jcla) + eros)**0.2594 )&
     &             / real( dtdays, dp )
!
!             Compute eroded shell height (mm), i.e. the cumulated erosion:
              hero(i, j, jcla) = hero(i, j, jcla)                       &
     &                         + eros * real( dtdays, dp )
!
!             Update the shell height (mm):
              hshe(i, j, jcla) = hshe(i, j, jcla)                       &
     &                         + (exte(i, jcla) - eros)                 &
     &                         * real( dtdays, dp )
              hshe(i, j, jcla) = max( hshe(i, j, jcla), 0._dp )
            end do
          end do ! i loop.
        end do ITER_LOOP
      end do J_LOOP
!
!     Transfer wshe,wtis,wgon from the scratch automatic arrays to NDBio2d.
      jcla = 1 ! Diploids.
      do j = Jstr, Jend
        do i = Istr, Iend
          DiaBio2d(i,j,iwshd) = real(wshe(i, j, jcla), r8) * rmask(i,j)
          DiaBio2d(i,j,iwtid) = real(wtis(i, j, jcla), r8) * rmask(i,j)
          DiaBio2d(i,j,iwgon) = real(wgon(i, j, jcla), r8) * rmask(i,j)
          DiaBio2d(i,j,ihshd) = real(hshe(i, j, jcla), r8) * rmask(i,j)
          DiaBio2d(i,j,iherd) = real(hero(i, j, jcla), r8) * rmask(i,j)
        end do
      end do
      jcla = 2 ! Triploids.
      do j = Jstr, Jend
        do i = Istr, Iend
          DiaBio2d(i,j,iwsht) = real(wshe(i, j, jcla), r8) * rmask(i,j)
          DiaBio2d(i,j,iwtit) = real(wtis(i, j, jcla), r8) * rmask(i,j)
          DiaBio2d(i,j,ihsht) = real(hshe(i, j, jcla), r8) * rmask(i,j)
          DiaBio2d(i,j,ihert) = real(hero(i, j, jcla), r8) * rmask(i,j)
        end do
      end do
      end subroutine ecooyster_tile
!
      subroutine calcite_saturation_state(                              &
     &       Istr, Iend, LBi, UBi, LBj, UBj, IminS, ImaxS, j,           &
     &       rmask,                                                     &
     &       hbot, TempC, Sal, TIC, pH, cass )
      implicit none
!     This subroutine is based on CO2SYS.m version 2.0 (20 Dec.2016). The
!       implementation uses the dissociation constants of Cai & Wang 1998 and
!       pH as input (assumed to be on the seawater scale). The subroutine also
!       takes DIC (in millimoles-C/m3) as input and assumes for simplicity that
!       1 mmol-C/m3 <=> 1 umol-C/kg-seawater (i.e., we accept the <3% error).
!     If used with double precision arithmetic (r8=dp), the subroutine produces
!       the same calcite saturation state as the aforementioned CO2SYS code.
!     Check value: DIC   = 1474.29 mmol-C/m3 <=> 1474.29 umol-C/kgSW,
!                  pH_sws= 8.62941, salinity= 19.9825, temp= 6.8575dec.C,
!                  press = 5 dbars, Omega_calcite = 5.443654.
!
      integer,    intent(in   ) ::                                      &
     &  Istr, Iend, LBi, UBi, LBj, UBj, IminS, ImaxS, j
!
      real( r8 ), intent(in   ) ::                                      &
     &  rmask( LBi : UBi, LBj : UBj ),                                  &
     &  hbot(  LBi : UBi, LBj : UBj ),                                  &
     &  TempC( IminS : ImaxS        ),                                  &
     &  Sal(   IminS : ImaxS        ),                                  &
     &  TIC(   IminS : ImaxS        ),                                  &
     &  pH(    LBi : UBi, LBj : UBj )
!
      real( r8 ), intent(  out) ::                                      &
     &  cass(  IminS : ImaxS        )
!
      real( r8 ) ::                                                     &
     &  Ca, TempK, logTempK, sqrSal, logKCa, KCa, deltaVKCa,            &
     &  KappaKCa, lnKCafac, H, CO3, Pbar, RT, F1, F2, K1, K2,           &
     &  pK1, pK2, fH, TC, deltaV, Kappa, lnK1fac, lnK2fac, K1fac, K2fac
      integer    :: i
!
      I_LOOP: do i = Istr, Iend
        if ( rmask(i, j) .le. tiny(0._r8) ) then ! Dry cell.
          cass( i ) = 0._r8
          continue
        end if
        Ca       = 0.02128_r8 / 40.087_r8  * (Sal( i ) / 1.80655_r8)
!
        TempK    = tempC( i ) + 273.15_r8
        logTempK = log( TempK )
        logKCa   = - 171.9065_r8 - 0.077993_r8 * TempK                  &
     &                           + 2839.319_r8 / TempK
        sqrSal   = sqrt( Sal( i ) )
        logKCa   = logKCa + 71.595_r8 * logTempK / log( 10._r8 )
        logKCa   = logKCa + ( - 0.77712_r8 + 0.0028426_r8 * TempK       &
     &                                     + 178.34_r8    / TempK )     &
     &           * sqrSal
        logKCa   = logKCa - 0.07711_r8   * Sal( i )                     &
     &           + 0.0041249_r8 * sqrSal * Sal( i )
        KCa      = 10._r8**logKCa ! (mol/kg-SW)^2
!
!       CO2SYS uses as input the hydrostatic pressure (zero at sea surface).
!         Here we neglect differences between the sea floor depth and the
!         position of the rho point, we assume an undisturbed sea surface, and
!        we assume that 1meter of seawater equates to a pressure of 1decibar.
!       Note that hbot is bounded to zero in case that WET_DRY is activated.
        Pbar      = max(hbot(i, j), 0._r8) * 1.e-1_r8 ! In bars.
!       Assume an ideal gas constant `R' of 83.1451 mL bar-1 K-1 mol-1.
        RT        = 83.1451_r8 * TempK
        deltaVKCa =  - 48.76_r8 + 0.5304_r8 * TempC( i )
        KappaKCa  = (- 11.76_r8 + 0.3692_r8 * TempC( i )) * 1.e-3_r8
        lnKCafac  = (- deltaVKCa + 0.5 * KappaKCa * Pbar) * Pbar / RT
        KCa       = KCa * exp( lnKCafac )
!
        F1      = 200.1_r8   / TempK + 0.3220_r8
        pK1     = 3404.71_r8 / TempK + 0.032786_r8 * TempK -14.8435_r8  &
     &          - 0.071692_r8  * F1 * sqrt(Sal( i ))                    &
     &          + 0.0021487_r8 * Sal( i )
        fH      = 1.2948_r8 - 0.002036_r8 * TempK                       &
     &          + (0.0004607_r8 - 0.000001475_r8 * TempK) * Sal( i )**2
        K1      = 10._r8**( - pK1 ) / fH
        F2      = - 129.24_r8 / TempK + 1.4381_r8
        pK2     = 2902.39_r8  / TempK + 0.02379_r8 * TempK - 6.4980_r8  &
     &          - 0.3191_r8 * F2 * sqrt(Sal(i)) + 0.0198_r8 * Sal(i)
        K2      = 10._r8**( - pK2      ) / fH
!
        deltaV  =  - 25.5_r8 + 0.1271_r8 * TempC(i)
        Kappa   = (- 3.08_r8 + 0.0877_r8 * TempC(i)) * 1.e-3_r8
        lnK1fac = (- deltaV + 0.5_r8 * Kappa * Pbar) * Pbar / RT
!
        deltaV  = - 15.82_r8 - 0.0219_r8 * TempC(i)
        Kappa   = (  1.13_r8 - 0.1475_r8 * TempC(i)) * 1.e-3_r8
        lnK2fac = (- deltaV + 0.5_r8 * Kappa * Pbar) * Pbar / RT
        K1fac   = exp( lnK1fac )
        K2fac   = exp( lnK2fac )
        K1      = K1 * K1fac
        K2      = K2 * K2fac
!
        H       = 10._r8**( - pH(i, j) )
!       Convert DIC from mmol-C/m3 to mol-C/kg-seawater (approximately):
        TC      = TIC( i ) * 1.e-6_r8
        CO3     = TC * K1 * K2 / (K1 * H + H**2 + K1 * K2)
!       Compute CAlcite Saturation State (cass):
        cass(i) = CO3 * Ca / KCa
      end do I_LOOP
      end subroutine calcite_saturation_state
!
#ifndef ECOOYSTER
      subroutine run_stand_alone()
      implicit none
      integer :: iter
!
      call initialize_variables()
!
      ng   =   1
      tile = - 1
      do iter = 1, ntimes
        call update_env_cond(     iter )
        call ecooyster( ng, tile, iter )
        call write_results(       iter )
      end do
      end subroutine run_stand_alone
#endif
!
#ifndef ECOOYSTER
      subroutine initialize_variables()
      implicit none
!     This subroutine is only relevant if the code is used in a stand-alone
!       context (i.e. not used within ROMS).
!
!     Initial oyster somatic tissue weight (g dry weight):
      DiaBio2d(:, :, iwtid) = 0.1_r8 ! Oyster class #1 (diploids ).
      DiaBio2d(:, :, iwtit) = 0.1_r8 ! Oyster class #2 (triploids).
!
!     Initial oyster gonadal tissue weight (g dry weight):
      DiaBio2d(:, :, iwgon) = 0._r8  ! (Exists only for diploids.)
!
!     Initial oyster shell weight (g dry weight):
      DiaBio2d(:, :, iwshd) = 58.05_r8                                  &
     &                      * ( DiaBio2d(:,:,iwtid)                     &
     &                        + DiaBio2d(:,:,iwgon) )**0.8681_r8
!     Catherine 20230717: Use the same initial shell weight for diploids and
!       triploids because of concerns in the triploids' allometric relation
!       for "shell weight as a function of tissue+gonadal weight".
      DiaBio2d(:, :, iwsht) = 58.05_r8                                  &
     &                      *   DiaBio2d(:,:,iwtit)**0.8681_r8
!     DiaBio2d(:, :, iwsht) = 29.70_r8                                  &
!    &                      *   DiaBio2d(:,:,iwtit)**1.6140_r8
!
!     Initial oyster shell height (mm):
      DiaBio2d(:, :, ihshd) = 73.85_r8                                  &
     &                      * ( DiaBio2d(:,:,iwtid)                     &
     &                        + DiaBio2d(:,:,iwgon) )**0.268_r8
!     Catherine 20230717: Use the same initial shell height for diploids and
!       triploids (we assume that the oysters are too young for reproduction
!       to have made much of a difference between diploid/triploid).
      DiaBio2d(:, :, ihsht) = 73.85_r8                                  &
     &                      *   DiaBio2d(:,:,iwtit)**0.268_r8
!     DiaBio2d(:, :, ihsht) = 63.243_r8                                 &
!    &                      *   DiaBio2d(:,:,iwtit)**0.4187_r8
!
!     Initial value for eroded shell height (mm):
      DiaBio2d(:, :, iherd) = 0._r8
      DiaBio2d(:, :, ihert) = 0._r8
!
!     Molar Carbon:Nitrogen ratios for phytoplankton and zooplankton:
      PhyCN = 106._r8 / 16._r8 ! (Redfield).
      ZooCN = 106._r8 / 16._r8
!
      t     = 0._r8 ! Begin with zero for everything, then refine.
!     Initial bottom salinity (practical salinity units):
      t(:,:,:,:, isalt               ) = 30._r8
!     Initial bottom water temperature (degrees C):
      t(:,:,:,:, itemp               ) = 25._r8
!     Initial bottom dissolved O2  (mmol-O2/m3):
      t(:,:,:,:, iOxyg               ) = 8._r8   / 32._r8 * 1.e3_r8
!     Initial bottom phytoplankton (mmol-N /m3):
      t(:,:,:,:, iPhyt               ) = 0.75_r8 / 12._r8 * 1.e3_r8     &
     &                                 / PhyCN(1)
!     Initial bottom sand conc.    (kg/m3):
      t(:,:,:,:, NAT + NPT + NCS + 1 ) = 12.064e-3_r8
!     Initial bottom TIC           (mmol-C /m3):
      t(:,:,:,:, iTIC_               ) = 1469.5_r8
!     Initial bottom total alkalin.(meq    /m3):
      t(:,:,:,:, iTAlk               ) = 1500._r8
!
!     Additional ROMS variables that need to be set if the code is used in a
!       stand-alone context:
      ntimes   = 1048320    ! Number of timesteps to complete (Jan1-Dec31,2014).
      time_ref = 721355._dp ! (days) <=> Jan.1,1975.
      dstart   =  14245._dp ! (days) <=> Jan.1,2014 if time_ref=721355.
      nstp     = 1
      nnew     = 1
      dt       = 30._r8     ! Baroclinic time-step (s).
!
      h               = 5._r8 ! Sea floor depth (m, >0) relative to geoid.
      rmask           = 0._r8
      rmask_wet       = 0._r8
      rmask(    1, 1) = 1._r8
      rmask_wet(1, 1) = 1._r8
      end subroutine initialize_variables
#endif
!
#ifndef ECOOYSTER
      subroutine update_env_cond( iter )
      implicit none
      integer,    intent(in)        :: iter
      integer,    save              :: unum, rec1, rec2
      integer                       :: iost, nlin
      real( dp ), save, allocatable :: envc(:, :)
      real( dp )                    :: line( 12 ), ctim
      real( r8 )                    :: wei2
!
! Overwrite the constant environmental conditions (originally set in
!   `initialize_variables') by time-varying environmental conditions from e.g.
!   a continuous monitoring station. We implicitly assume that the input file
!   is regular and continuous in time (no gap, no repetition) and that all
!   values are physically plausible (no NaN, no -999.).
! There are a total of 12 variables (corresponding to the length of vector
!   `line'):
! 1-Date in Matlab's DATENUM format, 2-temp[deg.C], 3-salt[PSU],
! 4-ph[total scale], 5-alkalinity [umol/kg], 6-TCO2 [umol/kg], 7-pCO2 [uatm],
! 8-Omega_ar, 9-Diss.O2[mg/L], 10-Omega_ca, 11-TSS[mg/L], 12-POC[g-C/m3].
! N.B. We make the assumption that 1umol/kg <=> 1millimol/m3.
!
      if ( iter == 1 ) then
        call get_unit( unum )
!
        open( unit = unum,  form = 'formatted', status = 'old',         &
     &    action = 'read',                                              &
     &    file = '/media/pierre/data/projects/rva_oysters2020/' //      &
     &           'ecooyster/CBP4.3_2014_ecooyster_sara20230713.txt' )
!    &           'ecooyster/gloucester_pt_conmon.txt' )
!    &    file = '/home/pierre/CBP_GB.txt' )
!    &    file = '/media/pierre/data/projects/rva_oysters2020/ecooyster/annual_cycle_le4p3_sara_v20220214.txt' )
        nlin = 0
        do
          read( unit = unum, fmt = *, iostat = iost ) line
          if ( iost == 0 ) then
            nlin = nlin + 1
          else
            exit
          end if
        end do
!
        allocate( envc( nlin, size( line ) ) )
        rewind( unit = unum )
        nlin = 0
        do
          read( unit = unum, fmt = *, iostat = iost ) line
          if ( iost == 0 ) then
            nlin = nlin + 1
            envc( nlin, : ) = line
          else
            exit
          end if
        end do
        close( unit = unum )
!
        write(*,*) 'Completed reading the input file.'
!
        if ( (dstart + time_ref) < envc(1, 1) ) then
          write( 0, * ) 'Error: The initial time (dstart) is earlier '  &
     &      // 'than the first date in the observational input file.'
          stop
        end if
!
!       Compute time at the last time-step:
        ctim = dstart + time_ref + real(ntimes, dp) * real(dt(1), dp)   &
     &       / 60._dp / 60._dp / 24._dp
!
        if ( ctim > envc(size(envc, 1), 1) ) then
          write( 0, * ) 'Error: The number of time-steps (ntimes) is '  &
     &      // 'such that it will exceed the last date in the '         &
     &      // 'observational input file.'
          stop
        end if
!
        rec1 = 1        ! Initialization.
        rec2 = rec1 + 1
      end if
!
!     Compute current time and find the two time records that bracket it:
      ctim = dstart + time_ref + real(iter - 1, dp) * real(dt(1), dp)   &
     &     / 60._dp / 60._dp / 24._dp
!
      do
        if ( envc( rec2, 1 ) < ctim ) then
          rec1 = rec1 + 1
          rec2 = rec2 + 1
        else
          exit
        end if
      end do
!
!     Interpolate linearly the environmental conditions at time `ctim' using a
!       weighted average of the two time records.
      wei2 = real( (ctim - envc( rec1, 1 ))                             &
     &           / (envc( rec2, 1 ) - envc( rec1, 1 )), r8 )
!
      t(:,:,:,:, itemp ) = (1._r8 - wei2) * real( envc( rec1,  2 ), r8 )&
     &                   +          wei2  * real( envc( rec2,  2 ), r8 )
      t(:,:,:,:, isalt ) = (1._r8 - wei2) * real( envc( rec1,  3 ), r8 )&
     &                   +          wei2  * real( envc( rec2,  3 ), r8 )
      t(:,:,:,:, iTAlk ) = (1._r8 - wei2) * real( envc( rec1,  5 ), r8 )&
     &                   +          wei2  * real( envc( rec2,  5 ), r8 )
      t(:,:,:,:, iTIC_ ) = (1._r8 - wei2) * real( envc( rec1,  6 ), r8 )&
     &                   +          wei2  * real( envc( rec2,  6 ), r8 )
!
!     Convert dissolved O2 from g-O2/m3 to millimoles-O2/m3:
      t(:,:,:,:, iOxyg ) = (1._r8 - wei2) * real( envc( rec1,  9 ), r8 )&
     &                   +          wei2  * real( envc( rec2,  9 ), r8 )
      t(:,:,:,:, iOxyg ) = t(:,:,:,:, iOxyg ) / 32.e-3_r8
!
!     Convert POC[g-C/m3] to phytoplankton concentration in millimol-N/m3
!       assuming a Redfield C:N ratio = 106:16:
      t(:,:,:,:, iPhyt ) = (1._r8 - wei2) * real( envc( rec1, 12 ), r8 )&
     &                   +          wei2  * real( envc( rec2, 12 ), r8 )
      t(:,:,:,:, iPhyt ) = t(:,:,:,:, iPhyt ) / 12._r8 * 1.e3_r8        &
     &                   * 16._r8 / 106._r8
!
!     Convert TSS[mg/L] and POC[g-C/m3] to inorganic suspended solids (ISS, in
!       kg/m3) using the relations of Turner et al. 2021:
!         TSS[mg/L] = VSS[mg/L] + FSS[mg/L], VSS[mg/L] = POC[mg-C/L] * 2.9,
!         FSS[mg/L] = 0.35 x VSS[mg/L] + ISS[mg/L],
!     which leads to:
!         ISS[mg/L] = TSS[mg/L] - 1.35 VSS[mg/L]
!                   = TSS[mg/L] - 1.35 x 2.9 x POC[mg-C/L]
      t(:,:,:,:, NAT + NPT + NCS + 1 )                                  &
     &   = ( (1._r8 - wei2) * real( envc( rec1, 11 ), r8 )              &
     &     +          wei2  * real( envc( rec2, 11 ), r8 ) )            &
     &   - ( (1._r8 - wei2) * real( envc( rec1, 12 ), r8 )              &
     &     +          wei2  * real( envc( rec2, 12 ), r8 ) )            &
     &   * 1.35_r8 * 2.9_r8
!     Convert ISS from mg/L to kg/m3:
      t(:,:,:,:, NAT + NPT + NCS + 1 ) = t(:,:,:,:, NAT + NPT + NCS + 1)&
     &                                 * 1.e-3_r8
      end subroutine update_env_cond
#endif
!
#ifndef ECOOYSTER
      subroutine get_unit( unum )
      implicit none
      integer, intent(out) :: unum
      logical              :: isop
!
!     Units (0,5,6) are usually reserved on UNIX-like systems.
      unum = 7
      do
        inquire( unit = unum, opened = isop )
        if ( isop ) then
          unum = unum + 1
        else
          exit
        end if
      end do
      end subroutine get_unit
#endif
!
#ifndef ECOOYSTER
      subroutine write_results( iter )
      implicit none
      integer, intent(in) :: iter
      integer, save       :: unum, nhou
      integer             :: i, j
      real( dp )          :: ocean_time
!
      if ( iter == 1 ) then
        call get_unit( unum )
        open( unit = unum,  form = 'formatted', status = 'replace',     &
     &    action = 'write', file = '/tmp/results_ecooyster.txt' )
!
        nhou = nint( 60._r8 * 60._r8 / dt(1) )
      end if
!
      if (   iter == 1                                                  &
     &  .or. iter == ntimes                                             &
     &  .or. mod( iter, nhou ) == 0 ) then
!       Time corresponding to output, in seconds since `time_ref':
        ocean_time = real( iter, dp ) * real( dt(1), dp )               &
     &             + dstart * 24._dp * 60._dp * 60._dp
        i = 1
        j = 1
        write( unum, * ) ocean_time,                                    &
     &                   DiaBio2d( i, j, iwtid ),                       &
     &                   DiaBio2d( i, j, iwtit ),                       &
     &                   DiaBio2d( i, j, iwshd ),                       &
     &                   DiaBio2d( i, j, iwsht ),                       &
     &                   DiaBio2d( i, j, iwgon ),                       &
     &                   DiaBio2d( i, j, iwgon ) * 0.,                  &
     &                   DiaBio2d( i, j, ihshd ),                       &
     &                   DiaBio2d( i, j, ihsht ),                       &
     &                   DiaBio2d( i, j, iherd ),                       &
     &                   DiaBio2d( i, j, ihert )
      end if
!
      if ( iter == ntimes ) then
        close( unit = unum )
      end if
      end subroutine write_results
#endif
!
#ifndef ECOOYSTER
      subroutine pCO2_water_RZ( Istr, Iend,                             &
     &                    LBi, UBi, LBj, UBj, IminS, ImaxS,             &
     &                    j, DoNewton,                                  &
     &                    rmask,                                        &
     &                    T, S, TIC, TAlk, pH, pCO2)
!***********************************************************************
!                                                                      !
!  This routine computes equilibrium partial pressure of CO2 (pCO2)    !
!  in the surface seawater.                                            !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Istr       Starting tile index in the I-direction.               !
!     Iend       Ending   tile index in the I-direction.               !
!     LBi        I-dimension lower bound.                              !
!     UBi        I-dimension upper bound.                              !
!     LBj        J-dimension lower bound.                              !
!     UBj        J-dimension upper bound.                              !
!     IminS      I-dimension lower bound for private arrays.           !
!     ImaxS      I-dimension upper bound for private arrays.           !
!     j          j-pipelined index.                                    !
!     DoNewton   Iteration solver:                                     !
!                  [0] Bracket and bisection.                          !
!                  [1] Newton-Raphson method.                          !
!     rmask      Land/Sea masking.                                     !
!     T          Surface temperature (Celsius).                        !
!     S          Surface salinity (PSS).                               !
!     TIC        Total inorganic carbon (millimol/m3).                 !
!     TAlk       Total alkalinity (milli-equivalents/m3).              !
!     pH         Best pH guess.                                        !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     pCO2       partial pressure of CO2 (ppmv).                       !
!                                                                      !
!  Check Value:  (T=24, S=36.6, TIC=2040, TAlk=2390, PO4=0,            !
!                 SiO3=0, pH=8)                                        !
!                                                                      !
!                pcO2= ppmv  (DoNewton=0)                              !
!                pCO2= ppmv  (DoNewton=1)                              !
!                                                                      !
!  This subroutine was adapted by Katja Fennel (Nov 2005) from         !
!  Zeebe and Wolf-Gladrow (2001).                                      !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Zeebe, R.E. and D. Wolf-Gladrow,  2005:  CO2 in Seawater:         !
!      Equilibrium, kinetics, isotopes, Elsevier Oceanographic         !
!      Series, 65, pp 346.                                             !
!                                                                      !
!***********************************************************************
      implicit none
!
!  Imported variable declarations.
!
      integer,  intent(in) :: LBi, UBi, LBj, UBj, IminS, ImaxS
      integer,  intent(in) :: Istr, Iend, j, DoNewton
!
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: T(IminS:ImaxS)
      real(r8), intent(in) :: S(IminS:ImaxS)
      real(r8), intent(in) :: TIC(IminS:ImaxS)
      real(r8), intent(in) :: TAlk(IminS:ImaxS)
      real(r8), intent(inout) :: pH(LBi:UBi,LBj:UBj)
!
      real(r8), intent(out) :: pCO2(IminS:ImaxS)
!
!  Local variable declarations.
!
      integer, parameter :: InewtonMax = 10
      integer, parameter :: IbrackMax = 30
!
      integer :: Hstep, Ibrack, Inewton, i
!
      real(r8) :: Tk, centiTk, invTk, logTk
      real(r8) :: scl, sqrtS
      real(r8) :: borate, alk, dic
      real(r8) :: ff, K1, K2, K12, Kb, Kw
      real(r8) :: p5, p4, p3, p2, p1, p0
      real(r8) :: df, fn, fni(3), ftest
      real(r8) :: deltaX, invX, invX2, X, X2, X3
      real(r8) :: pH_guess, pH_hi, pH_lo
      real(r8) :: X_guess, X_hi, X_lo, X_mid
      real(r8) :: CO2star, Htotal, Htotal2
      real(r8) :: pkcw, fhcw
!
!=======================================================================
!  Determine coefficients for surface carbon chemisty.  If land/sea
!  masking, compute only on water points.
!=======================================================================
!
      I_LOOP: DO i=Istr,Iend
        IF (rmask(i,j).gt.0.0_r8) THEN
          Tk=T(i)+273.15_r8
          centiTk=0.01_r8*Tk
          invTk=1.0_r8/Tk
          logTk=LOG(Tk)
          sqrtS=SQRT(S(i))
          scl=S(i)/1.80655_r8
!
          alk= TAlk(i)*0.000001_r8
          dic = TIC(i)*0.000001_r8
!
!-----------------------------------------------------------------------
!  Correction term for non-ideality, ff=k0*(1-pH2O). Equation 13 with
!  table 6 values from Weiss and Price (1980, Mar. Chem., 8, 347-359).
!-----------------------------------------------------------------------
!
          ff=EXP(-162.8301_r8+                                          &
     &         218.2968_r8/centiTk+                                     &
     &         LOG(centiTk)*90.9241_r8-                                 &
     &         centiTk*centiTk*1.47696_r8+                              &
     &         S(i)*(0.025695_r8-                                       &
     &               centiTk*(0.025225_r8-                              &
     &                        centiTk*0.0049867_r8)))
!
!-----------------------------------------------------------------------
!  Compute first (K1) and second (K2) dissociation constant of carbonic
!  acid:
!
!           K1 = [H][HCO3]/[H2CO3]
!           K2 = [H][CO3]/[HCO3]
!
!  Use K1,K2 of Cai and Wang, 1998, Limnology and Oceanography, p. 661.
!  The coefficients are on NBS scale, need to use "fhcw" to convert to SWS.
!  "fhcw" calculation is copied from CO2SYS.m (Takahashi et al, 1982)
!-----------------------------------------------------------------------
!
          fhcw =  1.2948_r8   - 2.036e-3_r8 * Tk                        &
     &         + (4.607e-4_r8 - 1.475e-6_r8 * Tk) * S(i)**2
!
          pkcw = 200.1_r8     * invTk + 0.322_r8
          pkcw = 3404.71_r8   * invTk + 0.032786_r8 * Tk - 14.8435_r8   &
     &         - 7.1692e-2_r8 * pkcw  * sqrtS       + 2.1487e-3_r8 *S(i)
          K1   = 10._r8**(- pkcw) ! This is on NBS scale
          K1   = K1 / fhcw        ! Convert to SWS scale
!
          pkcw = - 129.24_r8 * invTk + 1.4381_r8
          pkcw = 2902.39_r8  * invTk + 0.02379_r8 * Tk - 6.498_r8       &
     &         - 0.3191_r8   * pkcw  * sqrtS      + 0.0198_r8 * S(i)
          K2   = 10._r8**(- pkcw) ! This is on NBS scale
          K2   = K2 / fhcw        ! Convert to SWS scale
!
!-----------------------------------------------------------------------
!  Compute dissociation constant of boric acid, Kb=[H][BO2]/[HBO2].
!  From Millero (1995; page 669) using data from Dickson (1990).
!-----------------------------------------------------------------------
!
          Kb=EXP(-invTk*(8966.90_r8+                                    &
     &                 sqrtS*(2890.53_r8+                               &
     &                        sqrtS*(77.942_r8-                         &
     &                               sqrtS*(1.728_r8-                   &
     &                                      sqrtS*0.0996_r8))))-        &
     &         logTk*(24.4344_r8+                                       &
     &                sqrtS*(25.085_r8+                                 &
     &                       sqrtS*0.2474_r8))+                         &
     &         Tk*(sqrtS*0.053105_r8)+                                  &
     &         148.0248_r8+                                             &
     &         sqrtS*(137.1942_r8+                                      &
     &                sqrtS*1.62142_r8))
!
!-----------------------------------------------------------------------
!  Compute ion product of whater, Kw = [H][OH].
!  From Millero (1995; page 670) using composite data.
!-----------------------------------------------------------------------
!
          Kw=EXP(148.9652_r8-                                           &
     &         invTk*13847.26_r8-                                       &
     &         logTk*23.6521_r8-                                        &
     &         sqrtS*(5.977_r8-                                         &
     &                invTk*118.67_r8-                                  &
     &                logTk*1.0495_r8)-                                 &
     &         S(i)*0.01615_r8)
!
!-----------------------------------------------------------------------
! Calculate concentrations for borate (Uppstrom, 1974).
!-----------------------------------------------------------------------
!
          borate=0.000232_r8*scl/10.811_r8
!
!=======================================================================
!  Iteratively solver for computing hydrogen ions [H+] using either:
!
!    (1) Newton-Raphson method with fixed number of iterations,
!        use previous [H+] as first guess, or
!    (2) bracket and bisection
!=======================================================================
!
!  Solve for h in fifth-order polynomial. First calculate
!  polynomial coefficients.
!
          K12 = K1*K2
!
          p5 = -1.0_r8;
          p4 = -alk-Kb-K1;
          p3 = dic*K1-alk*(Kb+K1)+Kb*borate+Kw-Kb*K1-K12
          p2 = dic*(Kb*K1+2*K12)-alk*(Kb*K1+K12)+Kb*borate*K1           &
     &       +(Kw*Kb+Kw*K1-Kb*K12)
          p1 = 2.0_r8*dic*Kb*K12-alk*Kb*K12+Kb*borate*K12               &
     &       +Kw*Kb*K1+Kw*K12
          p0 = Kw*Kb*K12;
!
!  Set first guess and brackets for [H+] solvers.
!
          pH_guess=pH(i,j)         ! Newton-Raphson
          pH_hi=10.0_r8            ! high bracket/bisection
          pH_lo=5.0_r8             ! low bracket/bisection
!
!  Convert to [H+].
!
          X_guess=10.0_r8**(-pH_guess)
          X_lo=10.0_r8**(-pH_hi)
          X_hi=10.0_r8**(-pH_lo)
          X_mid=0.5_r8*(X_lo+X_hi)
!
!-----------------------------------------------------------------------
!  Newton-Raphson method.
!-----------------------------------------------------------------------
!
          IF (DoNewton.eq.1) THEN
            X=X_guess
!
            DO Inewton=1,InewtonMax
!
!  Evaluate f([H+]) = p5*x^5+...+p1*x+p0
!
              fn=((((p5*X+p4)*X+p3)*X+p2)*X+p1)*X+p0
!
!  Evaluate derivative, df([H+])/dx:
!
!     df= d(fn)/d(X)
!
              df=(((5*p5*X+4*p4)*X+3*p3)*X+2*p2)*X+p1
!
!  Evaluate increment in [H+].
!
              deltaX=-fn/df
!
!  Update estimate of [H+].
!
              X=X+deltaX
            END DO
!
!-----------------------------------------------------------------------
!  Bracket and bisection method.
!-----------------------------------------------------------------------
!
          ELSE
!
!  If first step, use Bracket and Bisection method with fixed, large
!  number of iterations
!
            BRACK_IT: DO Ibrack=1,IbrackMax
              DO Hstep=1,3
                IF (Hstep.eq.1) X=X_hi
                IF (Hstep.eq.2) X=X_lo
                IF (Hstep.eq.3) X=X_mid
!
!  Evaluate f([H+]) for bracketing and mid-value cases.
!
                fni(Hstep)=((((p5*X+p4)*X+p3)*X+p2)*X+p1)*X+p0
              END DO
!
!  Now, bracket solution within two of three.
!
              IF ( abs( fni(3) ) < 1.e-37_r8 ) THEN
!             IF (fni(3).eq.0) THEN
                EXIT BRACK_IT
              ELSE
                ftest=fni(1)/fni(3)
                IF (ftest.gt.0) THEN
                  X_hi=X_mid
                ELSE
                  X_lo=X_mid
                END IF
                X_mid=0.5_r8*(X_lo+X_hi)
              END IF
            END DO BRACK_IT
!
! Last iteration gives value.
!
            X=X_mid
          END IF
!
!-----------------------------------------------------------------------
!  Determine pCO2.
!-----------------------------------------------------------------------
!
!  Total Hydrogen ion concentration, Htotal = [H+].
!
          Htotal=X
          Htotal2=Htotal*Htotal
!
!  Calculate [CO2*] (mole/m3) as defined in DOE Methods Handbook 1994
!  Version 2, ORNL/CDIAC-74, Dickson and Goyet, Eds. (Chapter 2,
!  page 10, Eq A.49).
!
          CO2star=dic*Htotal2/(Htotal2+K1*Htotal+K1*K2)
!
!  Save pH is used again outside this routine.
!
          pH(i,j)=-LOG10(Htotal)
!
!  Add two output arguments for storing pCO2surf.
!
          pCO2(i)=CO2star*1000000.0_r8/ff
        ELSE
          pH(i,j)=0.0_r8
          pCO2(i)=0.0_r8
        END IF
      END DO I_LOOP
      end subroutine pCO2_water_RZ
#endif
!
#ifndef ECOOYSTER
      subroutine caldate( CurrentTime, mm_i )
      implicit none
      real( dp ), intent(in   ) :: CurrentTime
      integer,    intent(  out) :: mm_i
      integer                   :: z, a, c, y
      real( dp )                :: b
!
! Move day 0 from midnight -0001-12-31 to midnight 0000-3-1
      z = floor( CurrentTime + time_ref ) - 60
!
! Calculate number of centuries.
!   K1 = 0.25 is to avoid rounding problems.
      a = floor( (real(z, dp) - 0.25_dp) / 36524.25_dp )
!
! Days within century.
!   K2 = 0.25 is to avoid rounding problems.
      b = real(z, dp) - 0.25_dp + real( a )                             &
     &  - real(floor(real(a, dp) / 4._dp), dp)
!
! Calculate the year (year starts on March 1).
      y = floor(b / 365.25_dp)
!
! Calculate day in year.
      c = int( b - real( floor( 365.25_dp * real(y, dp) ), dp ) ) + 1
!
! Calculate month in year.
      mm_i = int( real(5 * c + 456, dp) / 153._dp )
!
! Move to Jan 1 as start of year.
      if ( mm_i > 12 ) then
        y    = y + 1
        mm_i = mm_i - 12
      end if
      end subroutine caldate
#endif
!
      end module mod_ecooyster
