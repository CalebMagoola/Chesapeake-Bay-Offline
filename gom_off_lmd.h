/*
**
** Options for the northern Gulf of Mexico Application.
**
*/

# define ROMS_MODEL
# define AVERAGES

# define MASKING
# define CURVGRID

/* OFFLINE */
# define OFFLINE
# define OFFLINE_BIOLOGY
# define OUT_DOUBLE
# define PERFECT_RESTART
# define FENNEL

# define OCLIMATOLOGY
# define ATCLIMATOLOGY
# undef AKXCLIMATOLOGY
# undef MIXCLIMATOLOGY

/* BASIC NUMERIC NONLINEAR */
# define UV_ADV
# define DJ_GRADPS
# define UV_COR
# define UV_LOGDRAG
# define UV_VIS2
# define MIX_S_UV
# define DIFF_GRID
# define VISC_GRID

# define SPLINES_VDIFF
# define SPLINES_VVISC

# define NONLIN_EOS
# define SOLVE3D

/* MIXING */
# define TS_MPDATA
# define TS_DIF2
# define MIX_GEO_TS

# undef MY25_MIXING
#ifdef MY25_MIXING
# define N2S2_HORAVG
# define RI_SPLINES
# define KANTHA_CLAYSON
#endif

# undef  GLS_MIXING
#ifdef GLS_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# define RI_SPLINES
#endif

# define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
# define RI_SPLINES
#endif

/* SURFACE AND BOTTOM */
# undef SOLAR_SOURCE
# undef DIURNAL_SRFLUX
# undef QCORRECTION
# undef BULK_FLUXES
# undef LONGWAVE

# define SALINITY
# define ANA_SSFLUX 

# undef ANA_CLOUD
# define ANA_SPFLUX
# define ANA_BPFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# undef ANA_M2OBC
# undef ANA_FSOBC

/* BIOLOGY */
# define BIOLOGY
# define BIO_FENNEL
# define OXYGEN
# define RW14_OXYGEN_SC
# define BIO_SEDIMENT
# define PO4
# define DENITRIFICATION
# define RIVER_BIOLOGY
# define RIVER_DON
