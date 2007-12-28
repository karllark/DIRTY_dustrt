#ifndef _DIRTY_DEBUG_
#define _DIRTY_DEBUG_

//**********************************************************************

// debug starts with photon this photon number if set
/* #ifndef DEBUG_PHOTON_NUM */
/* #define DEBUG_PHOTON_NUM 12282 */
/* #endif */

#define OUTNUM -1
//#define OUTNUM 4628428

//  overall debug statment
//#define DEBUG
// if DEBUG set, set all the debug statements on
#ifdef DEBUG
#define DEBUG_DIRTY
#define DEBUG_RT
#define DEBUG_INITOUT
#define DEBUG_SDG
#define DEBUG_NP1S
#define DEBUG_NPDS
#define DEBUG_NPGS
#define DEBUG_FFS
#define DEBUG_NS
#define DEBUG_CPT
#define DEBUG_CDD
#define DEBUG_DPPI
#define DEBUG_CSP
#define DEBUG_CSCP
#define DEBUG_OUTR
#define DEBUG_SWTO
#endif

// debug for main dirty.cc code
#ifndef DEBUG_DIRTY
//#define DEBUG_DIRTY
#endif

// debug for main radiative_transfer.cc code
#ifndef DEBUG_RT
//#define DEBUG_RT
#endif

// debug for initialize_outut.cc code
#ifndef DEBUG_INITOUT
//#define DEBUG_INITOUT
#endif

// debug for setup_dust_grid.cc code
#ifndef DEBUG_SDG
//#define DEBUG_SDG
#endif

// debug for new_photon_1star
#ifndef DEBUG_NP1S
//#define DEBUG_NP1S
#endif

// debug for forced_first_scatter.cc code
#ifndef DEBUG_FFS
//#define DEBUG_FFS
#endif

// debug for classify_stellar_photon.cc code
#ifndef DEBUG_CSP
//#define DEBUG_CSP
#endif

// debug for classify_scattered_photon.cc code
#ifndef DEBUG_CSCP
//#define DEBUG_CSCP
#endif

// debug for next_scatter.cc code
#ifndef DEBUG_NS
//#define DEBUG_NS
#endif

// debug for next_scatter.cc code
#ifndef DEBUG_OUTR
//#define DEBUG_OUTR
#endif

// debug for calc_photon_trajectory.cc code
#ifndef DEBUG_CPT
//#define DEBUG_CPT
#endif

// debug for calc_delta_dist.cc code
#ifndef DEBUG_CDD
//#define DEBUG_CDD
#endif

// debug for determine_photon_position_index.cc code
#ifndef DEBUG_DPPI
//#define DEBUG_DPPI
#endif

// debug for scattered_weight_towards_observer.cc code
#ifndef DEBUG_SWTO
//#define DEBUG_SWTO
#endif

// debug for scattered_weight_towards_observer.cc code
#ifndef DEBUG_RZFO
//#define DEBUG_RZFO
#endif

#endif

