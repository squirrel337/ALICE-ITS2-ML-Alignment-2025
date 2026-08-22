#ifndef ROOT_YO2Compat
#define ROOT_YO2Compat

// The O2 header-only helpers and constants that YImpactParameter and the
// curvature expressions use, in three modes. Two independent guards select them:
//
//   YGEOM_USE_O2         geometry backend: O2's GeometryTGeo, or the ROOT cache
//   YO2_LOCAL_CONSTANTS  take the constants and the curvature arithmetic from
//                        this header instead of from O2
//
// which the run console spells as one setting, GEOM_BACKEND:
//
//   o2        YGEOM_USE_O2                        O2 geometry, O2's own constants,
//                                                 and the original single-precision
//                                                 impact-parameter code. This mode
//                                                 preprocesses to the original-negB
//                                                 module byte for byte -- it is the
//                                                 reference, and nothing here is
//                                                 reachable from it.
//
//   o2-local  YGEOM_USE_O2 + YO2_LOCAL_CONSTANTS  O2 geometry, the constants below,
//                                                 and the impact parameter computed
//                                                 in double. Differs from o2 in the
//                                                 constants and the precision ONLY,
//                                                 with the geometry held fixed, so
//                                                 the pair isolates exactly what the
//                                                 single-precision path costs.
//
//   cache     YO2_LOCAL_CONSTANTS                 no O2 at all: cache geometry, the
//                                                 constants below, double arithmetic.
//
// ---------------------------------------------------------------------------
// WHY THE LOCAL CONSTANTS ARE DOUBLE
//
// O2 declares its math constants float, and routes some helpers through
// GPUCommonMath, which is a single-precision interface even when the caller's type
// is double -- detail::atan2<double> narrows to float, calls atan2f and widens the
// result back. Mode o2 reproduces all of that exactly, because it calls O2.
//
// These do not. They are the same quantities carried at full double precision:
// B2C to its exact decimal rather than the float rounding of it (4.09e-08 apart),
// PI and the sector span likewise, and ATan2 as std::atan2 rather than atan2f.
//
// That is a deliberate difference from O2, not an approximation of it. It exists
// because the impact-parameter path is where single precision bites: getImpactParams
// takes float arguments, computes every intermediate in float and returns float
// ip[2], and those feed the hard cuts |ip_r| > RANGE_IMPACTPARAMS_R and
// |ip_z| > RANGE_IMPACTPARAMS_Z. Float there quantises the result at ~1.2e-07
// relative, so an arbitrarily small upstream difference lands on a different side
// of a cut, a track leaves the sum, and with QUALITY_VERTEXING and
// Num_Of_Bad_Tracks > 2 above it a whole event can drop. Widening the path removes
// that amplifier; it does not, by itself, remove whatever feeds it.
//
// Use o2 when the answer must match the original module. Use o2-local to measure
// what the single-precision path is worth, with everything else identical.
// ---------------------------------------------------------------------------

#include <cmath>
// o2::gpu::gpustd::array<T, N> is "using array = std::array<T, N>" on the host
// (GPU/Common/GPUCommonArray.h:38), so the call sites spell it std::array in every
// mode -- that is the same type under O2, not a substitution for it.
#include <array>

#if defined(YGEOM_USE_O2) && !defined(YO2_LOCAL_CONSTANTS)


// ===========================================================================
//  O2 present: forward to O2. Same two headers the original module included.
// ===========================================================================
#include "CommonConstants/MathConstants.h"
#include "MathUtils/Utils.h"

namespace YO2 {

constexpr float kPI            = o2::constants::math::PI;
constexpr float kTwoPI         = o2::constants::math::TwoPI;
constexpr float kPIHalf        = o2::constants::math::PIHalf;
constexpr float kRad2Deg       = o2::constants::math::Rad2Deg;
constexpr float kDeg2Rad       = o2::constants::math::Deg2Rad;
constexpr float kB2C           = o2::constants::math::B2C;
constexpr int   kNSectors      = o2::constants::math::NSectors;
constexpr float kSectorSpanDeg = o2::constants::math::SectorSpanDeg;
constexpr float kSectorSpanRad = o2::constants::math::SectorSpanRad;

inline double ATan2(double y, double x)
{
   return o2::math_utils::detail::atan2<double>(y, x);
}

inline double Angle2Alpha(double phi)
{
   return o2::math_utils::detail::angle2Alpha<double>(phi);
}

inline double Abs(double x)
{
   return o2::math_utils::detail::abs<double>(x);
}

inline void SinCos(double a, double& s, double& c)
{
   o2::math_utils::detail::sincos(a, s, c);
}

inline void RotateZ(std::array<double, 3>& v, double alpha)
{
   o2::math_utils::detail::rotateZ<double>(v, alpha);
}

} // namespace YO2
#else

// ===========================================================================
//  Local constants, in double. Modes o2-local and cache.
// ===========================================================================
namespace YO2 {

// The same quantities O2 names in CommonConstants/MathConstants.h, carried at
// double precision instead of float. Every one of these is within 4.1e-08 of the
// value O2 uses; see the note at the top of this file for why that is on purpose.
constexpr double kPI           = 3.14159265358979323846;
constexpr double kTwoPI        = 2.0 * kPI;
constexpr double kPIHalf       = 0.5 * kPI;
constexpr double kRad2Deg      = 180.0 / kPI;
constexpr double kDeg2Rad      = kPI / 180.0;

// converts (kG, GeV/c, cm) to a curvature
constexpr double kB2C          = -0.299792458e-3;

// the barrel is divided into 18 sectors of 20 degrees
constexpr int    kNSectors     = 18;
constexpr double kSectorSpanDeg = 360.0 / kNSectors;
constexpr double kSectorSpanRad = kSectorSpanDeg * kDeg2Rad;

// o2::math_utils::detail::atan2, at full precision. O2's forwards to
// GPUCommonMath::ATan2(float, float) and so rounds to float; this does not.
inline double ATan2(double y, double x)
{
   return std::atan2(y, x);
}

// o2::math_utils::detail::abs
inline double Abs(double x)
{
   return std::abs(x);
}

// o2::math_utils::detail::sincos
inline void SinCos(double a, double& s, double& c)
{
#if defined(_GNU_SOURCE) || defined(__GNU_SOURCE__)
   ::sincos(a, &s, &c);
#else
   s = std::sin(a);
   c = std::cos(a);
#endif
}

// o2::math_utils::detail::angle2Alpha = sector2Angle(angle2Sector(phi)) -- O2's
// algorithm, evaluated in double. O2 computes SectorSpanRad * (0.5f + sect) in
// float and only then widens; here the whole expression stays double.
inline double Angle2Alpha(double phi)
{
   int sect = phi * kRad2Deg / kSectorSpanDeg;
   if (phi < 0.0) {
      sect += kNSectors - 1;
   }
   double ang = kSectorSpanRad * (0.5 + sect);
   if (ang > kPI) {
      ang -= kTwoPI;
   } else if (ang < -kPI) {
      ang += kTwoPI;
   }
   return ang;
}

// o2::math_utils::detail::rotateZ, acting in place on a 3-vector
inline void RotateZ(std::array<double, 3>& v, double alpha)
{
   double s, c;
   SinCos(alpha, s, c);
   const double x = v[0], y = v[1];
   v[0] = x * c - y * s;
   v[1] = x * s + y * c;
}

} // namespace YO2

#endif // O2 constants vs local constants

// ---------------------------------------------------------------------------
// Without O2 these three names have no other source. With O2 present -- modes o2
// and o2-local alike -- they come from its header chain exactly as they did in the
// original module, so the guards below are no-ops there.
// ---------------------------------------------------------------------------
#ifndef YGEOM_USE_O2

namespace YO2 {
// A sink for O2's LOG(level) macro, used once in YImpactParameter. O2 routes
// LOG(debug) to FairLogger, where debug level is suppressed by default, so
// discarding matches the behaviour the module already had.
struct NullStream {
   template <typename T> NullStream& operator<<(const T&) { return *this; }
};
} // namespace YO2

#ifndef LOG
 #define LOG(level) YO2::NullStream()
#endif

// DetectorConstant.h spells VERTEXFIT and VERTEX_DERIVATIVES as TRUE and FALSE,
// which are not ROOT's -- ROOT has kTRUE and kFALSE and defines neither. Spelling
// them as ROOT's own constants keeps them 0 and 1, and putting them here rather
// than in DetectorConstant.h leaves that file byte-identical to the original: a
// macro body is only expanded where it is used, and every use is in the .cxx.
#ifndef FALSE
 #define FALSE kFALSE
#endif
#ifndef TRUE
 #define TRUE kTRUE
#endif

#endif // !YGEOM_USE_O2

#endif
