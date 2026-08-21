#ifndef ROOT_YO2Compat
#define ROOT_YO2Compat

// The handful of O2 header-only helpers that YImpactParameter uses, in the same two
// backends as YDetectorGeometry and switched by the same guard.
//
//   YGEOM_USE_O2 defined  -- O2 is present. YO2:: is a set of inline one-line
//                            forwards to O2's own helpers, and the two O2 headers
//                            the original module included are included from here,
//                            at the point the original included them. Nothing in
//                            this file computes anything: the arithmetic executed
//                            is O2's, so this backend cannot drift from the
//                            original module however carefully or carelessly the
//                            other backend is written.
//
//   YGEOM_USE_O2 undefined -- the cache backend. No O2 anywhere, so the helpers are
//                            transcribed from the O2 source at tag nightly-20230501
//                            and stand on their own.
//
// The point of the split is that the transcriptions are unreachable when O2 is
// there. A run that must reproduce the original module runs the O2 backend and
// therefore runs the original code; the cache backend is what the transcriptions
// are for, and what they are validated against.
//
// ---------------------------------------------------------------------------
// WHAT THE CACHE BACKEND HAS TO GET RIGHT
//
// Transcribed from the O2 source, not reconstructed from the conventions --
// constants from Common/Constants/include/CommonConstants/MathConstants.h, the
// helpers from Common/MathUtils/include/MathUtils/detail/{trigonometric,basicMath}.h
// and GPU/Common/GPUCommonMath.h. The ARITHMETIC TYPES are part of the answer: O2
// declares its math constants float, and routes some of these helpers through
// GPUCommonMath, which is a single-precision interface even when the caller's type
// is double. Substituting the obvious double-precision equivalent is not neutral.
//
// The sharpest case is ATan2. O2's detail::atan2<T> forwards to
// GPUCommonMath::ATan2(float, float), which has no double overload, so
// atan2<double>(y, x) narrows both arguments to float, calls atan2f, and widens the
// result back. Using std::atan2 instead shifts the track azimuth by 0.04 urad on
// average and up to 0.26 urad -- 0.02 to 0.10 um at r = 40 cm.
//
// This is not cosmetic. YImpactParameter feeds the selection cut at
// YMultiLayerPerceptron.cxx:5666 (|ip_r| > RANGE_IMPACTPARAMS_R, |ip_z| >
// RANGE_IMPACTPARAMS_Z) and the skip at :6640, so getting it wrong changes which
// tracks enter the cost.
//
// Verified bit for bit against the O2 sources above: all nine constants, and
// Angle2Alpha, ATan2, SinCos, RotateZ and Abs over 4,000,001 samples each.
//
// 2025 NOTE - kB2C reaches further here than it did in the 2024 module.
// Besides YImpactParameter it sets the seed radius R = pT/q/(kB2C*B) used by
// circle3Dfit_useTPCMomentum, which TrackerFit fixes rather than fits. That radius
// feeds UpdateVertexByAlignment and the p*|d|/40um track-vertex quality gate, so a
// wrong kB2C moves the estimated vertex and changes which tracks and events survive
// -- not merely which impact parameter is reported. Seven call sites in
// YMultiLayerPerceptron.cxx depend on it.
// ---------------------------------------------------------------------------

#include <cmath>
// o2::gpu::gpustd::array<T, N> is "using array = std::array<T, N>" on the host
// (GPU/Common/GPUCommonArray.h:38), so the call sites spell it std::array in both
// backends -- that is the same type under O2, not a substitution for it.
#include <array>

#ifdef YGEOM_USE_O2

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
//  Cache backend: no O2. Transcribed from the O2 source, see the note above.
// ===========================================================================
namespace YO2 {

// Copied from o2::constants::math, tag nightly-20230501, INCLUDING the types, and
// checked constant by constant against that header: same expressions, same order of
// evaluation, identical bit patterns.
//
// These were previously double, reconstructed from the conventions rather than
// taken from the source, and every one of them was wrong in the eighth digit.
// O2 declares them float and derives PIHalf and the sector span from a float PI
// that is itself written out as an exact float literal, so a double rewrite of
// "the same" constant is a different number:
//
//   B2C           O2 -0.00029979244573041797   double -0.00029979245800000002   4.1e-08
//   PIHalf        O2  1.5707963705062866       double  1.5707963267948966       2.8e-08
//   SectorSpanRad O2  0.34906584024429321      double  0.34906585039886590      2.9e-08
//
// Small, but not nothing: kB2C sets the seed radius that TrackerFit fixes rather
// than fits, so it reaches the estimated vertex and the track-vertex quality gate,
// and a run that must reproduce another bit for bit cannot carry a different
// constant. Float it is -- matching O2 matters more here than looking precise.
constexpr float kPI = 3.14159274101257324e+00f;
constexpr float kTwoPI = 2.f * kPI;
constexpr float kPIHalf = 0.5f * kPI;
constexpr float kRad2Deg = 180.f / kPI;

// converts (kG, GeV/c, cm) to a curvature
constexpr float kB2C = -0.299792458e-3;

// the barrel is divided into 18 sectors of 20 degrees; derived exactly as O2 does
constexpr int   kNSectors = 18;
constexpr float kSectorSpanDeg = 360.f / kNSectors;
constexpr float kDeg2Rad = kPI / 180.f;
constexpr float kSectorSpanRad = kSectorSpanDeg * kDeg2Rad;

// o2::math_utils::detail::angle2Alpha — quantises an azimuth to the centre of the
// sector containing it.
inline double Angle2Alpha(double phi)
{
   // detail::angle2Alpha<T> = sector2Angle<T>(angle2Sector<T>(phi)), inlined here
   // from MathUtils/detail/trigonometric.h. The arithmetic types are part of the
   // answer: SectorSpanRad and (0.5f + sect) are both float, so the product is
   // evaluated in FLOAT and only then widened to T -- computing it in double gives a
   // different number in the eighth digit for 17 sectors out of 18.
   //
   // The negative-phi branch is not a mirror of the positive one: angle2Sector adds
   // NSectors-1 and lets sector2Angle's bringToPMPi wrap the result, where
   // subtracting one sector would differ by 18*SectorSpanRad - TwoPI, not zero in
   // float. The wrap itself is detail::toPMPi: one conditional step, > PI tested
   // before < -PI, not a loop.
   int sect = phi * kRad2Deg / kSectorSpanDeg;
   if (phi < 0.f) {
      sect += kNSectors - 1;
   }
   double ang = kSectorSpanRad * (0.5f + sect);
   if (ang > kPI) {
      ang -= kTwoPI;
   } else if (ang < -kPI) {
      ang += kTwoPI;
   }
   return ang;
}

// o2::math_utils::detail::atan2<double>.
//
// NOT std::atan2. O2's detail::atan2<T> is a one-line forward to
// o2::gpu::GPUCommonMath::ATan2, declared
//
//    GPUhdi() static float ATan2(float y, float x);
//    ... { return CHOICE(atan2f(y, x), atan2f(y, x), atan2(y, x)); }   // host: atan2f
//
// with no double overload anywhere in GPUCommonMath.h. Instantiating it at T = double
// therefore narrows y and x to float, evaluates atan2f in single precision, and widens
// the float result back to double through the template's return type. Reproducing that
// is the whole point of this function: see the note at the top of the file for what
// using the double-precision routine instead costs.
inline double ATan2(double y, double x)
{
   return (double)atan2f((float)y, (float)x);
}

// o2::math_utils::detail::abs<double>, which specialises to GPUCommonMath::Abs<double>
// and thence to fabs. Unlike ATan2 there IS a double overload, so this one is the
// obvious equivalent after all -- spelled as a function so that both backends offer
// the same name and the call sites need no #ifdef.
inline double Abs(double x)
{
   return std::abs(x);
}

// o2::math_utils::detail::sincos<double>, which specialises to
// GPUCommonMath::SinCosd -- glibc's ::sincos on Linux, separate sin and cos elsewhere.
// The two agree bit for bit on this platform (checked over 2,000,001 angles), but the
// point of this header is to transcribe rather than to rely on that.
inline void SinCos(double a, double& s, double& c)
{
#if defined(_GNU_SOURCE) || defined(__GNU_SOURCE__)
   ::sincos(a, &s, &c);
#else
   s = std::sin(a);
   c = std::cos(a);
#endif
}

// o2::math_utils::detail::rotateZ, acting in place on a 3-vector. O2 takes its sine
// and cosine from the same sincos<T> as above, hence the call rather than std::sin and
// std::cos here.
inline void RotateZ(std::array<double, 3>& v, double alpha)
{
   double s, c;
   SinCos(alpha, s, c);
   const double x = v[0], y = v[1];
   v[0] = x * c - y * s;
   v[1] = x * s + y * c;
}

// A sink for O2's LOG(level) macro, used once in YImpactParameter. O2 routes
// LOG(debug) to FairLogger, where debug level is suppressed by default, so
// discarding matches the behaviour the module already had.
struct NullStream {
   template <typename T> NullStream& operator<<(const T&) { return *this; }
};

} // namespace YO2

// ---------------------------------------------------------------------------
//  The names the module actually writes.
//
// YMultiLayerPerceptron.{h,cxx} spell these o2::..., exactly as the original
// module spells them, and are left that way: with O2 those names are O2's own, and
// without O2 they are these. So neither file needs an #ifdef or an edit at any of
// the arithmetic call sites -- and with YGEOM_USE_O2 defined, both preprocess to
// the original module byte for byte.
//
// Signatures follow the O2 declarations exactly, including the explicit template
// argument each call site passes, so that the same source binds under both
// backends. Everything forwards to the YO2:: definitions above.
// ---------------------------------------------------------------------------
namespace o2 {

namespace constants { namespace math {
constexpr float B2C    = YO2::kB2C;
constexpr float PI     = YO2::kPI;
constexpr float TwoPI  = YO2::kTwoPI;
constexpr float PIHalf = YO2::kPIHalf;
} }

// GPU/Common/GPUCommonArray.h:38 -- on the host this is exactly what O2 does.
namespace gpu { namespace gpustd {
template <class T, std::size_t N> using array = std::array<T, N>;
} }

namespace math_utils { namespace detail {

template <class T> T atan2(T y, T x);
template <> inline double atan2<double>(double y, double x) { return YO2::ATan2(y, x); }

template <class T> T angle2Alpha(T phi);
template <> inline double angle2Alpha<double>(double phi) { return YO2::Angle2Alpha(phi); }

template <class T> T abs(T x);
template <> inline double abs<double>(double x) { return YO2::Abs(x); }

// Called without an explicit template argument, and deduced as double at both sites.
inline void sincos(double a, double& s, double& c) { YO2::SinCos(a, s, c); }

template <class T> void rotateZ(gpu::gpustd::array<T, 3>& v, T alpha);
template <> inline void rotateZ<double>(gpu::gpustd::array<double, 3>& v, double alpha)
{
   YO2::RotateZ(v, alpha);
}

} } // namespace math_utils::detail

} // namespace o2

#endif // YGEOM_USE_O2

#ifndef LOG
 #define LOG(level) YO2::NullStream()
#endif

// DetectorConstant.h spells VERTEXFIT and VERTEX_DERIVATIVES as TRUE and FALSE, which
// are not ROOT's -- ROOT has kTRUE and kFALSE and defines neither of these. They reach
// the O2 build through its own header chain, and without O2 nothing defines them, so
// if(VERTEX_DERIVATIVES) at YMultiLayerPerceptron.cxx:6223 would not compile.
//
// Guarded, and here rather than in DetectorConstant.h, so that with O2 present this is
// a no-op and DetectorConstant.h stays byte-identical to the original module: a macro
// body is only expanded where it is used, and every use is in the .cxx, long after this
// header. Spelling them as ROOT's own constants keeps them 0 and 1 either way.
#ifndef FALSE
 #define FALSE kFALSE
#endif
#ifndef TRUE
 #define TRUE kTRUE
#endif

#endif
