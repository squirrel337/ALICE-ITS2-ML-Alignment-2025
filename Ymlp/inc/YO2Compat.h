#ifndef ROOT_YO2Compat
#define ROOT_YO2Compat

// Local stand-ins for the handful of O2 header-only helpers that YImpactParameter
// uses. Everything else in the module reaches O2 only through YDetectorGeometry,
// which has its own cache backend, so this header plus that backend is what lets
// the module build and run without O2.
//
// ---------------------------------------------------------------------------
// STATUS
//
// The constants are now copied from the O2 source at tag nightly-20230501, types
// included, and check out to every digit -- see the note on them below. They were
// previously reconstructed from the conventions, and were wrong in the eighth digit.
//
// Angle2Alpha's LOGIC is still reconstructed: it quantises an azimuth to the centre
// of the 20-degree sector containing it, which is what o2::math_utils::detail::
// angle2Alpha (MathUtils/detail/basicMath.h) does, but the branch structure here was
// written from the description rather than transcribed. The constants it uses are
// now O2's, so any remaining disagreement is in the wrapping, not the arithmetic.
//
// This is not cosmetic. YImpactParameter feeds the selection cut at
// YMultiLayerPerceptron.cxx:5666 (|ip_r| > RANGE_IMPACTPARAMS_R, |ip_z| >
// RANGE_IMPACTPARAMS_Z) and the skip at :6640, so getting it wrong changes which
// tracks enter the cost.
//
// 2025 NOTE — kB2C reaches further here than it did in the 2024 module.
// Besides YImpactParameter it sets the seed radius R = pT/q/(kB2C*B) used by
// circle3Dfit_useTPCMomentum, which TrackerFit fixes rather than fits. That radius
// feeds UpdateVertexByAlignment and the p*|d|/40um track-vertex quality gate, so a
// wrong kB2C moves the estimated vertex and changes which tracks and events survive
// -- not merely which impact parameter is reported. Seven call sites in
// YMultiLayerPerceptron.cxx depend on it.
// ---------------------------------------------------------------------------

#include <cmath>
#include <array>

namespace YO2 {

// Copied from o2::constants::math, tag nightly-20230501, INCLUDING the types.
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
constexpr float kPIHalf = 0.5f * kPI;

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
   const int sect = int(phi / kSectorSpanRad + (phi < 0 ? -0.0 : 0.0));
   // sector index counted from -pi so that negative phi lands in the right sector
   const int s = (phi >= 0) ? sect : sect - 1;
   return kSectorSpanRad * (double(s) + 0.5);
}

// o2::math_utils::detail::rotateZ, acting in place on a 3-vector
inline void RotateZ(std::array<double, 3>& v, double alpha)
{
   const double s = std::sin(alpha), c = std::cos(alpha);
   const double x = v[0], y = v[1];
   v[0] = x * c - y * s;
   v[1] = x * s + y * c;
}

inline void SinCos(double a, double& s, double& c) { s = std::sin(a); c = std::cos(a); }

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

#endif
