// ---------------------------------------------------------------------------
// Run this in an O2 environment. It writes the one thing needed to validate the
// ROOT-only geometry cache: the local-to-global transform of every ITS chip.
//
//   alienv enter O2/latest
//   root -l -b -q 'dump_L2G_from_O2.C'      # -> o2_L2G.txt
//
// Point it at the SAME o2sim_geometry.root and ITSAlignment.root the module uses,
// so the comparison isolates how the matrix is extracted, not which files fed it.
// ---------------------------------------------------------------------------
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "ITSBase/GeometryTGeo.h"
#include <TFile.h>
#include <cstdio>
#include <vector>

void dump_L2G_from_O2(const char* geomFile  = "o2sim_geometry.root",
                      const char* alignFile = "ITSAlignment.root",
                      const char* outFile   = "o2_L2G.txt")
{
   o2::base::GeometryManager::loadGeometry("", false, false);

   TFile fa(alignFile);
   std::vector<o2::detectors::AlignParam>* ali = nullptr;
   fa.GetObject("ccdb_object", ali);
   if (!ali) { ::Error("dump", "no ccdb_object in %s", alignFile); return; }
   o2::base::GeometryManager::applyAlignment(*ali);
   printf("applied %zu AlignParam\n", ali->size());

   auto* geom = o2::its::GeometryTGeo::Instance();
   geom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L,
                                                  o2::math_utils::TransformType::T2GRot,
                                                  o2::math_utils::TransformType::L2G));

   FILE* out = fopen(outFile, "w");
   fprintf(out, "# ITS2 L2G from O2  geometry=%s alignment=%s\n", geomFile, alignFile);
   fprintf(out, "# chipID Rxx Rxy Rxz Ryx Ryy Ryz Rzx Rzy Rzz Tx Ty Tz\n");
   for (int id = 0; id < 24120; ++id) {
      double Rxx,Rxy,Rxz,Tdx, Ryx,Ryy,Ryz,Tdy, Rzx,Rzy,Rzz,Tdz;
      geom->getMatrixL2G(id).GetComponents(Rxx,Rxy,Rxz,Tdx, Ryx,Ryy,Ryz,Tdy, Rzx,Rzy,Rzz,Tdz);
      fprintf(out, "%5d %+.12e %+.12e %+.12e %+.12e %+.12e %+.12e %+.12e %+.12e %+.12e %+.12e %+.12e %+.12e\n",
              id, Rxx,Rxy,Rxz, Ryx,Ryy,Ryz, Rzx,Rzy,Rzz, Tdx,Tdy,Tdz);
   }
   fclose(out);
   printf("wrote %s (24120 chips)\n", outFile);
}
