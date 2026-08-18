// Same 12 numbers per chip, taken from the ROOT-only cache, in the identical
// format. Run in this repository (no O2 needed):
//   root -l -b -q 'dump_L2G_from_cache.C'   # -> cache_L2G.txt
void dump_L2G_from_cache(const char* cache="geometry/its2_geom.root",
                         const char* outFile="cache_L2G.txt"){
   TFile* f=TFile::Open(cache); TTree* t=(TTree*)f->Get("geom");
   Int_t id; Double_t R[9],T[3];
   t->SetBranchAddress("chipID",&id); t->SetBranchAddress("R",R); t->SetBranchAddress("T",T);
   std::vector<std::array<double,12>> v(24120); std::vector<bool> seen(24120,false);
   for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
      for(int k=0;k<9;k++) v[id][k]=R[k];
      for(int k=0;k<3;k++) v[id][9+k]=T[k];
      seen[id]=true; }
   FILE* out=fopen(outFile,"w");
   fprintf(out,"# ITS2 L2G from ROOT-only cache %s\n",cache);
   fprintf(out,"# chipID Rxx Rxy Rxz Ryx Ryy Ryz Rzx Rzy Rzz Tx Ty Tz\n");
   for(int i=0;i<24120;i++){ if(!seen[i]) continue;
      fprintf(out,"%5d",i);
      for(int k=0;k<12;k++) fprintf(out," %+.12e",v[i][k]);
      fprintf(out,"\n"); }
   fclose(out); printf("wrote %s\n",outFile);
}
