#!/bin/bash
# Run the supplied monitoring macros on one epoch of a run.
#
#   tools/monitoring/make_monitoring_plots.sh <epoch> [outdir] [residual-dir]
#
# The macros address the tree by the global names TrkVtxer and ResMonitor. Both are pointed at
# the TrkVtxer tree here: that is the tree TrackVertexQualityEstimator writes from TrackerFit,
# the fit that uses the clusters alone. ResMonitor -- the tree from the refit that carries the
# primary vertex as an eighth point -- is deliberately not used, for either the DCA or the
# residuals.
set -e
M=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)   # this directory
R=$(cd "$M/../.." && pwd)                         # repository root
EP=${1:--1}
OUT=${2:-$M/plots}
DIR=${3:-$R/MLPTrain_Step901/Residual}
M=$M/user_macros
SCALE=2.0                     # the +-2 sigma pass of the supplied fitter
mkdir -p "$OUT"

root -l -b <<EOF 2>&1 | grep -vE "warning|note:|In file included|declared here|^ *\^|VLA|variable length|input_line|^Info in <TCanvas::Print>: png"
gStyle->SetOptStat(0);
gStyle->SetPadLeftMargin(0.15); gStyle->SetPadBottomMargin(0.13);
gStyle->SetTitleOffset(1.3,"Y");
TChain* TrkVtxer = new TChain("TrkVtxer");
TrkVtxer->Add("$DIR/Residual_Monitor_Epoch_At_${EP}_Train.root");
TrkVtxer->Add("$DIR/Residual_Monitor_Epoch_At_${EP}_Test.root");
TChain* ResMonitor = TrkVtxer;          // residuals also come from the cluster-only fit
printf(">>> epoch ${EP}: %lld tracks\n", TrkVtxer->GetEntries());

.L $M/check_vertex_pT_plots_trkvtxer_color.C
.L $M/check_vertex_eta_plots_trkvtxer_color.C
.L $M/check_vertex_z_plots_trkvtxer_color.C
.L $M/check_vertex_phi_plots_trkvtxer_color.C
.L $M/check_track_residualsHB_FitGaus.C

auto dump=[](const char* what){
   TH1D* wy=(TH1D*)gROOT->FindObject("h_width_dcaY"); TH1D* wz=(TH1D*)gROOT->FindObject("h_width_dcaZ");
   TH1D* my=(TH1D*)gROOT->FindObject("h_mean_dcaY");  TH1D* mz=(TH1D*)gROOT->FindObject("h_mean_dcaZ");
   if(!wy||!wz||!my||!mz) return;
   printf("\n### DCA vs %s  [um]\n%10s %10s %10s %10s %10s\n",what,"x","mean_y","width_y","mean_z","width_z");
   for(int i=1;i<=wy->GetNbinsX();++i){
      if(wy->GetBinContent(i)==0 && wz->GetBinContent(i)==0) continue;
      printf("%10.3f %10.2f %10.2f %10.2f %10.2f\n",wy->GetBinCenter(i),
             my->GetBinContent(i),wy->GetBinContent(i),mz->GetBinContent(i),wz->GetBinContent(i));
   }
};
auto save=[](const char* cname, const char* fname){
   TCanvas* c=(TCanvas*)gROOT->GetListOfCanvases()->FindObject(cname);
   if(!c){ printf("!! canvas %s missing\n",cname); return; }
   c->SaveAs(fname);
};
// The supplied fitter fits every x bin, including empty ones past the detector edge, and those
// come back with errors larger than the pad. Blank any bin whose parent projection is too thin.
auto prune=[](const char* mname, const char* hname, double minN){
   TH1D* m=(TH1D*)gROOT->FindObject(mname); TH2D* h=(TH2D*)gROOT->FindObject(hname);
   if(!m||!h) return;
   for(int i=1;i<=m->GetNbinsX();++i)
      if(h->Integral(i,i,1,h->GetNbinsY()) < minN){ m->SetBinContent(i,0); m->SetBinError(i,0); }
};
// the residual macro overlays all/HB0/HB1 in one colour; separate them after the fact
auto colourHB=[&](const char* stem){
   const char* sfx[3]={"","_HB0","_HB1"};
   int col[3]={kBlack,kAzure+2,kOrange+7}, mst[3]={20,21,22};
   for(int l=0;l<7;l++) for(int k=0;k<3;k++){
      prune(Form("m%s_%d%s",stem,l,sfx[k]), Form("%s_%d%s",stem,l,sfx[k]), 200);
      TH1D* a=(TH1D*)gROOT->FindObject(Form("m%s_%d%s",stem,l,sfx[k]));
      if(!a) continue;
      a->SetLineColor(col[k]); a->SetMarkerColor(col[k]);
      a->SetLineWidth(2); a->SetMarkerStyle(mst[k]); a->SetMarkerSize(0.8);
   }
};
auto colourPhi=[&](const char* stem){
   for(int l=0;l<7;l++){
      prune(Form("m%s_%d",stem,l), Form("%s_%d",stem,l), 200);
      TH1D* a=(TH1D*)gROOT->FindObject(Form("m%s_%d",stem,l));
      if(!a) continue;
      a->SetLineColor(kAzure+2); a->SetMarkerColor(kAzure+2);
      a->SetLineWidth(2); a->SetMarkerStyle(20); a->SetMarkerSize(0.8);
   }
};
auto keyHB=[](const char* cname){
   TCanvas* c=(TCanvas*)gROOT->GetListOfCanvases()->FindObject(cname); if(!c) return;
   c->cd(4);
   TLegend* l=new TLegend(0.05,0.35,0.95,0.65);
   TH1D* a=new TH1D(Form("k0_%s",cname),"",1,0,1); a->SetLineColor(kBlack);    a->SetMarkerColor(kBlack);    a->SetMarkerStyle(20); a->SetLineWidth(2);
   TH1D* b=new TH1D(Form("k1_%s",cname),"",1,0,1); b->SetLineColor(kAzure+2);  b->SetMarkerColor(kAzure+2);  b->SetMarkerStyle(21); b->SetLineWidth(2);
   TH1D* d=new TH1D(Form("k2_%s",cname),"",1,0,1); d->SetLineColor(kOrange+7); d->SetMarkerColor(kOrange+7); d->SetMarkerStyle(22); d->SetLineWidth(2);
   l->AddEntry(a,"both half-barrels","lp"); l->AddEntry(b,"HB0  (#phi > 0)","lp"); l->AddEntry(d,"HB1  (#phi < 0)","lp");
   l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.075); l->Draw();
};

// ---- DCA_y and DCA_z against pT, eta, z, phi --------------------------------------
check_vertex_pT_plots_trkvtxer_color(kAzure+2,21,$SCALE);
dump("pT [GeV/c]"); save("canvas_dca_pT","$OUT/dca_vs_pT_ep${EP}.png");
check_vertex_eta_plots_trkvtxer_color(kAzure+2,21,$SCALE);
dump("eta"); save("canvas_dca_eta","$OUT/dca_vs_eta_ep${EP}.png");
check_vertex_z_plots_trkvtxer_color(kAzure+2,21,$SCALE);
dump("z [cm]"); save("canvas_dca_vtxevtZ","$OUT/dca_vs_z_ep${EP}.png");
check_vertex_phi_plots_trkvtxer_color(kAzure+2,21,$SCALE);
dump("phi [rad]"); save("canvas_dca_phi","$OUT/dca_vs_phi_ep${EP}.png");

// ---- residuals per layer against phi and z -----------------------------------------
check_track_residualsHB_FitGaus(1,$SCALE);
colourHB("hds1_vs_z");    colourHB("hds2_vs_z");
colourPhi("hds1_vs_phi"); colourPhi("hds2_vs_phi");
keyHB("canvas_ds1_vs_z"); keyHB("canvas_ds2_vs_z");
for(auto o : *gROOT->GetListOfCanvases()){ ((TCanvas*)o)->Modified(); ((TCanvas*)o)->Update(); }
save("canvas_ds2_vs_phi","$OUT/residualXY_vs_phi_ep${EP}.png");
save("canvas_ds2_vs_z",  "$OUT/residualXY_vs_z_ep${EP}.png");
save("canvas_ds1_vs_phi","$OUT/residualZ_vs_phi_ep${EP}.png");
save("canvas_ds1_vs_z",  "$OUT/residualZ_vs_z_ep${EP}.png");
.q
EOF
