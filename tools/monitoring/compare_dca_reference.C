#include <fstream>
// ITS2 alignment monitoring.
//  - DCA_y / DCA_z profiled against pT, eta, z, phi, against BOTH vertex references
//  - residuals and DCA summarised by an iterative Gaussian fit
// Widths come from a fit, not from GetRMS: the initial window is the full width at
// half maximum, then the fit is repeated inside +-2 sigma of its own result until
// the sigma stops moving.
//
// Two references are kept apart on purpose.
//   v_est  = fvertex_TRKF, the adaptive vertex (tree: vtxevt*). This is what the module
//            itself measures the DCA against, and it is fitted from the very tracks being
//            measured -- and re-enters the circle fit as an 8th point carrying
//            Sigma_MEAS[7] = 4.74 um. The width therefore collapses onto that sigma and
//            loses its 1/p dependence: it is a fit residual, not an impact parameter.
//   v_reco = tv1,tv2,tv3 from the input file (tree: vtx*). Independent of the track, so
//            the width keeps the multiple-scattering 1/p behaviour.
// The transverse DCA is recomputed here from the fitted circle (curvX,curvY,curvR) so
// the reference can be swapped without re-running the module; the longitudinal one
// shifts rigidly with the reference z.
#include "../../Ymlp/inc/YMultiLayerPerceptron.h"

struct Fit { double mean, sigma, emean, esigma; bool ok; };

Fit gfit(TH1* h, int maxIter = 8, double tol = 1e-3)
{
   Fit r{0,0,0,0,false};
   if (!h || h->GetEntries() < 50) return r;
   // FWHM window around the peak
   int pb = h->GetMaximumBin();
   double pk = h->GetBinContent(pb), half = 0.5*pk;
   int lo = pb, hi = pb;
   while (lo > 1              && h->GetBinContent(lo) > half) --lo;
   while (hi < h->GetNbinsX() && h->GetBinContent(hi) > half) ++hi;
   double x1 = h->GetBinCenter(lo), x2 = h->GetBinCenter(hi);
   if (!(x2 > x1)) { x1 = h->GetMean()-h->GetRMS(); x2 = h->GetMean()+h->GetRMS(); }
   TF1 f("f","gaus",x1,x2);
   f.SetParameters(pk, h->GetBinCenter(pb), 0.5*(x2-x1)/1.177);  // FWHM = 2.355 sigma
   if (h->Fit(&f,"QNR") != 0) return r;
   double prev = f.GetParameter(2);
   for (int it = 0; it < maxIter; ++it) {
      double m = f.GetParameter(1), s = std::fabs(f.GetParameter(2));
      if (s <= 0) break;
      f.SetRange(m-2*s, m+2*s);
      if (h->Fit(&f,"QNR") != 0) break;
      double now = std::fabs(f.GetParameter(2));
      if (std::fabs(now-prev) < tol*std::max(now,1e-12)) { prev = now; break; }
      prev = now;
   }
   r.mean=f.GetParameter(1); r.sigma=std::fabs(f.GetParameter(2));
   r.emean=f.GetParError(1); r.esigma=f.GetParError(2); r.ok=true;
   return r;
}

// transverse distance of closest approach of the fitted circle to an arbitrary point
static double d0xy(double cx, double cy, double R, double vx, double vy)
{
   return std::fabs(R) - std::sqrt((cx-vx)*(cx-vx) + (cy-vy)*(cy-vy));
}

// one differential set: a variable binned, each bin holding a DCA distribution
struct Diff {
   TH1D** hy; TH1D** hz; int n; double lo, hi; TString name, unit;
   void book(const char* nm, const char* un, const char* tag, int nb, double a, double b,
             double rngY, double zctr, double rngZ){
      name=nm; unit=un; n=nb; lo=a; hi=b;
      hy=new TH1D*[n]; hz=new TH1D*[n];
      for(int i=0;i<n;i++){
         hy[i]=new TH1D(Form("dy_%s%s_%d",nm,tag,i),"",120,-rngY,rngY);
         hz[i]=new TH1D(Form("dz_%s%s_%d",nm,tag,i),"",120,zctr-rngZ,zctr+rngZ);
      }
   }
   int bin(double v) const { if(v<lo||v>=hi) return -1; return int((v-lo)/(hi-lo)*n); }
   void fill(double v,double dy,double dz){ int b=bin(v); if(b<0) return; hy[b]->Fill(dy); hz[b]->Fill(dz); }
};

void draw(Diff& d, const char* out, const char* tag, const char* ref)
{
   std::vector<double> x, my, sy, mz, sz, ey, ez, esy, esz;
   for(int i=0;i<d.n;i++){
      Fit fy=gfit(d.hy[i]), fz=gfit(d.hz[i]);
      if(!fy.ok||!fz.ok) continue;
      x.push_back(d.lo+(i+0.5)*(d.hi-d.lo)/d.n);
      my.push_back(1e4*fy.mean); ey.push_back(1e4*fy.emean); sy.push_back(1e4*fy.sigma); esy.push_back(1e4*fy.esigma);
      mz.push_back(1e4*fz.mean); ez.push_back(1e4*fz.emean); sz.push_back(1e4*fz.sigma); esz.push_back(1e4*fz.esigma);
   }
   if(x.size()<2) return;
   TCanvas* c=new TCanvas(Form("c_%s%s",d.name.Data(),tag),"",1200,520); c->Divide(2,1);
   auto mk=[&](std::vector<double>& v,std::vector<double>& e,int col,int mst){
      TGraphErrors* g=new TGraphErrors(x.size(),&x[0],&v[0],nullptr,&e[0]);
      g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(mst);
      g->SetMarkerSize(1.2); g->SetLineWidth(2); return g; };
   c->cd(1); gPad->SetGridy(); gPad->SetLeftMargin(0.15);
   TGraphErrors* g1=mk(my,ey,kAzure+2,20); TGraphErrors* g2=mk(mz,ez,kOrange+7,21);
   TMultiGraph* m1=new TMultiGraph(); m1->Add(g1,"lp"); m1->Add(g2,"lp"); m1->Draw("A");
   m1->GetXaxis()->SetTitle(Form("%s %s",d.name.Data(),d.unit.Data()));
   m1->GetYaxis()->SetTitle(Form("DCA mean vs %s [#mum]",ref)); m1->GetYaxis()->SetTitleOffset(1.5);
   TLegend* l1=new TLegend(0.62,0.76,0.9,0.9); l1->AddEntry(g1,"DCA_{y}","lp");
   l1->AddEntry(g2,"DCA_{z}","lp"); l1->SetBorderSize(0); l1->SetFillStyle(0); l1->Draw();
   c->cd(2); gPad->SetGridy(); gPad->SetLeftMargin(0.15);
   TGraphErrors* g3=mk(sy,esy,kAzure+2,20); TGraphErrors* g4=mk(sz,esz,kOrange+7,21);
   TMultiGraph* m2=new TMultiGraph(); m2->Add(g3,"lp"); m2->Add(g4,"lp"); m2->Draw("A");
   m2->GetXaxis()->SetTitle(Form("%s %s",d.name.Data(),d.unit.Data()));
   m2->GetYaxis()->SetTitle(Form("DCA width vs %s [#mum]",ref)); m2->GetYaxis()->SetTitleOffset(1.5);
   TLegend* l2=new TLegend(0.62,0.76,0.9,0.9); l2->AddEntry(g3,"DCA_{y}","lp");
   l2->AddEntry(g4,"DCA_{z}","lp"); l2->SetBorderSize(0); l2->SetFillStyle(0); l2->Draw();
   c->SaveAs(Form("%s/dca_vs_%s%s.png",out,d.name.Data(),tag));
   printf("\n  [%s] %-6s %10s %10s %10s %10s\n",ref,d.name.Data(),"mean_y","width_y","mean_z","width_z");
   for(size_t i=0;i<x.size();i++)
      printf("  %6.2f %10.2f %10.2f %10.2f %10.2f\n",x[i],my[i],sy[i],mz[i],sz[i]);
}

void compare_dca_reference(const char* dir="MLPTrain_Step901/Residual", const char* out="tools/monitoring/plots", int FINAL=4)
{
   gSystem->mkdir(out, kTRUE);
   gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   gStyle->SetPadLeftMargin(0.15); gStyle->SetPadBottomMargin(0.14);

   // vs v_est (what the module reports) and vs v_reco (independent of the track)
   Diff ept, eet, ez, eph;      // v_est
   Diff rpt, ret, rz, rph;      // v_reco
   ept.book("pT",  "[GeV/c]","_vest", 12, 0.0,  6.0,          0.02, 0.0, 0.02);
   eet.book("eta", "",       "_vest", 20, -1.0, 1.0,          0.02, 0.0, 0.02);
   ez .book("zvtx","[cm]",   "_vest", 60, -15.0,15.0,         0.02, 0.0, 0.02);
   eph.book("phi", "[rad]",  "_vest", 12, 0.0, 2*TMath::Pi(), 0.02, 0.0, 0.02);
   // v_reco: dca_z carries the known ~ -1.8 mm vertex-z discrepancy, so centre its window
   rpt.book("pT",  "[GeV/c]","_vreco",12, 0.0,  6.0,          0.05, -0.18, 0.08);
   ret.book("eta", "",       "_vreco",20, -1.0, 1.0,          0.05, -0.18, 0.08);
   rz .book("zvtx","[cm]",   "_vreco",60, -15.0,15.0,         0.05, -0.18, 0.08);
   rph.book("phi", "[rad]",  "_vreco",12, 0.0, 2*TMath::Pi(), 0.05, -0.18, 0.08);

   // residuals: coarse binning, one histogram per layer
   TH1D *r1[8], *r2[8], *r1f[8], *r2f[8];
   for(int l=0;l<8;l++){
      double rng = (l<3||l==7) ? 100 : 500;              // um
      r1[l] =new TH1D(Form("r1_%d",l), "",60,-rng,rng);
      r2[l] =new TH1D(Form("r2_%d",l), "",60,-rng,rng);
      r1f[l]=new TH1D(Form("r1f_%d",l),"",60,-rng,rng);
      r2f[l]=new TH1D(Form("r2f_%d",l),"",60,-rng,rng);
   }
   TH1D* hdvz=new TH1D("hdvz","",120,-0.4,0.1);
   TH1D* hipy=new TH1D("hipy","",120,-0.02,0.02),  *hipz=new TH1D("hipz","",120,-0.02,0.02);
   TH1D* hrpy=new TH1D("hrpy","",120,-0.05,0.05),  *hrpz=new TH1D("hrpz","",120,-0.26,-0.10);
   TProfile* pdvz_eta=new TProfile("pdvz_eta","",20,-1,1);
   TProfile* pdvz_z  =new TProfile("pdvz_z","",30,-15,15);
   TProfile* pdvz_pt =new TProfile("pdvz_pt","",12,0,6);
   TProfile* pdvz_ph =new TProfile("pdvz_ph","",12,0,2*TMath::Pi());

   auto scan=[&](int epoch, bool base){
      for(const char* sfx : {"_Train","_Test"}){
         TFile* f=TFile::Open(Form("%s/Residual_Monitor_Epoch_At_%d%s.root",dir,epoch,sfx));
         if(!f||f->IsZombie()) continue;
         TTree* t=(TTree*)f->Get("ResMonitor"); if(!t){f->Close();continue;}
         YResidualMonitor* m=new YResidualMonitor(); t->SetBranchAddress("monitor",&m);
         for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
            for(int l=0;l<8;l++){
               if(l<7 && m->fchipID[l]<0) continue;
               double a=1e4*m->fds1[l], b=1e4*m->fds2[l];
               if(std::abs(a)<1e3){ (base?r1[l]:r1f[l])->Fill(a); }
               if(std::abs(b)<1e3){ (base?r2[l]:r2f[l])->Fill(b); }
            }
            if(!base) continue;
            double dvz = m->vtxevtZ - m->vtxZ;
            // vs v_est: as reported by the module
            double ey = m->fip[0], ezz = m->fip[1];
            // vs v_reco: transverse recomputed from the fitted circle, longitudinal
            // shifted rigidly by the reference-z difference
            double ry = d0xy(m->curvX, m->curvY, m->curvR, m->vtxX, m->vtxY);
            double rzz = ezz + dvz;
            hipy->Fill(ey); hipz->Fill(ezz);
            hrpy->Fill(ry); hrpz->Fill(rzz);
            double phi=m->phi; if(phi<0) phi+=2*TMath::Pi();
            ept.fill(m->pT,ey,ezz);  eet.fill(m->eta,ey,ezz);
            ez .fill(m->vtxZ,ey,ezz);eph.fill(phi,ey,ezz);
            rpt.fill(m->pT,ry,rzz);  ret.fill(m->eta,ry,rzz);
            rz .fill(m->vtxZ,ry,rzz);rph.fill(phi,ry,rzz);
            hdvz->Fill(dvz);
            pdvz_eta->Fill(m->eta,1e4*dvz); pdvz_z->Fill(m->vtxZ,1e4*dvz);
            pdvz_pt->Fill(m->pT,1e4*dvz);   pdvz_ph->Fill(phi,1e4*dvz);
         }
         f->Close();
      }
   };
   scan(-1,true); scan(FINAL,false);

   const char* LAY[8]={"L0","L1","L2","L3","L4","L5","L6","VTX"};
   printf("\n=== RESIDUALS : iterative +-2sigma gaussian fit [um] ===\n");
   printf("%-5s | %9s %9s %9s %9s | %9s %9s %9s %9s\n","lay",
          "m1(-1)",Form("m1(%d)",FINAL),"w1(-1)",Form("w1(%d)",FINAL),
          "m2(-1)",Form("m2(%d)",FINAL),"w2(-1)",Form("w2(%d)",FINAL));
   for(int l=0;l<8;l++){
      Fit a=gfit(r1[l]),b=gfit(r1f[l]),c=gfit(r2[l]),d=gfit(r2f[l]);
      printf("%-5s | %9.2f %9.2f %9.2f %9.2f | %9.2f %9.2f %9.2f %9.2f\n",LAY[l],
             a.mean,b.mean,a.sigma,b.sigma,c.mean,d.mean,c.sigma,d.sigma);
   }
   Fit fy=gfit(hipy), fz=gfit(hipz), gy=gfit(hrpy), gz=gfit(hrpz);
   printf("\n=== DCA (fit) : reference matters ===\n");
   printf("  %-26s %12s %12s\n","","mean [um]","width [um]");
   printf("  %-26s %12.2f %12.2f\n","DCA_y vs v_est  (module)",1e4*fy.mean,1e4*fy.sigma);
   printf("  %-26s %12.2f %12.2f\n","DCA_y vs v_reco (tv1,tv2)",1e4*gy.mean,1e4*gy.sigma);
   printf("  %-26s %12.2f %12.2f\n","DCA_z vs v_est  (module)",1e4*fz.mean,1e4*fz.sigma);
   printf("  %-26s %12.2f %12.2f\n","DCA_z vs v_reco (tv3)",1e4*gz.mean,1e4*gz.sigma);
   Fit fv=gfit(hdvz);
   printf("\n=== dvz = v_est,z - v_reco,z ===\n  fit mean %+.1f um   width %.1f um   (hist mean %+.1f um)\n",
          1e4*fv.mean,1e4*fv.sigma,1e4*hdvz->GetMean());

   // residual canvases: 8 pads, coarse bins
   auto rescan=[&](TH1D** A, TH1D** B, const char* nm, const char* ttl){
      TCanvas* c=new TCanvas(nm,"",1300,620); c->Divide(4,2);
      for(int l=0;l<8;l++){ c->cd(l+1); gPad->SetLeftMargin(0.17);
         A[l]->SetLineColor(kAzure+2); A[l]->SetLineWidth(2); A[l]->SetStats(0);
         A[l]->GetXaxis()->SetTitle(Form("%s  %s [#mum]",LAY[l],ttl));
         A[l]->GetYaxis()->SetTitle("tracks"); A[l]->Draw("hist");
         B[l]->SetLineColor(kRed+1); B[l]->SetLineStyle(2); B[l]->SetLineWidth(2); B[l]->Draw("hist same");
         Fit ff=gfit(A[l]);
         if(ff.ok){ TF1* g=new TF1(Form("g_%s_%d",nm,l),"gaus",ff.mean-2*ff.sigma,ff.mean+2*ff.sigma);
            g->SetParameters(A[l]->GetMaximum(),ff.mean,ff.sigma); g->SetLineColor(kBlack); g->SetLineWidth(1); g->Draw("same"); }
      }
      c->SaveAs(Form("%s/%s.png",out,nm));
   };
   rescan(r1,r1f,"residual_s1","#Deltas1");
   rescan(r2,r2f,"residual_s2","#Deltas2");

   draw(ept,out,"_vest","v_{est}");  draw(eet,out,"_vest","v_{est}");
   draw(ez ,out,"_vest","v_{est}");  draw(eph,out,"_vest","v_{est}");
   draw(rpt,out,"_vreco","v_{reco}");draw(ret,out,"_vreco","v_{reco}");
   draw(rz ,out,"_vreco","v_{reco}");draw(rph,out,"_vreco","v_{reco}");

   // the headline comparison: DCA_y width vs pT for both references
   {
      std::vector<double> x,we,wr,ee,er;
      for(int i=0;i<ept.n;i++){
         Fit a=gfit(ept.hy[i]), b=gfit(rpt.hy[i]);
         if(!a.ok||!b.ok) continue;
         x.push_back(ept.lo+(i+0.5)*(ept.hi-ept.lo)/ept.n);
         we.push_back(1e4*a.sigma); ee.push_back(1e4*a.esigma);
         wr.push_back(1e4*b.sigma); er.push_back(1e4*b.esigma);
      }
      if(x.size()>1){
         TCanvas* cd=new TCanvas("cdref","",900,620); cd->SetGridy(); gPad->SetLeftMargin(0.16);
         TGraphErrors* ga=new TGraphErrors(x.size(),&x[0],&we[0],nullptr,&ee[0]);
         TGraphErrors* gb=new TGraphErrors(x.size(),&x[0],&wr[0],nullptr,&er[0]);
         ga->SetLineColor(kGray+2); ga->SetMarkerColor(kGray+2); ga->SetMarkerStyle(24); ga->SetLineWidth(2); ga->SetMarkerSize(1.3);
         gb->SetLineColor(kAzure+2); gb->SetMarkerColor(kAzure+2); gb->SetMarkerStyle(20); gb->SetLineWidth(3); gb->SetMarkerSize(1.4);
         TMultiGraph* mg=new TMultiGraph(); mg->Add(ga,"lp"); mg->Add(gb,"lp"); mg->Draw("A");
         mg->GetXaxis()->SetTitle("p_{T} [GeV/c]");
         mg->GetYaxis()->SetTitle("DCA_{y} width [#mum]"); mg->GetYaxis()->SetTitleOffset(1.5);
         mg->SetMinimum(0.);
         TLegend* lg=new TLegend(0.42,0.72,0.9,0.88);
         lg->AddEntry(ga,"vs v_{est} (module, circular)","lp");
         lg->AddEntry(gb,"vs v_{reco} = tv1,tv2 (independent)","lp");
         lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036); lg->Draw();
         cd->SaveAs(Form("%s/dca_reference.png",out));
      }
   }

   // dvz diagnostics
   TCanvas* cv=new TCanvas("cv","",1300,620); cv->Divide(3,2);
   auto pp=[&](int p,TProfile* h,const char* xt){ cv->cd(p); gPad->SetGridy(); gPad->SetLeftMargin(0.17);
      h->SetLineColor(kMagenta+2); h->SetMarkerColor(kMagenta+2); h->SetMarkerStyle(20); h->SetStats(0);
      h->GetXaxis()->SetTitle(xt); h->GetYaxis()->SetTitle("dv_{z} [#mum]");
      h->GetYaxis()->SetTitleOffset(1.6); h->Draw("e1"); };
   pp(1,pdvz_pt,"p_{T} [GeV/c]"); pp(2,pdvz_eta,"#eta");
   pp(3,pdvz_z,"v_{z}^{reco} [cm]"); pp(4,pdvz_ph,"#phi [rad]");
   cv->cd(5); gPad->SetLeftMargin(0.17); hdvz->SetLineColor(kMagenta+2); hdvz->SetLineWidth(2);
   hdvz->GetXaxis()->SetTitle("dv_{z} [cm]"); hdvz->GetYaxis()->SetTitle("tracks"); hdvz->Draw("hist");
   cv->SaveAs(Form("%s/dvz.png",out));
   {  // cost vs epoch
      std::vector<double> ep,lr,te; std::ifstream in("tools/monitoring/cost_epoch.txt");
      double a,b,c; while(in>>a>>b>>c){ep.push_back(a);lr.push_back(b);te.push_back(c);}
      if(ep.size()>1){
         TCanvas* cc=new TCanvas("cc","",900,620); cc->SetGridy(); gPad->SetLeftMargin(0.16);
         TGraph* g1=new TGraph(ep.size(),&ep[0],&lr[0]); TGraph* g2=new TGraph(ep.size(),&ep[0],&te[0]);
         double lo=1e9,hi=-1e9; for(double v:lr){lo=std::min(lo,v);hi=std::max(hi,v);}
         for(double v:te){lo=std::min(lo,v);hi=std::max(hi,v);} double pd=0.15*(hi-lo);
         g1->SetTitle(""); g1->GetXaxis()->SetTitle("epoch"); g1->GetYaxis()->SetTitle("cost / ndf");
         g1->GetYaxis()->SetRangeUser(lo-pd,hi+pd); g1->GetYaxis()->SetTitleOffset(1.7);
         g1->GetXaxis()->SetTitleSize(0.05); g1->GetYaxis()->SetTitleSize(0.05);
         g1->GetXaxis()->SetLabelSize(0.045); g1->GetYaxis()->SetLabelSize(0.045); g1->GetYaxis()->SetMaxDigits(5);
         g1->SetLineColor(kAzure+2); g1->SetLineWidth(3); g1->SetMarkerColor(kAzure+2); g1->SetMarkerStyle(20); g1->SetMarkerSize(1.4);
         g2->SetLineColor(kOrange+7); g2->SetLineWidth(3); g2->SetMarkerColor(kOrange+7); g2->SetMarkerStyle(21); g2->SetMarkerSize(1.4);
         g1->Draw("ALP"); g2->Draw("LP same");
         TLegend* lg=new TLegend(0.58,0.74,0.9,0.88); lg->AddEntry(g1,"training","lp"); lg->AddEntry(g2,"test","lp");
         lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042); lg->Draw();
         cc->SaveAs(Form("%s/cost_epoch.png",out));
         printf("\n=== Cost vs Epoch ===\n%8s %14s %14s\n","epoch","training","test");
         for(size_t i=0;i<ep.size();i++) printf("%8.0f %14.6f %14.6f\n",ep[i],lr[i],te[i]);
         printf("  training %+.4f%%   test %+.4f%%\n",100*(lr.back()/lr.front()-1),100*(te.back()/te.front()-1));
      }
   }
   printf("\nplots -> %s/{cost_epoch,residual_s1,residual_s2,dca_reference,dvz}.png\n",out);
   printf("         %s/dca_vs_{pT,eta,zvtx,phi}_{vest,vreco}.png\n",out);
}
