#include "TH1F.h"
#include "TRandom2.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"

double low_tgpol = 0.78;
double high_tgpol = 0.85;

int target_pol = 1; // define target polarization

double mean_beampol = 0.6;
double sigma_beampol = 0.1;
double beampol_err = 0.10; // 10% error for polarization

double G_v = 0.5;
double Sigma_v = 0.3;

int tot_events = 1000;

TRandom *fRandom;

double GetTG_pol(int event,int pol) {

  double tg_pol = double(pol) * fRandom->Uniform(low_tgpol,high_tgpol);
  
  return tg_pol;
}


double GetBeam_pol(int event,int pol) {
  double beam_pol;
  double at_mean_pol = mean_beampol + double(pol)*0.05;
  if ( fRandom->Rndm() < 0.5 ) {
    sigma_beampol = 0.1;
    beam_pol = fRandom->Gaus(0.0, sigma_beampol);
    beam_pol = double(pol) * (at_mean_pol + TMath::Abs(beam_pol));
  }
  else {
    sigma_beampol = 0.025;
    beam_pol = fRandom->Gaus(0.0, sigma_beampol);
    beam_pol =double(pol) *( at_mean_pol - TMath::Abs(beam_pol) );
 
  }
  return beam_pol;
}

Double_t dsigma_domega_pol(Double_t *x, Double_t *par)
{
  Float_t xx =x[0];
  Double_t f = 0.0;
  if(((-180<xx/TMath::Pi()*180.) && (xx/TMath::Pi()*180.<-155))||((-145<xx/TMath::Pi()*180.)&&(xx/TMath::Pi()*180.<-95))||((-85<xx/TMath::Pi()*180.)&&(xx/TMath::Pi()*180.<-35))||((-25<xx/TMath::Pi()*180.)&&(xx/TMath::Pi()*180.<25))||((35<xx/TMath::Pi()*180.)&&(xx/TMath::Pi()*180.<85))||((95<xx/TMath::Pi()*180.)&&(xx/TMath::Pi()*180.<145))||((155<xx/TMath::Pi()*180.)&&(xx/TMath::Pi()*180.<180))){
    f = par[0]*(1+par[1]*par[2] * cos(2*xx) + par[3] * par[2] *par[4] * sin(2*xx));  
    // par[0] = cross section unpolarized
    // par[1] = Sigma
    // par[2] = Beam Polarization
    // par[3] = G
    // par[4] = Target Polarization
  }
   return f;
}


Double_t Get_flat_amo() {

  Double_t xx = 0.0;
  Double_t yy ;
  while (xx == 0.0) {
    yy = fRandom->Uniform(-TMath::Pi(),TMath::Pi());
    if(((-180<yy/TMath::Pi()*180.) && (yy/TMath::Pi()*180.<-155))||((-145<yy/TMath::Pi()*180.)&&(yy/TMath::Pi()*180.<-95))||((-85<yy/TMath::Pi()*180.)&&(yy/TMath::Pi()*180.<-35))||((-25<yy/TMath::Pi()*180.)&&(yy/TMath::Pi()*180.<25))||((35<yy/TMath::Pi()*180.)&&(yy/TMath::Pi()*180.<85))||((95<yy/TMath::Pi()*180.)&&(yy/TMath::Pi()*180.<145))||((155<yy/TMath::Pi()*180.)&&(yy/TMath::Pi()*180.<180))) xx = yy; 
  }
  return xx;

}

void Fill_histo(TH1F *histo , float  val, float pol) {
  int at_bin = histo->FindBin(val);
  float  at_val = histo->GetBinContent(at_bin);
  float at_err = histo->GetBinError(at_bin);
  float pol_err = TMath::Abs(pol) * beampol_err;
  histo->Fill(val,TMath::Abs(1./pol));
  float val_err = at_err;
  if (pol != 0.0 && pol_err != 0.0 ) val_err = pow(1./pow(pol,2) + pow(pol_err,2)/pow(pol,4),0.5);
  histo->SetBinError(at_bin,pow(pow(at_err,2)+pow(val_err,2),0.5));

}

void Run(){
  fRandom = new TRandom2(0);
  double beam_pol, tg_pol;
  TFile *file_out = new TFile("output_test_polasym5_1000.root","RECREATE");

  TF1 *para_func = new TF1("para_func",dsigma_domega_pol,-TMath::Pi(),TMath::Pi(),5);
  TF1 *perp_func = new TF1("para_func",dsigma_domega_pol,-TMath::Pi(),TMath::Pi(),5);

  TH1F *para_h = new TH1F("para_h","PARA histogram; #phi",100,-TMath::Pi(),TMath::Pi());
  TH1F *perp_h = new TH1F("perp_h","PERP histogram; #phi",100,-TMath::Pi(),TMath::Pi());

  TH1F *para2_h = new TH1F("para2_h","PARA histogram weighted pol; #phi",100,-TMath::Pi(),TMath::Pi());
  TH1F *perp2_h = new TH1F("perp2_h","PERP histogram weighted pol; #phi",100,-TMath::Pi(),TMath::Pi());

  TH1F *amo_h = new TH1F("amo_h","AMO histogram; #phi",100,-TMath::Pi(),TMath::Pi());

  TH1F *beampol_0h = new TH1F("beampol_0h","Beam Polarization -1 ; Pol",100,0.,1.0);
  TH1F *beampol_1h = new TH1F("beampol_1h","Beam Polarization 1; Pol",100,0.,1.0);

  double val_phi;
  double meas_beampol;

  for (int j=0; j<2; j++) { // beam polarization
    for (int i =0 ; i< tot_events; i++) {
      if (i % 100 == 0) printf("Polarization %d at event %d \n",(j*2-1),i);

      if (i<int(tot_events * 0.8) ) {
	amo_h->Fill(Get_flat_amo());
      }
      beam_pol = GetBeam_pol(i,(j*2)-1);
      meas_beampol = fRandom->Gaus(beam_pol,sigma_beampol);
      if (j==0) beampol_0h->Fill(TMath::Abs(meas_beampol));
      else beampol_1h->Fill(TMath::Abs(meas_beampol));
   
      //      printf("Got Beam Pol \n");
      tg_pol = GetTG_pol(i,target_pol);
      //     printf("Got TG pol \n");
      if (j==0) {
	para_func->SetParameters(10.,Sigma_v,beam_pol,G_v,tg_pol);
	val_phi = para_func->GetRandom();
	para_h->Fill(val_phi);
	// para2_h->Fill(val_phi,TMath::Abs(1/beam_pol));
	Fill_histo(para2_h,val_phi,meas_beampol);
      }

      if (j==1) {
	perp_func->SetParameters(10.,Sigma_v,beam_pol,G_v,tg_pol);
	val_phi = perp_func->GetRandom();
	perp_h->Fill(val_phi);
	//	perp2_h->Fill(val_phi,TMath::Abs(1/beam_pol));
	Fill_histo(perp2_h,val_phi,meas_beampol);
      }

    }
    
  }
  para_h->Write();
  perp_h->Write();
  para2_h->Write();
  perp2_h->Write();
  amo_h->Write();
  beampol_0h->Write();
  beampol_1h->Write();

  file_out->Close();


}


