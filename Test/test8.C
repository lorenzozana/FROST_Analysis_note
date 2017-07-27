#include "TH1F.h"
#include "TRandom2.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>
using namespace std;

double low_tgpol = 0.78;
double high_tgpol = 0.85;

int target_pol = 1; // define target polarization

double mean_beampol = 0.6;
double sigma_beampol = 0.1;
double beampol_err = 0.10; // 10% error for polarization

double G_val[21] = {-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
double Sigma_val[21] = {-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};

int events_sim[8]={300,600,800,1000,2000,3000,4000,5000};
  
double G_v = 0.5;
double Sigma_v = 0.3;

double G_calc = 0.0;
double Sigma_calc = 0.0;
double G_calc_e = 0.0;
double Sigma_calc_e = 0.0;

double avg_tgpol;
double mean_beampol_e,avg_tgpol_e;
int tot_events = 1000;

int binE = 1;
int bintheta = 2;
double count_histo[3];

TRandom *fRandom;

TH1F *h_tgpol_para;
TH1F *h_bpol_para;
TH1F *h_tgpol_perp;
TH1F *h_bpol_perp;

std::ofstream ofs("output_model.txt", std::ofstream::out);
char what_txt[100];



void Get_h_pol(char *filename, int tg_id,int E_id, int theta_id) {
  
  char name_histo[60];
  TH1F *histo_count[3];
  TFile *file = TFile::Open(filename);
  double m_bpol_para,m_bpol_perp;
  double m_bpol_para_e,m_bpol_perp_e;
  double m_tgpol_para,m_tgpol_perp;
  double m_tgpol_para_e,m_tgpol_perp_e;
  for (int pol=0; pol<2; pol++) {
    if (pol == 0) {
      sprintf(name_histo,"tg_pol_PARA_Target%d_E%d_Th%d",tg_id,E_id,theta_id);
      h_tgpol_para = (TH1F*)file->Get(name_histo);
      sprintf(name_histo,"beam_pol_PARA_Target%d_E%d_Th%d",tg_id,E_id,theta_id);
      h_bpol_para = (TH1F*)file->Get(name_histo);
      sprintf(name_histo,"phi_pip_PARA_Target%d_E%d_Th%d",tg_id,E_id,theta_id);
      histo_count[pol] = (TH1F*)file->Get(name_histo);
      count_histo[pol] = histo_count[pol]->GetEntries();
      m_bpol_para = h_bpol_para->GetMean();
      m_bpol_para_e = h_bpol_para->GetRMS();
      m_tgpol_para = h_tgpol_para->GetMean();
      m_tgpol_para_e = h_tgpol_para->GetRMS();
    }
    else if (pol == 1) {
      sprintf(name_histo,"tg_pol_PERP_Target%d_E%d_Th%d",tg_id,E_id,theta_id);
      h_tgpol_perp = (TH1F*)file->Get(name_histo);
      sprintf(name_histo,"beam_pol_PERP_Target%d_E%d_Th%d",tg_id,E_id,theta_id);
      h_bpol_perp = (TH1F*)file->Get(name_histo);
      sprintf(name_histo,"phi_pip_PERP_Target%d_E%d_Th%d",tg_id,E_id,theta_id);
      histo_count[pol] = (TH1F*)file->Get(name_histo);
      count_histo[pol] = histo_count[pol]->GetEntries();
      m_bpol_perp = h_bpol_perp->GetMean();
      m_bpol_perp_e = h_bpol_perp->GetRMS();
      m_tgpol_perp = h_tgpol_perp->GetMean();
      m_tgpol_perp_e = h_tgpol_perp->GetRMS();
    }
  }
  mean_beampol = (m_bpol_para + m_bpol_perp)/2.;
  mean_beampol_e = pow(pow(m_bpol_para_e,2) + pow(m_bpol_perp,2),0.5);
  avg_tgpol = (m_tgpol_para + m_tgpol_perp)/2.;
  avg_tgpol_e = pow(pow(m_tgpol_para_e,2) + pow(m_tgpol_perp,2),0.5);
  
  sprintf(name_histo,"phi_pip_AMO_Target%d_E%d_Th%d",tg_id,E_id,theta_id);
  histo_count[2] = (TH1F*)file->Get(name_histo);
  count_histo[2] = histo_count[2]->GetEntries();

  
}

double GetTG_pol(int event,int pol) {

  double tg_pol = double(pol) * fRandom->Uniform(low_tgpol,high_tgpol);
  
  return tg_pol;
}

double GetTG_pol_h(int event,int pol) {
  double tg_pol =  0.0;

  if (pol == 0) tg_pol =  h_tgpol_para->GetRandom();
  else if (pol == 1) tg_pol =  h_tgpol_perp->GetRandom();
  
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

double GetBeam_pol_h(int event,int pol) {
  double beam_pol =  0.0;

  if (pol == 0) beam_pol =  h_bpol_para->GetRandom();
  else if (pol == 1) beam_pol =  h_bpol_perp->GetRandom();
  
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

void calculate_G_Sigma(TH1F *para_h,TH1F *perp_h, TH1F *amo_h){

  TH1F *para3_h = (TH1F*)para_h->Clone();
  //  para3_h->Sumw2();

  TH1F *PARA_asym = (TH1F*)para_h->Clone();
  //  PARA_asym->Sumw2();
  PARA_asym->Divide(amo_h);
  TH1F *PERP_asym = (TH1F*)perp_h->Clone();
  //  PERP_asym->Sumw2();
  PERP_asym->Divide(amo_h);

  TF1 *cos2phi_sin2phi = new TF1("cos2phi_sin2phi","([0]+([1]*cos((2*(x-[2]))))+([3]*sin((2*(x-[2])))))",-TMath::Pi(),TMath::Pi());
  cos2phi_sin2phi->FixParameter(2,0);
  PARA_asym->Fit("cos2phi_sin2phi","QB");
  double PARA_asym_ratio = cos2phi_sin2phi->GetParameter(0);
  PERP_asym->Fit("cos2phi_sin2phi","QB");
  double PERP_asym_ratio = cos2phi_sin2phi->GetParameter(0);
  double val_flux_ratio = PERP_asym_ratio/PARA_asym_ratio;
  para3_h->Scale(val_flux_ratio);
  TH1F *asym = (TH1F*)perp_h->GetAsymmetry(para3_h);

  // TCanvas* c1 = new TCanvas();
  asym->Fit("cos2phi_sin2phi","QB");
  //0  asym->Draw("E");
  G_calc=((cos2phi_sin2phi->GetParameter(3))/(avg_tgpol * mean_beampol));
  G_calc_e = ((cos2phi_sin2phi->GetParError(3))/(avg_tgpol * mean_beampol));
  Sigma_calc=((cos2phi_sin2phi->GetParameter(1))/(mean_beampol));
  Sigma_calc_e = ((cos2phi_sin2phi->GetParError(1))/(mean_beampol));
  printf ("G= %.3f +/- %.3f \n",G_calc,G_calc_e);
  printf ("Sigma= %.3f +/- %.3f \n",Sigma_calc,Sigma_calc_e);




}





void Run(){
  fRandom = new TRandom2(0);
  double beam_pol, tg_pol;
  Get_h_pol("output_pol05_1700_pos.root",1,binE,bintheta);
  TFile *file_out = new TFile("output_test_polasym_file.root","RECREATE");

  TF1 *para_func = new TF1("para_func",dsigma_domega_pol,-TMath::Pi(),TMath::Pi(),5);
  TF1 *perp_func = new TF1("para_func",dsigma_domega_pol,-TMath::Pi(),TMath::Pi(),5);

  TH1F *para_h = new TH1F("para_h","PARA histogram; #phi",50,-TMath::Pi(),TMath::Pi());
  TH1F *perp_h = new TH1F("perp_h","PERP histogram; #phi",50,-TMath::Pi(),TMath::Pi());

  TH1F *para2_h = new TH1F("para2_h","PARA histogram weighted pol; #phi",50,-TMath::Pi(),TMath::Pi());
  TH1F *perp2_h = new TH1F("perp2_h","PERP histogram weighted pol; #phi",50,-TMath::Pi(),TMath::Pi());

  TH1F *amo_h = new TH1F("amo_h","AMO histogram; #phi",50,-TMath::Pi(),TMath::Pi());

  TH1F *beampol_0h = new TH1F("beampol_0h","Beam Polarization -1 ; Pol",100,0.,1.0);
  TH1F *beampol_1h = new TH1F("beampol_1h","Beam Polarization 1; Pol",100,0.,1.0);

  double val_phi;
  double meas_beampol;

  //first guess for G is the starting value
  double G_gues = G_v;
  double Sigma_gues = Sigma_v;

  

  
  G_calc = 1000.;
  Sigma_calc = 1000.;

  printf("Starting DeltaG=%.3f ; DeltaS=%.3f \n",TMath::Abs(G_calc-G_v),TMath::Abs(Sigma_calc-Sigma_v));
  
  while (TMath::Abs(G_calc-G_v)>0.01) {
    printf("G_calc = %.3f ; G_v = %.3f ; G_gues = %.3f ; Sigma_calc = %.3f ; Sigma_v = %.3f ; Sigma_gues = %.3f \n",G_calc,G_v,G_gues,Sigma_calc,Sigma_v,Sigma_gues);
    for (int j=0; j<3; j++) { // beam polarization PARA=0 PERP=1 AMO=2
      for (int i =0 ; i< int(count_histo[j]); i++) {
	//    for (int i =0 ; i< 200; i++) {

	if (i % 100 == 0) printf("Polarization %d at event %d \n",(j*2-1),i);

	if (j<2) {
	  beam_pol = GetBeam_pol_h(i,j);
	  //	  printf("Beam pol=%.3f \n",beam_pol);
	  //	meas_beampol = fRandom->Gaus(beam_pol,sigma_beampol); // Could be in the future to consider that the measured polarization is not the real one
	  if (j==0) beampol_0h->Fill(TMath::Abs(beam_pol));
	  else beampol_1h->Fill(TMath::Abs(beam_pol));
   
	  //      printf("Got Beam Pol \n");
	  tg_pol = GetTG_pol_h(i,target_pol);
	  //     printf("Got TG pol \n");
	}
	if (j==0) {
	  para_func->SetParameters(10.,Sigma_gues,beam_pol,G_gues,tg_pol);
	  val_phi = para_func->GetRandom();
	  para_h->Fill(val_phi);
	  // para2_h->Fill(val_phi,TMath::Abs(1/beam_pol));
	  Fill_histo(para2_h,val_phi,beam_pol);
	}

	if (j==1) {
	  perp_func->SetParameters(10.,Sigma_gues,beam_pol,G_gues,tg_pol);
	  val_phi = perp_func->GetRandom();
	  perp_h->Fill(val_phi);
	  //	perp2_h->Fill(val_phi,TMath::Abs(1/beam_pol));
	  Fill_histo(perp2_h,val_phi,beam_pol);
	}
	if (j==2 ) {
	  amo_h->Fill(Get_flat_amo());
	}


      }
    }
    //Calculate G and Sigma
    calculate_G_Sigma(para2_h,perp2_h, amo_h);
    G_gues = G_gues + G_v - G_calc;
    Sigma_gues =Sigma_gues + Sigma_v - Sigma_calc;
    para_h->Reset();
    perp_h->Reset();
    amo_h->Reset();
    para2_h->Reset();
    perp2_h->Reset();
    beampol_0h->Reset();
    beampol_1h->Reset();
    
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


void Run_grid(){
  fRandom = new TRandom2(0);
  double beam_pol, tg_pol;
  Get_h_pol("output_pol05_1700_pos.root",1,binE,bintheta);
  TFile *file_out = new TFile("output_test_polasym_file.root","RECREATE");

  TF1 *para_func = new TF1("para_func",dsigma_domega_pol,-TMath::Pi(),TMath::Pi(),5);
  TF1 *perp_func = new TF1("para_func",dsigma_domega_pol,-TMath::Pi(),TMath::Pi(),5);

  TH1F *para_h = new TH1F("para_h","PARA histogram; #phi",50,-TMath::Pi(),TMath::Pi());
  TH1F *perp_h = new TH1F("perp_h","PERP histogram; #phi",50,-TMath::Pi(),TMath::Pi());

  TH1F *para2_h = new TH1F("para2_h","PARA histogram weighted pol; #phi",50,-TMath::Pi(),TMath::Pi());
  TH1F *perp2_h = new TH1F("perp2_h","PERP histogram weighted pol; #phi",50,-TMath::Pi(),TMath::Pi());

  TH1F *amo_h = new TH1F("amo_h","AMO histogram; #phi",50,-TMath::Pi(),TMath::Pi());

  TH1F *beampol_0h = new TH1F("beampol_0h","Beam Polarization -1 ; Pol",100,0.,1.0);
  TH1F *beampol_1h = new TH1F("beampol_1h","Beam Polarization 1; Pol",100,0.,1.0);

  double val_phi;
  double meas_beampol;

  //first guess for G is the starting value
  double G_gues = G_v;
  double Sigma_gues = Sigma_v;

  

  
  G_calc = 1000.;
  Sigma_calc = 1000.;

  //  printf("Starting DeltaG=%.3f ; DeltaS=%.3f \n",TMath::Abs(G_calc-G_v),TMath::Abs(Sigma_calc-Sigma_v));
  
  for (int ii=0; ii<21 ; ii++) { // G loop
    for (int jj=0; jj<21 ; jj++) { // Sigma loop
      for (int kk=0; kk<8 ; kk++) { // number of events loop

	// printf("G_calc = %.3f ; G_v = %.3f ; G_gues = %.3f ; Sigma_calc = %.3f ; Sigma_v = %.3f ; Sigma_gues = %.3f \n",G_calc,G_v,G_gues,Sigma_calc,Sigma_v,Sigma_gues);
	for (int j=0; j<3; j++) { // beam polarization PARA=0 PERP=1 AMO=2
	  for (int i =0 ; i< int(events_sim[kk]); i++) {
	    //    for (int i =0 ; i< 200; i++) {

	    if (i % 100 == 0) printf("Polarization %d at event %d \n",(j*2-1),i);

	    if (j<2) {
	      beam_pol = GetBeam_pol_h(i,j);
	      //	  printf("Beam pol=%.3f \n",beam_pol);
	      //	meas_beampol = fRandom->Gaus(beam_pol,sigma_beampol); // Could be in the future to consider that the measured polarization is not the real one
	      if (j==0) beampol_0h->Fill(TMath::Abs(beam_pol));
	      else beampol_1h->Fill(TMath::Abs(beam_pol));
   
	      //      printf("Got Beam Pol \n");
	      tg_pol = GetTG_pol_h(i,target_pol);
	      //     printf("Got TG pol \n");
	    }
	    if (j==0) {
	      para_func->SetParameters(10.,Sigma_val[jj],beam_pol,G_val[ii],tg_pol);
	      val_phi = para_func->GetRandom();
	      para_h->Fill(val_phi);
	      // para2_h->Fill(val_phi,TMath::Abs(1/beam_pol));
	      Fill_histo(para2_h,val_phi,beam_pol);
	    }

	    if (j==1) {
	      perp_func->SetParameters(10.,Sigma_val[jj],beam_pol,G_val[ii],tg_pol);
	      val_phi = perp_func->GetRandom();
	      perp_h->Fill(val_phi);
	      //	perp2_h->Fill(val_phi,TMath::Abs(1/beam_pol));
	      Fill_histo(perp2_h,val_phi,beam_pol);
	    }
	    if (j==2 ) {
	      amo_h->Fill(Get_flat_amo());
	    }


	  }
	}
	//Calculate G and Sigma
	calculate_G_Sigma(para2_h,perp2_h, amo_h);
	G_gues = G_gues + G_v - G_calc;
	Sigma_gues =Sigma_gues + Sigma_v - Sigma_calc;
	sprintf(what_txt,"%d \t %.3f  \t %.3f \t %.3f \t %.3f \n",events_sim[kk],G_val[ii],G_calc,Sigma_val[jj],Sigma_calc);
	ofs << what_txt;
	sprintf(what_txt,"para_%d_%.1f_%.1f",events_sim[kk],G_val[ii],Sigma_val[jj]);
	para2_h->SetName(what_txt);
	para2_h->Write();
	sprintf(what_txt,"perp_%d_%.1f_%.1f",events_sim[kk],G_val[ii],Sigma_val[jj]);
	perp2_h->SetName(what_txt);
	perp2_h->Write();
	sprintf(what_txt,"amo_%d_%.1f_%.1f",events_sim[kk],G_val[ii],Sigma_val[jj]);
	amo_h->SetName(what_txt);
	amo_h->Write();
	sprintf(what_txt,"beampol0_%d_%.1f_%.1f",events_sim[kk],G_val[ii],Sigma_val[jj]);
	beampol_0h->SetName(what_txt);
	beampol_0h->Write();
	sprintf(what_txt,"beampol1_%d_%.1f_%.1f",events_sim[kk],G_val[ii],Sigma_val[jj]);
	beampol_1h->SetName(what_txt);
	beampol_1h->Write();
	printf("events = %d ; G_calc = %.3f ; G_v = %.3f ; Sigma_calc = %.3f ; Sigma_v = %.3f \n",events_sim[kk],G_calc,G_val[ii],Sigma_calc,Sigma_val[jj]);

	
	para_h->Reset();
	perp_h->Reset();
	amo_h->Reset();
	para2_h->Reset();
	perp2_h->Reset();
	beampol_0h->Reset();
	beampol_1h->Reset();
    
      }
    }
  }
  ofs.close();  

  file_out->Close();


}


