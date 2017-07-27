{
  double low_tgpol = 0.78;
  double high_tgpol = 0.85;
  double avg_tgpol = (low_tgpol + high_tgpol) /2;

  int target_pol = 1; // define target polarization
 
  double mean_beampol = 0.6;
  double sigma_beampol = 0.1;


  TFile *file_in = TFile::Open("output_test_pol_1000.root");
  
  TH1F *para_h = file_in->Get("para_h");
  TH1F *perp_h = file_in->Get("perp_h");
  
  TH1F *para2_h = file_in->Get("para2_h");
  TH1F *perp2_h = file_in->Get("perp2_h");

  TH1F *amo_h = file_in->Get("amo_h");
  
  para_h->Rebin(2);
  perp_h->Rebin(2);
  para2_h->Rebin(2);
  perp2_h->Rebin(2);
  amo_h->Rebin(2);

  TH1F *1_para_h = (TH1F*)para_h->Clone();
  TH1F *1_perp_h = (TH1F*)perp_h->Clone();
  TH1F *1_para2_h = (TH1F*)para2_h->Clone();
  TH1F *1_perp2_h = (TH1F*)perp2_h->Clone();

  TH1F *PARA_asym2 = (TH1F*)para2_h->Clone();
  PARA_asym2->Sumw2();
  PARA_asym2->Divide(amo_h);
  TH1F *PERP_asym2 = (TH1F*)perp2_h->Clone();
  PERP_asym2->Sumw2();
  PERP_asym2->Divide(amo_h);

  TH1F *para3_h = (TH1F*)para2_h->Clone();
  para3_h->Sumw2();

  TH1F *PARA_asym = (TH1*)para_h->Clone();
  PARA_asym->Sumw2();
  PARA_asym->Divide(amo_h);
  TH1F *PERP_asym = (TH1*)perp_h->Clone();
  PERP_asym->Sumw2();
  PERP_asym->Divide(amo_h);

  TH1F *para4_h = para_h->Clone();
  para4_h->Sumw2();

  TCanvas *c1 = new TCanvas();
  TH1F *asym = 1_perp_h->GetAsymmetry(1_para_h);
  TH1F *asym2 = 1_perp2_h->GetAsymmetry(1_para2_h);
  TF1 *cos2phi_sin2phi = new TF1("cos2phi_sin2phi","([0]+([1]*cos((2*(x-[2]))))+([3]*sin((2*(x-[2])))))",-TMath::Pi(),TMath::Pi());
  cos2phi_sin2phi->FixParameter(2,0);
  asym->Fit("cos2phi_sin2phi","QB");
  asym->Draw("E");
  double G_pol=((cos2phi_sin2phi->GetParameter(3))/(avg_tgpol * mean_beampol));
  double G_pol_e = ((cos2phi_sin2phi->GetParError(3))/(avg_tgpol * mean_beampol));
  printf ("G= %.3f +/- %.3f \n",G_pol,G_pol_e);

  TCanvas *c2 = new TCanvas();

  asym2->Fit("cos2phi_sin2phi","QB");
  asym2->Draw("E");
  double G_pol2=((cos2phi_sin2phi->GetParameter(3))/(avg_tgpol));
  double G_pol2_e = ((cos2phi_sin2phi->GetParError(3))/(avg_tgpol));
  printf ("weighted polarization G= %.3f +/- %.3f \n",G_pol2,G_pol2_e);
  


  PARA_asym2->Fit("cos2phi_sin2phi","QB");
  double PARA_asym2_ratio = cos2phi_sin2phi->GetParameter(0);
  PERP_asym2->Fit("cos2phi_sin2phi","QB");
  double PERP_asym2_ratio = cos2phi_sin2phi->GetParameter(0);
  double val_flux_ratio2 = PERP_asym2_ratio/PARA_asym2_ratio; // I have defined the simulation in the other way for para and perp
  
  para3_h->Scale(val_flux_ratio2);

  TH1F *asym3 = perp2_h->GetAsymmetry(para3_h);
  
  TCanvas *c3 = new TCanvas();

  asym3->Fit("cos2phi_sin2phi","QB");
  asym3->Draw("E");
  double G_pol3=((cos2phi_sin2phi->GetParameter(3))/(avg_tgpol));
  double G_pol3_e = ((cos2phi_sin2phi->GetParError(3))/(avg_tgpol));
  printf ("weighted polarization G= %.3f +/- %.3f \n",G_pol3,G_pol3_e);



  PARA_asym->Fit("cos2phi_sin2phi","QB");
  double PARA_asym_ratio = cos2phi_sin2phi->GetParameter(0);
  PERP_asym->Fit("cos2phi_sin2phi","QB");
  double PERP_asym_ratio = cos2phi_sin2phi->GetParameter(0);
  double val_flux_ratio = PERP_asym_ratio/PARA_asym_ratio; // I have defined the simulation in the other way for para and perp
  
  para4_h->Scale(val_flux_ratio);

  TH1F *asym4 = perp_h->GetAsymmetry(para4_h);
  
  TCanvas *c4 = new TCanvas();

  asym4->Fit("cos2phi_sin2phi","QB");
  asym4->Draw("E");
  double G_pol4=((cos2phi_sin2phi->GetParameter(3))/(avg_tgpol* mean_beampol));
  double G_pol4_e = ((cos2phi_sin2phi->GetParError(3))/(avg_tgpol* mean_beampol));
  printf ("weighted polarization G= %.3f +/- %.3f \n",G_pol4,G_pol4_e);


}
