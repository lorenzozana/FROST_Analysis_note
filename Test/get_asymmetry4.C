{
  double low_tgpol = 0.78;
  double high_tgpol = 0.85;
  double avg_tgpol = (low_tgpol + high_tgpol) /2;

  int target_pol = 1; // define target polarization
 
  double mean_beampol = 0.63;
  double sigma_beampol = 0.1;


  //  TFile *file_in = TFile::Open("output_test_pol_10000.root");
  TFile *file_in = TFile::Open("output_test_polasym4_1000.root");

  TH1F *para_h = file_in->Get("para2_h");
  TH1F *perp_h = file_in->Get("perp2_h");
  // TH1F *para_h = file_in->Get("para_h");
  // TH1F *perp_h = file_in->Get("perp_h");
  TH1F *amo_h = file_in->Get("amo_h");
  
  para_h->Rebin(2);
  perp_h->Rebin(2);
  amo_h->Rebin(2);
  para_h->Sumw2();
  perp_h->Sumw2();
  amo_h->Sumw2();

  TH1F *para3_h = (TH1F*)para_h->Clone();
  para3_h->Sumw2();

  TH1F *PARA_asym = (TH1*)para_h->Clone();
  PARA_asym->Sumw2();
  PARA_asym->Divide(amo_h);
  TH1F *PERP_asym = (TH1*)perp_h->Clone();
  PERP_asym->Sumw2();
  PERP_asym->Divide(amo_h);

  TF1 *cos2phi_sin2phi = new TF1("cos2phi_sin2phi","([0]+([1]*cos((2*(x-[2]))))+([3]*sin((2*(x-[2])))))",-TMath::Pi(),TMath::Pi());
  cos2phi_sin2phi->FixParameter(2,0);
  PARA_asym->Fit("cos2phi_sin2phi","QB");
  double PARA_asym_ratio = cos2phi_sin2phi->GetParameter(0);
  PERP_asym->Fit("cos2phi_sin2phi","QB");
  double PERP_asym_ratio = cos2phi_sin2phi->GetParameter(0);
  double val_flux_ratio = PERP_asym_ratio/PARA_asym_ratio;
  para3_h->Scale(val_flux_ratio);
  TH1F *asym = perp_h->GetAsymmetry(para3_h);

 
  asym->Fit("cos2phi_sin2phi","QB");
  asym->Draw("E");
  double G_pol=((cos2phi_sin2phi->GetParameter(3))/(avg_tgpol * mean_beampol));
  double G_pol_e = ((cos2phi_sin2phi->GetParError(3))/(avg_tgpol * mean_beampol));
  printf ("G= %.3f +/- %.3f \n",G_pol,G_pol_e);
}
