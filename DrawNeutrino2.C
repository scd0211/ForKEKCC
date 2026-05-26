{
	TFile* inFile = new TFile("neu_20m_all.root", "READ");
	TTree* neutrinoTree = (TTree*)inFile->Get("neutrinoTree");

	int numNu = 0;
	int PDGCode = 0;
	int PDGCodeM = 0;
	int volumeID = 0;

	std::vector<int>* nuPDGCode = 0;
	std::vector<double>* nuEnergy = 0;
	std::vector<double>* xDirNu = 0;
	std::vector<double>* yDirNu = 0;
	std::vector<double>* zDirNu = 0;
	std::vector<double>* xPosNu = 0;
	std::vector<double>* yPosNu = 0;
	std::vector<double>* zPosNu = 0;
	std::vector<double>* weightPL = 0;
	std::string* volumeName = nullptr;

	neutrinoTree->SetBranchAddress("numNu", &numNu);
	neutrinoTree->SetBranchAddress("PDGCodeM", &PDGCodeM);
	neutrinoTree->SetBranchAddress("PDGCode", &nuPDGCode);
	neutrinoTree->SetBranchAddress("energy", &nuEnergy);
	neutrinoTree->SetBranchAddress("volumeName", &volumeName);

	int nEntries = neutrinoTree->GetEntries();


	double bwidth=2.5;
	double emax=1000.;
	double nbins=emax/bwidth;

	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);
	gStyle->SetTitleFont(132,"x");
	gStyle->SetTitleFont(132,"y");
	gStyle->SetTitleFont(132,"z");
	gStyle->SetLabelFont(132,"x");
	gStyle->SetLabelFont(132,"y");
	gStyle->SetLabelFont(132,"z");

	TH1D* nuESpectrum_4pi = new TH1D("nuESpectrum_4pi", "#nu_{e};E_{#nu} (MeV);",nbins, 0.0, emax);
	TH1D* nuEBarSpectrum_4pi = new TH1D("nuEBarSpectrum_4pi", "#bar{#nu}_{e};E_{#nu} (MeV);", nbins, 0.0, emax);
	TH1D* nuMuSpectrum_4pi = new TH1D("nuMuSpectrum_4pi", "#nu_{#mu};E_{#nu} (MeV);", nbins, 0.0, emax);
	TH1D* nuMuBarSpectrum_4pi = new TH1D("nuMuBarSpectrum_4pi", "#bar{#nu}_{#mu};E_{#nu} (MeV);", nbins, 0.0, emax);



	cout<<nEntries<<endl;

	int n_NuMu=0;
	int n_NuMu_pi=0;
	int n_NuMu_k=0;

	for (int ientry = 0; ientry < nEntries; ientry++) {
		neutrinoTree->GetEntry(ientry);

		for (int inu = 0; inu < nuPDGCode->size(); inu++) {
			int PDGCode=nuPDGCode->at(inu);
			double energy = nuEnergy->at(inu);

			if (PDGCode == -14) {
				nuMuBarSpectrum_4pi->Fill(energy);
			}
			else if (PDGCode == 14) {
				nuMuSpectrum_4pi->Fill(energy);
				if(energy>250){
					n_NuMu++;
					if(PDGCodeM==211) n_NuMu_pi++;
					else if(PDGCodeM==321) n_NuMu_k++;
					else cout<<PDGCodeM<<" "<<volumeName<<endl;
				}

			}
			else if (PDGCode == 12) {
				nuESpectrum_4pi->Fill(energy);
			}
			else if (PDGCode == -12) {
				nuEBarSpectrum_4pi->Fill(energy);
			}

		}
	}

	TCanvas * c1 = new TCanvas();
	c1->SetLogy();

	nuMuSpectrum_4pi->Draw("hist");
	nuMuBarSpectrum_4pi->Draw("histsame");
	nuESpectrum_4pi->Draw("histsame");
	nuEBarSpectrum_4pi->Draw("histsame");

	nuMuSpectrum_4pi->SetLineWidth(2);
	nuMuBarSpectrum_4pi->SetLineWidth(2);
	nuESpectrum_4pi->SetLineWidth(2);
	nuEBarSpectrum_4pi->SetLineWidth(2);

	nuMuSpectrum_4pi->SetLineColor(2);
	nuMuBarSpectrum_4pi->SetLineColor(1);
	nuESpectrum_4pi->SetLineColor(6);
	nuEBarSpectrum_4pi->SetLineColor(4);

	nuMuSpectrum_4pi->SetTitle("");
	nuMuSpectrum_4pi->GetYaxis()->SetTitle("#nu/2x10^{7} POT/2.5 MeV");

	cout<<nuMuSpectrum_4pi->GetBinWidth(1)<<endl;

	TLegend * leg = new TLegend(0.6,0.6,0.89,0.89);
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextFont(132);
	leg->AddEntry(nuESpectrum_4pi,"#nu_{e}");
	leg->AddEntry(nuMuSpectrum_4pi,"#nu_{#mu}");
	leg->AddEntry(nuEBarSpectrum_4pi,"#bar{#nu}_{e}");
	leg->AddEntry(nuMuBarSpectrum_4pi,"#bar{#nu}_{#mu}");
	leg->Draw();

	cout<<"#nu_{e} : "<<nuESpectrum_4pi->Integral()<<endl;
	cout<<"#nu_{#mu} : "<<nuMuSpectrum_4pi->Integral()<<endl;
	cout<<"#bar{#nu}_{e} : "<<nuEBarSpectrum_4pi->Integral()<<endl;
	cout<<"#bar{#nu}_{#mu} : "<<nuMuBarSpectrum_4pi->Integral()<<endl;

	cout<<"numu with E>250MeV : "<<n_NuMu<<endl;
	cout<<"from pion : "<<n_NuMu_pi<<" ("<<100.*double(n_NuMu_pi)/double(n_NuMu)<<"%)"<<endl;
	cout<<"from kaon : "<<n_NuMu_k<<" ("<<100.*double(n_NuMu_k)/double(n_NuMu)<<"%)"<<endl;


}
