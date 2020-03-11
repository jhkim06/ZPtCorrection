void GetZptReweight(bool isElectron = true, int massBin = 2, TString histFile = "hists/ISR_detector_plots_electron.root")
{
    gROOT->SetBatch();
    TH1::AddDirectory(kFALSE);

    // Input histograms for Data and Backgrounds at reconstruction level
    if(!isElectron)
    {
        histFile = "hists/ISR_detector_plots_muon.root";
    }
    TFile fhist(histFile);
    // Output file
    TString output_postfix = "electron";
    TString massBinStr = "m80to100";
    if(massBin ==0)
    {
        massBinStr = "m50to65";
        if(!isElectron)
            massBinStr = "m40to60";
    }
    if(massBin ==1)
    {
        massBinStr = "m65to80";
        if(!isElectron)
            massBinStr = "m60to80";
    }
    if(massBin ==3)
    {
        massBinStr = "m100to200";
    }
    if(massBin ==4)
    {
        massBinStr = "m200to350";
    }

    if(!isElectron)
    {
        output_postfix = "muon";
    }

    TFile outfile("ZptWeight_"+output_postfix+"_"+massBinStr+".root","UPDATE");


    TH1* hDYPtReweighted = NULL;

    int stop = 0;
    int iter = 0;

    TH1* hDataPt = NULL;
    TH1* hBkgPt = NULL;
    TH1* hZptWeight = NULL;
    TH1* hZptWeight_temp = NULL;
    TH1* hDYPtBfReweighted = NULL;
    TH1* hDYPtReweighted_previous = NULL;

    int previous_iter = 0;
    TString previous_iter_;
    do
    {
        previous_iter_.Form("%d", previous_iter + 1);
        hZptWeight_temp = (TH1*)outfile.Get("ZptWeight_iter"+previous_iter_);
        if(hZptWeight_temp != NULL) cout << previous_iter_ << " th iteration found!" << endl;
        previous_iter++;
    }
    while(hZptWeight_temp != NULL);
    iter = previous_iter - 1;
    previous_iter_.Form("%d", iter);

    if(iter != 0)
    {
        TString current_iter_;
        current_iter_.Form("%d", iter + 1);
        hZptWeight_temp = (TH1*)outfile.Get("ZptWeightInputFor_iter"+current_iter_);
    }

    TString dir = "detector_level_"+massBinStr;
    TString dataName = "DoubleEG";
    TString channelName = "EE";
    if(isElectron == false)
    {
        dataName = "DoubleMuon";
        channelName = "MuMu";
    }

    // Get Data and Background histograms
    hDataPt=(TH1*)fhist.Get(dir+"/ptll_variable_bin/histo_"+dataName);
    hBkgPt=(TH1*)fhist.Get(dir+"/ptll_variable_bin/histo_ZZ_pythia");
    hBkgPt->Add((TH1*)fhist.Get(dir+"/ptll_variable_bin/histo_WZ_pythia"));
    hBkgPt->Add((TH1*)fhist.Get(dir+"/ptll_variable_bin/histo_WW_pythia"));
    hBkgPt->Add((TH1*)fhist.Get(dir+"/ptll_variable_bin/histo_TTLL_powheg"));
    hBkgPt->Add((TH1*)fhist.Get(dir+"/ptll_variable_bin/histo_DYJets10to50ToTauTau"));
    hBkgPt->Add((TH1*)fhist.Get(dir+"/ptll_variable_bin/histo_DYJetsToTauTau"));

    hDYPtBfReweighted=(TH1*)fhist.Get(dir+"/ptll_variable_bin/histo_DYJetsTo"+channelName);
    hDYPtBfReweighted->Add((TH1*)fhist.Get(dir+"/ptll_variable_bin/histo_DYJets10to50To"+channelName));
    hDYPtReweighted=(TH1*)fhist.Get(dir+"/ptll_variable_bin/histo_DYJetsTo"+channelName);
    hDYPtReweighted->Add((TH1*)fhist.Get(dir+"/ptll_variable_bin/histo_DYJets10to50To"+channelName));

    if(hZptWeight_temp != NULL)
    {
        cout << "Get " + previous_iter_ + " th reweighted DY histogram..." << endl;
        hDYPtReweighted_previous=(TH1*)outfile.Get("hDYReweighted_iter" + previous_iter_);
    }
    else
    {
        hDYPtReweighted_previous=(TH1*)fhist.Get(dir+"/ptll_variable_bin/histo_DYJetsTo"+channelName);
        hDYPtReweighted_previous->Add((TH1*)fhist.Get(dir+"/ptll_variable_bin/histo_DYJets10to50To"+channelName));
    }

    if(hZptWeight_temp == NULL)
    {
        cout << "No previous ZptWeight exits..." << endl;
        hZptWeight_temp=(TH1*)fhist.Get(dir+"/ptll_variable_bin/histo_"+dataName);
        hZptWeight_temp->Add(hBkgPt, -1);                   // Subtract bkg from data
        hZptWeight_temp->Divide(hDYPtBfReweighted);
    }

    TCanvas* c1=new TCanvas("c1", "c1", 50, 50, 750, 550);
    c1->cd();
    hZptWeight_temp->Draw();
    hZptWeight_temp->SetMinimum(0.8);
    hZptWeight_temp->SetMaximum(1.2);
    c1->SaveAs("comparison.png");
    delete c1;

    // DY MC tree
    // Reconstruction level variables
    bool evt_tag_analysisevnt_sel_rec;
    bool evt_tag_dielectron_rec;
    bool evt_tag_dimuon_rec;

    Double_t evt_weight_total_rec;
    Double_t evt_weight_recoSF_rec;
    Double_t evt_weight_idSF_rec;
    Double_t evt_weight_isoSF_rec;
    Double_t evt_weight_trigSF_rec;
    Double_t evt_weight_trigSFDZ_rec;

    Double_t dilep_pt_rec;
    Double_t dilep_mass_rec;

    // Generator level variables
    bool evt_tag_ditau_hardprocess;
    bool evt_tag_dielectron_hardprocess;
    bool evt_tag_dimuon_hardprocess;
    bool evt_tag_dimuon_promptfinal;
    bool evt_tag_dielectron_promptfinal;

    Double_t dilep_pt_FSRgamma_gen_ispromptfinal;
    Double_t dilep_mass_FSRgamma_gen_ispromptfinal;
    Double_t evt_weight_total_gen;

    Long64_t nentries;

    TTree* tsignal;

    TChain * chain = new TChain("recoTree/SKFlat","ZptCorrection");
    string line;
    ifstream myfile("/home/jhkim/ISR_Run2/ZptCorrection/2016ISR/ntuples/list.txt");
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
            if( line.length() != 0 )
            {
                cout << line << '\n';
                char file_path[300];
                strcpy(file_path, line.c_str());
                chain->Add(file_path);
            }
        }
        myfile.close();
    }
    tsignal = chain;

    tsignal->SetBranchAddress("evt_tag_analysisevnt_sel_rec_Nominal",&evt_tag_analysisevnt_sel_rec);
    tsignal->SetBranchAddress("evt_tag_dielectron_rec_Nominal",&evt_tag_dielectron_rec);
    tsignal->SetBranchAddress("evt_tag_dimuon_rec_Nominal",&evt_tag_dimuon_rec);
    tsignal->SetBranchAddress("evt_weight_total_rec",&evt_weight_total_rec);
    tsignal->SetBranchAddress("evt_weight_recoSF_rec_Nominal",&evt_weight_recoSF_rec);
    tsignal->SetBranchAddress("evt_weight_idSF_rec_Nominal",&evt_weight_idSF_rec);
    tsignal->SetBranchAddress("evt_weight_isoSF_rec_Nominal",&evt_weight_isoSF_rec);
    tsignal->SetBranchAddress("evt_weight_trigSF_rec_Nominal",&evt_weight_trigSF_rec);
    tsignal->SetBranchAddress("evt_weight_trigSFDZ_rec_Nominal",&evt_weight_trigSFDZ_rec);
    tsignal->SetBranchAddress("dilep_pt_rec_Nominal",&dilep_pt_rec);
    tsignal->SetBranchAddress("dilep_mass_rec_Nominal",&dilep_mass_rec);

    tsignal->SetBranchAddress("evt_tag_ditau_hardprocess",&evt_tag_ditau_hardprocess);
    tsignal->SetBranchAddress("evt_tag_dielectron_hardprocess",&evt_tag_dielectron_hardprocess);
    tsignal->SetBranchAddress("evt_tag_dimuon_hardprocess",&evt_tag_dimuon_hardprocess);
    tsignal->SetBranchAddress("evt_tag_dimuon_promptfinal",&evt_tag_dimuon_promptfinal);
    tsignal->SetBranchAddress("evt_tag_dielectron_promptfinal",&evt_tag_dielectron_promptfinal);
    tsignal->SetBranchAddress("dilep_pt_FSRgamma_gen_ispromptfinal",&dilep_pt_FSRgamma_gen_ispromptfinal);
    tsignal->SetBranchAddress("dilep_mass_FSRgamma_gen_ispromptfinal",&dilep_mass_FSRgamma_gen_ispromptfinal);
    tsignal->SetBranchAddress("evt_weight_total_gen",&evt_weight_total_gen);

    nentries=tsignal->GetEntries();

    do
    {
        cout << "Start " << iter+1 << " th iteration..." << endl;

        TH1* hDYPtReweighted_temp;
        hDYPtReweighted_temp = (TH1*)hDYPtReweighted->Clone("hDYReweighted_temp");
        hDYPtReweighted_temp->Reset(); // Reset before apply the current Z pt reweight in this iteration

        cout<<"start signal loop"<<endl;
        for(Long64_t i=0;i<nentries;i++)
        //for(int i=0;i<100000;i++)
        {
            tsignal->GetEntry(i);
            if(i%10000000==0) cout<<i<< " /" << nentries << " (" << (double)(100.*i/nentries) << "%)" << endl;

            if(evt_tag_ditau_hardprocess == false)
            {
                Double_t GenZpt;          // FIXME electron or muon
                if(isElectron)
                    if(evt_tag_dimuon_promptfinal) continue;
                if(!isElectron)
                    if(evt_tag_dielectron_promptfinal) continue;

                GenZpt = dilep_pt_FSRgamma_gen_ispromptfinal;

                bool recChannel = evt_tag_dielectron_rec;
                Double_t totWeight = evt_weight_total_gen*evt_weight_total_rec*evt_weight_recoSF_rec*evt_weight_idSF_rec*evt_weight_trigSF_rec*evt_weight_trigSFDZ_rec;
                if(!isElectron)
                {
                    recChannel = evt_tag_dimuon_rec;
                    totWeight = evt_weight_total_gen*evt_weight_total_rec*evt_weight_isoSF_rec*evt_weight_idSF_rec*evt_weight_trigSF_rec;
                }

                if(evt_tag_analysisevnt_sel_rec && recChannel)      // TODO option for muon
                {
                    double low_mass = 80.;
                    double high_mass = 100.;
                    if(massBin ==0)
                    {
                        low_mass = 50.;
                        high_mass = 65.;
                        if(!isElectron)
                        {
                            low_mass = 40.;
                            high_mass = 60.;
                        }
                    }
                    if(massBin ==1)
                    {
                        low_mass = 65.;
                        high_mass = 80.;
                        if(!isElectron)
                        {
                            low_mass = 60.;
                            high_mass = 80.;
                        }
                    }
                    if(massBin ==3)
                    {
                        low_mass = 100.;
                        high_mass = 200.;
                    }
                    if(massBin ==4)
                    {
                        low_mass = 200.;
                        high_mass = 350.;
                    }
                    if(dilep_mass_rec > low_mass && dilep_mass_rec < high_mass && dilep_pt_rec < 100.)
                    {
                        if( GenZpt > 100.) GenZpt = 99.5;
                        int bin = hZptWeight_temp->FindBin(GenZpt);
                        Double_t zpt_weight = hZptWeight_temp->GetBinContent(bin);
                        hDYPtReweighted_temp->Fill(dilep_pt_rec, totWeight*zpt_weight);
                    }// Z peak
                }
            }
        }// DY loop
        iter++;
        TString iter_; iter_.Form("%d",iter);

        double norm_;
        norm_ = hDYPtBfReweighted->Integral()/ hDYPtReweighted_temp->Integral();
        cout << "norm: " << norm_ << endl;
        hDYPtReweighted_temp->Scale(norm_);
        hDYPtReweighted_temp->SetName("hDYReweighted_iter" + iter_);
        hZptWeight_temp->Scale(norm_);

        // Update reweight histogram
        TH1* hDataPt_temp;
        TH1* hBkgPt_temp;

        hDataPt_temp=(TH1*)hDataPt->Clone("hDataPt_temp");
        hBkgPt_temp=(TH1*)hBkgPt->Clone("hBkgPt_temp");

        hDataPt_temp->Add(hBkgPt_temp, -1);
        hDataPt_temp->Divide(hDYPtReweighted_temp); // Next weight before NORMALISATION applied

        // Save weight and normalisation information, this can be used for reweighting
        hZptWeight = (TH1*)hZptWeight_temp->Clone("ZptWeight_iter" + iter_);

        for(int bin = 1; bin < hZptWeight_temp->GetNbinsX()+1; bin++)
        {
            float current_weight = hZptWeight_temp->GetBinContent(bin);
            hZptWeight_temp->SetBinContent(bin, current_weight*(hDataPt_temp->GetBinContent(bin))); // for next iteration
            cout << bin << " th bin.." << " next weight: " << hZptWeight_temp->GetBinContent(bin) << " current weight: " << current_weight << endl;
        }
        TString next_iter_; next_iter_.Form("%d",iter + 1);
        TH1* hZptWeightInputForNextIter = (TH1*)hZptWeight_temp->Clone("ZptWeightInputFor_iter" + next_iter_);

        // TODO save this comparison
        TH1* hDYPt_comparison = (TH1*)hDYPtReweighted_previous->Clone("hDYcomparison_iter"+ iter_);
        hDYPt_comparison->Divide(hDYPtReweighted_temp);
        hDYPtReweighted_previous->Reset(); // Reset for next iteration
        hDYPtReweighted_previous->Add(hDYPtReweighted_temp); // Save current ZPt reweighted histogram

        cout << iter << " th iteration finished " << endl;

        outfile.cd();
        hZptWeight->Write();
        hZptWeightInputForNextIter->Write();
        hDYPtReweighted_temp->Write();
        hDYPt_comparison->Write();

        delete hDataPt_temp;
        delete hBkgPt_temp;
        delete hZptWeight;
        delete hDYPtReweighted_temp;
        delete hDYPt_comparison;
        delete hZptWeightInputForNextIter;
        //delete tsignal;
    }
    while(iter!=10); // TODO add condition to exit
    delete chain;

    std::cout << "initial N: " << hDYPtBfReweighted->Integral() << " final N: " << hDYPtReweighted->Integral() << std::endl;

    outfile.cd();
    outfile.Close();
}

