// Standard Library Includes
#include <iostream>
#include <string>
#include <math.h>

// ROOT Includes
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TParameter.h"

// Geant4 Includes
#include "G4RunManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4DecayTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DynamicParticle.hh"
#include "G4Track.hh"
#include "G4DecayProducts.hh"
#include "G4SystemOfUnits.hh"
#include "G4MuMinusCapturePrecompound.hh"
#include "G4PhysicalVolumeStore.hh"

// Action Initialization
#include "MLFBeamGeometry.hh"
#include "MLFBeamActionInit.hh"
#include "QGSP_BERT.hh"

// Detector location relative to the target
const G4double DET_Y = 11600.0*mm;
const G4double DET_Z = 21000.0*mm;

// Detector size parameters (inner volume)
const G4double DET_RADIUS = 1600.0*mm;
const G4double DET_HALF_HEIGHT = 1250.0*mm;
const G4double MAX_PATH_LENGTH = 2.0*sqrt(DET_RADIUS*DET_RADIUS + DET_HALF_HEIGHT*DET_HALF_HEIGHT);

int main(int argc, char** argv) {

	if (argc != 4) { 
		std::cout << "Takes three arguments. Usage: " << std::endl;
		std::cout << "Usage: ./decayParents <input filename> <output filename> <# of decays/parent>" << std::endl;
		return 1;
	}
	std::string inFilename = argv[1];
	std::string outFilename = argv[2];

	// Number of times to decay each parent
	std::stringstream ss(argv[3]);
	int numDecays;
	if (!(ss >> numDecays)) {
		std::cout << "Invalid number!" << std::endl;
		return 1;
	}

	// Turn decay-at-rest (DAR) on and off in case we are uninterested in it
	bool useDAR = true;

	// Load the required user initialization classes; we only need the
	// physics list so GenericIon is loaded for the muon capture process,
	// but you are required to load all 3 initialization classes before
	// calling initialize or Geant4 will crash
	G4RunManager* runManager = new G4RunManager;
	runManager->SetUserInitialization(new QGSP_BERT());
	runManager->SetUserInitialization(new MLFBeamGeometry());
	runManager->SetUserInitialization(new MLFBeamActionInit());
	runManager->Initialize();

	// Get some geo objects that we'll need later
	MLFBeamGeometry* geo = (MLFBeamGeometry*)runManager->GetUserDetectorConstruction();
	G4PhysicalVolumeStore* vStore = G4PhysicalVolumeStore::GetInstance();

	// Create the muon capture process to invoke it later
	G4MuMinusCapturePrecompound muCapture;

	// Load the other particle definitions we need so they get
	// filled into the appropriate particle tables (not just parents but also daughters)
	G4ParticleDefinition* muPlus = G4MuonPlus::MuonPlusDefinition();
	G4ParticleDefinition* muMinus = G4MuonMinus::MuonMinusDefinition();
	G4ParticleDefinition* muNeutrino = G4NeutrinoMu::NeutrinoMuDefinition();
	G4ParticleDefinition* muAntiNeutrino = G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
	G4ParticleDefinition* positron = G4Positron::PositronDefinition();
	G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
	G4ParticleDefinition* elecNeutrino = G4NeutrinoE::NeutrinoEDefinition();
	G4ParticleDefinition* elecAntiNeutrino = G4AntiNeutrinoE::AntiNeutrinoEDefinition();
	G4ParticleDefinition* piPlus = G4PionPlus::PionPlusDefinition();
	G4ParticleDefinition* piMinus = G4PionMinus::PionMinusDefinition();
	G4ParticleDefinition* piZero = G4PionZero::PionZeroDefinition();
	G4ParticleDefinition* kaonPlus = G4KaonPlus::KaonPlusDefinition();
	G4ParticleDefinition* kaonMinus = G4KaonMinus::KaonMinusDefinition();
	G4ParticleDefinition* kaonLong = G4KaonZeroLong::KaonZeroLongDefinition();
	G4ParticleDefinition* kaonShort = G4KaonZeroShort::KaonZeroShortDefinition();

	// Count the number of parents of each type
	int numMuPlus = 0, numMuMinus = 0;
	int numPiPlus = 0, numPiMinus = 0;
	int numKPlus = 0, numKMinus = 0, numKShort = 0, numKLong = 0;

	// Add the helicity suppressed pion decays
	G4DecayTable* piPlusTable = piPlus->GetDecayTable();
	G4DecayTable* piMinusTable = piMinus->GetDecayTable();
	// Reduce the pi -> mu nu modes to the PDG values
	G4VDecayChannel* piPlusMuonMode = piPlusTable->GetDecayChannel(0);
	piPlusMuonMode->SetBR(1.0 - 1.23e-4);
	G4VDecayChannel* piMinusMuonMode = piMinusTable->GetDecayChannel(0);
	piMinusMuonMode->SetBR(1.0 - 1.23e-4);
	// Add the helicity suppressed modes
	G4VDecayChannel* piPlusPositronMode = new G4PhaseSpaceDecayChannel("pi+", 1.23e-4, 2, "e+", "nu_e");
	piPlusTable->Insert(piPlusPositronMode);
	G4VDecayChannel* piMinusElectronMode = new G4PhaseSpaceDecayChannel("pi-", 1.23e-4, 2, "e-", "anti_nu_e");
	piMinusTable->Insert(piMinusElectronMode);

	// Print the decay tables so we know we're getting it right
	muPlus->GetDecayTable()->DumpInfo();
	muMinus->GetDecayTable()->DumpInfo();
	piPlus->GetDecayTable()->DumpInfo();
	piMinus->GetDecayTable()->DumpInfo();
	kaonPlus->GetDecayTable()->DumpInfo();
	kaonMinus->GetDecayTable()->DumpInfo();
	kaonShort->GetDecayTable()->DumpInfo();
	kaonLong->GetDecayTable()->DumpInfo();

	// Setup the branches to get the parent tree variables
	TFile* inFile = new TFile(inFilename.c_str());	
	TTree* parentTree = (TTree*)inFile->Get("parentTree");
	int PDGCode, volumeID, numNeutrinos;
	double xMom, yMom, zMom, kinE, xPos, yPos, zPos;
	parentTree->SetBranchAddress("PDGCode", &PDGCode);
	parentTree->SetBranchAddress("xMom", &xMom);
	parentTree->SetBranchAddress("yMom", &yMom);
	parentTree->SetBranchAddress("zMom", &zMom);
	parentTree->SetBranchAddress("kinE", &kinE);
	parentTree->SetBranchAddress("xPos", &xPos);
	parentTree->SetBranchAddress("yPos", &yPos);
	parentTree->SetBranchAddress("zPos", &zPos);
	parentTree->SetBranchAddress("volumeID", &volumeID);
	parentTree->SetBranchAddress("numNeutrinos", &numNeutrinos);

	int numEntries = parentTree->GetEntries();
	std::cout << "Loading parent from file " << inFilename << std::endl;
	std::cout << numEntries << " parents will be decayed." << std::endl;
	std::cout << "Each parent will be decayed " << numDecays << " times." << std::endl;

	// Output neutrino property tree (saved to a separate file)
	TFile* outFile = new TFile(outFilename.c_str(), "RECREATE");
	int numNu;
	std::vector<int>* nuPDGCode = 0;
	std::vector<double>* nuEnergy = 0;
	std::vector<double>* xDirNu = 0;
	std::vector<double>* yDirNu = 0;
	std::vector<double>* zDirNu = 0;
	std::vector<double>* xPosNu = 0;
	std::vector<double>* yPosNu = 0;
	std::vector<double>* zPosNu = 0;
	std::vector<double>* weightPL = 0;

	TTree* neutrinoTree = new TTree("neutrinoTree", "");
	neutrinoTree->Branch("numNu", &numNu);
	neutrinoTree->Branch("PDGCode", &nuPDGCode);
	neutrinoTree->Branch("energy", &nuEnergy);
	neutrinoTree->Branch("xDir", &xDirNu);
	neutrinoTree->Branch("yDir", &yDirNu);
	neutrinoTree->Branch("zDir", &zDirNu);
	neutrinoTree->Branch("xPos", &xPosNu);
	neutrinoTree->Branch("yPos", &yPosNu);
	neutrinoTree->Branch("zPos", &zPosNu);
	neutrinoTree->Branch("weightPL", &weightPL);
	neutrinoTree->Branch("volumeID", &volumeID);

	std::cout << "Output neutrinoTree + plots will be written to " << outFilename << "\n";

	// Save the number of decays to the output file as a TParameter
	TParameter<int> numDecaysParam("NUM_DECAYS", numDecays);
	numDecaysParam.Write("", TObject::kOverwrite);

	int nbins=250;

	// Plots
	TH2D* xyIntersections = new TH2D("xyIntersections", ";x (cm);y (cm)", 230, -1.1*DET_RADIUS/cm, 1.1*DET_RADIUS/cm, 230, -1.1*DET_RADIUS/cm, 1.1*DET_RADIUS/cm);
	TH2D* xzIntersections = new TH2D("xzIntersections", ";x (cm);z (cm)", 230, -1.1*DET_RADIUS/cm, 1.1*DET_RADIUS/cm, 175, -1.1*DET_HALF_HEIGHT/cm, 1.1*DET_HALF_HEIGHT/cm);
	TH1D* cosThetaBeam = new TH1D("cosThetaBeam", ";cos#theta;", 100, 0.0, 1.0);
	TH1D* nuESpectrum_4pi = new TH1D("nuESpectrum_4pi", "#nu_{e};E_{#nu} (MeV);",nbins, 0.0, 500.0);
	TH1D* nuEBarSpectrum_4pi = new TH1D("nuEBarSpectrum_4pi", "#bar{#nu}_{e};E_{#nu} (MeV);", nbins, 0.0, 500.0);
	TH1D* nuMuSpectrum_4pi = new TH1D("nuMuSpectrum_4pi", "#nu_{#mu};E_{#nu} (MeV);", nbins, 0.0, 500.0);
	TH1D* nuMuBarSpectrum_4pi = new TH1D("nuMuBarSpectrum_4pi", "#bar{#nu}_{#mu};E_{#nu} (MeV);", nbins, 0.0, 500.0);
	TH1D* nuESpectrum_JSNS2 = new TH1D("nuESpectrum_JSNS2", "#nu_{e};E_{#nu} (MeV);", nbins, 0.0, 500.0);
	TH1D* nuEBarSpectrum_JSNS2 = new TH1D("nuEBarSpectrum_JSNS2", "#bar{#nu}_{e};E_{#nu} (MeV);", nbins, 0.0, 500.0);
	TH1D* nuMuSpectrum_JSNS2 = new TH1D("nuMuSpectrum_JSNS2", "#nu_{#mu};E_{#nu} (MeV);", nbins, 0.0, 500.0);
	TH1D* nuMuBarSpectrum_JSNS2 = new TH1D("nuMuBarSpectrum_JSNS2", "#bar{#nu}_{#mu};E_{#nu} (MeV);", nbins, 0.0, 500.0);

	TH1D* nuEVolume_4pi = new TH1D("nuEVolume_4pi", "nuEVolume_4pi",25,0,25);
	TH1D* nuEBarVolume_4pi = new TH1D("nuEBarVolume_4pi", "nuEBarVolume_4pi",25,0,25);
	TH1D* nuMuVolume_4pi = new TH1D("nuMuVolume_4pi", "nuMuVolume_4pi",25,0,25);
	TH1D* nuMuBarVolume_4pi = new TH1D("nuMuBarVolume_4pi", "nuMuBarVolume_4pi",25,0,25);

	static bool labelsInitialized = false;

	for (int iParent = 0; iParent < numEntries; iParent++) {
		if (iParent % 1000 == 0) std::cout << (double)iParent / (double)numEntries * 100.0 << "%" << std::endl;
		parentTree->GetEntry(iParent);

		// Clear the neutrino branch variables
		numNu = 0;
		nuPDGCode->clear();
		nuEnergy->clear();
		xDirNu->clear();
		yDirNu->clear();
		zDirNu->clear();
		xPosNu->clear();
		yPosNu->clear();
		zPosNu->clear();
		weightPL->clear();

		// Skip DAR if useDAR == false
		if (kinE == 0 && !useDAR) continue;

		// Set the parent properties and build the parent particle
		// Need to make sure the normalized direction is set correctly
		double momNorm = sqrt(xMom*xMom + yMom*yMom + zMom*zMom);
		G4ThreeVector parentDirection;
		if (momNorm == 0) {
			parentDirection.setX(0.0);
			parentDirection.setY(0.0);
			parentDirection.setZ(0.0);
		}
		else {
			parentDirection.setX(xMom/momNorm);
			parentDirection.setY(yMom/momNorm);
			parentDirection.setZ(zMom/momNorm);
		}

		G4DynamicParticle parent;
		if (PDGCode == -13) {
			parent = G4DynamicParticle(muPlus, parentDirection, kinE); 
			numMuPlus++;
		}
		else if (PDGCode == 13) {
			parent = G4DynamicParticle(muMinus, parentDirection, kinE); 
			numMuMinus++;
		}
		else if (PDGCode == 211) {
			parent = G4DynamicParticle(piPlus, parentDirection, kinE); 
			numPiPlus++;
		}
		else if (PDGCode == -211) {
			parent = G4DynamicParticle(piMinus, parentDirection, kinE); 
			numPiMinus++;
		}
		else if (PDGCode == 321) {
			parent = G4DynamicParticle(kaonPlus, parentDirection, kinE); 
			numKPlus++;
		}
		else if (PDGCode == -321) {
			parent = G4DynamicParticle(kaonMinus, parentDirection, kinE); 
			numKMinus++;
		}
		else if (PDGCode == 130) {
			parent = G4DynamicParticle(kaonLong, parentDirection, kinE); 
			numKLong++;
		}
		else if (PDGCode == 310) {
			parent = G4DynamicParticle(kaonShort, parentDirection, kinE); 
			numKShort++;
		}
		else {
			std::cout << "Unknown parent particle; PDG = " << PDGCode << std::endl;
			continue;
		}

		G4double parentEnergy = parent.GetTotalEnergy();
		G4ThreeVector pos = G4ThreeVector(xPos, yPos, zPos);

		// Transform to detector coordinates from accelerator coordinates
		pos.rotateX(-90.0*deg);
		pos.rotateY(180.0*deg);
		pos -= G4ThreeVector(0.0, DET_Y, DET_Z);

		// Get the decay table for this parent particle
		G4DecayTable* decayTable = decayTable = parent.GetParticleDefinition()->GetDecayTable();

		// Repeatedly decay the parent and fill the output tree with the properties
		for (int iDecay = 0; iDecay < numDecays; iDecay++) {
			// Neutrino variables
			std::vector<int> nuPDG;
			std::vector<G4ThreeVector> nuDir;
			std::vector<double> energyNu;	

			// Check if this was a mu- that captured instead of decaying
			// In that case, we do the capture process
			if (PDGCode == 13 && numNeutrinos == 1) {
				const G4HadProjectile hadProj(parent);
				// Lookup the material for the volume this muon captured in
				G4Material* material = 0;	

				if (!labelsInitialized) {

					std::set<int> usedID;
					std::vector<G4VPhysicalVolume*>::iterator iter;

					for (iter = vStore->begin(); iter != vStore->end(); ++iter) {

						int id = geo->GetVolumeID((*iter)->GetName());

						if (usedID.count(id)) continue;

						usedID.insert(id);

						const char* name = (*iter)->GetName().c_str();

						nuEVolume_4pi    ->GetXaxis()->SetBinLabel(id+1, name);
						nuEBarVolume_4pi ->GetXaxis()->SetBinLabel(id+1, name);
						nuMuVolume_4pi   ->GetXaxis()->SetBinLabel(id+1, name);
						nuMuBarVolume_4pi->GetXaxis()->SetBinLabel(id+1, name);
					}

					labelsInitialized = true;
				}

				std::vector<G4VPhysicalVolume*>::iterator iter;
				for (iter = vStore->begin(); iter != vStore->end(); ++iter) {
					if (geo->GetVolumeID((*iter)->GetName()) == volumeID) {
						material = (*iter)->GetLogicalVolume()->GetMaterial();
						//	std::cout<<(*iter)->GetName()<<std::endl;
					}
				}
				// Do the capture on a nucleus from the volume
				G4Nucleus nuc = G4Nucleus(material);
				G4HadFinalState* finalState = muCapture.ApplyYourself(hadProj, nuc);

				// Look for the muon neutrino in the final state
				for (int i = 0; i < finalState->GetNumberOfSecondaries(); ++i) {
					const G4HadSecondary* hadSec = finalState->GetSecondary(i);
					const G4DynamicParticle* part = hadSec->GetParticle();
					int PDG = part->GetPDGcode();
					if (PDG == 14) {
						nuPDG.push_back(PDG);
						energyNu.push_back(part->GetTotalEnergy());
						nuDir.push_back(part->GetMomentumDirection());
						nuMuSpectrum_4pi->Fill(part->GetTotalEnergy());
						nuMuVolume_4pi->Fill(volumeID);
					}
				}
			}

			// Otherwise, we're going to do a decay
			// Choose a decay channel and decay the parent
			else {
				G4VDecayChannel* decayChannel = decayTable->SelectADecayChannel();
				G4DecayProducts* products = decayChannel->DecayIt(parent.GetMass());

				// Boost the decay products into the lab frame
				products->Boost(parentEnergy, parentDirection);

				// Loop over the daughters and find the neutrino(s)
				G4int numDaughters = products->entries();
				for (int iDaughter = 0; iDaughter < numDaughters; iDaughter++) {
					G4DynamicParticle* daughter = products->PopProducts();
					int PDG = daughter->GetPDGcode();

					// Look for neutrino daughters
					if (PDG == 12 || PDG == -12 || PDG == 14 || PDG == -14) {	
						nuPDG.push_back(PDG);
						energyNu.push_back(daughter->GetTotalEnergy());
						nuDir.push_back(daughter->GetMomentumDirection());

						if (PDG == 12){
							nuESpectrum_4pi->Fill(daughter->GetTotalEnergy());
							nuEVolume_4pi->Fill(volumeID);							
						}
						else if (PDG == -12){
							nuEBarSpectrum_4pi->Fill(daughter->GetTotalEnergy());
							nuEBarVolume_4pi->Fill(volumeID);							
						}
						else if (PDG == 14){
							nuMuSpectrum_4pi->Fill(daughter->GetTotalEnergy());
							nuMuVolume_4pi->Fill(volumeID);							
						}
						else if (PDG == -14){
							nuMuBarSpectrum_4pi->Fill(daughter->GetTotalEnergy());
							nuMuBarVolume_4pi->Fill(volumeID);							
						}
						delete daughter; // Only you can prevent memory leaks
					}
				}
				delete products; // Only you can prevent memory leaks
			}

			// If there wasn't a neutrino, continue
			if (nuPDG.size() == 0) continue;

			for (int iNu = 0; iNu < nuPDG.size(); ++iNu) {
				int PDG = nuPDG[iNu];
				double energy = energyNu[iNu];
				G4ThreeVector dir = nuDir[iNu];

				// Transform the neutrino direction into the detector coordinates	
				dir.rotateX(-90.0*deg);
				dir.rotateY(180.0*deg);

				// Skip neutrinos which are going down or backwards
				if (dir.z() < 0.0 || dir.y() < 0.0) continue;	

				// Check if the neutrino will pass through the detector
				// First, look to see if the neutrino has intersections in the xy-plane
				double temp1 = pos.x()*dir.x() + pos.y()*dir.y();
				double temp2 = dir.x()*dir.x() + dir.y()*dir.y();
				double temp3 = pos.x()*pos.x() + pos.y()*pos.y() - DET_RADIUS*DET_RADIUS;
				double discrim = temp1*temp1 - temp2*temp3;
				if (discrim <= 0.0) continue;

				// Next, check the z-coordinates of those intersections with the xy-plane
				double lambda1 = (-temp1 + sqrt(discrim))/temp2;
				double lambda2 = (-temp1 - sqrt(discrim))/temp2;
				double z1 = pos.z() + lambda1*dir.z();
				double z2 = pos.z() + lambda2*dir.z();

				// Skip neutrinos which intersect the cylinder completely above or below the detector
				if ( (z1 < -DET_HALF_HEIGHT && z2 < -DET_HALF_HEIGHT) ||
						(z1 > DET_HALF_HEIGHT && z2 > DET_HALF_HEIGHT) ) {
					continue;
				}

				// If we've made it here, we have a valid set of intersection points. We just
				// need to find them. We can correct intersections above or below the detector
				double x1, y1 = 0.0;

				// Case 1: z-coordinate is above the detector
				if (z1 > DET_HALF_HEIGHT) {
					z1 = DET_HALF_HEIGHT;
					lambda1 = (DET_HALF_HEIGHT - pos.z())/dir.z();
					x1 = pos.x() + lambda1*dir.x();
					y1 = pos.y() + lambda1*dir.y();
				}
				// Case 2: z-coordinate is below the detector
				else if (z1 < -DET_HALF_HEIGHT) {
					z1 = -DET_HALF_HEIGHT;
					lambda1 = (-DET_HALF_HEIGHT - pos.z())/dir.z();
					x1 = pos.x() + lambda1*dir.x();
					y1 = pos.y() + lambda1*dir.y();
				}
				// Case 3: z-coordinate is within the detector bounds
				else {
					x1 = pos.x() + lambda1*dir.x();
					y1 = pos.y() + lambda1*dir.y();
				}

				// Repeated for the second intersection point
				double x2, y2 = 0.0;

				// Case 1: z-coordinate is above the detector
				if (z2 > DET_HALF_HEIGHT) {
					z2 = DET_HALF_HEIGHT;
					lambda2 = (DET_HALF_HEIGHT - pos.z())/dir.z();
					x2 = pos.x() + lambda2*dir.x();
					y2 = pos.y() + lambda2*dir.y();
				}
				// Case 2: z-coordinate is below the detector
				else if (z2 < -DET_HALF_HEIGHT) {
					z2 = -DET_HALF_HEIGHT;
					lambda2 = (-DET_HALF_HEIGHT - pos.z())/dir.z();
					x2 = pos.x() + lambda2*dir.x();
					y2 = pos.y() + lambda2*dir.y();
				}
				// Case 3: z-coordinate is within the detector bounds
				else {
					x2 = pos.x() + lambda2*dir.x();
					y2 = pos.y() + lambda2*dir.y();
				}

				numNu++;
				nuPDGCode->push_back(PDG);
				nuEnergy->push_back(energy);
				xDirNu->push_back(dir.x());
				yDirNu->push_back(dir.y());
				zDirNu->push_back(dir.z());

				double pathLength = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
				weightPL->push_back(pathLength/MAX_PATH_LENGTH);

				// Save the lower (smaller z) intersection point
				if (z1 < z2) {
					xPosNu->push_back(x1);
					yPosNu->push_back(y1);
					zPosNu->push_back(z1);
				}
				else {
					xPosNu->push_back(x2);
					yPosNu->push_back(y2);
					zPosNu->push_back(z2);
				}

				// Fill the histograms
				xyIntersections->Fill(x1/cm, y1/cm);
				xyIntersections->Fill(x2/cm, y2/cm);
				xzIntersections->Fill(x1/cm, z1/cm);
				xzIntersections->Fill(x2/cm, z2/cm);
				cosThetaBeam->Fill(dir.y()/sqrt(dir.y()*dir.y() + dir.z()*dir.z()));
				if (PDG == 12) nuESpectrum_JSNS2->Fill(energy);
				else if (PDG == -12) nuEBarSpectrum_JSNS2->Fill(energy);
				else if (PDG == 14) nuMuSpectrum_JSNS2->Fill(energy);
				else if (PDG == -14) nuMuBarSpectrum_JSNS2->Fill(energy);
			} // End of loop over neutrinos	

		} // End loop over decays
		neutrinoTree->Fill();

	} // End loop over parents

	std::cout << "===========================" << std::endl;
	std::cout << "     Number of Parents     " << std::endl;
	std::cout << "===========================" << std::endl;
	std::cout << "mu+: " << numMuPlus << std::endl;
	std::cout << "mu-: " << numMuMinus << std::endl;
	std::cout << "pi+: " << numPiPlus << std::endl;
	std::cout << "pi-: " << numPiMinus << std::endl;
	std::cout << "K+: " << numKPlus << std::endl;
	std::cout << "K-: " << numKMinus << std::endl;
	std::cout << "K_S: " << numKShort << std::endl;
	std::cout << "K_L: " << numKLong << std::endl;
	std::cout << "===========================" << std::endl;

	// Close up shop; Overwrite all old objects in case we've run
	// the decay code on this file befor and they're already in the file
	std::cout << "Writing output.\n";
	outFile->cd();
	xyIntersections->Write("", TObject::kOverwrite);	
	xzIntersections->Write("", TObject::kOverwrite);
	cosThetaBeam->Write("", TObject::kOverwrite);
	nuESpectrum_4pi->Write("", TObject::kOverwrite);	
	nuEBarSpectrum_4pi->Write("", TObject::kOverwrite);	
	nuMuSpectrum_4pi->Write("", TObject::kOverwrite);	
	nuMuBarSpectrum_4pi->Write("", TObject::kOverwrite);	
	nuESpectrum_JSNS2->Write("", TObject::kOverwrite);	
	nuEBarSpectrum_JSNS2->Write("", TObject::kOverwrite);	
	nuMuSpectrum_JSNS2->Write("", TObject::kOverwrite);	
	nuMuBarSpectrum_JSNS2->Write("", TObject::kOverwrite);	
	nuEVolume_4pi->Write("", TObject::kOverwrite);
	nuEBarVolume_4pi->Write("", TObject::kOverwrite);
	nuMuVolume_4pi->Write("", TObject::kOverwrite);
	nuMuBarVolume_4pi->Write("", TObject::kOverwrite);

	neutrinoTree->Write("", TObject::kOverwrite);
	outFile->Close();
	inFile->Close();	

	return 0;
}
