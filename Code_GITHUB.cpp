#include "TFile.h"
#include "TTree.h"
#include "TBox.h"
#include <TF1.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <iostream>
#include <limits>
#include <vector>
#include <fstream>
#include <string>
//#include <sstream>
//#include <map>
#include <algorithm>
//#include <TGraph.h>
//#include <TGraphErrors.h>
//#include <climits>
#include <cmath>
#include <iomanip>

//#include <locale>
//#include <codecvt>

using namespace std;

#include "sndisplay.cc"

// TREE 1

TFile* run_file = nullptr;
TTree* run_tree = nullptr;

// BRANCHES

vector<int>* digitracker_side = new vector<int>;
vector<int>* digitracker_column = new vector<int>;
vector<int>* digitracker_layer = new vector<int>;
vector<double>* digitracker_anodetimestampR0 = new vector<double>; // long !
vector<long>* digitracker_bottomcathodetimestamp = new vector<long>;
vector<long>* digitracker_topcathodetimestamp = new vector<long>;

vector<int>* digicalo_type = new vector<int>;
vector<int>* digicalo_side = new vector<int>;
vector<int>* digicalo_column = new vector<int>;
vector<int>* digicalo_row = new vector<int>;
vector<int>* digicalo_wall = new vector<int>;
vector<double>* digicalo_timestamp = new vector<double>; // long !
vector<bool>* digicalo_lowthresholdonly = new vector<bool>;
vector<bool>* digicalo_highthreshold = new vector<bool>;
vector<short>* digicalo_peakamplitude = new vector<short>;
vector<int>* digicalo_charge = new vector<int>;
vector<int>* digicalo_falling_cell = new vector<int>;

// TREE 2

TFile* class_file = nullptr;
TTree* class_tree = nullptr;

// BRANCHES

vector<double>* dt_good_corr = new vector<double>;
vector<int>* calo_corr_ij = new vector<int>;
vector<int>* index_un_calo = new vector<int>;
int nb_total_elec = -1;
int nb_total_alpha = -1;
int nb_total_weirdies = -1;
bool candidate_2e = false;
bool pure_event = false;
vector<int>* calo_tablo = new vector<int>;

// structure for cluster storage

struct digicluster
{
	// indexes of tracker hit in UDD tree
	vector<int> tracker_indexes;

	// indexes of calorimeter hit in UDD tree
	vector<int> calo_indexes;
	vector<int> associated_calo_indexes;   // electron-like hit 
	vector<int> unassociated_calo_indexes; // gamma-like hit

	// extra tracker information
	int index_tracker_layer_min;
	int column_at_flag;
	int tracker_side;
	int tracker_column_min;
	int tracker_column_max;
	int tracker_layer_min;
	int tracker_layer_max;
	int test_flag;
	int tracker_flag;
	// tracker_flag: bit-wise storage of boolean information
	//  0001 = 1  ->  (column_min,layer_min) corner populated
	//  0010 = 2  ->  (column_min,layer_max) corner populated
	//  0100 = 4  ->  (column_max,layer_min) corner populated
	//  1000 = 8  ->  (column_max,layer_max) corner populated

	double tracker_first_anode_time;
};

vector<digicluster>* digiclusters = new vector<digicluster>;

void open_run_file(const char* run_filename)
{
	// delete the TFile if one was already opened
	if (run_file != nullptr)
		delete run_file;

	// open the file
	run_file = TFile::Open(run_filename, "READ");

	// quit if the file does not exist
	if (!run_file->IsOpen())
	{
		cout << "BOUH" << endl;
		return;
	}

	// retrieve the TTree
	run_tree = (TTree*)(run_file->Get("SimData"));

	// quit if 'SimData' was not found
	if (run_tree == nullptr)
		return;

	// enable loading of only wished branches (for I/O speed performances)
	run_tree->SetBranchStatus("*", 0);
	run_tree->SetBranchStatus("digitracker.side", 1);
	run_tree->SetBranchStatus("digitracker.column", 1);
	run_tree->SetBranchStatus("digitracker.layer", 1);
	run_tree->SetBranchStatus("digitracker.anodetimestampR0", 1);
	run_tree->SetBranchStatus("digitracker.bottomcathodetimestamp", 1);
	run_tree->SetBranchStatus("digitracker.topcathodetimestamp", 1);
	run_tree->SetBranchStatus("digicalo.type", 1);
	run_tree->SetBranchStatus("digicalo.side", 1);
	run_tree->SetBranchStatus("digicalo.column", 1);
	run_tree->SetBranchStatus("digicalo.row", 1);
	run_tree->SetBranchStatus("digicalo.wall", 1);
	run_tree->SetBranchStatus("digicalo.timestamp", 1);
	run_tree->SetBranchStatus("digicalo.lowthresholdonly", 1);
	run_tree->SetBranchStatus("digicalo.highthreshold", 1);
	run_tree->SetBranchStatus("digicalo.peakamplitude", 1);
	run_tree->SetBranchStatus("digicalo.charge", 1);
	run_tree->SetBranchStatus("digicalo.falling_cell", 1);
	run_tree->SetBranchStatus("digicalo.highthreshold", 1);

	// setup branch addresses
	run_tree->SetBranchAddress("digitracker.side", &digitracker_side);
	run_tree->SetBranchAddress("digitracker.column", &digitracker_column);
	run_tree->SetBranchAddress("digitracker.layer", &digitracker_layer);
	run_tree->SetBranchAddress("digitracker.anodetimestampR0", &digitracker_anodetimestampR0);
	run_tree->SetBranchAddress("digitracker.bottomcathodetimestamp", &digitracker_bottomcathodetimestamp);
	run_tree->SetBranchAddress("digitracker.topcathodetimestamp", &digitracker_topcathodetimestamp);
	run_tree->SetBranchAddress("digicalo.type", &digicalo_type);
	run_tree->SetBranchAddress("digicalo.side", &digicalo_side);
	run_tree->SetBranchAddress("digicalo.column", &digicalo_column);
	run_tree->SetBranchAddress("digicalo.row", &digicalo_row);
	run_tree->SetBranchAddress("digicalo.wall", &digicalo_wall);
	run_tree->SetBranchAddress("digicalo.timestamp", &digicalo_timestamp);
	run_tree->SetBranchAddress("digicalo.lowthresholdonly", &digicalo_lowthresholdonly);
	run_tree->SetBranchAddress("digicalo.highthreshold", &digicalo_highthreshold);
	run_tree->SetBranchAddress("digicalo.peakamplitude", &digicalo_peakamplitude);
	run_tree->SetBranchAddress("digicalo.charge", &digicalo_charge);
	run_tree->SetBranchAddress("digicalo.falling_cell", &digicalo_falling_cell);
	run_tree->SetBranchAddress("digicalo.highthreshold", &digicalo_highthreshold);

	// print a message
	cout << "=> run file \"" << run_filename << "\" opened (" << run_tree->GetEntries() << " events)" << endl;
}

void open_class_file(const char* run_filename) // S'ASSURER QUE L'ARBRE A BIEN ETE CREE : PARCE QUE C'EST UNE FONCTION PLUS BAS QUI LE CREE...
{
	// delete the TFile if one was already opened
	if (class_file != nullptr)
		delete class_file;

	// open the file
	class_file = TFile::Open(run_filename, "READ");

	// quit if the file does not exist
	if (!class_file->IsOpen())
	{
		cout << "BOUH" << endl;
		return;
	}

	// retrieve the TTree
	class_tree = (TTree*)(class_file->Get("classification"));

	// quit if 'SimData' was not found
	if (class_tree == nullptr)
		return;

	// setup branch addresses
	class_tree->SetBranchAddress("Nb.of.electrons", &nb_total_elec);
	class_tree->SetBranchAddress("Nb.of.alpha", &nb_total_alpha);
	class_tree->SetBranchAddress("Nb.of.weirdies", &nb_total_weirdies);
	class_tree->SetBranchAddress("Calo.gamma", &index_un_calo);
	class_tree->SetBranchAddress("candidates.2e", &candidate_2e);
	class_tree->SetBranchAddress("dt.correlation", &dt_good_corr);
	class_tree->SetBranchAddress("calo.correlation", &calo_corr_ij);
	class_tree->SetBranchAddress("pure.event", &pure_event);

	// print a message
	cout << "=> run file \"" << run_filename << "\" opened (" << class_tree->GetEntries() << " events)" << endl;
}

// compute cellnum for tracker hit of index i
int digitracker_cellnum(int i)
{
	return digitracker_side->at(i) * 9 * 113 + digitracker_column->at(i) * 9 + digitracker_layer->at(i);
}

// compute omnum for calo hit of index i
int digicalo_omnum(int i)
{
	if (digicalo_type->at(i) == 0) // main-wall case
		return 260 * digicalo_side->at(i) + 13 * digicalo_column->at(i) + digicalo_row->at(i);

	if (digicalo_type->at(i) == 1) // x-wall case
		return 520 + 64 * digicalo_side->at(i) + 32 * digicalo_wall->at(i) + 16 * digicalo_column->at(i) + digicalo_row->at(i);

	if (digicalo_type->at(i) == 2) // g-veto case
		return 648 + 32 * digicalo_side->at(i) + 16 * digicalo_wall->at(i) + digicalo_column->at(i);

	// we should not reach this part
	return -1;
}

// Fonction pour lire les corrections à partir d'un fichier et créer un vecteur de corrections
vector<double> readAndCreateCorrectionVector(const string& filename) {
	ifstream rfile(filename);
	if (!rfile.is_open()) {
		cerr << "Unable to open file" << endl;
		return {};
	}

	vector<pair<int, double>> corrections;
	int maxIndex = 0;

	string line;
	while (getline(rfile, line)) {
		istringstream iss(line);
		int number;
		double value;
		if (iss >> number >> value) {
			corrections.emplace_back(number, value);
			maxIndex = max(maxIndex, number);
		}
		else {
			cerr << "Error parsing line: " << line << endl;
		}
	}

	rfile.close();

	vector<double> correctionVector(maxIndex + 1, 0.0);
	for (const auto& correction : corrections) {
		correctionVector[correction.first] = correction.second;
	}

	return correctionVector;
}

vector<double> readAndCreateEnergyVector(const string& filename) {
	ifstream rfile(filename);
	if (!rfile.is_open()) {
		cerr << "Unable to open file" << endl;
		return {};
	}

	vector<pair<int, double>> energies;
	int maxIndex = 0;

	string line;
	while (getline(rfile, line)) {
		istringstream iss(line);
		int number;
		double value;
		if (iss >> number >> value) {
			energies.emplace_back(number, value);
			maxIndex = max(maxIndex, number);
		}
		else {
			cerr << "Error parsing line: " << line << endl;
		}
	}

	rfile.close();

	vector<double> energyVector(maxIndex + 1, 0.0);
	for (const auto& correction : energies) {
		energyVector[correction.first] = correction.second;
	}

	return energyVector;
}

// spatial clusterisation of tracker hit within a radius of sqrt(6)
// and time clusterisation with delta_anode_time within 15 us
//
//  X X X X X X X
//  X X O O O X X  [#] hit tracker_i
//  X O O O O O X
//  X O O # O O X  [O] hit tracker_j clusterised with tracker_i
//  X O O O O O X
//  X X O O O X X
//  X X X X X X X
//
const double tracker_cluster_r2_threshold = 6.0;   // in cell^2 unit
const double tracker_cluster_dt_threshold = 15E-6; // in seconds

const double calo_amplitude_threshold = 50.; // in mV

// clusterisation of tracker hit with calo hit is based on deltat
// between the first anode time of the cluste and calo time within:
const double cluster_deltat_anode_calo_min = -0.22E-6; // in seconds
const double cluster_deltat_anode_calo_max = +5.00E-6; // in seconds

// spatial association of calo hit with tracker hit :
//
//   effective tracker row centered
//     in front of the calo hit:
//                |
//                v
// ---+-------+-------+-------+---
//    |       |   '   |       |
//    |       |   '   |       |
//    |       |   '   |       |
// ---+-------+-------+-------+---
//    X X X'O O O O O O O'X X X X     8
//    X X X'O O O O O O O'X X X X     7
//    X X X'X X X X X X X'X X X X     6
//         '      '      '          layer
//         '<-----'----->'
//                '
// associate the calo hit to the tracker hit if the cell is
// within the two last layer and within a range of +/-5 of rows
const float cluster_delta_row_tracker_calo = 5;

sndisplay::demonstrator* event_display = new sndisplay::demonstrator;

vector<TBox*> cluster_boxes;

void event_clusterisation()
{
	// we assume run_tree->GetEntry(x) was already called

	digiclusters->clear();

	// temporary storage of UDD's indexes of clusterised tracker hit:
	// clusters_tracker_indexes[CLUSTER][CELL]
	vector<vector<int>> clusters_tracker_indexes;

	const int nb_calo_hits = digicalo_side->size();
	const int nb_tracker_hits = digitracker_side->size();

	// main tracker hit loop
	for (int tracker_i = 0; tracker_i < nb_tracker_hits; tracker_i++)
	{
		// the tracker hit must have an anode R0 time
		if (digitracker_anodetimestampR0->at(tracker_i) <= 0)
			continue;

		// storage of cluster indexes for which the hit tracker_i matches
		// space/time criteria for at least 1 cell of thise cluster
		vector<int> matching_clusters_id;

		const int nb_clusters = clusters_tracker_indexes.size();

		// iterate over all cluster already existing
		for (int cluster_id = 0; cluster_id < nb_clusters; cluster_id++)
		{
			// iterate over all cells in cluster_id
			for (const int& tracker_j : clusters_tracker_indexes.at(cluster_id))
			{
				// retrieve cell sides skip if they are different
				const bool digitracker_side_i = digitracker_side->at(tracker_i);
				const bool digitracker_side_j = digitracker_side->at(tracker_j);

				if (digitracker_side_i != digitracker_side_j)
					continue;

				// compute delta_radius^2 and skip if the space criteria failed
				const int digitracker_column_i = digitracker_column->at(tracker_i);
				const int digitracker_column_j = digitracker_column->at(tracker_j);
				const int digitracker_delta_column = digitracker_column_j - digitracker_column_i;

				const int digitracker_layer_i = digitracker_layer->at(tracker_i);
				const int digitracker_layer_j = digitracker_layer->at(tracker_j);
				const int digitracker_delta_layer = digitracker_layer_j - digitracker_layer_i;

				const double digitracker_delta_r2 = pow(digitracker_delta_column, 2) + pow(digitracker_delta_layer, 2);

				if (digitracker_delta_r2 > tracker_cluster_r2_threshold)
					continue;

				// compute delta_anode_time and skip if the time criteria failed
				const double digitracker_anodetime_i = digitracker_anodetimestampR0->at(tracker_i) * 12.5E-9;
				const double digitracker_anodetime_j = digitracker_anodetimestampR0->at(tracker_j) * 12.5E-9;
				const double digitracker_delta_anodetime = digitracker_anodetime_j - digitracker_anodetime_i;

				if (TMath::Abs(digitracker_delta_anodetime) > tracker_cluster_dt_threshold)
					continue;

				// the hit #tracker_i succedded space+time correlation with the #tracker_j (taken
				// from the cluster #cluster_id), so we add this cluster in the matching list:
				matching_clusters_id.push_back(cluster_id);

				// no need to check further cells from this cluster, so we stop the tracker_j loop 
				break;

			} // for (tracker_j)

		} // for (cluster_id)


		if (matching_clusters_id.size() == 0)
		{
			// the hit tracker_i does not match with any existing clusters:
			// => we create a new cluster and fill it with its tracker index
			vector<int> new_cluster_tracker_indexes;
			new_cluster_tracker_indexes.push_back(tracker_i);

			// add this new cluster in the cluster storage
			clusters_tracker_indexes.push_back(new_cluster_tracker_indexes);
		}

		else if (matching_clusters_id.size() == 1)
		{
			// the hit tracker_i does match an UNIQUE existing cluster
			// => we append it into this cluster
			clusters_tracker_indexes[matching_clusters_id.front()].push_back(tracker_i);
		}

		else
		{
			// the hit tracker_i is matching several existing clusters
			// => we append it into the first matching cluster
			clusters_tracker_indexes[matching_clusters_id.front()].push_back(tracker_i);

			// => we merge all other clusters into the first one
			for (int c = matching_clusters_id.size() - 1; c > 0; --c)
			{
				const int cluster_id = matching_clusters_id.at(c);

				// append all tracker indexes from the cluster into first one
				for (const int& tracker_index : clusters_tracker_indexes[cluster_id])
					clusters_tracker_indexes[matching_clusters_id.front()].push_back(tracker_index);

				// delete the cluster from the storage
				clusters_tracker_indexes.erase(clusters_tracker_indexes.begin() + cluster_id);
			}
		}

	} // for (tracker_i)

  // now we iterate over cluster identified to fill the vector<digicluster>

	const int nb_clusters = clusters_tracker_indexes.size();

	for (int cluster_id = 0; cluster_id < nb_clusters; cluster_id++)
	{
		digicluster new_digicluster;

		const vector<int>& cluster_tracker_indexes = clusters_tracker_indexes.at(cluster_id);
		new_digicluster.tracker_indexes = cluster_tracker_indexes;

		// initialise cluster's extra informations
		new_digicluster.tracker_side = digitracker_side->at(new_digicluster.tracker_indexes.front());
		new_digicluster.tracker_column_min = 113;
		new_digicluster.tracker_column_max = -1;
		new_digicluster.tracker_layer_min = 9;
		new_digicluster.tracker_layer_max = -1;
		new_digicluster.column_at_flag = -1;
		new_digicluster.tracker_flag = 0;
		new_digicluster.test_flag = 0;
		new_digicluster.tracker_first_anode_time = numeric_limits<double>::max();

		for (const int& tracker_index : new_digicluster.tracker_indexes)
		{
			const int& tracker_column = digitracker_column->at(tracker_index);
			if (tracker_column < new_digicluster.tracker_column_min)
				new_digicluster.tracker_column_min = tracker_column;
			if (tracker_column > new_digicluster.tracker_column_max)
				new_digicluster.tracker_column_max = tracker_column;

			const int& tracker_layer = digitracker_layer->at(tracker_index);
			if (tracker_layer < new_digicluster.tracker_layer_min)
				new_digicluster.tracker_layer_min = tracker_layer;
			if (tracker_layer > new_digicluster.tracker_layer_max)
				new_digicluster.tracker_layer_max = tracker_layer;

			const double digitracker_anodetime = digitracker_anodetimestampR0->at(tracker_index) * 12.5E-9;
			if (digitracker_anodetime < new_digicluster.tracker_first_anode_time)
				new_digicluster.tracker_first_anode_time = digitracker_anodetime;

		} // for (tracker_index)

		  // new iteration on cluster's tracker hits to analyse population
		  // of box corners (now that we know min/max of column and layer)

		for (const int& tracker_index : new_digicluster.tracker_indexes)
		{
			const int& tracker_column = digitracker_column->at(tracker_index);
			const int& tracker_layer = digitracker_layer->at(tracker_index);

			if ((tracker_column == new_digicluster.tracker_column_min) && (tracker_layer == new_digicluster.tracker_layer_min))
			{
				new_digicluster.tracker_flag |= (1 << 0); // (column_min, layer_min) populated
				new_digicluster.column_at_flag = tracker_column; /////////////////////////////// ajout !
			}

			if ((tracker_column == new_digicluster.tracker_column_min) && (tracker_layer == new_digicluster.tracker_layer_max))
				new_digicluster.tracker_flag |= (1 << 1); // (column_min, layer_max) populated

			if ((tracker_column == new_digicluster.tracker_column_max) && (tracker_layer == new_digicluster.tracker_layer_min))
			{
				new_digicluster.tracker_flag |= (1 << 2); // (column_max, layer_min) populated
				new_digicluster.column_at_flag = tracker_column; /////////////////////////////// ajout !
			}

			if ((tracker_column == new_digicluster.tracker_column_max) && (tracker_layer == new_digicluster.tracker_layer_max))
				new_digicluster.tracker_flag |= (1 << 3); // (column_max, layer_max) populated

		}

		const int corner_flag = new_digicluster.tracker_flag & 0xF;
		/*if ((corner_flag == 0x6 || corner_flag == 0x9))
			new_digicluster.test_flag = 1;*/

		if (new_digicluster.tracker_column_max - new_digicluster.tracker_column_min > 1 && (corner_flag == 0x6 || corner_flag == 0x9))
			new_digicluster.test_flag = 1;

		if (new_digicluster.tracker_column_max - new_digicluster.tracker_column_min == 1 && (corner_flag == 0x6 || corner_flag == 0x9 || corner_flag == 0xE || corner_flag == 0xD || corner_flag == 0xB || corner_flag == 0x7))
		{
			new_digicluster.test_flag = 1;
			//cout << new_digicluster.tracker_flag << " tracker !" << endl;
		}

		if (new_digicluster.tracker_column_max - new_digicluster.tracker_column_min == 0)
		{
			new_digicluster.test_flag = 1;
		}

		///////////////////////////////////// CONDITION 1

		// now we iterate on calorimeter hits to look for time/space correlation with the tracker cluster

		for (int calo_i = 0; calo_i < nb_calo_hits; calo_i++)
		{
			const double calo_amplitude = -(2500. / 4096.) * (1. / 8.) * digicalo_peakamplitude->at(calo_i);

			// skip calo hit if bellow amplitude threshold
			if (calo_amplitude < calo_amplitude_threshold)
				continue;

			const double calo_time = digicalo_timestamp->at(calo_i) * 6.25E-9 - 250E-9 + digicalo_falling_cell->at(calo_i) * (1. / 256.) * (400.E-9 / 1024.);

			// perform calo/tracker time correlation
			const double deltat_anode_calo = new_digicluster.tracker_first_anode_time - calo_time;
			if (deltat_anode_calo > cluster_deltat_anode_calo_max) continue;
			if (deltat_anode_calo < cluster_deltat_anode_calo_min) continue;

			//////////////////////////////////////////////////////
			int tracker_max_time = 0;
			for (const int& tracker_index : new_digicluster.tracker_indexes)
			{
				if (digitracker_anodetimestampR0->at(tracker_index) > tracker_max_time)
				{
					tracker_max_time = digitracker_anodetimestampR0->at(tracker_index);
					//cout << tracker_max_time << " at " << tracker_index << endl;
				}
			}

			// the calo hit is intime with the tracker cluster
			new_digicluster.calo_indexes.push_back(calo_i);

			// perform calo/tracker spatial correlation
			bool calo_is_associated = false;

			//////////////////// consider main-wall calo hit !!!!!!!!!!!
			if (digicalo_type->at(calo_i) == 0)
			{

				// search tracker association on same side only
				if (digicalo_side->at(calo_i) == new_digicluster.tracker_side)
				{
					// convert the (middle of) calorimeter column into tracker column
					const float calo2tracker_column = 0.667 + 5.8771579 * digicalo_column->at(calo_i); // column association for main-wall
					// x-wall requires a layer association !

					for (const int& tracker_i : new_digicluster.tracker_indexes)
					{
						// search tracker association within 2 last layers [7-8]
						if (digitracker_layer->at(tracker_i) < 7)
							continue;

						const int& tracker_column = digitracker_column->at(tracker_i);
						if (TMath::Abs(tracker_column - calo2tracker_column) > cluster_delta_row_tracker_calo)
							continue;

						calo_is_associated = true;
						break;
					}

				} // (if calo/tracker same side)

				if (calo_is_associated)
					new_digicluster.associated_calo_indexes.push_back(calo_i);
				else
					new_digicluster.unassociated_calo_indexes.push_back(calo_i);
			}

		} // for (calo_i)

		//if ( digiclusters->at(c).associated_calo_indexes.size() == 0)

	  //////////////////////////////////
	  // CONDITIOOOOOOOOONS ///////
	  //////////////////////////////////

		//cout << new_digicluster.tracker_indexes.size() << endl;

	/*for (int i = 0; i < new_digicluster.associated_calo_indexes.size(); i++)
		cout << new_digicluster.associated_calo_indexes.at(i) << endl;*/

		// append the cluster
		digiclusters->push_back(new_digicluster);

	} // for (cluster_id)

}

void tree_new_class_clusterised_event() // ASSOCIER LES FICHIERS VERSION GITHUB
{
	open_run_file("snemo_run-1307_event-2589_github_version.root");
	vector<double> corrVector = readAndCreateCorrectionVector("calo-t0.txt");
	TFile* Classification_events = new TFile("Classification_events_1307_github_version.root", "RECREATE");
	TTree classification("classification", "");

	vector<int> true_calo_test = {};

	classification.Branch("Nb.of.electrons", &nb_total_elec);
	classification.Branch("Nb.of.alpha", &nb_total_alpha);
	classification.Branch("Nb.of.weirdies", &nb_total_weirdies);
	classification.Branch("Calo.gamma", &index_un_calo);
	classification.Branch("dt.correlation", &dt_good_corr);
	classification.Branch("calo.correlation", &calo_corr_ij);
	classification.Branch("candidates.2e", &candidate_2e);
	classification.Branch("pure.event", &pure_event);

	int entries = run_tree->GetEntries();

	for (int entry = 0; entry < entries; entry++)
	{
		run_tree->GetEntry(entry);

		int asso_weirdies = 0;
		int nb_total_caloasso = 0;
		vector<int> indice_calo = {};
		vector<int> index_as_calo = {};
		index_un_calo->clear();
		nb_total_alpha = 0;
		nb_total_elec = 0;
		nb_total_weirdies = 0;
		pure_event = false;
		candidate_2e = false;
		dt_good_corr->clear();
		calo_corr_ij->clear();
		calo_tablo->clear();
		true_calo_test = {};

		event_clusterisation();

		const int nb_calo_hits = digicalo_side->size();
		const int nb_tracker_hits = digitracker_side->size();
		const int nb_clusters = digiclusters->size();

		for (int cluster_i = 0; cluster_i < nb_clusters; cluster_i++)
		{
			const digicluster& digicluster_i = digiclusters->at(cluster_i);

			if (digicluster_i.tracker_indexes.size() == 1)
				continue;

			int nb_caloasso = digicluster_i.associated_calo_indexes.size();

			nb_total_caloasso = nb_total_caloasso + nb_caloasso;

			if (nb_caloasso == 0)
			{
				if (digicluster_i.tracker_indexes.size() < 7) nb_total_alpha++;
				else nb_total_weirdies++;
			}

			for (int i = 0; i < nb_caloasso; i++)
			{
				int true_calo = digicluster_i.associated_calo_indexes.at(i);
				double calo_amplitude_i = -(2500. / 4096.) * (1. / 8.) * digicalo_peakamplitude->at(true_calo);
				if (calo_amplitude_i < 50)
					continue;

				//// GAMMA COUNT 1st part ////
				index_as_calo.push_back(true_calo);


				true_calo_test.push_back(true_calo);
				if (digicluster_i.tracker_indexes.size() > 7 && digicluster_i.tracker_layer_min < 2 && digicluster_i.tracker_layer_max > 7)
				{
					nb_total_elec++;
				}
				else
				{
					nb_total_weirdies++;
					asso_weirdies++;
				}

				if (digicluster_i.tracker_indexes.size() < 7 || digicluster_i.tracker_layer_min >= 2 || digicluster_i.tracker_layer_max <= 7 || nb_caloasso == 0)
					continue;
				for (int cluster_j = cluster_i + 1; cluster_j < nb_clusters; cluster_j++)
				{
					const digicluster& digicluster_j = digiclusters->at(cluster_j);
					int nb_caloasso_j = digicluster_j.associated_calo_indexes.size();
					if (digicluster_j.tracker_indexes.size() <= 7 || digicluster_j.tracker_layer_min >= 2 || digicluster_j.tracker_layer_max <= 7 || nb_caloasso_j == 0)
						continue;

					for (int j = 0; j < nb_caloasso_j; j++)
					{
						int true_calo_j = digicluster_j.associated_calo_indexes.at(j);
						double calo_amplitude_j = -(2500. / 4096.) * (1. / 8.) * digicalo_peakamplitude->at(true_calo_j);
						if (calo_amplitude_j < 50)
							continue;

						const double calo_time = digicalo_timestamp->at(true_calo) * 6.25 - 250 + digicalo_falling_cell->at(true_calo) * (1. / 256.) * (400 / 1024.) - corrVector[digicalo_omnum(true_calo)];
						const double calo_time_j = digicalo_timestamp->at(true_calo_j) * 6.25 - 250 + digicalo_falling_cell->at(true_calo_j) * (1. / 256.) * (400 / 1024.) - corrVector[digicalo_omnum(true_calo_j)];

						if (fabs(calo_time_j - calo_time) < 2)
						{
							calo_corr_ij->push_back(true_calo);
							calo_corr_ij->push_back(true_calo_j);
							dt_good_corr->push_back(calo_time_j - calo_time); // calo time correlation
						}
					}
				}

			} // end caloasso with criteria amp
		} // end cluster i

		//////////////////////////////////////// SET VALID CALO COUNT //////////////////////////////////////// i.e. GAMMA COUNT 2nd part

		for (int k = 0; k < nb_calo_hits; k++)
		{
			double calo_amplitude_i = -(2500. / 4096.) * (1. / 8.) * digicalo_peakamplitude->at(k);
			if (calo_amplitude_i < 50)
				continue;
			int asso_unasso = 0;

			calo_tablo->push_back(digicalo_omnum(k)); // needed to draw !!

			for (int l : index_as_calo)
			{
				if (l == k)
				{
					asso_unasso = 1;
				}
			}
			if (asso_unasso == 0)
				index_un_calo->push_back(k);
		}

		///////////////////// Classification

		if (nb_total_elec != 0 && nb_total_alpha == 0 && nb_total_weirdies == 0)
		{
			pure_event = true;
		}

		classification.Fill();
	}
	Classification_events->cd();
	classification.Write();
	Classification_events->Close();
}

void class_3_analysis() //	NEED GITHUB VERSION
{
	open_run_file("snemo_run-1307_event-2589_github_version.root");
	open_class_file("Classification_events_1307_github_version.root");

	TH1F* histo_corr = new TH1F("hist_corr", "Histo corr", 100, -2, 2);
	vector<double> corrVector = readAndCreateCorrectionVector("calo-t0.txt");
	vector<double> enerVector = readAndCreateEnergyVector("run-1307_digicalo-charge-to-energy.txt");

	vector<int> good_entries = {};
	const int entries = class_tree->GetEntries();

	/////////////////////
	// ENERGY ANALYSIS //
	/////////////////////

	run_tree->GetEntry(0);
	class_tree->GetEntry(0);

	int nb_calo_hits = digicalo_charge->size();

	for (int i = 0; i < nb_calo_hits; i++)
	{
		double calo_amplitude_i = -(2500. / 4096.) * (1. / 8.) * digicalo_peakamplitude->at(i);
		if (calo_amplitude_i < 50)
			continue;
		cout << "Energy (MeV) : " << digicalo_charge->at(i) * enerVector[digicalo_omnum(i)] << " at calo " << digicalo_omnum(i) << endl;
	}

	////////////////////
	// FIRST CRITERIA //
	////////////////////

	event_clusterisation();
	if (nb_total_elec == 2 && pure_event == true && dt_good_corr->size() != 0)
		//if (digitracker_side->size() >= 50 && digiclusters->size() == 2)
	{
		good_entries.push_back(0);
		/*for (int j = 0; j < dt_good_corr->size(); j++)
		{
			histo_corr->Fill(dt_good_corr->at(j));
		}*/
	}

	/////////////////////
	// SECOND CRITERIA //
	/////////////////////

	//for (int entry = 0; entry < 10000; entry++)
	//{
	//	class_tree->GetEntry(entry);
	//	if (nb_total_elec == 2 && pure_event == true && dt_good_corr->size() != 0)
	//	{
	//		//good_entries.push_back(entry);
	//		cout << entry << " : " << dt_good_corr->at(0) << endl;
	//	}
	//}

	for (int i = 0; i < good_entries.size(); i++)
	{
		cout << "Entry " << good_entries.at(i) << " : ";
		class_tree->GetEntry(good_entries.at(i));
		run_tree->GetEntry(good_entries.at(i));
		cout << "dt = " << dt_good_corr->at(0) << " between calo " << calo_corr_ij->at(0) << " and " << calo_corr_ij->at(1) << endl;

		for (int j = 0; j < digicalo_side->size() - 1; j++)
		{
			double calo_amplitude_j = -(2500. / 4096.) * (1. / 8.) * digicalo_peakamplitude->at(j);
			if (calo_amplitude_j < 50)
				continue;

			const double calo_time_j = (digicalo_timestamp->at(j) * 6.25 - 250 + digicalo_falling_cell->at(j) * (1. / 256.) * (400 / 1024.) - corrVector[digicalo_omnum(j)]);

			for (int k = j + 1; k < digicalo_side->size(); k++)
			{
				double calo_amplitude_k = -(2500. / 4096.) * (1. / 8.) * digicalo_peakamplitude->at(k);
				if (calo_amplitude_k < 50)
					continue;

				const double calo_time_k = (digicalo_timestamp->at(k) * 6.25 - 250 + digicalo_falling_cell->at(k) * (1. / 256.) * (400 / 1024.) - corrVector[digicalo_omnum(k)]);
				cout << "Absolute time : " << digicalo_omnum(j) << " - " << digicalo_omnum(k) << " has a dt = " << calo_time_j - calo_time_k << endl;
			}
		}
	}

	//TCanvas* C = new TCanvas();
	//C->SetLogy();
	//histo_corr->SetLineColor(kBlue);
	//histo_corr->Draw();

	/*auto legend = new TLegend(0.1, 0.7, 0.48, 0.9);
	legend->AddEntry(histo_corr, "Histogram", "l");
	legend->Draw();

	gPad->BuildLegend();*/

} // Quels seraient les critères same-sided ? 1e avec 2 caloasso de delta column ? Eventuellement 1 alpha ou weirdy ? !!!!!!!!!!!!!!! récupérer les delta_t < 2ns same-sided !!!!!!!!!!!!!!!

void show_clusterised_event() // NEED GITHUB VERSION
{
	open_class_file("Classification_events_1307_github_version.root");
	open_run_file("snemo_run-1307_event-2589_github_version.root");
	event_display->reset();

	for (TBox* box : cluster_boxes)
		delete box;
	cluster_boxes.clear();

	vector<int> true_calo_test = {};
	double first_hit_time = numeric_limits<double>::max();
	double last_hit_time = numeric_limits<double>::min();

	run_tree->GetEntry(0);
	class_tree->GetEntry(0);
	event_clusterisation();

	const int nb_calo_hits = digicalo_side->size();
	const int nb_tracker_hits = digitracker_side->size();
	const int nb_clusters = digiclusters->size();
	calo_tablo->clear();

	for (int calo_i = 0; calo_i < nb_calo_hits; calo_i++)
	{
		double calo_amplitude_i = -(2500. / 4096.) * (1. / 8.) * digicalo_peakamplitude->at(calo_i);
		if (calo_amplitude_i < 50)
			continue;

		calo_tablo->push_back(digicalo_omnum(calo_i)); // needed to draw all calo

		const int omnum = digicalo_omnum(calo_i);
		const double calo_hit_time = digicalo_timestamp->at(calo_i) * 6.25E-9 - 250E-9 + digicalo_falling_cell->at(calo_i) * (1. / 256.) * (400.E-9 / 1024.);
		event_display->setomcontent(omnum, calo_hit_time); // we use the time as content

		if (calo_hit_time < first_hit_time)
			first_hit_time = calo_hit_time;

		if (calo_hit_time > last_hit_time)
			last_hit_time = calo_hit_time;
	}

	for (int tracker_i = 0; tracker_i < nb_tracker_hits; tracker_i++)
	{
		const int cellnum = digitracker_cellnum(tracker_i);
		const double anode_hit_time = digitracker_anodetimestampR0->at(tracker_i) * 12.5E-9;

		// the tracker hit must have an anode R0 time
		if (anode_hit_time <= 0)
			continue;

		event_display->setggcontent(cellnum, anode_hit_time); // we use the time as content

		if (anode_hit_time < first_hit_time)
			first_hit_time = anode_hit_time;

		if (anode_hit_time > last_hit_time)
			last_hit_time = anode_hit_time;
	}

	event_display->setrange(first_hit_time, last_hit_time);
	event_display->settitle(Form("Event %d", 2589));
	event_display->draw_top();

	// add cluster's boxes
	for (const digicluster& cluster : *digiclusters)
	{
		double cluster_box_x[2] = { 1, 0 };
		double cluster_box_y[2] = { 1, 0 };
		float mean_cluster_color = 0;

		for (const int& tracker_i : cluster.tracker_indexes)
		{
			if (cluster.tracker_indexes.size() == 1) continue;
			const int cellnum = digitracker_cellnum(tracker_i);

			// retrieve dimensions of the cell (TEllipse object) in the event_display
			const double cx = event_display->top_gg_ellipse[cellnum]->GetX1();
			const double cy = event_display->top_gg_ellipse[cellnum]->GetY1();
			const double r1 = event_display->top_gg_ellipse[cellnum]->GetR1();
			const double r2 = event_display->top_gg_ellipse[cellnum]->GetR2();
			const double x1 = cx - r1; const double y1 = cy - r2;
			const double x2 = cx + r1; const double y2 = cy + r2;

			if (x1 < cluster_box_x[0]) cluster_box_x[0] = x1;
			if (x2 > cluster_box_x[1]) cluster_box_x[1] = x2;
			if (y1 < cluster_box_y[0]) cluster_box_y[0] = y1;
			if (y2 > cluster_box_y[1]) cluster_box_y[1] = y2;

			// retrieve the color of the cell
			mean_cluster_color += event_display->getggcolorindex(cellnum);

		} // for (tracker_index)

		TBox* cluster_box = new TBox(cluster_box_x[0], cluster_box_y[0], cluster_box_x[1], cluster_box_y[1]);
		cluster_box->SetFillStyle(0);
		cluster_box->SetLineWidth(5);

		mean_cluster_color /= cluster.tracker_indexes.size();
		cluster_box->SetLineColor((int)(0.5 + mean_cluster_color));
		cluster_box->Draw();

		cluster_boxes.push_back(cluster_box);

	} // for (cluster)

	sncalo = new sndisplay::calorimeter("sndisplay_dis");
	for (int i : *calo_tablo) {
		sncalo->setcontent(i, 1);
	}
	sncalo->draw1();
	sncalo->draw1();
	cout << "Number of : e - " << nb_total_elec << " - g - " << index_un_calo->size() << " - alpha - " << nb_total_alpha << " - weirdies - " << nb_total_weirdies << "." << endl;
	cout << "Test : " << dt_good_corr->at(0) << endl;
}