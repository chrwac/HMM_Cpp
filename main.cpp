#include <iostream>
#include <string>
#include <vector>

#include "HMM_oo.h"
#include "HMM_binner.h"
using namespace std;

CHMMData<double> hmm_data;
CHMM<double> hmm;
CBinner<double> hmm_binner;

int main (int argc, char * const argv[]) 
{
	string filename("C:\\Studium\\Promotion\\C-Codes\\HMM_object_oriented\\HMM_object_oriented\\sequence_part.txt");
	hmm_data.ReadHMMDataFromTxt(filename,HMM_DATA_FP);
	
	int loc_num_states=3;
	int loc_num_bins=300;
	
	vector<double> loc_init_probs;
	loc_init_probs.push_back(0.3f);
	loc_init_probs.push_back(0.6f);
	loc_init_probs.push_back(0.1f);
	
	vector<vector<double> > loc_trans_probs;
	
	
	vector<double> row_probs;
	row_probs.push_back(0.1f);
	row_probs.push_back(0.8f);
	row_probs.push_back(0.1f);
	
	
	loc_trans_probs.push_back(row_probs);
	loc_trans_probs.push_back(row_probs);
	loc_trans_probs.push_back(row_probs);
	
	
	hmm.SetNumberOfStates(loc_num_states);
	hmm.SetOrderOfMarkovChain(1);
	hmm.SetNumberOfBins(loc_num_bins);
	hmm.SetPointerToHMMData(&hmm_data);
	
	hmm.SetInitialProbs(loc_num_states,loc_init_probs);
	hmm.SetTransitionProbs(loc_num_states, loc_trans_probs);
	
	hmm.DoBinningOfObservations();
	
	vector<vector<double> > gauss_parameters;
	vector<double> gauss_parameters_row;
	gauss_parameters_row.push_back(10.0f);
	gauss_parameters_row.push_back(1.0f);
	gauss_parameters.push_back(gauss_parameters_row);
	
	gauss_parameters_row.clear();
	gauss_parameters_row.push_back(15.0f);
	gauss_parameters_row.push_back(1.0f);
	gauss_parameters.push_back(gauss_parameters_row);
	
	gauss_parameters_row.clear();
	gauss_parameters_row.push_back(20.0f);
	gauss_parameters_row.push_back(1.0f);
	gauss_parameters.push_back(gauss_parameters_row);
	
	hmm.InitializeEmissionProbabilities(EMISSION_TYPE_GAUSSIAN,gauss_parameters);
	hmm.PrintOutCurrentValues();
	
	for(int i=0;i<20;i++)
		hmm.Iterate();
	
	string filename_output("C:\\Studium\\Promotion\\C-Codes\\HMM_object_oriented\\HMM_object_oriented\\sequence_part_assigned.txt");
	cout << "assigning max-state-per-time values: " << endl;
	hmm.MaxAssignment();
	cout << "assigning Viterbi states: " << endl;
	hmm.ViterbiAssignment();

	cout << "writing state assignment to file: " << endl;
	hmm.WriteAssignmentToFile(filename_output);

	cout << "writing Viterbi state assignment to file: " << endl;
	string filename_output2("C:\\Studium\\Promotion\\C-Codes\\HMM_object_oriented\\HMM_object_oriented\\sequence_part_assigned_Vit.txt");
	hmm.WriteViterbiAssignmentToFile(filename_output2);
	
	cout << "Eingabe zum Beenden: " << endl;
	int d;
	cin >> d;
    return 0;
}
