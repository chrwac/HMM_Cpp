/*
 *  HMM_oo.h
 *  HMM_object_oriented
 *
 *  Created by Christian on 8/19/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef HMM_OO_H
#define HMM_OO_H

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <deque>

#include "HMM_binner.h"
#include "EmissionProb.h"

using namespace std;

const int HMM_DATA_INT = 1;
const int HMM_DATA_FP = 2;
const double MIN_EMISSION_PROBABILITY=1e-6; // 1e-6 as minimal emission probability...

template <typename T>
class CHMMData
{
private:
	vector<T> m_data;
	vector<short int> m_classifications;
	int m_num_data_points;
public:
	CHMMData();
	void ReadHMMDataFromTxt(string filename,const int data_type);
	void WriteHMMDataToTxt(string filename);			// writes the (classified!!) HMM-Data to a file..
	vector<T>& ReturnRefToData() {return m_data;}
	vector<short int>& RetRefToClassifications(){return m_classifications;}
	int ReturnNumberOfDataPoints(){return m_num_data_points;}
};

template <typename T>
class CHMM
{
private:
	CHMMData<T>* m_p_to_data;
	CBinner<T> m_binner;
	CEmissionPObject<T> m_emission_p_object; 
	
	unsigned short int m_num_bins; // set binning of emission-values
	int m_num_states; // number of states with which to model the data
	int m_order_of_markov_chain; // oder of Markov-Chain, for now 1, might be extended at some point !!!
	vector<double> m_initial_state_probs;
	vector<vector<double> > m_state_transition_probs;
	vector<unsigned short int> m_binned_values;
	vector<vector<double> > m_state_emission_probabilities;
	
	vector<vector<double> > m_alphas; // return values of the forward algorithm
	vector<vector<double> > m_betas; // return values of the backward algorithm
	vector<vector<double> > m_gammas;
	vector<vector<vector<double> > > m_zetas;

	vector<double> m_scaling_factors;	// scaling factors obtained from the Forward algorithm...
	//vector<unsigned short int> m_viterbi_states;
	unsigned short int* m_viterbi_states;
	vector<unsigned short int> m_max_states;

	void ForwardAlgorithm();
	void BackwardAlgorithm();
	void CalculateGammas();
	void CalculateZetas();
	void UpdateParameters();
	void BackTrack(double** ppVs,unsigned short int** ppPsis);
public:
	//CHMM();
	~CHMM();
	void SetPointerToHMMData(CHMMData<T>* p_to_data){ m_p_to_data=p_to_data;}
	void SetNumberOfStates(int nstates){ m_num_states=nstates;}
	void SetOrderOfMarkovChain(int order){ m_order_of_markov_chain=order;}
	void SetNumberOfBins(int nbins){m_num_bins = nbins;}
	void SetInitialProbs(int nstates,vector<T>& initial_probs);
	void SetTransitionProbs(int nstates,vector<vector<T> >& transition_probs);
	void PrintOutCurrentValues();	// prints out the current values (number of states, initial state probabilities, transition probabilities) to the screen....
	void Iterate();
	void DoBinningOfObservations();	// do the binning of the values...
	void InitializeEmissionProbabilities(const int type, vector<vector<T> >& metadata);
	void ViterbiAssignment();
	void MaxAssignment();
	void WriteAssignmentToFile(string filename);
	void WriteViterbiAssignmentToFile(string filename);
};

template<typename T>
void CHMM<T>::WriteViterbiAssignmentToFile(string filename)
{
	int num_points = m_p_to_data->ReturnNumberOfDataPoints();

	FILE* fp_output = fopen(filename.c_str(),"wt");
	for(int i=0;i<num_points;++i)
	{
		fprintf(fp_output,"%d\n",this->m_viterbi_states[i]);
	}
	fclose(fp_output);
}

template <typename T>
CHMM<T>::~CHMM()
{
	if(m_viterbi_states!=NULL)
		delete[] m_viterbi_states;
}

template <typename T>
void CHMM<T>::WriteAssignmentToFile(string filename)
{
	// oldschool C, change it at some point..
	int num_points = m_p_to_data->ReturnNumberOfDataPoints();
	FILE* fp_outputfile = fopen(filename.c_str(),"wt");
	for(int i=0;i<num_points;++i)
		fprintf(fp_outputfile,"%d\n",m_max_states[i]);
	fclose(fp_outputfile);
}

template <typename T>
void CHMM<T>::BackTrack(double** ppVs,unsigned short int** ppPsis)
{
	// fill m_viterbi_states with the most
	// probable single state sequence (as obtained by Viterbi-algorithm)
	int num_points = m_p_to_data->ReturnNumberOfDataPoints();

	m_viterbi_states = new unsigned short int[num_points];
	int state_max_index=0;
	double max_V = ppVs[num_points-1][0];
	for(int i=1;i<m_num_states;++i)
	{
		if(ppVs[num_points-1][i]>max_V)
		{
			max_V = ppVs[num_points-1][i];
			state_max_index=i;
		}
	}
	m_viterbi_states[num_points-1]=state_max_index;
	for(int t=num_points-2;t>=0;t--)
		m_viterbi_states[t] = ppPsis[t+1][m_viterbi_states[t+1]];
}

template <typename T>
void CHMM<T>::ViterbiAssignment()
{
	// calculate the log-probabilities
	int num_points = m_p_to_data->ReturnNumberOfDataPoints();

	double* log_state_initial_probs = new double[m_num_states];
	double** log_state_trans_probs = new double*[m_num_states];
	double** log_state_em_probs = new double*[m_num_states];
	double** Vs = new double*[num_points];
	unsigned short int** psis = new unsigned short int*[num_points];

	for(int i=0;i<num_points;++i)
	{
		psis[i] = new unsigned short int[m_num_states];
		Vs[i] = new double[m_num_states];
	}
//	int num_points = m_p_to_data->ReturnNumberOfDataPoints();
	for(int i=0;i<m_num_states;++i)
	{
		log_state_trans_probs[i] = new double[m_num_states];
		log_state_em_probs[i] = new double[m_num_bins];
		
		log_state_initial_probs[i] = log10(m_initial_state_probs[i]);
		for(int j=0;j<m_num_states;j++)
			log_state_trans_probs[i][j] = log10(m_state_transition_probs[i][j]);
		for(int j=0;j<m_num_bins;j++)
			log_state_em_probs[i][j] = log10(m_state_emission_probabilities[i][j]);
	}

	// initialize values for Viterbi-algorithm:
	for(int i=0;i<m_num_states;++i)
	{
		Vs[0][i] = log_state_initial_probs[i] + log_state_em_probs[i][m_binned_values[0]];
		psis[0][i] = 0;
	}
		//

	double log_prob_max,log_prob;
	int state_max_index;
	for(int i=1;i<num_points;++i)
		for(int j=0;j<m_num_states;++j)
		{
			state_max_index=0;
			log_prob_max = Vs[i-1][0]+log_state_trans_probs[0][j]+log_state_em_probs[j][m_binned_values[i]];
			for(int k=1;k<m_num_states;++k)
			{
				log_prob = Vs[i-1][k]+log_state_trans_probs[k][j]+log_state_em_probs[j][m_binned_values[i]];
				if(log_prob>log_prob_max)
				{
					state_max_index=k;
					log_prob_max=log_prob;
				}
			}
			psis[i][j]=state_max_index;
			Vs[i][j]=log_prob_max;
		}

	// compute the actual Viterbi-sequence via
	// back-tracking...
	BackTrack(Vs,psis);

	for(int i=0;i<m_num_states;++i)
	{
		delete[] log_state_em_probs[i];
		delete[] log_state_trans_probs[i];
	}	
	
	delete[] log_state_em_probs;
	delete[] log_state_trans_probs;
	delete[] log_state_initial_probs;

	for(int i=0;i<num_points;++i)
	{
		delete[] Vs[i];
		delete[] psis[i];
	}
	delete[] Vs;
	delete[] psis;
}

template <typename T>
void CHMM<T>::MaxAssignment()
{
	m_max_states.clear();
	int num_points = m_p_to_data->ReturnNumberOfDataPoints();
	for(int i=0;i<num_points;++i)
	{
		double max_value=m_gammas[i][0];
		unsigned short int max_index = 0;
		for(int j=1;j<m_num_states;++j)
		{
			if(m_gammas[i][j]>max_value)
			{
				max_value=m_gammas[i][j];
				max_index=j;	
			}
		}
		m_max_states.push_back(max_index);
	}
}


template <typename T>
void CHMM<T>::UpdateParameters()
{
	
	int num_points = m_p_to_data->ReturnNumberOfDataPoints();
	double* temp_denominators = new double[m_num_states];
	for(int i=0;i<m_num_states;++i)
	{
		// for each state, initial state probability is given
		// by gamma_i_0
		m_initial_state_probs[i]=m_gammas[0][i];
		// set state-transition probabilities
		// and state emission probabilities to zero...
		for(int j=0;j<m_num_states;++j)
			m_state_transition_probs[i][j]=0.0f;
		for(int j=0;j<m_num_bins;++j)
			m_state_emission_probabilities[i][j]=0.0f;
		
		temp_denominators[i]=0.0f;

		// update the state transition probabilities
		// using the precalculated zetas
		for(int t=0;t<(num_points-1);t++)
		{
			for(int j=0;j<m_num_states;++j)
			{
				m_state_transition_probs[i][j]+=m_zetas[t][i][j];
				temp_denominators[i]+=m_zetas[t][i][j];
			}
		}
		// normalize the state transition probabilities
		for(int j=0;j<m_num_states;++j)
			m_state_transition_probs[i][j]/=temp_denominators[i];
	}
	delete[] temp_denominators;

	double* sums = new double[m_num_states];
	for(int i=0;i<m_num_states;++i)
	{
		sums[i]=0.0f;
		for(int t=0;t<num_points;++t)
		{
			sums[i]+=m_gammas[t][i];
			int curr_index = m_binned_values[t];
			m_state_emission_probabilities[i][curr_index]+=m_gammas[t][i];
		}
	}
	
	// first round of normalization
	for(int i=0;i<m_num_states;++i)
		for(int j=0;j<m_num_bins;++j)
			m_state_emission_probabilities[i][j]/=sums[i];
	
	// second round of normalization
	for(int i=0;i<m_num_states;++i)
	{
		double sum_ems=0.0f;
		for(int j=0;j<m_num_states;++j)
		{
			if(m_state_emission_probabilities[i][j]<MIN_EMISSION_PROBABILITY)
				m_state_emission_probabilities[i][j]=MIN_EMISSION_PROBABILITY;
			sum_ems+=m_state_emission_probabilities[i][j];
		}
		for(int j=0;j<m_num_states;++j)
			m_state_emission_probabilities[i][j]/=sum_ems;
	}
}

template <typename T>
void CHMM<T>::CalculateZetas()
{
	int num_points = m_p_to_data->ReturnNumberOfDataPoints();
	vector<vector<double> > temp_matrix;
	vector<double> temp_vector;

	for(int i=0;i<m_num_states;++i)
		temp_vector.push_back(0.0f);
	
	for(int j=0;j<m_num_states;++j)
		temp_matrix.push_back(temp_vector);

	
	for(int t=0;t<num_points-1;++t)
	{
		m_zetas.push_back(temp_matrix);	
		for(int i=0;i<m_num_states;++i)
			for(int j=0;j<m_num_states;++j)
		m_zetas[t][i][j] = m_gammas[t][i]*m_state_transition_probs[i][j]*m_state_emission_probabilities[j][m_binned_values[t+1]]*m_scaling_factors[t+1]*m_betas[t+1][j]/(m_betas[t][i]);
	}
}

template <typename T>
void CHMM<T>::Iterate()
{
	cout << "ForwardAlgorithm" << endl;
	ForwardAlgorithm();
	cout << "BackwardAlgorithm" << endl;
	BackwardAlgorithm();
	cout << "CalculateGammas" << endl;
	CalculateGammas();
	cout << "CalculateZetas" << endl;
	CalculateZetas();
	cout << "UpdateParameters" << endl;
	UpdateParameters();
	PrintOutCurrentValues();

}
template <typename T>
void CHMM<T>::CalculateGammas()
{
	m_gammas.clear();
	vector<double> temp_gammas;
	
	double denominator;
	
	int num_points=m_p_to_data->ReturnNumberOfDataPoints();
	for(int index_point=0;index_point<num_points;++index_point)
	{
		denominator=0.0f;
		temp_gammas.clear();
		for(int index_state=0;index_state<m_num_states;index_state++)
			denominator+=m_alphas[index_point][index_state]*m_betas[index_point][index_state];
		for(int index_state=0;index_state<m_num_states;index_state++)
		{
			temp_gammas.push_back(0.0f);
			temp_gammas[index_state] = m_alphas[index_point][index_state]*m_betas[index_point][index_state]/denominator;
		}
		m_gammas.push_back(temp_gammas);
	}
}

template <typename T>
void CHMM<T>::ForwardAlgorithm() 
{
	// clear alphas and scaling factors, which will be calculated again from scratch...
	m_alphas.clear(); 
	m_scaling_factors.clear();
	
	int num_points = m_p_to_data->ReturnNumberOfDataPoints();
	m_alphas.reserve(num_points);
	
	vector<double> alphas_temp;
	m_scaling_factors.push_back(0.0f);
	
	for(int index_state=0;index_state<m_num_states;index_state++)
	{
		alphas_temp.push_back(0.0f);
		alphas_temp[index_state]=m_initial_state_probs[index_state]*m_state_emission_probabilities[index_state][m_binned_values[0]];
		m_scaling_factors[0]+=alphas_temp[index_state];
	}
	m_scaling_factors[0]=1.0f/m_scaling_factors[0];
	// now scale the temporary alphas and push them back...
	for(int index_state=0;index_state<m_num_states;index_state++)
		alphas_temp[index_state]*=m_scaling_factors[0];
	m_alphas.push_back(alphas_temp);
	
	
	// calculate the rest of the alphas and scaling factors..
	for(int index_point=1;index_point<num_points;index_point++)
	{
		alphas_temp.clear();
		m_scaling_factors.push_back(0.0f);
		for(int index_state=0;index_state<m_num_states;index_state++)
		{
			alphas_temp.push_back(0.0f);
			for(int k=0;k<m_num_states;k++)
				alphas_temp[index_state]+=m_alphas[index_point-1][k]*m_state_transition_probs[k][index_state]*m_state_emission_probabilities[index_state][m_binned_values[index_point]];
		
			m_scaling_factors[index_point]+=alphas_temp[index_state];
		}
		m_scaling_factors[index_point] = 1.0f/m_scaling_factors[index_point];
		for(int index_state=0;index_state<m_num_states;index_state++)
			alphas_temp[index_state]*=m_scaling_factors[index_point];
		m_alphas.push_back(alphas_temp);
	}
}

template <typename T>
void CHMM<T>::BackwardAlgorithm()
{
	
	m_betas.clear();
	int num_points = m_p_to_data->ReturnNumberOfDataPoints();
	vector<double> temp_vector;
	for(int i=0;i<m_num_states;++i)
		temp_vector.push_back(0.0f);
	
	for(int i=0;i<num_points;i++)
	{
		m_betas.push_back(temp_vector);
	}
	
	for(int index_state=0;index_state<m_num_states;index_state++)
	{
		m_betas[num_points-1][index_state]=1.0f;
	}
	for(int index_point=num_points-2;index_point>=0;index_point--)
		for(int index_state=0;index_state<m_num_states;index_state++)
		{
			m_betas[index_point][index_state]=0.0f;
			for(int k=0;k<m_num_states;k++)
				m_betas[index_point][index_state]+=m_betas[index_point+1][k]*m_state_transition_probs[index_state][k]*m_state_emission_probabilities[k][m_binned_values[index_point+1]];
			m_betas[index_point][index_state]*=m_scaling_factors[index_point+1];																																	 
		}
}



template <typename T>
void CHMM<T>::InitializeEmissionProbabilities(const int type, vector<vector<T> >& metadata)
{
	m_emission_p_object.SetNumberOfStates(m_num_states);
	m_emission_p_object.SetType(type);
	if(type==EMISSION_TYPE_GAUSSIAN)
	{
		cout << "Gaussian Emission Probabilities: " << endl;
		m_emission_p_object.SetNumberOfMetadataValues(2);
		m_emission_p_object.SetMetaData(metadata);
		
		// now calculate for each value of the binned values
		// 1. The unbinned value (using m_binner)
		// 2. the emission value for each state (using m_emission_p_object)
		for(int index_state=0;index_state<m_num_states;index_state++)
		{
			vector<double> temp_vec;
			temp_vec.clear();
			double sum=0.0f;
			for(int index_bin=0;index_bin<m_num_bins;index_bin++)
			{
				// get value from index (using m_binner)
				T value = m_binner.ReturnValueFromIndex((unsigned short int)index_bin);
				
				// get emission probability (UNNORMALIZED!!!!)
				// (using m_emission_p_object)
				double prob = m_emission_p_object.ReturnEmissionProbability(index_state,value);
				sum+=prob;
				temp_vec.push_back(prob);
			}
			double sum2=0.0f;
			// normalize and account for minimal emission probability
			for(int index_bin=0;index_bin<m_num_bins;index_bin++)
			{
				temp_vec[index_bin]/=sum;
				// NOTE: A minimal emission probability is applied
				// No emisssion is supposed to be smaller than this value....
				if(temp_vec[index_bin]<MIN_EMISSION_PROBABILITY)
					temp_vec[index_bin]=MIN_EMISSION_PROBABILITY;
				
				sum2+=temp_vec[index_bin];
			}
			// normalize again and push back
			for(int index_bin=0;index_bin<m_num_bins;index_bin++)
				temp_vec[index_bin]/=sum2;
			
			m_state_emission_probabilities.push_back(temp_vec);
		}
	}
}

template <typename T>
void CHMM<T>::DoBinningOfObservations()
{
	// create "binning-object:"
	m_binner.SetNumBins(m_num_bins);
	m_binner.DoBinning(m_p_to_data->ReturnRefToData(),m_binned_values);
	m_binner.PrintOutBinnerParameters();
}

template <typename T>
void CHMM<T>::SetInitialProbs(int nstates,vector<T>& initial_probs)
{
	// clear all existing elements
	m_initial_state_probs.clear();
	for(int i=0;i<nstates;++i)
		m_initial_state_probs.push_back(initial_probs[i]);
}

template <typename T>
void CHMM<T>::SetTransitionProbs(int nstates,vector<vector<T> >& transition_probs)
{
	m_state_transition_probs.clear();
	vector<T> row;
	
	for(int i=0;i<nstates;++i)
	{
		row.clear();
		for(int j=0;j<nstates;++j)
			row.push_back(transition_probs[i][j]);
		
		m_state_transition_probs.push_back(row);
	}
}

template <typename T>
void CHMM<T>::PrintOutCurrentValues()
{
	// print out basic information, like
	// number of states, order of Markov Chain,
	// the initial state probabilities and
	// the state transition probabilities....
	cout << "number of states: " << m_num_states << endl;
	cout << "order of Markov-Chain: " << m_order_of_markov_chain << endl;
	
	cout << "initial state probabilities: " << endl;
	for(int i=0;i<m_num_states;++i)
		cout << "state" << i << ": " << m_initial_state_probs[i] << endl;

	cout << "state transition probabilities: " << endl;
	for(int i=0;i<m_num_states;++i)
	{
		for(int j=0;j<m_num_states;++j)
			cout << m_state_transition_probs[i][j] << "\t";
		
		cout << endl;
	}	 
	/*
	cout << "emission values: " << endl;
	for(int index_state=0;index_state<m_num_states;index_state++)
	{
		cout << "state: " << index_state << endl;
		for(int index_bin=0;index_bin<m_num_bins;++index_bin)
		{
			cout << "emission at: " << m_binner.ReturnValueFromIndex((unsigned short int)index_bin) << " " << m_state_emission_probabilities[index_state][index_bin] << endl;
		}
	}
	*/
}

template <typename T> 
CHMMData<T>::CHMMData()
{
	m_num_data_points=0;
}

template <typename T> 
void CHMMData<T>::ReadHMMDataFromTxt(string filename,const int data_type)
{
	string line;
	ifstream data_file(filename.c_str(),ios::in);
	int num_data_points_read=0;
	char* pEnd;
	T temp_value;
	if(data_file.is_open())
	{
		while (getline(data_file,line))
		{
		//	cout << line << endl;
			switch(data_type)
			{
				case HMM_DATA_INT:
					temp_value = strtol(line.c_str(),&pEnd,10);
					break;
				case HMM_DATA_FP: 
					temp_value = strtod(line.c_str(),&pEnd);
					break;
			}
			m_data.push_back(temp_value);
			num_data_points_read++;
		}
	}
	else
		cout << "unable to open file --> check filename" << endl;
	
	m_num_data_points=num_data_points_read;
	cout << num_data_points_read << " data points read from file: " << filename << endl;
	data_file.close();
}

template <typename T> 
void CHMMData<T>::WriteHMMDataToTxt(string filename)
{
	ofstream data_file(filename.c_str(),ios::out);
	if(data_file.is_open())
	{
		for(int i=0;i<m_num_data_points;++i)
		{
			data_file << m_data[i] << "\t" << m_classifications[i] << endl;
		}
	}
	else 
		cout << "unable to open file --> check filename" << endl;
}

#endif