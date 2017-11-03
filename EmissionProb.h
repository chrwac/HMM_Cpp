/*
 *  EmissionProb.h
 *  HMM_object_oriented
 *
 *  Created by Christian on 8/19/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef HMM_EMISSION_PROB_OBJ_H
#define HMM_EMISSION_PROB_OBJ_H

#include <vector>
#include <math.h>
#include "EmissionProb.h"

const int EMISSION_TYPE_GAUSSIAN=1;

template<typename T>
T unnormalized_normal_pdf(T x,T m,T s)
{
	T a = (x-m)/s;
	return exp(-0.5f*a*a);
}

template <typename T> 
class CEmissionPObject
{
private:
	int m_num_states;
	int m_num_of_metadata_values; // the number of "metadata"-values (e.g. 2 for a Gaussian (mean and stdev))
	vector<vector<T> > m_metadata; // for example, means and stdevs of Gaussians...
	int m_type; // for example Gaussian emission probabilities....
public:
	CEmissionPObject();
	void SetNumberOfStates(int n_states){m_num_states=n_states;}
	void SetNumberOfMetadataValues(int n_values){m_num_of_metadata_values=n_values;}
	void SetType(const int type){m_type=type;}
	
	void SetMetaData(vector<vector<T> > &metadata);
	vector<vector<T> >& ReturnRefToMetaData(){return m_metadata;}
	int ReturnType(){return m_type;}
	int ReturnNumberOfStates(){return m_num_states;}
	int ReturnNumberOfMetadataValues(){return m_num_of_metadata_values;}
	double ReturnEmissionProbability(int state_index,T value);
};

template <typename T> 
CEmissionPObject<T>::CEmissionPObject()
{
	m_num_states=2;
	m_type=EMISSION_TYPE_GAUSSIAN;
	m_num_of_metadata_values=2;
	m_type=EMISSION_TYPE_GAUSSIAN;
}

template <typename T>
double CEmissionPObject<T>::ReturnEmissionProbability(int state_index,T value)
{
	if(m_type==EMISSION_TYPE_GAUSSIAN)
	{
		double res=0.0f;
		T mean = m_metadata[state_index][0];
		T stdev = m_metadata[state_index][1];
		return unnormalized_normal_pdf(value,mean,stdev);
	}
}

template <typename T>
void CEmissionPObject<T>::SetMetaData(vector<vector<T> > &metadata)
{
	m_metadata.clear();
	vector<T> md_row;
	for(int index_state=0;index_state<m_num_states;index_state++)
	{
		md_row.clear();
		for(int index_md_value=0;index_md_value<m_num_of_metadata_values;index_md_value++)
			md_row.push_back(metadata[index_state][index_md_value]);
		m_metadata.push_back(md_row);
	}
}

#endif


