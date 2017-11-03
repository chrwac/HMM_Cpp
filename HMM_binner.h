/*
 *  HMM_binner.h
 *  HMM_object_oriented
 *
 *  Created by Christian on 8/19/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef HMM_BINNER_H
#define HMM_BINNER_H

#include <vector>

#define HMM_BINNER_TESTMODE
using namespace std;

template <typename T>
class CBinner
{
private:
	unsigned short int m_num_bins;
	T m_min_value;
	T m_max_value;
	T m_increment; 
public:
	CBinner();
	CBinner(int num_bins); 
	void SetNumBins(int num_bins);
	void DoBinning(vector<T>& original_values,vector<unsigned short int>& binned_values);
	void PrintOutBinnerParameters();
	T ReturnValueFromIndex(unsigned short int bin_index);
	unsigned short int ReturnIndexFromValue(T value);
};

template <typename T>
CBinner<T>::CBinner()
{m_num_bins=100;}

template <typename T>
CBinner<T>::CBinner(int num_bins)
{m_num_bins=num_bins;}

template <typename T>
void CBinner<T>::SetNumBins(int num_bins)
{m_num_bins=num_bins;}

template <typename T>
void CBinner<T>::DoBinning(vector<T>& original_values,vector<unsigned short int>& binned_values)
{
	binned_values.clear();
	// find minimum and maximum....
	m_min_value = *min_element(original_values.begin(), original_values.end());
	m_max_value = *max_element(original_values.begin(), original_values.end());
	
	// careful, the actual "maximum" is very slightly shifted.....
	m_max_value+=0.0001*(m_max_value-m_min_value);
	
	m_increment=(m_max_value-m_min_value)/(m_num_bins);

	int curr_index=0;
	for(int i=0;i<original_values.size();++i)
	{
		curr_index=(original_values[i]-m_min_value)/m_increment;
		binned_values.push_back(curr_index);
	}
}

template <typename T>
void CBinner<T>::PrintOutBinnerParameters()
{
	cout << "number of bins: " << m_num_bins << endl;
	cout << "minimum: " << m_min_value << endl;
	cout << "maximum: " << m_max_value << endl;
	cout << "increment: " << m_increment << endl;
}

template <typename T>
T CBinner<T>::ReturnValueFromIndex(unsigned short int bin_index)
{
	T return_value = m_min_value + bin_index*m_increment+(0.5f*m_increment);
	return return_value;
}

template <typename T>
unsigned short int CBinner<T>::ReturnIndexFromValue(T value)
{
	unsigned short int index_to_return = (value-m_min_value)/m_increment;
	return index_to_return;
}

#endif