/*
	Written by Jelle Mulyadi, 2021
*/

#include "PivotalSearch.h"
#include "../General/Tools.h"
#include <iostream>
#include <algorithm>
#include <functional>
#include <ctime>

//Algorithms and code based on:
// - Deng, Dong, Guoliang Li, and Jianhua Feng. "A pivotal prefix based filtering algorithm for string similarity search." Proceedings of the 2014 ACM SIGMOD international conference on Management of data. 2014.
// - Qin, Jianbin, and Chuan Xiao. "Pigeonring: A principle for faster thresholded similarity search." arXiv preprint arXiv:1804.01614 (2018).

PivotalSearch::PivotalSearch(Database* db, int q, int threshold, int chainLength)
{
	this->db = db;
	this->q = q;
	this->threshold = threshold;
	this->chainLength = chainLength;

	//Indexing
	std::clock_t start = std::clock();
	Indexing();
	double duration = (std::clock() - start) / (CLOCKS_PER_SEC / 1000);
	indexTime = (int)duration;
	std::cout << "Time elapsed indexing: " << duration << "\n";
}

//Sorting orders
bool PositionOrdering(Qgram i, Qgram j)
{
	return (i.pos < j.pos);
}

//Indexing
void PivotalSearch::CountFrequency()
{
	//For every qgram count frequency in the database
	for (int i = 0; i < db->db.size(); i++)
	{
		std::string entry = db->db[i].sequence;
		for (int j = 0; j < entry.length() - q + 1; j++)
			qgramFrequency[entry.substr(j, q)]++;
	}
}

void PivotalSearch::GeneratePrefixPivotal(std::string entry, std::vector<Qgram>& prefix, std::vector<Qgram>& pivotal)
{
	//Generate q-grams & sort by frequency
	std::vector<Qgram> qgrams;
	for (int i = 0; i < entry.length() - q + 1; i++)
	{
		std::string qgram = entry.substr(i, q);
		qgrams.push_back(Qgram(qgram, qgramFrequency[qgram], i));
	}
	std::sort(qgrams.begin(), qgrams.end());

	//Generate prefixes & sort by position
	int max = q * threshold + 1;
	if (max > qgrams.size())
		max = qgrams.size();
	for (int i = 0; i < max; i++)
		prefix.push_back(qgrams[i]);
	std::sort(prefix.begin(), prefix.end(), PositionOrdering);

	//Generate pivotals
	max = threshold + 1;
	if (max > prefix.size())
		max = prefix.size();
	int pos = -1;
	int size = 0;
	for (int i = 0; i < prefix.size(); i++)
	{
		if (prefix[i].pos < pos)
			continue;
		pivotal.push_back(prefix[i]);
		pos = prefix[i].pos + q;
		size++;
		if (size == max)
			break;
	}
}

void PivotalSearch::Indexing()
{
	//Count q-gram frequency (frequency counting and division into q-grams has been split to save memory)
	std::sort(db->db.begin(), db->db.end(), CompareLength);
	CountFrequency();

	//Iterate over db, generate prefixes and pivotals and create index
	for (int i = 0; i < db->db.size(); i++)
	{
		std::vector<Qgram> prefix;
		std::vector<Qgram> pivotal;
		GeneratePrefixPivotal(db->db[i].sequence, prefix, pivotal);

		for (int j = 0; j < prefix.size(); j++)
		{
			std::string qgram = prefix[j].s;
			indexPrefixes[db->db[i].sequence.length()][qgram].push_back(PrefixEntry(i, prefix[j].pos));
		}
		lastPrefixFrequency.push_back(prefix[prefix.size() - 1].frequency);

		for (int j = 0; j < pivotal.size(); j++)
		{
			std::string qgram = pivotal[j].s;
			indexPivotals[db->db[i].sequence.length()][qgram].push_back(PivotalEntry(i, pivotal[j].pos, j));
		}
		pivotals.push_back(pivotal);
	}
}

//Searching
bool PivotalSearch::PigeonRing(std::string candidate, int pivotalNr, std::string query, std::vector<Qgram>& pivotal)
{
	if (chainLength == 0) //Chain length 0 = off, chain length threshold + 1 = equal to alignment filter
		return true;

	int minLength = std::min(chainLength, (int)pivotal.size());

	//Pre calculate thresholds
	double single = (double)threshold / minLength;
	int nrOfSegments = pivotal.size();
	int nr = pivotalNr + 1;
	//Find prefix-viable chain
	double errors = 0; //double because hammingdistance returns double
	for (int j = 1; j < minLength; j++)
	{
		int index = nr % nrOfSegments;
		Qgram piv = pivotal[index];
		int startQ = std::max(0, piv.pos - threshold);
		int lengthQ = std::min((int)query.length(), (piv.pos + q + threshold)) - startQ;
		int startC = piv.pos;
		int lengthC = q;
		errors += SubstringEditDistance(query, candidate, startQ, lengthQ, startC, lengthC);
		if (errors > ((double)(j + 1) * single))
			return false;
		nr++;
	}
	//Prefix-viable chain found -> pigeonring filter passed
	return true;
}

bool PivotalSearch::AlignmentFilter(std::string query, std::string candidate, std::vector<Qgram>& pivotal)
{
	//Alignment filter is worse than pigeonring for Pivotal
	int errors = 0;
	for (int i = 0; i < pivotal.size(); i++)
	{
		Qgram piv = pivotal[i];
		int startQ = std::max(0, piv.pos - threshold);
		int lengthQ = std::min((int)query.length(), (piv.pos + q + threshold)) - startQ;
		int startC = piv.pos;
		int lengthC = q;
		errors += SubstringEditDistance(query, candidate, startQ, lengthQ, startC, lengthC);
		if (errors > threshold)
			return false;
	}
	//Alignment filter passed
	return true;
}

std::vector<Result> PivotalSearch::SearchSequence(std::string query)
{
#pragma region Initialization
	std::clock_t startC = std::clock();
	std::vector <Result> result;
	std::vector<bool> checked = std::vector(db->db.size(), false);
	double duration = (std::clock() - startC) / (CLOCKS_PER_SEC / 1000);
	std::cout << duration << ";"; //Initialization
#pragma endregion

#pragma region Searching
	startC = std::clock();
	std::vector<Qgram> prefix;
	std::vector<Qgram> pivotal;
	GeneratePrefixPivotal(query, prefix, pivotal);

	for (int i = 0; i < prefix.size(); i++)
	{
		std::string pre = prefix[i].s;
		auto startLength = indexPivotals.lower_bound((int)query.length() - threshold);
		auto endLength = indexPivotals.upper_bound((int)query.length() + threshold);
		for (auto lengthIT = startLength; lengthIT != endLength; lengthIT++)
		{
			std::vector<PivotalEntry> entryList = (*lengthIT).second[pre];
			for (int j = 0; j < entryList.size(); j++)
			{
				PivotalEntry entry = entryList[j];																					//Pigeonhole: piv(entry) intersect pre(query) = non empty
				if (!checked[entry.index]																							//Check if already checked
					&& (prefix[prefix.size() - 1].frequency > lastPrefixFrequency[entry.index])										//last(pre(query)) > last(pre(entry))
					&& abs(entry.pos - prefix[i].pos) <= threshold)																	//Q-gram position filter	
				{
					if (PigeonRing(db->db[entry.index].sequence, entry.pivotalNr, query, pivotals[entry.index]))					//Pigeonring
					{
						int score = LengthAwareED(query, db->db[entry.index].sequence, 0, query.length(), 0, db->db[entry.index].sequence.length(), threshold);
						if (score <= threshold)
							result.push_back(Result(db->db[entry.index].index, score, true));
						checked[entry.index] = true;
					}
				}
			}
		}
	}
	for (int i = 0; i < pivotal.size(); i++)
	{
		std::string piv = pivotal[i].s;
		auto startLength = indexPrefixes.lower_bound((int)query.length() - threshold);
		auto endLength = indexPrefixes.upper_bound((int)query.length() + threshold);
		for (auto lengthIT = startLength; lengthIT != endLength; lengthIT++)
		{
			std::vector<PrefixEntry> entryList = (*lengthIT).second[piv];
			for (int j = 0; j < entryList.size(); j++)
			{
				PrefixEntry entry = entryList[j];																//Pigeonhole: piv(query) intersect pre(entry) = non empty
				if (!checked[entry.index]																		//Check if already checked
					&& (prefix[prefix.size() - 1].frequency <= lastPrefixFrequency[entry.index])				//last(pre(query)) <= last(pre(entry))
					&& abs(entry.pos - pivotal[i].pos) <= threshold)											//Q-gram position filter 
				{
					if (PigeonRing(query, i, db->db[entry.index].sequence, pivotal))							//Pigeonring
					{
						int score = LengthAwareED(query, db->db[entry.index].sequence, 0, query.length(), 0, db->db[entry.index].sequence.length(), threshold);
						if (score <= threshold)
							result.push_back(Result(db->db[entry.index].index, score, true));
						checked[entry.index] = true;
					}
				}
			}
		}
	}
	duration = (std::clock() - startC) / (CLOCKS_PER_SEC / 1000);
	std::cout << ";" << duration << ";"; //Initialization
#pragma endregion

#pragma region Sorting results
	startC = std::clock();
	std::sort(result.rbegin(), result.rend());
	duration = (std::clock() - startC) / (CLOCKS_PER_SEC / 1000);
	std::cout << duration << ";\n"; //Sorting
#pragma endregion
	return result;
}
