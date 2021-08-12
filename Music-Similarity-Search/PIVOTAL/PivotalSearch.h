/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
#include "../General/SimilaritySearch.h"
#include "Qgram.h"
#include "IndexEntry.h"
#include <map>
#include <unordered_map>

class PivotalSearch : public SimilaritySearch
{
private:
	//Variables
	int q;
	int threshold;
	int chainLength;
	std::unordered_map<std::string, int> qgramFrequency;
	std::vector<int> lastPrefixFrequency;
	std::map<int, std::unordered_map<std::string, std::vector<PrefixEntry>>> indexPrefixes;
	std::map<int, std::unordered_map<std::string, std::vector<PivotalEntry>>> indexPivotals;
	std::vector<std::vector<Qgram>> pivotals;
	//Methods
	void CountFrequency();
	void GeneratePrefixPivotal(std::string entry, std::vector<Qgram>& prefix, std::vector<Qgram>& pivotal);
	void Indexing();
	bool PigeonRing(std::string candidate, int pivotalNr, std::string query, std::vector<Qgram>& pivotal);
	bool AlignmentFilter(std::string query, std::string candidate, std::vector<Qgram>& pivotal);
public:
	//Methods
	PivotalSearch() {};
	PivotalSearch(Database* db, int q, int threshold, int chainLength);
	std::vector<Result> SearchSequence(std::string query) override;
	int indexTime;
};
