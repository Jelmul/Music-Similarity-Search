/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
#include "../General/SimilaritySearch.h"
#include <map>
#include <unordered_map>

class PassJoin : public SimilaritySearch
{
private:
	//Variables
	int threshold;
	int chainLength;
	std::map<int, std::vector<std::unordered_map<std::string, std::vector<int>>>> index;
	//Methods
	void Indexing();
	void SubstringSelection(std::string query, int pos, int segment, int segLength, int entryLength, int& start, int& end);
	void Verification(std::string query, std::vector<int>& candidates, int entryPos, int entrySeg, int segLength, std::vector<double>& thresholds, std::vector<Result>& result, std::vector<bool>& checked);
	bool PigeonRing(std::string candidate, int segmentNr, int segmentPos, std::string query, std::vector<double>& thresholds);
	bool AlignmentFilter(std::string candidate, std::string query);
public:
	PassJoin() {};
	PassJoin(Database* db, int threshold, int chainLength);
	std::vector<Result> SearchSequence(std::string query) override;
	int indexTime;
};