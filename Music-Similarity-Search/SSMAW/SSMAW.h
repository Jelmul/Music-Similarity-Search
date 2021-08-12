/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
#include "../General/SimilaritySearch.h"
#include "Trie.h"
#include <set>

class SSMAW : public SimilaritySearch
{
private:
	//Variables
	int min;
	int max;
	int nrOfResults;
	std::vector<int> MAWcounts;
	Trie MAWsTrie = Trie();
	//Methods
	void CalculateMAWs(std::string entry, std::vector<int>& SA, std::vector<int>& LCP,
		std::vector<std::bitset<53>>& B1, std::vector<std::bitset<53>>& B2, std::set<std::string>& MAWs);
	void Indexing();
public:
	SSMAW() {};
	SSMAW(Database* db, int min, int max, int nrOfResults);
	std::vector<Result> SearchSequence(std::string query) override;
	int indexTime;
};