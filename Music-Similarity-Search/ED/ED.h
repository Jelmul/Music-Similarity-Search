/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
#include "../General/SimilaritySearch.h"

class ED : public SimilaritySearch
{
private:
	int k;
	int WagnerFischer(std::string query, std::string candidate);
public:
	ED() {}
	ED(Database* db, int k);
	std::vector<Result> SearchSequence(std::string query) override;
};