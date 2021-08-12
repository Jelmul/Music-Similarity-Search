/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
#include "../General/SimilaritySearch.h"

class GA : public SimilaritySearch
{
private:
	int k;
	int substitutionMatrix[53][53];
	int Gotoh(std::string query, std::string candidate);
public:
	GA() {}
	GA(Database* db, int k);
	std::vector<Result> SearchSequence(std::string query) override;
};