/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
#include "../General/SimilaritySearch.h"

class LA : public SimilaritySearch
{
private:
	int k;
	int substitutionMatrix[53][53];
	int SmithWaterman(std::string query, std::string candidate);
public:
	LA() {}
	LA(Database* db, int k);
	std::vector<Result> SearchSequence(std::string query) override;
};