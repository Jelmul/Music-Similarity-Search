/*
	Written by Jelle Mulyadi, 2021
*/

#include "SimilaritySearch.h"

std::vector<Result> SimilaritySearch::SearchSequenceID(std::string queryID)
{
	std::string query = db->index[queryID];
	return SearchSequence(query);
}