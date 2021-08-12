/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
#include "Database.h"
#include "Result.h"
#include <string>
#include <vector>

class SimilaritySearch
{
public:
	//Methods
	Database* db;
	SimilaritySearch() {};
	virtual std::vector<Result> SearchSequence(std::string query) = 0;
	std::vector<Result> SearchSequenceID(std::string queryID);
};