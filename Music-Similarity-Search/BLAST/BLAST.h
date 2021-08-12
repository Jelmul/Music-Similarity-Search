/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
#include "../General/SimilaritySearch.h"
#include "HSW.h"

class BLAST : public SimilaritySearch
{
private:
	int T;
	int A;
	int X;
	int q;
	int substitutionMatrix[53][53];
	std::unordered_map<std::string, std::vector<HSW>> GenerateHSWIndex(std::string query);
	void UngappedExtension(int qpos, int epos, std::string query, std::string entry, int& score, int& ql, int& qr, int& el, int& er);
public:
	BLAST() {}
	BLAST(Database* db, int T, int A, int X, int q);
	std::vector<Result> SearchSequence(std::string query) override;
};