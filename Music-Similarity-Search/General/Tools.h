/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
#include<string>
#include "Result.h"
#include "Entry.h"

int LengthAwareED(std::string query, std::string candidate, int startQ, int lengthQ, int startC, int lengthC, int threshold);
int SubstringEditDistance(std::string query, std::string candidate, int startQ, int lengthQ, int startC, int lengthC);
double SubstringHammingDistance(std::string query, std::string candidate, int startQ, int lengthQ, int startC, int lengthC);
bool CompareLength(Entry i, Entry j);