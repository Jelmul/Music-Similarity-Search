/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
#include "Entry.h"
#include <string>
#include <vector>
#include <unordered_map>

class Database
{
public:
	//Variables
	std::unordered_map<std::string, std::string> index;
	std::vector<std::string> rIndex;
	std::vector<Entry> db;
	int maxLength = 0;
	int minLength = INT_MAX;
	int sumLengths = 0;
	std::string alphabet;
	//Methods
	Database() {};
	Database(std::string name);
};