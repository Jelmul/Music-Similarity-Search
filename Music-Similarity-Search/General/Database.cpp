/*
	Written by Jelle Mulyadi, 2021
*/

#include "Database.h"
#include "Tools.h"
#include <fstream>
#include <algorithm>
#include <iostream>

Database::Database(std::string name)
{
	std::fstream file;
	file.open(name, std::ios::in);
	if (file.is_open())
	{
		std::string line;
		while (std::getline(file, line))
		{
			int j = 0;
			while (line[j] != ' ')
				j++;
			std::string id = line.substr(0, j);
			std::string sub = line.substr(j + 1);
			index[id] = sub;
			rIndex.push_back(id);
			int length = sub.length();
			if (length > maxLength)
				maxLength = length;
			if (length < minLength)
				minLength = length;
			sumLengths += length;
			db.push_back(Entry(id, sub));
		}
		alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz-";
	}
	file.close();
}