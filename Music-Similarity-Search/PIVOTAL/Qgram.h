/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
#include <string>

struct Qgram
{
	//Variables
	std::string s;
	int pos;
	int frequency;
	//Methods
	Qgram() {}
	Qgram(std::string s, int frequency, int pos)
	{
		this->s = s;
		this->frequency = frequency;
		this->pos = pos;
	}
	bool operator<(const Qgram& rhs) const noexcept
	{
		return this->frequency < rhs.frequency;
	}
};