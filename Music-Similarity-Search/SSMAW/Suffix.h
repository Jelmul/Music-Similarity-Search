/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
#include <string>

struct Suffix
{
	std::string suffix;
	int pos;
	Suffix() {};
	Suffix(std::string suffix, int pos)
	{
		this->suffix = suffix;
		this->pos = pos;
	}

	bool operator<(const Suffix& rhs) const noexcept
	{
		return this->suffix < rhs.suffix;
	}
};