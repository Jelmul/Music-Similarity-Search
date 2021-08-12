/*
	Written by Jelle Mulyadi, 2021
*/

#pragma once
struct HSW
{
	int score;
	int pos;
	HSW() {};
	HSW(int score, int pos)
	{
		this->score = score;
		this->pos = pos;
	}
};