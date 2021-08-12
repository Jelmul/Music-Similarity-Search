/*
	Written by Jelle Mulyadi, 2021
*/

#include "Scoring.h"
#include <algorithm>
int CharacterScore(char a)
{
	if ('A' <= a && a <= 'Z')
		return a - 64;
	else if ('a' <= a && a <= 'z')
		return -(a - 96);
	else // '-'
		return 0;
}

int CharacterIndex(char a)
{
	if ('A' <= a && a <= 'Z')
		return a - 'A';
	else if ('a' <= a && a <= 'z')
		return (a - 'a') + 26;
	else // '-'
		return 52;
}

int IntervalScore(char a, char b)
{
	//Score ranging from 1 to -2 and back to 1, the further away from a note the lesser the score
	int distance = (abs(CharacterScore(a) - CharacterScore(b)) % 7);
	if (distance == 0)
		return 1;
	else if (distance == 1)
		return 0;
	else if (distance == 2)
		return -1;
	else if (distance == 3)
		return -2;
	else if (distance == 4)
		return -2;
	else if (distance == 5)
		return -1;
	else if (distance == 6)
		return 0;
}

int clamp(int val)
{
	if (val < 0)
		return 0;
	else
		return val;
}