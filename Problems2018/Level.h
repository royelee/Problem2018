//
//  Level.h
//  XProject
//
//  Created by Roye Li on 12/22/15.
//  Copyright (c) 2015 Roye Li. All rights reserved.
//

#ifndef __XProject__Level__
#define __XProject__Level__
#include <vector>
#include <string>
#include <iostream>

struct TreeNode {
	int val;
	TreeNode *left;
	TreeNode *right;
	TreeNode(int x) : val(x), left(nullptr), right(nullptr) {}
};

struct ListNode {
	int val;
	ListNode *next;
	ListNode(int x) : val(x), next(nullptr) {}
};

struct TreeLinkNode {
	int val;
	TreeLinkNode *left, *right, *next;
	TreeLinkNode(int x) : val(x), left(nullptr), right(nullptr), next(nullptr) {}
};

struct UndirectedGraphNode {
	int label;
	std::vector<UndirectedGraphNode *> neighbors;
	UndirectedGraphNode(int x) : label(x) {};
};

struct Point {
	int x;
	int y;
	Point() : x(0), y(0) {}
	Point(int a, int b) : x(a), y(b) {}
};

inline std::string BoolToStr(bool a) { return a ? "true" : "false"; }

inline void Verify( bool condition, const std::string& message )
{
	if( condition )
	{
		std::cout << "Test " << message << " success." << std::endl;
	}
	else
	{
		std::cout << "Test " << message << " failed. " << std::endl;
	}
}

template< typename T >
void Verify( const T& actual, const T& expect, const std::string& message )
{
	std::cout << "==========================" << std::endl;
	if( expect == actual )
	{
		std::cout << "Test " << message << " success." << std::endl;
	}
	else
	{
		std::cout << "Test " << message << " failed. " << std::endl;
		std::cout << "Expect:" << expect << "  "
				  << "Acutal:" << actual << std::endl;
	}
}

inline void Verify( bool condition )
{
	Verify( condition, "Message" );
}

#define D_L(LEVELNAME) class Level##LEVELNAME\
{\
public:\
	void Run();\
};

void preOrderTrasveral( TreeNode* node );

D_L(1);
D_L(2);
D_L(3);
D_L(4);
D_L(5);
D_L(6);
D_L(7);
D_L(8);
D_L(9);
D_L(10);
D_L(11);
D_L(12);
D_L(13);
D_L(14);
D_L(15);
D_L(16);
D_L(17);
D_L(18);
D_L(19);
D_L(20);
D_L(21);
D_L(22);
D_L(23);
D_L(24);
D_L(25);


#endif /* defined(__XProject__Level__) */
