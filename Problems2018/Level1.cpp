//
//  Level1.cpp
//  XProject
//
//  Created by Roye Li on 6/17/15.
//  Copyright (c) 2015 Roye Li. All rights reserved.
//
#include "Level.h"
#include <algorithm>
#include <map>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>

extern "C"
{
#include <ctype.h>
}

using namespace std;

#pragma mark -

static int myAtoi(string str)
{
	int64_t returnV = 0;
	size_t i = 0;
	while( i < str.size() && str[i] == ' ')
		i++;
	
	int8_t sign = 1;
	if( i < str.size() && ( str[i] == '-' || str[i] == '+' ) )
	{
		sign = str[i] == '-' ? -1 : 1;
		i++;
	}
	
	while( i < str.size() )
	{
		if( str[i] >= '0' && str[i] <= '9' )
		{
			returnV = returnV * 10 + str[i] - '0';
			if( returnV * sign <= std::numeric_limits<int32_t>::min() )
				return std::numeric_limits<int32_t>::min();
			else if( returnV * sign >= std::numeric_limits<int32_t>::max() )
				return std::numeric_limits<int32_t>::max();
		}
		else
		{
			break;
		}
		
		i++;
	}
	return returnV * sign;
}

static void testMyAtoi()
{
	Verify( myAtoi( "42" ) == 42, "42" );
	Verify( myAtoi( "   -42" ) == -42, "-42" );
	Verify( myAtoi( "4193 with words" ) == 4193, "4193" );
	Verify( myAtoi( "words and 987" ) == 0, "987" );
	Verify( myAtoi( "-91283472332" ) == -2147483648, "-2147483648" );
	Verify( myAtoi( "-2147483648" ) == -2147483648, "-2147483648" );
	Verify( myAtoi( "2147483647" ) == 2147483647, "2147483647" );
	Verify( myAtoi( "0" ) == 0, "0" );
	Verify( myAtoi( "+49" ) == 49, "+49" );
	Verify( myAtoi( "+-49" ) == 0, "+-49" );
	Verify( myAtoi( "-+49" ) == 0, "-+49" );
}

#pragma mark -

static int gcd( int a, int b )
{
	if( b == 0 )    return a;
	else            return gcd(b, a%b );
}

static int maxPoints(vector<Point>& points)
{
	int maxValue = 0;
	for( size_t i = 0; i < points.size(); i++ )
	{
		int vertical = 0, samePoint = 0;
		map< pair< int, int >, int > m;
		for( size_t j = i + 1; j < points.size(); j++ )
		{
			auto& pointA = points[ i ];
			auto& pointB = points[ j ];
			
			if( pointA.x == pointB.x && pointA.y == pointB.y)
			{
				samePoint++;
			}
			else if( pointA.y == pointB.y )
			{
				vertical++;
			}
			else
			{
				auto diffX = pointA.x - pointB.x;
				auto diffY = pointA.y - pointB.y;
				auto g = gcd( diffX, diffY );
				if( g != 0 )
				{
					diffX /= g;
					diffY /= g;
				}
				
				auto p = make_pair( diffX, diffY );
				if( m.find( p ) != m.end() )
				{
					m[p]++;
				}
				else
				{
					m[p] = 1;
				}
			}
		}
		
		int maxInMap = 0;
		for( const auto& v : m )
		{
			if(v.second > maxInMap)
			{
				maxInMap = v.second;
			}
		}
		int localMax = max( vertical, maxInMap );
		maxValue = max( localMax + samePoint + 1, maxValue );
	}
	
	return maxValue;
}

static void testMaxPoints()
{
	vector<pair<vector<Point>, int>> tests = {
		make_pair( vector<Point>(), 0 ),
		make_pair( vector<Point>{ {-1, -1} } , 1 ),
		make_pair( vector<Point>{ {1,1}, {2,2}, {3,3} }, 3 ),
		make_pair( vector<Point>{ {1,1},{3,2},{5,3},{4,1},{2,3},{1,4} }, 4 ),
		make_pair( vector<Point>{ {0,0},{1,1},{0,0} }, 3 )
	};
	size_t index = 0;
	for( auto& test : tests )
	{
		cout << "Testing ... " << index++ << endl;
		Verify( maxPoints( test.first ) == test.second, to_string( test.second ) );
	}
}

#pragma mark -

static int divide(int dividend, int divisor) {
	if( divisor == 0 || ( dividend == INT_MIN && divisor == -1 ) )
		return INT_MAX;
	
	int sign = ( dividend > 0 ) ^ ( divisor > 0 ) ? -1 : 1;
	int64_t divisor_a = abs( static_cast< int64_t >( divisor ) );
	int64_t remainder = abs( static_cast< int64_t >( dividend ) );
	int result = 0;
	while( remainder >= divisor_a )
	{
		int offset = 0;
		while( divisor_a << ( offset + 1 ) <= remainder )
		{
			offset++;
		}
		
		result += ( 1 << offset );
		remainder -= ( divisor_a << offset );
	}
	
	return result * sign;
}

static void testDivide()
{
	int32_t i = numeric_limits<int32_t>::min();
	int32_t j = -1;
	std::cout << std::hex << i << endl;
	std::cout << std::hex << j << endl;
	
	Verify( divide( 10, 3 ), 10 / 3, "10/3" );
	Verify( divide( 10, 2 ), 10 / 2, "10/2" );
	Verify( divide( 10, 1 ), 10 / 1, "10/1" );
	Verify( divide( 10, -3 ), -10 / 3, "10/-3" );
	Verify( divide( 3, 10 ), 3 / 10, "3/10" );
	Verify( divide( -10, 3 ), -10 / 3, "-10/3" );
	Verify( divide( -10, -3 ), -10 / -3, "-10/-3" );
	Verify( divide( -10, 0 ), INT_MAX, "10/0" );
	Verify( divide( 2, 2 ), 1, "2/2" );
	Verify( divide( INT_MAX, INT_MIN ), INT_MAX / INT_MIN, "max/min" );
	Verify( divide( INT_MIN, -1 ), INT_MAX, "min/-1" );
}

#pragma mark -

static string fractionToDecimal(int numerator, int denominator) {
	if( numerator == 0 )
		return "0";
	if( denominator == 0 )
		return "";
	
	int sign = ( numerator > 0 ) ^ ( denominator > 0 ) ? -1 : 1;
	int64_t num = abs(static_cast<int64_t>(numerator));
	int64_t den = abs(static_cast<int64_t>(denominator));
	
	unordered_map< int, int > m;
	string rtn;
	rtn += to_string( num / den );
	const size_t intPartSize = rtn.size();
	num = num % den;
	int index = rtn.size();
	while( num != 0 )
	{
		num = num * 10;
		
		if( index == intPartSize )
		{
			rtn += ".";
			index++;
		}
		
		// Found, it's recursive.
		if( m.find( num ) != m.end() )
		{
			rtn += ")";
			rtn.insert( m[num], "(" );
			break;
		}
		else
		{
			int64_t div = num / den;
			rtn += to_string( div );
			m[num] = index;
			index++;
		}
		
		num = num % den;
	}
	
	return ( sign == -1 ? "-" : "" ) + rtn;
}

static void testFractionToDecimal()
{
	cout << abs(INT_MIN) << endl;
	Verify<string>( fractionToDecimal( INT_MIN, -1 ), "2147483648", "INT_MIN/-1" );
	Verify<string>( fractionToDecimal( 1, 2 ), "0.5", "1/2" );
	Verify<string>( fractionToDecimal( 2, 1 ), "2", "2" );
	Verify<string>( fractionToDecimal( 2, 3 ), "0.(6)", "2/3" );
	Verify<string>( fractionToDecimal( 20, 7 ), "2.(857142)", "20/7" );
	Verify<string>( fractionToDecimal( 2, 0 ), "", "2/0" );
	Verify<string>( fractionToDecimal( 0, -1 ), "0", "0/-1" );
}

#pragma mark -

static void helperDFS(vector<vector<char>>& board, int i, int j)
{
	if( i >= 0 && i < board.size() && j >= 0 && j < board[0].size() )
	{
		if( board[i][j] == 'O' )
		{
			board[i][j] = 'F';	// flip
			helperDFS(board, i+1, j );
			helperDFS(board, i-1, j );
			helperDFS(board, i, j+1 );
			helperDFS(board, i, j-1 );
		}
	}
}

static void solve(vector<vector<char>>& board) {
	if( board.size() == 0 )
		return;
	if( board[0].size() == 0 )
		return;
	
	size_t w = board.size();
	size_t h = board[0].size();
	// Four borders
	for( int i = 0; i < w; i++ )
	{
		for( int j = 0; j < h; j++ )
		{
			if( i == 0 || i == w - 1 || j == 0 || j == h -1 )
				helperDFS( board, i, j );
		}
	}
	
	// flip all 'O' -> 'X', 'F' -> 'O'
	for( int i = 0; i < w; i++ )
	{
		for( int j = 0; j < h; j++ )
		{
			if( board[i][j] == 'O' )
				board[i][j] = 'X';
			else if( board[i][j] == 'F' )
				board[i][j] = 'O';
		}
	}
}

#pragma mark -

static int ladderLength(string beginWord, string endWord, vector<string>& wordList)
{
	unordered_set< string > wordSet;
	for( const auto& word : wordList )
		wordSet.insert( word );
	if( wordSet.find( endWord) == wordSet.end() )
		return 0;
	
	int step = 1;
	unordered_set< string > nbsFront, nbsEnd;
	unordered_set< string >* pFront = &nbsFront;
	unordered_set< string >* pEnd = &nbsEnd;
	nbsFront.insert( beginWord );
	nbsEnd.insert( endWord );
	while( !pFront->empty() && !pEnd->empty() )
	{
		if( pFront->size() > pEnd->size() )
		{
			swap( pFront, pEnd );
		}
		
		unordered_set<string> newFront;
		for( const auto& word : *pFront )
		{
			for( int i = 0; i < word.size(); i++ )
			{
				string currentWord = word;
				for( int ch = 0; ch < 26; ch++ )
				{
					char c = 'a' + ch;
					currentWord[i] = c;
					if( pEnd->find( currentWord ) != pEnd->end() )
					{
						return step + 1;
					}
					else if( wordSet.find( currentWord ) != wordSet.end() )
					{
						newFront.insert( currentWord );
						wordSet.erase( currentWord );
					}
				}
			}
		}
		
		pFront->swap(newFront);
		step++;
	}
	
	return 0;
}

static void testLadderLength()
{
	vector<string> dict = {"hot","dot","dog","lot","log","cog"};
	Verify( ladderLength( "hit","cog", dict ), 5, "default test case" );
	vector<string> dict1 = {"cog","hot","dot","dog","lot","log"};
	Verify( ladderLength( "hit","cog", dict1 ), 5, "default test case" );
	vector<string> dict2 = {"si","go","se","cm","so","ph","mt","db","mb","sb","kr","ln","tm","le","av","sm","ar","ci","ca","br","ti","ba","to","ra","fa","yo","ow","sn","ya","cr","po","fe","ho","ma","re","or","rn","au","ur","rh","sr","tc","lt","lo","as","fr","nb","yb","if","pb","ge","th","pm","rb","sh","co","ga","li","ha","hz","no","bi","di","hi","qa","pi","os","uh","wm","an","me","mo","na","la","st","er","sc","ne","mn","mi","am","ex","pt","io","be","fm","ta","tb","ni","mr","pa","he","lr","sq","ye"};
	Verify( ladderLength( "qa","sq", dict2 ), 6, "default test case" );
}

#pragma mark -

struct DoubleLinkListNode
{
	DoubleLinkListNode* pre{ nullptr };
	DoubleLinkListNode* next{ nullptr };
	int                 value{ -1 };
	int                 key{ -1 };
};

class DoubleLinkList
{
public:
	DoubleLinkList()
	{
		m_pHead = new DoubleLinkListNode;
		m_pTail = new DoubleLinkListNode;
		m_pHead->next = m_pTail;
		m_pTail->pre = m_pHead;
	}
	
	void PushBack( DoubleLinkListNode* pNode )
	{
		pNode->next = m_pTail;
		pNode->pre = m_pTail->pre;
		m_pTail->pre->next = pNode;
		m_pTail->pre = pNode;
	}
	
	DoubleLinkListNode* Head() const
	{
		return m_pHead->next;
	}
	
	void PopHead()
	{
		DoubleLinkListNode* pHeadNext = m_pHead->next;
		pHeadNext->next->pre = m_pHead;
		m_pHead->next = pHeadNext->next;
	}
	
	void MoveToEnd( DoubleLinkListNode* pNode )
	{
		DoubleLinkListNode* pNodeNext = pNode->next;
		DoubleLinkListNode* pNodePre = pNode->pre;
		pNodePre->next = pNodeNext;
		pNodeNext->pre = pNodePre;
		PushBack( pNode );
	}
	
private:
	DoubleLinkListNode* m_pTail{ nullptr };
	DoubleLinkListNode* m_pHead{ nullptr };
};

class LRUCache {
public:
	LRUCache(int capacity) {
		m_capacity = capacity;
	}
	
	int get(int key) {
		auto pNode = getNode( key );
		return pNode ? pNode->value : -1;
	}
	
	DoubleLinkListNode* getNode(int key) {
		if( m_map.find( key ) != m_map.end() )
		{
			auto pNode = m_map[key];
			m_doubleLinkList.MoveToEnd( pNode );
			return pNode;
		}
		
		return nullptr;
	}
	
	void put(int key, int value) {
		DoubleLinkListNode* pNode = getNode( key );
		if( pNode )
		{
			pNode->value = value;
		}
		else
		{
			if( m_map.size() >= m_capacity )
			{
				auto pHead = m_doubleLinkList.Head();
				auto pHeadKey = pHead->key;
				m_doubleLinkList.PopHead();
				m_map.erase( pHeadKey );
			}
			
			DoubleLinkListNode* pNode = new DoubleLinkListNode;
			pNode->key = key;
			pNode->value = value;
			m_doubleLinkList.PushBack( pNode );
			m_map[ key ] = pNode;
		}
	}
	
private:
	unordered_map< int, DoubleLinkListNode* > m_map;
	DoubleLinkList                            m_doubleLinkList;
	int 					                  m_capacity;
};

void testLRUCache()
{
	LRUCache cache( 2 /* capacity */ );
	cache.put(1, 1);
	cache.put(2, 2);
	Verify( cache.get(1), 1, "default test case" );
	cache.put(3, 3);    // evicts key 2
	Verify( cache.get(2), -1, "default test case" );
	cache.put(4, 4);    // evicts key 1
	Verify( cache.get(1), -1, "default test case" );
	Verify( cache.get(3), 3, "default test case" );
	Verify( cache.get(4), 4, "default test case" );
	
	LRUCache cache1( 1 /* capacity */ );
	cache1.put(2, 1);
	Verify( cache1.get(2), 1, "default test case" );
	cache1.put(3, 2);
	Verify( cache1.get(2), -1, "default test case" );
	Verify( cache1.get(3), 2, "default test case" );
}

#pragma mark -

static int numDecodings(string s)
{
	// "231" -> 2 + 31, 23 + 1
	// State change ( 231 ) = 1 + state(31) + 2 + state(1)
	// c < '3' && c+1 < '7' '? (numDecodings( s - 1 ) + numDecodings( s - 2) ) : numDecodings( s - 1 )
	// if c == '0' || c > '3' || c+1 > '6'
	// 	numDecodings( s - 1 )
	// else
	//	numDecodings( s - 1 ) + numDecodings( s - 2 )
	if(s.size() == 0 )
		return 0;
	
	vector<int> dp( s.size() + 2 );
	dp[s.size()] = 1;
	dp[s.size() + 1] = 1;
	for( int i = s.size() - 1; i >= 0; i-- )
	{
		char c = s[i];
		if( c == '0' )
			continue;
		else if( i + 1 < s.size() &&
				( ( c == '1' && s[i+1] > '0' ) || ( c == '2' && s[i+1] > '0' && s[i+1] <= '6' ) ) )
			dp[i] = dp[i+1] + dp[i+2];
		else if( i + 1 < s.size() && ( ( c == '1' || c == '2' ) && s[i+1] == '0' ) )
			dp[i] = dp[i+2];
		else
			dp[i] = dp[i+1];
	}
	
	return dp[0];
}

void testNumDecoding()
{
	int i = 0;
	Verify( numDecodings("0"), 0, "Test case " + to_string(++i) );
	Verify( numDecodings("17494"), 2, "Test case " + to_string(++i) );
	Verify( numDecodings("226"), 3, "Test case " + to_string(++i) );
	Verify( numDecodings("12"), 2, "Test case " + to_string(++i) );
	Verify( numDecodings("10"), 1, "Test case " + to_string(++i) );
	Verify( numDecodings("11"), 2, "Test case " + to_string(++i) );
	Verify( numDecodings("27"), 1, "Test case " + to_string(++i) );
	Verify( numDecodings("4757562545844617494555774581341211511296816786586787755257741178599337186486723247528324612117156948"), 589824, "Test case " + to_string(++i) );
}

#pragma mark -

static bool isMatch(string s, string p)
{
	// if s[1..n], p[1..n]
	// s[n] match p[n] so isMatch( s[1..n-1], p[1..n-1] )
	// p[n] could be c, ?, *
	// p[n] - c or ?, if c match s[n] <- isMatch( s[1..n-1], p[1..n-1] )
	// p[n] - *, isMatch( s[1..n], p[1..n-1] ) or isMatch( [i..n-1], p[1..n])
	
	// Table
	// p/s 0 1 2 3 4 5 ... n
	// 0
	// 1
	// 2
	// .
	// n
	
	vector< vector< bool > > dp( p.size() + 1, vector< bool >( s.size() + 1 ) );
	// Initialization
	dp[0][0] = true;
	for( int i = 1; i < dp.size(); i++ )
	{
		if( p[i-1] == '*' )
		{
			dp[i][0] = dp[i-1][0];
		}
	}
	
	// Bottom up DP
	for( int i = 1; i < dp.size(); i++ )
	{
		for( int j = 1; j < dp[i].size(); j++ )
		{
			char pc = p[i-1];
			char sc = s[j-1];
			if( pc == '*' )
				dp[i][j] = ( dp[i-1][j] | dp[i][j-1] );
			else if( pc == '?' || pc == sc )
				dp[i][j] = dp[i-1][j-1];
		}
	}
	
	return dp[p.size()][s.size()];
}

static bool isMatchGreedy(string s, string p)
{
	int si = 0, pi = 0, starIdx = -1, matchS = 0;
	while( si < s.size() )
	{
		if( pi < p.length() && ( p[pi] == '?' || p[pi] == s[si] ) )
		{
			pi++;
			si++;
		}
		else if( pi < p.length() && p[pi] == '*' )
		{
			starIdx = pi;
			matchS = si;
			pi++;
		}
		else if( starIdx != -1 )
		{
			pi = starIdx + 1;
			matchS++;
			si = matchS;
		}
		else return false;
	}
	
	while( pi < p.size() && p[pi] == '*' )
		pi++;
		
	return pi == p.size();
}

void testIsMatch()
{
	auto* pFun = isMatchGreedy;
	vector< pair< pair< string, string >, bool > > tests =
	{
		make_pair( make_pair("aa", "*"), true ),
		make_pair( make_pair("aa", "a"), false ),
		make_pair( make_pair("cb", "?a"), false ),
		make_pair( make_pair("adceb", "*a*b"), true ),
		make_pair( make_pair("acdcb", "a*c?b"), false )
	};
	int i = 0;
	for( const auto& t : tests )
	{
		Verify( (*pFun)(t.first.first, t.first.second), t.second, "Test case " + to_string(++i) );
	}
}

#pragma mark -

static vector<vector<int>> threeSum(vector<int>& nums)
{
	vector<vector<int>> results;
	sort( nums.begin(), nums.end() );
	for( int i = 0; i < (int)nums.size() - 2; i++ )
	{
		int left = i + 1;
		int right = (int)nums.size() - 1;
		int target = -nums[i];
		while( left < right )
		{
			if( nums[left] + nums[right] == target )
			{
				results.push_back( { nums[i], nums[left], nums[right] } );
				while( left + 1 < nums.size() && nums[left] == nums[left+1] )	left++;
				left++;
				while( right - 1 >= 0 && nums[right] == nums[right-1] )	        right--;
				right--;
			}
			else if( nums[left] + nums[right] < target )
			{
				left++;
			}
			else
			{
				right--;
			}
		}
		
		while( i + 1 < nums.size() && nums[i+1]==nums[i] )
			i++;
	}
	
	return results;
}

void testThreeSum()
{
	vector<vector<int>> tests = {{-1,0,1,2,-1,-4}, {}};
	for( auto& t : tests )
	{
		cout << " =======   " << endl;
		auto r = threeSum( t );
		for( auto p : r )
		{
			for( auto i : p )
			{
				cout << i << " ";
			}
			cout << endl;
		}
	}
}

#pragma mark -

static string getRange( int64_t start, int64_t end )
{
	return end == start + 2 ? to_string( start + 1 ) : to_string( start + 1 ) + "->" + to_string( end - 1);
}
/*
 * @param nums: a sorted integer array
 * @param lower: An integer
 * @param upper: An integer
 * @return: a list of its missing ranges
 */
static vector<string> findMissingRanges(vector<int> &nums, int lower, int upper)
{
	vector<int64_t> nums64( nums.size() );
	std::copy ( nums.begin(), nums.end(), nums64.begin() );
	nums64.insert( nums64.end(), int64_t(upper) + 1 );
	nums64.insert( nums64.begin(), int64_t(lower) - 1 );
	
	vector<string> results;
	for( int i = 0; i < nums64.size() - 1; i++ )
	{
		int64_t cur = nums64[i];
		int64_t next = nums64[i+1];
		if( cur == next || cur + 1 == next )
			continue;
		else
			results.push_back( getRange( cur, next ) );
	}
	return results;
}

static void testFindMissingRanges()
{
	vector<int> t2 = {};
	auto rs = findMissingRanges( t2, -2147483648, INT_MAX );
	for( auto& r : rs )
		cout << r << " ";
	cout << endl << "*******" << endl;
	
	vector<int> t1 = {0, 1, 3, 50, 75};
	rs = findMissingRanges( t1, 0, 99 );
	for( auto& r : rs )
		cout << r << " ";
}

#pragma mark -

double findMedianSortedArrays(vector<int>& nums1, vector<int>& nums2) {
	int m = nums1.size(), n = nums2.size();
	if (m > n) {
		return findMedianSortedArrays(nums2, nums1);
	}
	int l = 0, r = m, mid1, mid2, maxl, minr;
	while (l <= r) {
		mid1 = (l+r) / 2;
		mid2 = ((m + n + 1) >> 1) - mid1;
		if (mid1 < m && nums1[mid1] < nums2[mid2 - 1]) {
			l = mid1 + 1;
		}
		else if (mid1 > 0 && nums1[mid1-1] > nums2[mid2]) {
			r = mid1 - 1;
		}
		else {
			if (mid1 == 0) {
				maxl = nums2[mid2 - 1];
			} else if (mid2 == 0) {
				maxl = nums1[mid1 - 1];
			} else {
				maxl = max(nums1[mid1 - 1], nums2[mid2 - 1]);
			}
			if ((m + n) % 2 == 1) {
				return maxl;
			}
			if (mid1 == m) {
				minr = nums2[mid2];
			} else if (mid2 == n) {
				minr = nums1[mid1];
			} else {
				minr = min(nums1[mid1], nums2[mid2]);
			}
			return (maxl + minr) / 2.0;
		}
	}
	return -1;
}

#pragma mark - 

static string largestNumber(vector<int>& nums) {
    vector<string> results; 
	for( auto& i : nums )
		results.push_back( to_string( i ) );
		
	sort( results.begin(), results.end(), []( const string& left, const string& right ){
		return left+right > right+left;
	});
	
	if( !results.empty() && results[0] == "0" )
		return "0";
	
	string result;
	for( auto& r : results )
		result += r;
	return result;
}

#pragma mark -

static bool isValidBSTHelper( TreeNode* node, int64_t min, int64_t max )
{
	if( node )
	{
		if( node->val > min && node->val < max )
		{
			return isValidBSTHelper( node->left, min, node->val ) && 
				   isValidBSTHelper( node->right, node->val, max );
		}
		
		return false;
	}
	
	return true;
}

static bool isValidBST(TreeNode* root) {
    return isValidBSTHelper( root, numeric_limits<int64_t>::min(), numeric_limits<int64_t>::max() );
}

static bool isValidBSTInorderHelper(TreeNode* node, TreeNode*& pre )
{
	if( node == nullptr )	return true;
	if( !isValidBSTInorderHelper( node->left, pre ) )	return false;
	if( pre != nullptr && pre->val >= node->val ) 		return false;
	pre = node;
	return isValidBSTInorderHelper( node->right, pre );
}

static bool isValidBSTInorder(TreeNode* root)
{
	TreeNode* pre = nullptr;
	return isValidBSTInorderHelper( root, pre );
}

#pragma mark -

static bool isMatchReg(string s, string p) 
{
    // dp[s.size()+1][p.size() + 1] = false
    // need to init.
    // p[i] == * then it's true, otherwise it's false.
    // 3 case
    // cur -> 'a-z' 
    // cur -> '.'
    // cur -> '*' -> 3 options
    //               match axxx mulitple a
    //               match a single a
    //               match empty
    vector<vector<bool>> dp(p.size() + 1, vector<bool>(s.size() + 1));
    dp[0][0] = true;
    for( int i = 2; i < dp.size(); i++ )
    {
    	if( p[i-1] == '*' && ( p[i-2] == '.' || ( p[i-2] >= 'a' && p[i-2] <= 'z' ) ) )
    		dp[i][0] = dp[i-2][0];
    }
    
    for( int i = 1; i < dp.size(); i++ )
    {
    	for( int j = 1; j < dp[i].size(); j++ )
    	{
    		char pi = p[i - 1];
    		char si = s[j - 1];
    		if( pi == si || pi == '.' )
    		{
    			dp[i][j] = dp[i-1][j-1];
    		}
    		else if( pi == '*' && i - 2 >= 0 )
    		{
    			if( p[i-2] == si || p[i-2] == '.' )
    			{
					dp[i][j] = dp[i-1][j]   | 	// match one
						   	   dp[i][j-1]   |   // match multiple
						   	   dp[i-2][j];		// match zero.
    			}
    			else
    			{
    				// match zero.
    				dp[i][j] = dp[i-2][j];
    			}
    		}
    	}
    }
    
    return dp[p.size()][s.size()];
}

#pragma mark - 

static int reverse(int x) {
	if( x == INT_MIN )	return 0;
	
	int result = 0;
	int flag = x < 0 ? -1 : 1;
	x = abs(x);
	// int from -2^31 ~ 2^31 - 1
	// -2147483648 ~ 2147483647
	while( x != 0 )
	{
		int last = x % 10;
		if( result > INT_MAX / 10  )
			return 0;
		else if( result == INT_MAX / 10 )
		{
			// check last one
			if( flag == 1 && last > 7 )
				return 0;
			else if( flag == -1 && last > 8 )
				return 0;
		}
		
		result = result * 10 + last;
		x = x / 10;
	}
	
	return result * flag;
}

#pragma mark - 

static int lengthOfLongestSubstringV1(string s) {
    // Given "abcabcbb", the answer is "abc", which the length is 3.
	// Given "bbbbb", the answer is "b", with the length of 1.
	// Given "pwwkew", the answer is "wke", with the length of 3. Note that the answer must be a substring, "pwke" is a subsequence and not a substring.	
	int r = 0;
	unordered_set< char > myset;
	for( int i = 0, j = 0; i < s.size() && j < s.size(); )
	{
		if( myset.find( s[j] ) == myset.end() )
		{
			myset.insert( s[j] );
			r = max( r, j - i + 1 );
			j++;
		}
		else
		{
			myset.erase( s[i] );
			i++;
		}
	}
	
	return r;
}

static int lengthOfLongestSubstringV2(string s) {
	int r = 0;
	unordered_map< char, int > mymap;
	for( int i = 0, j = 0; j < s.size(); j++ )
	{
		if( mymap.find( s[j] ) != mymap.end() )
			i = max(mymap[ s[j] ], i);

		r = max( r, j - i + 1 );
		mymap[ s[j] ] = j + 1;
	}
	
	return r;
}

static int lengthOfLongestSubstring(string s) {
	return lengthOfLongestSubstringV2(s);
}

static void testLengthOfLongestSubstring()
{
	Verify( lengthOfLongestSubstring( "abcabcbb" ), 3, "" );
	Verify( lengthOfLongestSubstring( "bbbbb" ), 1, "" );
	Verify( lengthOfLongestSubstring( "pwwkew" ), 3, "" );
}

#pragma mark - 

namespace Working1 {
	static vector<string> findWordBreak( 	const string& s, unordered_set<string>& wordDict, unordered_map< string, vector<string> >& m )
	{
		if( m.find( s ) != m.end() )
		{
			return m[s];
		}
		
		vector<string> r;
		if( wordDict.count(s) )	 //a whole string is a word
		{
			r.push_back(s);
		}
		
		for( int i = 1; i < s.size(); i++ )
		{
			auto firstPartStr = s.substr( 0, i );
			auto secondPartStr =  s.substr( i, s.size() - i );
			if( wordDict.find( firstPartStr ) != wordDict.end() )
			{
				auto rV = findWordBreak( secondPartStr, wordDict, m );
				for( const auto& w : rV )
				{
					r.push_back( firstPartStr + " " + w );
				}
			}
		}
		m[s] = r;
		return r;
	}
	
	vector<string> wordBreak(string s, vector<string>& wordVec ) {
		unordered_map< string, vector<string> > m;
		unordered_set< string > wordDict;
		for( const auto& w : wordVec )
		{
			wordDict.insert( w );
		}
		
		return findWordBreak( s, wordDict, m );
	}
}

static void findWordBreak( 	const string& s, unordered_set<string>& wordDict, unordered_map< string, vector<string> >& m )
{
	if( m.find( s ) != m.end() )
	{
		return;
	}
	
	if( wordDict.find( s ) != wordDict.end() )
	{
		if( m[s].empty() )
			m[s].push_back( s );
	}
	
	for( int i = 1; i < s.size(); i++ )
	{
		auto firstPartStr = s.substr( 0, i );
			auto secondPartStr =  s.substr( i, s.size() - i );
			if( wordDict.find( firstPartStr ) != wordDict.end() )
			{
				findWordBreak( secondPartStr, wordDict, m );
				for( const auto& v2 : m[secondPartStr] )
				{
					m[s].push_back( firstPartStr + " " + v2 );
				}
			}
	}
}

vector<string> wordBreak(string s, vector<string>& wordVec ) {
	unordered_map< string, vector<string> > m;
	unordered_set< string > wordDict;
	for( const auto& w : wordVec )
	{
		wordDict.insert( w );
	}
	
	findWordBreak( s, wordDict, m );
	return m[s];
}

namespace Other
{
	unordered_map<string, vector<string>> m;

	vector<string> combine(string word, vector<string> prev){
		for(int i=0;i<prev.size();++i){
			prev[i]+=" "+word;
		}
		return prev;
	}

	vector<string> wordBreak(string s, unordered_set<string>& dict) {
		if(m.count(s)) return m[s]; //take from memory
		vector<string> result;
		if(dict.count(s)){ //a whole string is a word
			result.push_back(s);
		}
		for(int i=1;i<s.size();++i){
			string word=s.substr(i);
			if(dict.count(word)){
				string rem=s.substr(0,i);
				vector<string> prev=combine(word,wordBreak(rem,dict));
				result.insert(result.end(),prev.begin(), prev.end());
			}
		}
		m[s]=result; //memorize
		return result;
	}
}

static void testWordBreak()
{
 	// s = "catsanddog"
	// wordDict = ["cat", "cats", "and", "sand", "dog"]
	string s = "aaaaaaa";
	vector< string > wordDict = {"aaaa","aa", "a"};
	string s1 = "catsanddog";
	vector< string > wordDict1 = {"cat", "cats", "and", "sand", "dog"};
	auto r = wordBreak( s, wordDict );
	Verify( r.size(), Working1::wordBreak( s, wordDict ).size(), "" );
	for( const auto& str : r )
		cout << str << endl;
	
}

#pragma mark -
class TrieTree
{
public:
    TrieTree()
    {
        root = make_unique<TrieNode>();
    }
    
    void Insert( const string& word )
    {
        TrieNode* node = root.get();
        for( const auto& c : word )
        {
            int index = c - 'a';
            if( node->children[index] == nullptr )
            {
                node->children[index] = make_unique<TrieNode>();
            }
            
            node = node->children[index].get();
        }
        
        node->isEnd = true;
    }
    
    bool Search( const string& word )
    {
        TrieNode* node = root.get();
        for( const auto& c : word )
        {
            int index = c - 'a';
            if( node->children[index] == nullptr )
                return false;
            
            node = node->children[index].get();
        }
        
        return node && node->isEnd;
    }
    
    bool StartWith( const string& word )
    {
        TrieNode* node = root.get();
        for( const auto& c : word )
        {
            int index = c - 'a';
            if( node->children[index] == nullptr )
                return false;
            
            node = node->children[index].get();
        }
        
        return node;
    }
    
    TrieNode* GetRoot() const { return root.get(); }

private:
    unique_ptr< TrieNode > root;
};

static void dfsFindWord(vector<vector<char>>& board, int i, int j, TrieNode* node, string word, vector<string>& results )
{
    if( board.size() == 0 ) return;
    if( i < 0 || j < 0 || i >= board.size() || j >= board[0].size() )   return;
    
    char c = board[i][j];
    if( c == '*' )          return;

    word = word + c;
    int index = c - 'a';
    TrieNode* next = node->children[index].get();
    if( next != nullptr )
    {
        board[i][j] = '*';
        if( next->isEnd )
        {
            results.push_back( word );
            next->isEnd = false;    // decouple.
        }
        
        dfsFindWord(board, i-1, j, next, word, results);
        dfsFindWord(board, i+1, j, next, word, results);
        dfsFindWord(board, i, j-1, next, word, results);
        dfsFindWord(board, i, j+1, next, word, results);
        
        board[i][j] = c;
    }
}

static vector<string> findWords(vector<vector<char>>& board, vector<string>& words) {
    TrieTree tree;
    for( const auto& word : words )
        tree.Insert( word );

    vector<string> results;
    for( int i = 0; i < board.size(); i++ )
    {
        for( int j = 0; j < board[i].size(); j++ )
        {
            string word;
            dfsFindWord(board, i, j, tree.GetRoot(), word, results );
        }
    }
    
    return results;
}

static void testFindWord()
{
    vector<vector<char>> board ={{'o','a','a','n'},{'e','t','a','e'},{'i','h','k','r'},{'i','f','l','v'}};
    vector<string> words = {"oath","pea","eat","rain"};
    board = {{'a','b'},{'c','d'}};
    words = {"ab","cb","ad","bd","ac","ca","da","bc","db","adcb","dabc","abb","acb"};
    board = {{'a', 'a'}};
    words = {"a"};
    auto results = findWords( board, words );
    for( auto w : results )
        cout << w << " ";
}

#pragma mark - run


void Level1::Run()
{
	testFindWord();
}
