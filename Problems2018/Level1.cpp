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
    int                 key{ -1 };
    int                 value{ -1 };
};

class DoubleLinkList
{
public:
    DoubleLinkList()
    {
        m_pSentinel = new DoubleLinkListNode;
        m_pTail = m_pSentinel;
    }
    
    DoubleLinkListNode* Push( int key, int value )
    {
        auto pNode = new DoubleLinkListNode;
        pNode->key = key;
        pNode->value = value;
        pNode->next = nullptr;
        pNode->pre = m_pTail;
        
        m_pTail->next = pNode;
        m_pTail = pNode;
        if( m_pSentinel->next == nullptr )
            m_pSentinel->next = pNode;
        
        return pNode;
    }
    
    DoubleLinkListNode* Head() const
    {
        return m_pSentinel->next;
    }
    
    void Pop()
    {
        auto* pHead = Head();
        if( pHead )
        {
            m_pSentinel->next = pHead->next;
            pHead->pre = m_pSentinel;
            delete pHead;
        }
    }

    void MoveToEnd( DoubleLinkListNode* pNode )
    {
        if( pNode->next == nullptr )    // already tail
            return;

        auto* pNodeNext = pNode->next;
        auto* pNodePre = pNode->pre;
        
        pNodePre->next = pNodeNext;
        
        pNode->next = nullptr;
        pNode->pre = m_pTail;
        m_pTail->next = pNode;
        m_pTail = pNode;
    }

private:
    DoubleLinkListNode* m_pTail{ nullptr };
    DoubleLinkListNode* m_pSentinel{ nullptr };
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
                m_doubleLinkList.Pop();
                m_map.erase( pHeadKey );
            }
            
            auto* pNode = m_doubleLinkList.Push( key, value );
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


#pragma mark - run

    
void Level1::Run()
{
    testLRUCache();
}
