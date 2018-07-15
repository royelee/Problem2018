//
//  main.cpp
//  Problems2018
//
//  Created by Roye Li on 6/30/18.
//  Copyright © 2018 Roye Li. All rights reserved.
//

#include <iostream>
#include "Level.h"

template <typename T>
void runLevel() {
	T level;
	level.Run();
}

int main(int argc, const char * argv[]) {
	// insert code here...
	runLevel<Level1>();
	return 0;
}
