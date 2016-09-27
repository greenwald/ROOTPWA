
#include <iostream>
#include <string>
#include <stdio.h>

#include "TFhh.h"
#include "TJSS.h"

#include "primeNumbers.h"
#include "reportingUtils.hpp"

using namespace std;

int main(int narg, char* carg[])
{

    if (narg < 4) {
        cout << endl
             << "This program requires 3 input strings for the mother and " << endl
             << "the 2 decay particles, each of the form Jp," << endl
             << "where J is the spin of the particle and p = +/- its parity" << endl;

        return 0;
    }

    string fileName = (narg == 5) ? carg[5] : "primeNumberCache.root";
    if (not rpwa::primeNumbers::instance().readCacheFile(fileName)) {
        printErr << "could not read prime number cache file. Aborting..." << endl;
        return 1;
    }

    int j[3];
    int p[3];

    for (int i = 0; i < 3; ++i) {
        char c;
        sscanf(carg[i + 1], "%d%c", &j[i], &c);
        p[i] = (c == '+') ? 1 : -1;
    }

    cout << "Mother particle:   " << j[0] << parity_to_string(p[0]) << endl;
    cout << "1. decay particle: " << j[1] << parity_to_string(p[1]) << endl;
    cout << "2. decay particle: " << j[2] << parity_to_string(p[2]) << endl;

    cout << j[0] << "," << p[0] << ","
         << j[1] << "," << p[1] << ","
         << j[2] << "," << p[2] << endl;

    TJSS jss(j[0], p[0], j[1], p[1], j[2], p[2]);
    jss.CalcAmpl();


}

