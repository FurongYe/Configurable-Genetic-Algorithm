//
//  main.cpp
//  ConfigurableGA
//
//  Created by Furong Ye on 26/10/2020.
//  Copyright © 2020 Furong Ye. All rights reserved.
//

#include "ioha.hpp"
#include "ioh.hpp"

#include <vector>

using namespace std;


int main(int argc, const char *argv[])
{
    string algorithm_name = "fGA-beta1.5";
    string dir = "./";
    int runs = 100;
    int budget = 1000000;
    unsigned seed = 666;

    auto suite = ioh::suite::SuiteRegistry<ioh::problem::IntegerSingleObjective>::instance()
                     .create("PBO", {3}, {1}, {100,500,1000,5000});
    ioha::alg::FastGA fGA(1, 1);
    // fGA.set_beta_f(1.01);
    fGA.run(dir, algorithm_name, suite, budget, std::numeric_limits<int>::max(), runs, seed);

    return 0;
}