#include "../tools/getArgs.hpp"

#include "../graph/graph.hpp"

#include "../kClique/ccpath.hpp"
#include "../kClique/pivoterMsk.hpp"

#include <cassert>
#include <string>
#include <iostream>
using std::string;
#include <cassert>

int main(int argc, char * argv[])
{
    argsController * aC = new argsController(argc, argv);
    Graph * g = new Graph();
    g->load(aC->get("-edge"), aC->get("-idx"), atoi(aC->get("-v").c_str()));

    PivoterMsk * pt = new PivoterMsk(g);
    v_size deb = atoi(aC->get("-debug").c_str());
    // if(aC->get("-debug").c_str() != "") deb = ;
    pt->runV(atoi(aC->get("-k").c_str()), deb);

    delete pt;
    delete g;

    return 0;
}
