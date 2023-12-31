#include "../tools/getArgs.hpp"

#include "../graph/graph.hpp"

#include "../kclique/ccpath.hpp"
#include "../kclique/ccpathParallel.hpp"
#include "../kclique/pivoterMsk.hpp"
#include "../kclique/ccParallel.hpp"

#include <cassert>
#include <string>
#include <iostream>
using std::string;

#include <winbase.h>

int main(int argc, char * argv[])
{
    argsController * aC = new argsController(argc, argv);
    
    v_size deb = 0;
    if(aC->exist("-deb")) deb = atoi(aC->get("-deb").c_str());
    string filePath = "data/skitter/";
    if(aC->exist("-f")) filePath = aC->get("-f");

    printf("\n %s \n",filePath.c_str());

    e_size N = 5000000;

    if(aC->exist("-N")) {
        string tmp = aC->get("-N");
        N = 0;
    }


    double alpha = 1.0;
    if(aC->exist("-a")) alpha = atof(aC->get("-a").c_str());

    v_size k = 10;
    if(aC->exist("-k")) k = atoi(aC->get("-k").c_str());


    if(deb == 0) {
        Graph * g = new Graph();
        v_size n, maxK, tmp;
        double exCnt = 0.0;
        FILE * f = fopen((filePath+"s.txt").c_str(), "r");
        int err = fscanf(f, "%u", &n);
        if(err != 1) {
            printf("s.txt not exist\n");
            return  0;
        }
        if(~fscanf(f, "%u", &maxK)) {
            if(maxK >= k) {
                while(~fscanf(f, "%u-clique: %lf", &tmp, &exCnt)) {
                    if(tmp == k) break;
                }
            }
        }
        printf("%d \n",n);
        printf("\n %s \n",filePath.c_str());
        g->load(filePath+"edge.bin", filePath+"idx.bin", n);
        samplePlusExact * pt = new samplePlusExact(g);

        if(aC->exist("-cc")){
            pt->runCC(k, deb, exCnt, alpha, N);}
        else if(aC->exist("-cccpath")){
            pt->runCCPath(k, deb, exCnt, alpha, N);}
        else{
            pt->run(k, deb, exCnt, alpha, N);}

        delete g;
        delete pt;
    }
    else if(deb == 1) {

    	printf("bass 3 \n");
        int num = 9;
        string name[] = {
            "skitter/","google/","berkstan/","stanford/",
            "amazon0601/","gowalla/","comlj/", "orkut/", 
            "friender/"
        };

        bool ok[] = {
            false,false,false,false,
            false,false,false, false, 
            true
        };


        double ks[500];
        printf("Burasi 1 \n");
        for(int i = 0; i < num; i++) {
            if(!ok[i]) continue;
            Graph * g = new Graph();
            string filePath = "data/";
            filePath += name[i];

            v_size n, maxK, tmp;
            double exCnt = 0.0;

            FILE * f = fopen((filePath+"s.txt").c_str(), "r");
            int err = fscanf(f, "%u", &n);
            if(err != 1) {
                printf("s.txt not exist\n");
                return  0;
            }
            if(~fscanf(f, "%u", &maxK)) {
                while(~fscanf(f, "%u-clique: %lf", &tmp, &exCnt)) {
                    ks[tmp] = exCnt;
                }
            }

            g->load(filePath+"edge.bin", filePath+"idx.bin", n);

            samplePlusExact * pt = new samplePlusExact(g);
            pt->multiRunInit(15, 1.0);
            
            std::cout << filePath << std::endl;

            for(v_size k = 3; k <= 15; k++) {

                for(v_size tt = 5000000; tt <= 4000000*10; tt +=10000000)
                    pt->multiRun(k, 1, ks[k], tt);
            }
            
            delete g;
            delete pt;
        }
    }
    else if(deb == 3) {


    	Graph * g = new Graph();
        v_size n, maxK, tmp;
        double exCnt = 0.0;

        FILE * f = fopen((filePath+"s.txt").c_str(), "r");
        int err = fscanf(f, "%u", &n);
        if(err != 1) {
            printf("s.txt not exist\n");
            return  0;
        }
        
        if(~fscanf(f, "%u", &maxK)) {
            if(maxK >= k) {
                while(~fscanf(f, "%u-clique: %lf", &tmp, &exCnt)) {
                    if(tmp == k) break;
                }
            }
        }
        fclose(f);

        v_size threads = 10;
        if(aC->exist("-p")) threads = atoi(aC->get("-p").c_str());
        
        g->load(filePath+"edge.bin", filePath+"idx.bin", n);

        if(aC->exist("-cc")) {
            ccParallel * pt = new ccParallel(g, threads);
            pt->run(k, deb, exCnt, alpha, N);
            delete pt;
        }
        else {
            samplePlusExactParallel * pt = 
                new samplePlusExactParallel(g, threads);
            pt->run(k, deb, exCnt, alpha, N);
            delete pt;
        }
        delete g;
        
    }

    Sleep(100000);

    return 0;
}
