#include "PLSA.h"

void run_EM_PLSA(string mut_file, int experiment);

int main(int argc,char *argv[]){
    if(argc != 3){
        cout << "The number of argument is invalid." << endl;
        return(0);
    }
    string mut_file = argv[1];
    int experiment = atoi(argv[2]);
    run_EM_PLSA(mut_file, experiment);
}

void run_EM_PLSA(string mut_file, int experiment){
    PLSA plsa(mut_file, experiment);
    plsa.load_data();
    plsa.run_EM();
    plsa.write_data();
}
