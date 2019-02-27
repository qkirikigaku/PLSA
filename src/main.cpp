#include "PLSA.h"

void run_EM_PLSA(int num_topic, string mut_file, int experiment);

int main(int argc,char *argv[]){
    if(argc != 4){
        cout << "The number of argument is invalid." << endl;
        return(0);
    }
    int num_topic = atoi(argv[1]) + 1;
    string mut_file = argv[2];
    int experiment = atoi(argv[3]);
    run_EM_PLSA(num_topic, mut_file, experiment);
}

void run_EM_PLSA(int num_topic, string mut_file, int experiment){
    PLSA plsa(num_topic, mut_file, experiment);
    plsa.load_data();
    plsa.run_EM();
    plsa.write_data();
}
