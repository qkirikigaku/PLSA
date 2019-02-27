#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <random>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>

using namespace std;

class PLSA {
    private:
        int num_topic; int num_vocab; int num_doc;
        int experiment;
        vector<vector<int> > train_doc;
        int num_sum_word;
        string mut_file; // ex. sample1 | sample2 | real_sample1 | real_sample2
        vector<vector<vector<double> > > qz;
        vector<vector<double> > theta; vector<vector<double> > phi;
        double temp_ll;
        double old_ll;
    public:
        PLSA(int x, string y, int z);
        void run_EM();
        void initialize();
        void Update_qz();
        void Update_parameter();
        void calc_ll();
        void load_data();
        void write_data();
        void show_ll();
        void Normalize(vector<double> &vec);
        double log_sum_exp(vector<double> &vec);
};

PLSA::PLSA(int x, string y, int z){
    num_topic = x; mut_file = y; experiment = z;
}

void PLSA::run_EM(){
    initialize();
    old_ll = -10e15;
    Update_parameter();
    for (int i=0; i < 1000; i++){
        cout << "iter : " << i << endl;
        Update_qz();
        Update_parameter();
        show_ll();
        if(fabs(temp_ll - old_ll) < 1){
            break;
        }
        old_ll = temp_ll;
        cout << endl;
    }
}

void PLSA::initialize(){
    int d,i,k;
    random_device rnd;
    qz.resize(num_doc);
    for (d=0; d < num_doc; d++){
        qz[d].resize(train_doc[d].size());
        for (i=0; i < train_doc[d].size(); i++){
            qz[d][i].resize(num_topic, 0);
            for (k=0; k < num_topic; k++){
                qz[d][i][k] = (double) -rnd();
            }
            Normalize(qz[d][i]);
        }
    }
    theta.resize(num_doc);
    for (d=0; d < num_doc; d++){
        theta[d].resize(num_topic, 0);
    }
    phi.resize(num_topic);
    for (k=0; k < num_topic; k++){
        phi[k].resize(num_vocab, 0);
    }
    num_sum_word = 0;
    for (d=0; d < num_doc; d++){
        num_sum_word += train_doc[d].size();
    }
}

void PLSA::Update_qz(){
    // update qz: responsibility
    int d,i,k;
    for (d=0; d < num_doc; d++){
        for (i=0; i < train_doc[d].size(); i++){
            for (k=0; k < num_topic; k++){
                qz[d][i][k] = theta[d][k] + phi[k][train_doc[d][i]];
            }
            Normalize(qz[d][i]);
        }
    }
}

void PLSA::Update_parameter(){
    int d,i,k,v;
    // update theta: topic distriubtion
    for (d=0; d < num_doc; d++){
        for (k=0; k < num_topic; k++){
            vector<double> q_dk;
            q_dk.resize(train_doc[d].size());
            for (i=0; i < train_doc[d].size(); i++){
                q_dk[i] = qz[d][i][k];
            }
            theta[d][k] = log_sum_exp(q_dk);
        }
        Normalize(theta[d]);
    }
    // update phi: word distribution
    for (k=0; k < num_topic; k++){
        for (v=0; v < num_vocab; v++){
            vector<double> q_kv;
            q_kv.reserve(num_sum_word);
            for (d=0; d < num_doc; d++){
                for (i=0; i < train_doc[d].size(); i++){
                    if(train_doc[d][i] == v){
                        q_kv.push_back(qz[d][i][k]);
                    }
                }
            }
            phi[k][v] = log_sum_exp(q_kv);
        }
        Normalize(phi[k]);
    }
}

void PLSA::calc_ll(){
    int d,i,k,v;
    temp_ll = 0;
    for (d=0; d < num_doc; d++){
        for (i=0; i < train_doc[d].size(); i++){
            vector<double> sum_theta_phi;
            sum_theta_phi.resize(num_topic);
            for (k=0; k < num_topic; k++){
                sum_theta_phi[k] = theta[d][k] + phi[k][train_doc[d][i]];
            }
            temp_ll += log_sum_exp(sum_theta_phi);
        }
    }
}

void PLSA::load_data(){
    ifstream ifs;
    string input_file_name = "data/" + mut_file + ".txt";
    ifs.open(input_file_name.c_str(), ios::in);
    if(!ifs){
        cout << "Cannot open " + input_file_name << endl;
        exit(1);
    }
    char buf[1000000];
    char *temp;
    vector<vector<int> > raw_document;
    vector<int> words_number;
    ifs.getline(buf, 1000000);
    temp = strtok(buf, " ");
    num_doc = atoi(temp);
    raw_document.resize(num_doc);
    words_number.resize(num_doc, 0);
    train_doc.resize(num_doc);
    temp = strtok(NULL, " "); num_vocab = atoi(temp);
    int temp_word_number;
    for (int d=0; d < num_doc; d++){
        ifs.getline(buf, 1000000);
        for (int v=0; v < num_vocab; v++){
            if(v == 0) temp_word_number = atoi(strtok(buf, " "));
            else temp_word_number = atoi(strtok(NULL, " "));
            for (int i=0; i < temp_word_number; i++){
                raw_document[d].push_back(v);
                words_number[d]++;
            }
        }
    }
    for (int d=0; d < num_doc; d++){
        int count = 0;
        train_doc[d].resize(words_number[d]);
        for (int i=0; i < words_number[d]; i++){
            train_doc[d][i] = raw_document[d][i];
            count ++;
        }
    }
    ifs.close();
}

void PLSA::write_data(){
    ofstream ofs;
    string output_file_name = "result/" + mut_file + "_" +
        to_string(experiment) + "/result_k";
    if(num_topic < 10){
        output_file_name += "0" + to_string(num_topic) + ".txt";
    }
    else{
        output_file_name += to_string(num_topic) + ".txt";
    }
    ofs.open(output_file_name, ios::out);
    int num_param = num_doc * num_topic - 1 +
                num_topic * num_vocab - 1;
    double BIC = 2 * temp_ll - num_param * log(num_sum_word);
    ofs << to_string(BIC) << "\n";
    ofs << to_string(temp_ll) << "\n";

    int d, i, k, v;
	for (k=0; k < num_topic; k++){
        for (v=0; v < num_vocab; v++){
            ofs << to_string(exp(phi[k][v])) << " ";
        }
        ofs << "\n";
    }
    for (d=0; d < num_doc; d++){
        for (k=0; k < num_topic; k++){
            ofs << to_string(exp(theta[d][k])) << " ";
        }
        ofs << "\n";
    }
    ofs.close();
}

void PLSA::show_ll(){
    calc_ll();
    cout << "   Log Likelihood: " << temp_ll << endl;
    cout << "   Improvement point: " << temp_ll - old_ll << endl;
}

void PLSA::Normalize(vector<double> &vec){
    double sum = log_sum_exp(vec);
    int length = vec.size();
    for (int i=0; i < length; i++){
        vec[i] -= sum;
    }
}

double PLSA::log_sum_exp(vector<double> &vec){
    int length = vec.size();
    double max_iter = 0;
    for (int i=1; i < length; i++){
        if(vec[i] > vec[max_iter]) max_iter = i;
    }
    double sum = 0;
    for (int i=0; i < length; i++){
        sum += exp(vec[i] -vec[max_iter]);
    }
    double return_value = vec[max_iter] + log(sum);
    return(return_value);
}

void run_EM_PLSA(int num_topic, string cancer_type, int experiment);
