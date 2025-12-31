//
//  main.cpp
//  bacterial_speciation
//
//  Created by Haruna Nakamura on 2024/01/16.
//

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <random>
#include <algorithm>
#include <unistd.h>
#include <map>
#include <set>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <utility>
#include <sstream>

using namespace std;

struct ind_t
{
    int pop;
    gsl_vector *hap_inc; // incompatibility loci
    gsl_vector *hap_neu; // neutral loci
    double w; // fitness
    int hap_number; // haplotype number in pop_hap_t
    int new_hap_number; //
};

struct pop_hap_t
{
    int pop;
    string hap; //0000000000, 0000010000 etc.
    double freq;
    
    int birth_gen;
    vector <int> fix_start_gen;
    vector <int> fix_end_gen;
    int lost_gen;
    
};

ind_t *pop;
ind_t *oldpop;
pop_hap_t *incompatibility;
pop_hap_t *old_incompatibility;

random_device rd;
mt19937 mt(rd());
uniform_real_distribution<> dist(0.0, 1.0);

gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

/* default setting */
int gen = 1;
int final_gen = 1e4;

int pop_size = 1000;
int pop_number = 10;
map <int, int> pop_data;
int total_ind;
int maximum_ind;
double u = 1.0e-6; // mutation rate
map <int, float> u_data;
double m = 0.0001; // migration rate
double p = 0.2; // proportion of recombinants
double r = 1.0; // mean number of EM loci in recombinants
double s = 0.1; // effect of incompatibility loci
double q = 0.0; // effect of mutation
double thres = 0.95; // threshold frequency for fixed loci
int l_inc = 10; // number of loci in restrict-modification system
int l_neu = 0; // number of neutral sites
int maximum_haplotype_number = pow(2, l_inc);
int nHap; // number of haplotypes at current generation
vector <vector <double>> size_change;
vector <double> size_change_t;
vector <vector <double>> mutation_rate_change;
vector <double> mutation_rate_change_t;

gsl_vector *templete;
gsl_vector *diff;


double e_time;
double e_pop;
double e_value;
size_t nind;
size_t n;

string output_folder_path = "/Users/hnakamura/research_bacterial_speciation/result";
string prerun_result = "";
string hoge;
string pop_file_path = "";
string mut_file_path = "";
string haplo_file_path = "";

// life-cycle stage
extern void selection();
extern void migration();
extern void reproduction();
extern void recombinants(int, int, int);
extern void mutation();
extern void input_param();
extern void input_prerun();
extern void check_haplotype();
extern void record_haplotype();
extern void record_fixed_haplotype(int);
extern void record_final_gen();
extern void record_param();

// function
vector <string> commonElement(const vector<string>, const vector<string>);
vector <string> diffElement(const vector<string>, const vector<string>);

int main(int argc, const char * argv[]) {
    
    /* initial condition */
    
    int a=0;
    while (a<argc){
        if (strcmp(argv[a], "-N") == 0){ // population size
            pop_size = atoi(argv[a+1]);
        }
        if (strcmp(argv[a], "-mut") == 0){ // mutation rate
            u = atof(argv[a+1]);
        }
        if (strcmp(argv[a], "-mig") == 0){ // migration rate
            m = atof(argv[a+1]);
        }
        if (strcmp(argv[a], "-rec_p") == 0){ // proportion of recombinants
            p = atof(argv[a+1]);
        }
        if (strcmp(argv[a], "-rec_n") == 0){ // recombinant size
            r = atof(argv[a+1]);
        }
        if (strcmp(argv[a], "-s") == 0){ // effects of incompatibility between parents
            s = atof(argv[a+1]);
        }
        if (strcmp(argv[a], "-q") == 0){ // effects of mutations on survival rate
            q = atof(argv[a+1]);
        }
        if (strcmp(argv[a], "-nl_inc") == 0){
            l_inc = atoi(argv[a+1]);
        }
        if (strcmp(argv[a], "-nl_neu") == 0){
            l_neu = atoi(argv[a+1]);
        }
        if (strcmp(argv[a], "-o") == 0){
            output_folder_path = argv[a+1];
        }
        if (strcmp(argv[a], "-Fi") == 0){
            prerun_result = argv[a+1];
        }
        if (strcmp(argv[a], "-final_gen") == 0){
            final_gen = atoi(argv[a+1]);
        }
        if (strcmp(argv[a], "-thres") == 0){
            thres = atof(argv[a+1]);
        }
        if (strcmp(argv[a], "-eN") == 0){ // change in population size
            size_change_t.clear();
            
            if (string(argv[a+1])=="ALL"){
                e_pop = -1;
            }else{
                e_pop = atoi(argv[a+1]);
            }
            e_value = atof(argv[a+2]);
            e_time = atoi(argv[a+3]);
            
            
            
            size_change_t.push_back(e_pop);
            size_change_t.push_back(e_value);
            size_change_t.push_back(e_time);
            
            size_change.push_back(size_change_t);
            
        }
        
        
        if (strcmp(argv[a], "-eU") == 0){ // change in mutation rate
            mutation_rate_change_t.clear();
            
            if (string(argv[a+1])=="ALL"){
                e_pop = -1;
            }else{
                e_pop = atoi(argv[a+1]);
            }
            
            e_value = atof(argv[a+2]);
            e_time = atoi(argv[a+3]);
            
            
            mutation_rate_change_t.push_back(e_pop);
            mutation_rate_change_t.push_back(e_value);
            mutation_rate_change_t.push_back(e_time);
            
            mutation_rate_change.push_back(mutation_rate_change_t);
            
        }
        
        if (strcmp(argv[a], "-Fp") == 0){ // pop file
            pop_file_path = argv[a+1];
        }
        if (strcmp(argv[a], "-Fu") == 0){ // mutation rate file
            mut_file_path = argv[a+1];
        }
        if (strcmp(argv[a], "-Fh") == 0){ // haplotype file
            haplo_file_path = argv[a+1];
        }
        
        ++a;
    }
    
    total_ind = 0;
    
    maximum_ind = 100000;
    
    /* paramter change */
    if (size_change.size()>0){
        size_change_t.clear();
        for (size_t e=0, n=size_change.size(); e<n; ++e){
            size_change_t.push_back(size_change[e][2]);
        }
    }
    
    if (mutation_rate_change.size()>0){
        mutation_rate_change_t.clear();
        for (size_t e=0, n=mutation_rate_change.size(); e<n; ++e){
            mutation_rate_change_t.push_back(mutation_rate_change[e][2]);
        }
    }
    
    input_param();
    
    
    pop = new ind_t[maximum_ind];
    oldpop = new ind_t[maximum_ind];
    incompatibility = new pop_hap_t[maximum_haplotype_number*pop_number];
    old_incompatibility = new pop_hap_t[maximum_haplotype_number*pop_number];
    
    templete = gsl_vector_alloc(l_inc);
    diff = gsl_vector_alloc(l_inc);
    
    
    
    for (int i=0; i<maximum_ind; ++i){
        pop[i].hap_inc = gsl_vector_calloc(l_inc);
        pop[i].hap_neu = gsl_vector_calloc(l_neu);
        oldpop[i].hap_inc = gsl_vector_calloc(l_inc);
        oldpop[i].hap_neu = gsl_vector_calloc(l_neu);
    }
    
    
    if (prerun_result==""){
        total_ind = 0;
        
        
        for (int k=0; k<pop_number; ++k){
            for (int i=total_ind; i<total_ind+pop_data[p]; ++i){
                pop[i].pop = k;
                pop[i].w = 1.0;
                pop[i].hap_number = k;
            }
            total_ind += pop_data[p];
        }
        
        for (int k=0; k<pop_number; ++k){
            incompatibility[k].pop = k;
            
            hoge = "";
            for (int l=0; l<l_inc; ++l){
                hoge += "0";
            }
            incompatibility[k].hap = hoge;
            incompatibility[k].birth_gen = gen;
            incompatibility[k].fix_start_gen.push_back(gen);
            incompatibility[k].lost_gen = -1;
            
            ++nHap;
        }
        
        
    }else{
        
        input_prerun();
        
    }
    
    
    
    
    /* main execution */
    while (gen<final_gen+1){
        
        // change of population size
        auto it = find(size_change_t.begin(), size_change_t.end(), gen);
        if (it!=size_change_t.end()){
            //cout << "change population size to ";
            if (size_change[it-size_change_t.begin()][0]==-1){
                for (int k=0; k<pop_number; ++k){
                    pop_data[k] = size_change[it-size_change_t.begin()][1];
                }
                //cout << size_change[it-size_change_t.begin()][1] << " in all populations";
            
            }else{
                int k = size_change[it-size_change_t.begin()][0];
                pop_data[k] = size_change[it-size_change_t.begin()][1];
                //cout << pop_data[k] << " in pop " << k;
            }
            //cout << " at generation " << gen << endl;
        }
            
        // change of mutation rate
        auto it2 = find(mutation_rate_change_t.begin(), mutation_rate_change_t.end(), gen);
        if (it2!=mutation_rate_change_t.end()){ //
            if (mutation_rate_change[it2-mutation_rate_change_t.begin()][0]==-1){
                for (int k=0; k<pop_number; ++k){
                    u_data[k] = mutation_rate_change[it2-mutation_rate_change_t.begin()][1];
                }
            }else{
                int k = mutation_rate_change[it2-mutation_rate_change_t.begin()][0];
                u_data[k] = mutation_rate_change[it2-mutation_rate_change_t.begin()][1];
            }
            
        }
        
        
        
        migration();
        reproduction();
        mutation();
    
        
        check_haplotype();
        
        
        /*if (gen%int(pow(10, floor(log10(gen))))==0 && gen>=100){
            record_haplotype();
            
        }*/
        
        if (gen>=1000 && gen%pop_size==0){
            record_haplotype();
        }
        
        if (gen==final_gen){
            record_param();
            record_final_gen();
        }
        //cout << gen << endl;
        ++gen;
        
    }
    
    
    delete [] pop;
    delete [] incompatibility;
    delete [] old_incompatibility;
    
    return 0;
    
}



void migration(){
    
    /* stepping-stone model */
    int old_pop_number;
    int new_pop_number;
    vector <int> cand_pop_list;
    
    cand_pop_list.reserve(2);
    
    for (int i=0; i<total_ind; ++i){
        if (dist(mt)<m){
            old_pop_number = pop[i].pop;
            cand_pop_list.clear();
            
            if (old_pop_number==0){
                cand_pop_list.push_back(old_pop_number+1);
            }else if (old_pop_number==pop_number-1){
                cand_pop_list.push_back(old_pop_number-1);
            }else{
                cand_pop_list.push_back(old_pop_number-1);
                cand_pop_list.push_back(old_pop_number+1);
            }
            
            shuffle(cand_pop_list.begin(), cand_pop_list.end(), mt);
            new_pop_number = cand_pop_list.back();
            
            // migration
            pop[i].pop = new_pop_number;
        }
    }
}


void reproduction(){
    int recipient;
    int donor;
    vector <int> ind_list;
    
    ind_list.reserve(maximum_ind);
    
    for (int i=0; i<total_ind; ++i){
        oldpop[i].pop = pop[i].pop;
        oldpop[i].w = pop[i].w;
        gsl_vector_memcpy(oldpop[i].hap_inc, pop[i].hap_inc);
        gsl_vector_memcpy(oldpop[i].hap_neu, pop[i].hap_neu);
        oldpop[i].hap_number = pop[i].hap_number;
        oldpop[i].new_hap_number = pop[i].new_hap_number;
    }
    
    nind = 0;
    for (int k=0; k<pop_number; ++k){
        
        ind_list.clear();
        
        for (int j=0; j<total_ind; ++j){
            if (oldpop[j].pop==k){
                
                // death based on fitness
                if (dist(mt)<oldpop[j].w){
                    ind_list.push_back(j);
                }
            }
        }
        
        
        n = ind_list.size();
        for (size_t i=nind; i<nind+pop_data[k]; ++i){
            pop[i].pop = k;
            pop[i].w = 1.0;
            
            /* recombinant */
            if (dist(mt)<p){
                
                // recombination
                recipient = ind_list[rand() % n];
                donor = ind_list[rand() % n];
                
                recombinants(recipient, donor, (int)i);
                
            }else{
                recipient = ind_list[rand() % n];
                
                pop[i].hap_number = oldpop[recipient].hap_number;
                pop[i].new_hap_number = oldpop[recipient].new_hap_number;
                gsl_vector_memcpy(pop[i].hap_inc, oldpop[recipient].hap_inc);
                gsl_vector_memcpy(pop[i].hap_neu, oldpop[recipient].hap_neu);
                
            }
        }
        
        nind += pop_data[k];
    }
    
    total_ind = (int)nind;
    
}

// Mutations reduce survival by q. In recombinant individuals, additional mutations further decrease survival to x(1 − q).

void mutation(){
    double theta;
    int num_mutation;
    int dnm_loci;
    int dnm_ind;
    vector <int> ind_list;
    
    ind_list.reserve(maximum_ind);
    //uniform_int_distribution<int> mut_ind_dist(0, total_ind-1);
    uniform_int_distribution<int> mut_EM_dist(0, l_inc-1);
    
    for (int k=0; k<pop_number; ++k){
        ind_list.clear();
        
        theta = pop_data[k] * (double)u_data[k] * l_inc;
        
        poisson_distribution<int> mut_dist(theta);
        num_mutation = mut_dist(mt);
        
        for (int i=0; i<total_ind; ++i){
            if (pop[i].pop == k){
                ind_list.push_back(i);
            }
        }
        
        n = ind_list.size();
        if (num_mutation>0){
            //cout << k << "," << pop_data[k] << "," << u_data[k] << ":" << theta << "-" << num_mutation << endl;
            for (int l=0; l<num_mutation; ++l){
                dnm_loci = mut_EM_dist(mt);
                dnm_ind = ind_list[rand() % n];
                
                if (gsl_vector_get(pop[dnm_ind].hap_inc, dnm_loci)==0){
                    gsl_vector_set(pop[dnm_ind].hap_inc, dnm_loci, 1);
                }else{
                    gsl_vector_set(pop[dnm_ind].hap_inc, dnm_loci, 0);
                }
                
                pop[dnm_ind].w = pop[dnm_ind].w*(1.-q);
            }
        }
        
    }
    
}


void recombinants(int recipient, int donor, int new_ind_number){
    poisson_distribution<int> recombinant_EM(r);
    uniform_int_distribution<int> mut_EM_dist(0, 2*l_inc-1);
    
    int rec_n;
    int rr = 0;
    int rec_pos;
    double d; // number of different alleles between donor and recipient
    double survival_rate; // fitness of recombinants
    vector <int> recombinants_list;
    gsl_vector *r_recipient;
    gsl_vector *r_donor;
        
    gsl_vector_memcpy(pop[new_ind_number].hap_inc, oldpop[recipient].hap_inc);
    gsl_vector_memcpy(pop[new_ind_number].hap_neu, oldpop[recipient].hap_neu);
    
    while (rr==0){
        rr = recombinant_EM(mt);
    }
    rec_n = rr;
    rec_pos = mut_EM_dist(mt);
    
    recombinants_list.clear();
    if (rec_pos%2==0){
        for (int i=0; i<rec_n; ++i){
            if (rec_pos/2+i < l_inc){
                recombinants_list.push_back(rec_pos/2 + i);
            }else{
                recombinants_list.push_back(rec_pos/2 + i - l_inc);
            }
        }
    }else{
        for (int i=0; i<rec_n; ++i){
            if (rec_pos/2-i >= 0){
                recombinants_list.push_back(rec_pos/2-i);
            }else{
                recombinants_list.push_back(l_inc+rec_pos/2-i);
            }
        }
    }
        
    if (rec_n>l_inc){ //Genotypes are replaced by the donor when recombined loci outnumber incompatibility loci.
        gsl_vector_memcpy(pop[new_ind_number].hap_inc, oldpop[donor].hap_inc);
    }else if(rec_n==0){ //When zero loci are recombined, replacement defaults to the recipient (this case is already accounted for by the recombination rule).
        gsl_vector_memcpy(pop[new_ind_number].hap_inc, oldpop[recipient].hap_inc);
            
    }else{
        for (size_t i=0, n=recombinants_list.size(); i<n; ++i){
            gsl_vector_set(pop[new_ind_number].hap_inc, recombinants_list[(int)i], gsl_vector_get(oldpop[donor].hap_inc, recombinants_list[(int)i]));
        }
    }
    
    
    
}


void input_param(){
    int np;
    
    vector <vector <string>> data;
    string data_column;
    vector <string> data_row;
    data_row.reserve(pop_number);
    
    string line;
    int row;
    int column;
    
    // input of population size
    if (pop_file_path==""){
        
        for (int k=0; k<pop_number; ++k){
            pop_data[k] = pop_size;
        }
        
    }else{
        
        ifstream ifs(pop_file_path); // read file
        
        if (!ifs.is_open()){
            cerr << "Error! Unable to open the pop size file." << endl;
        }
        
        row = 0;
        while(getline(ifs, line)){
            istringstream pop_line(line);
            data_row.clear();
            
            if (row>0){
                while(getline(pop_line, data_column, '\t')){
                    data_row.push_back(data_column);
                
                }
                data.push_back(data_row);
            }
            ++row;
        }
        ifs.close();
        
        
        for (const auto & data_row: data){
            column = 0;
            
            for (const auto & data_column: data_row){
                if (column==1){
                    np = atoi(data_column.c_str());
                }else if (column==2){
                    pop_data[np] = atoi(data_column.c_str());
                }
                
                ++column;
            }
        }
    }
    
    
    // input mutation rate
    if (mut_file_path==""){
        
        for (int k=0; k<pop_number; ++k){
            u_data[k] = u;
        }
        
    }else{
        
        ifstream ifs2(mut_file_path); // read file
        
        if (!ifs2.is_open()){
            cerr << "Error! Unable to open the pop size file." << endl;
        }

        row = 0;
        while(getline(ifs2, line)){
            istringstream mut_line(line);
            data_row.clear();
            
            if (row>0){
                while(getline(mut_line, data_column, '\t')){
                    data_row.push_back(data_column);
                }
                data.push_back(data_row);
            }
            ++row;
        }
        ifs2.close();
        
        for (const auto & data_row: data){
            column = 0;
            np = -1;
            for (const auto & data_column: data_row){
                if (column==1){
                    np = atoi(data_column.c_str());
                }else if (column==2){
                    u_data[np] = atof(data_column.c_str());
                }
                
                ++column;
            }
        }
    }
    
    
    
}



void input_prerun(){
    ifstream ifs(prerun_result); // read file
    vector <vector <string>> data;
    string data_column;
    vector <string> data_row;
    data_row.reserve(maximum_ind);
    
    string line;
    int row;
    int column;
    string hap;
    int alle;
    
    if (!ifs.is_open()){
        cerr << "Error! Unable to open the prerun file." << endl;
    }
    
    
    while(getline(ifs, line)){
        istringstream d_line(line);
        data_row.clear();
        
        while(getline(d_line, data_column, '\t')){
            data_row.push_back(data_column);
        }
        data.push_back(data_row);
    }
    ifs.close();
    
    row = -1;
    for (const auto& data_row: data){
        column = 0;
        
        for (const auto& data_column : data_row){
            
            if(row>-1){
                if (column==0){ // gen
                    gen = atoi(data_column.c_str());
                }else if (column==1){ // pop number
                    pop[row].pop = atoi(data_column.c_str());
                }else if (column==2){ // haplotype
                    hap = data_column;
                    for (int k=0; k<l_inc; ++k){
                        alle = hap[k]-'0';
                        
                        gsl_vector_set(pop[row].hap_inc, k, alle);
                        
                    }
            
                }else if (column==3){ // fitness
                    pop[row].w = atof(data_column.c_str());
                }
            }
            ++column;
        }
        ++row;
        
    }
    total_ind = row;
    
    // allocation haplotypes
    vector <string> haplotype_list;
    vector <int> ind_list;
    
    haplotype_list.reserve(maximum_ind);
    ind_list.reserve(maximum_ind);
    map <string, size_t> hap_info;
    
    size_t c;
    double freq;
    
    for (int k=0; k<pop_number; ++k){
        ind_list.clear();
        
        for (int i=0; i<row; ++i){
            if (pop[i].pop==k){
                ind_list.push_back(i);
            }
        }
        nind = ind_list.size();
        
        for (int j : ind_list){
            hap = "";
            for (int l=0; l<l_inc; ++l){
                hap += to_string((int)gsl_vector_get(pop[j].hap_inc, l));
            }
            haplotype_list.push_back(hap);
        }
        
        string ind_hap;
        while (haplotype_list.size()>0){
            hap = haplotype_list[0];
            c = std::count(haplotype_list.begin(), haplotype_list.end(), hap);
            hap_info[hap] = c;
            haplotype_list.erase(remove(haplotype_list.begin(), haplotype_list.end(), hap), haplotype_list.end());
            
            freq = c/nind;
            
            for (int j:ind_list){
                ind_hap = "";
                for (int l=0; l<l_inc; ++l){
                    ind_hap += to_string((int)gsl_vector_get(pop[j].hap_inc, l));
                }
                
                if (ind_hap==hap){
                    pop[j].hap_number = nHap;
                }
            }
            
            incompatibility[nHap].pop = k;
            incompatibility[nHap].hap = hap;
            incompatibility[nHap].freq = freq;
            ++nHap;
        }
        
        
    }
    
    
    
    ifstream ifs2(haplo_file_path); // read file
    data.clear();
    data_row.clear();
    data_row.reserve(maximum_haplotype_number);
    
    if (!ifs2.is_open()){
        cerr << "Error! Unable to open the haplotype file." << endl;
    }
    
    int np;
    
    while(getline(ifs2, line)){
        istringstream d_line(line);
        data_row.clear();
        
        while(getline(d_line, data_column, '\t')){
            data_row.push_back(data_column);
        }
        data.push_back(data_row);
    }
    ifs2.close();
    
    int hap_number;
    string token;
    int start = 0;
    int end = 0;
    
    row = -1;
    for (const auto& data_row: data){
        column = 0;
        hap_number = 0;
        np = -1;
        
        for (const auto& data_column : data_row){
            if(row>-1){
                if (column==0){ // pop number
                    np = atoi(data_column.c_str());
                }else if (column==1){ // haplotype
                    hap = data_column.c_str();
                    
                    for (int l=0; l<nHap; ++l){
                        if (incompatibility[l].pop == np && incompatibility[l].hap== hap){
                            hap_number = l;
                            break;
                        }
                    }
                    
                }else if (column==2){ // birth_gen
                    incompatibility[hap_number].birth_gen = atoi(data_column.c_str());
                    
                }else if (column==3){ // fix_gen
                    incompatibility[hap_number].fix_start_gen.clear();
                    incompatibility[hap_number].fix_end_gen.clear();
                    
                    istringstream iss(data_column.c_str());
                    while (getline(iss, token, ';')){
                        size_t pos = token.find('-');
                        
                        if (pos != string::npos){
                            start = stoi(token.substr(0, pos));
                            
                            if (token.size() == pos+1){
                                ;
                            }else{
                                end = stoi(token.substr(pos+1));
                            }
                            
                            incompatibility[hap_number].fix_start_gen.push_back(start);
                            incompatibility[hap_number].fix_end_gen.push_back(end);
                                                        
                        }
                        
                    }
                    
                }else if (column==4){ // lost_gen
                    incompatibility[hap_number].lost_gen = atoi(data_column.c_str());
                    
                }
            }
            ++column;
        }
        ++row;
        
    }
    
}


void check_haplotype(){ // calculate haplotype frequency at end stage of a cycle. If a haplotype is lost within a population, remove the information from incompatibility loci data.
    
    int new_nHap=0;
    
    vector <int> polymorphic_loci_list;
    polymorphic_loci_list.reserve(nHap);
    polymorphic_loci_list.clear();
    
    for (int l=0; l<nHap; ++l){
        old_incompatibility[l].pop = incompatibility[l].pop;
        old_incompatibility[l].hap = incompatibility[l].hap;
        old_incompatibility[l].freq = incompatibility[l].freq;
        old_incompatibility[l].birth_gen = incompatibility[l].birth_gen;
        old_incompatibility[l].fix_start_gen = incompatibility[l].fix_start_gen;
        old_incompatibility[l].fix_end_gen = incompatibility[l].fix_end_gen;
        old_incompatibility[l].lost_gen = incompatibility[l].lost_gen;
    }
    
    int old_hap_number;
    int alle;
    string ind_hap;
    double freq;
    vector <string> haplotype_old_list;
    haplotype_old_list.reserve(maximum_haplotype_number);
    
    vector <string> haplotype_list;
    haplotype_list.reserve(maximum_haplotype_number);
    
    vector <string> denovo_list; // only new haplotype list
    vector <string> lost_list; // only lost haplotype list
    vector <string> polymorphic_list; //
    
    vector <int> ind_list;
    ind_list.reserve(maximum_ind);
    
    
    int score;
    
    for (int k=0; k<pop_number; ++k){
        // Haplotypes present in each subpopulation
        // 1) haplotypes in a former generation
        haplotype_old_list.clear();
        for (int l=0; l<nHap; ++l){
            if (old_incompatibility[l].pop == k){
                haplotype_old_list.push_back(old_incompatibility[l].hap);
            }
        }
        
        // 2) haplotype in a current generation
        haplotype_list.clear();
        ind_list.clear();
        
        for (int i=0; i<total_ind; ++i){
            if (pop[i].pop == k){
                ind_hap = "";
                for (int n=0; n<l_inc; ++n){
                    ind_hap += to_string((int)gsl_vector_get(pop[i].hap_inc, n));
                }
                haplotype_list.push_back(ind_hap);
                
                ind_list.push_back(i);
            }
        }
        
        set <string> unique_haplotype(haplotype_list.begin(), haplotype_list.end());
        vector <string> haplotype_new_list(unique_haplotype.begin(), unique_haplotype.end());
        
        
        denovo_list.clear();
        lost_list.clear();
        polymorphic_list.clear();
        
        denovo_list = diffElement(haplotype_new_list, haplotype_old_list);
        lost_list = diffElement(haplotype_old_list, haplotype_new_list);
        polymorphic_list = commonElement(haplotype_new_list, haplotype_old_list);
        
        nind = ind_list.size();
        
        for (string hap : polymorphic_list){
            old_hap_number=-1;
            
            for (int l=0; l<nHap; ++l){
                if (old_incompatibility[l].pop==k && old_incompatibility[l].hap==hap){
                    old_hap_number = l;
                    break;
                }
            }
            
            for (int n=0; n<l_inc; ++n){
                alle = hap[n]-'0';
                gsl_vector_set(templete, n, alle);
            }
            
            freq = 0;
            for (int i: ind_list){
                gsl_vector_memcpy(diff, pop[i].hap_inc);
                gsl_vector_sub(diff, templete);
                score = 0;
                
                for (int l=0; l<l_inc; ++l){
                    score += fabs(gsl_vector_get(diff, l));
                }
                
                if (score==0){
                    pop[i].new_hap_number = new_nHap;
                    ++freq;
                }
            }
            freq = freq/nind;
            
            //cout << "poly: " << hap << " " << k << " " << freq << endl;
            
            incompatibility[new_nHap].pop = k;
            incompatibility[new_nHap].hap = hap;
            incompatibility[new_nHap].freq = freq;
            
            incompatibility[new_nHap].birth_gen = old_incompatibility[old_hap_number].birth_gen;
            incompatibility[new_nHap].fix_start_gen = old_incompatibility[old_hap_number].fix_start_gen;
            incompatibility[new_nHap].fix_end_gen = old_incompatibility[old_hap_number].fix_end_gen;
            incompatibility[new_nHap].lost_gen = old_incompatibility[old_hap_number].lost_gen;
            
            if (freq>=thres){ // fixation
                if (old_incompatibility[old_hap_number].fix_start_gen.size()==old_incompatibility[old_hap_number].fix_end_gen.size()){
                    incompatibility[new_nHap].fix_start_gen.push_back(gen);
                }
            }else{
                if (old_incompatibility[old_hap_number].fix_start_gen.size()>old_incompatibility[old_hap_number].fix_end_gen.size()){
                    incompatibility[new_nHap].fix_end_gen.push_back(gen);
                }
            }
            ++new_nHap;
        }
            
        
        for (string hap : denovo_list){
            
            for (int n=0; n<l_inc; ++n){
                alle = hap[n]-'0';
                gsl_vector_set(templete, n, alle);
            }
            
            freq = 0;
            
            for (int i: ind_list){
                gsl_vector_memcpy(diff, pop[i].hap_inc);
                gsl_vector_sub(diff, templete);
                
                score = 0;
                for (int n=0; n<l_inc; ++n){
                    score += fabs(gsl_vector_get(diff, n));
                }
                    
                if (score==0){
                    pop[i].new_hap_number = new_nHap;
                    ++freq;
                }
            }
            freq = freq/nind;
            
            //cout << "new: " << hap << " " << p << " " << freq << endl;
            incompatibility[new_nHap].pop = k;
            incompatibility[new_nHap].hap = hap;
            incompatibility[new_nHap].freq = freq;
            
            incompatibility[new_nHap].birth_gen = gen;
            incompatibility[new_nHap].fix_start_gen.clear();
            incompatibility[new_nHap].fix_end_gen.clear();
            incompatibility[new_nHap].lost_gen = -1;
            
            
            ++new_nHap;
        }
        
        
        
        
        for (string hap : lost_list){
            old_hap_number=-1;
            
            for (int l=0; l<nHap; ++l){
                if (old_incompatibility[l].pop==k && old_incompatibility[l].hap==hap){
                    old_hap_number = l;
                    break;
                }
            }
            
            old_incompatibility[old_hap_number].lost_gen = gen;
            
            if (old_incompatibility[old_hap_number].fix_start_gen.size() > old_incompatibility[old_hap_number].fix_end_gen.size()){
                old_incompatibility[old_hap_number].fix_end_gen.push_back(gen);
            }
            
            // Haplotypes that reached fixation at least once
            if (old_incompatibility[old_hap_number].fix_start_gen.size()>0){
                record_fixed_haplotype(old_hap_number);
            }
            
            //cout << "lost: " << hap << " " << k << " " << endl;
        }
        
    }
    
    
    for (int i=0; i<total_ind; ++i){
        pop[i].hap_number = pop[i].new_hap_number;
    }
    
    
    nHap = new_nHap;
}



void record_haplotype(){
    
    // record
    const int bufsize=1000000;
    char buf[bufsize];
    ofstream ofs;
    ofs.rdbuf()->pubsetbuf(buf,bufsize);

    stringstream FileName;
    FileName << output_folder_path << "/haplotype.tsv";
    string filename(FileName.str());
    
    ifstream fin(filename);
    if (!fin){
        ofs.open(filename, ios_base::trunc | ios_base::in);
        // columns
        ofs << "gen\tpop\thaplotype\tcount\n";
    }else{
        ofs.open(filename, ios_base::app | ios_base::in);
    }
    
    vector <int> haplotype_list;
    haplotype_list.reserve(maximum_haplotype_number);
    
    vector <int> ind_list;
    ind_list.reserve(maximum_ind);
    
    string hap;
    double count;
    
    for (int k=0; k<pop_number; ++k){
        haplotype_list.clear();
        ind_list.clear();
        
        for (int l=0; l<nHap; ++l){
            if (incompatibility[l].pop==k){
                haplotype_list.push_back(l);
            }
        }
        
        for (int i=0; i<total_ind; ++i){
            if (pop[i].pop==k){
                ind_list.push_back(i);
            }
        }
        
        for (int l:haplotype_list){
            hap = incompatibility[l].hap;
            count = 0;
            
            for (int j:ind_list){
                if (pop[j].hap_number==l){
                    ++count;
                }
            }
            
            ofs << gen << "\t" << k << "\t" << hap << "\t" << count << endl;
        }
        
    }
    
    
}

void record_fixed_haplotype(int old_hap_number){
    
    // record
    const int bufsize=1000000;
    char buf[bufsize];
    ofstream ofs;
    ofs.rdbuf()->pubsetbuf(buf,bufsize);

    stringstream FileName;
    FileName << output_folder_path << "/fixed_haplotype.tsv";
    string filename(FileName.str());
    
    ifstream fin(filename);
    if (!fin){
        ofs.open(filename, ios_base::trunc | ios_base::in);
        // columns
        ofs << "pop\thaplotype\tbirth_gen\tfix_gen\tlost_gen\n";
    }else{
        ofs.open(filename, ios_base::app | ios_base::in);
    }
    
    ofs << old_incompatibility[old_hap_number].pop << "\t" << old_incompatibility[old_hap_number].hap << "\t" << old_incompatibility[old_hap_number].birth_gen;
    
    string fix_gen;
    for (size_t j=0, n=old_incompatibility[old_hap_number].fix_start_gen.size(); j<n; ++j){
        fix_gen = fix_gen + to_string(old_incompatibility[old_hap_number].fix_start_gen[j]) + "-" + to_string(old_incompatibility[old_hap_number].fix_end_gen[j]) + ";";
    }
    
    ofs << "\t" << fix_gen << "\t" << old_incompatibility[old_hap_number].lost_gen << endl;
    
}


void record_final_gen(){
    
    string hap;
    
    
    // record
    const int bufsize=1000000;
    char buf[bufsize];
    ofstream ofs;
    ofs.rdbuf()->pubsetbuf(buf,bufsize);

    stringstream FileName;
    FileName << output_folder_path << "/run_" << gen << ".tsv";
    string filename(FileName.str());
    
    ofs.open(filename, ios_base::trunc | ios_base::in);
    
    // columns
    ofs << "gen\tpop\thaplotype\tw\n";
    
    for (int i=0; i<total_ind; ++i){
        
        // haplotype information
        hap = "";
        for (int n=0; n<l_inc; ++n){
            hap += to_string((int)gsl_vector_get(pop[i].hap_inc, n));
        }
        
        ofs << gen << "\t" << pop[i].pop << "\t" << hap << "\t" << pop[i].w << endl;
    }
    
    
    // record
    ofstream ofs2;
    ofs2.rdbuf()->pubsetbuf(buf,bufsize);

    stringstream FileName2;
    FileName2 << output_folder_path << "/prerun_haplotype_" << gen << ".tsv";
    string filename2(FileName2.str());
    
    ofs2.open(filename2, ios_base::trunc | ios_base::in);
    // columns
    ofs2 << "pop\thaplotype\tbirth_gen\tfix_gen\tlost_gen\n";
    
    string fix_gen;
    for (int l=0; l<nHap; ++l){
        ofs2 << incompatibility[l].pop << "\t" << incompatibility[l].hap << "\t" << incompatibility[l].birth_gen;
        
        fix_gen = "";
        for (size_t j=0, n=incompatibility[l].fix_start_gen.size(); j<n; ++j){
            
            if (j==n-1){
                fix_gen += to_string(incompatibility[l].fix_start_gen[j]) + "-";
            }else{
                fix_gen += to_string(incompatibility[l].fix_start_gen[j]) + "-" + to_string(incompatibility[l].fix_end_gen[j]);
            }
            fix_gen += ";";
        }
        
        ofs2 << "\t" << fix_gen << "\t" << incompatibility[l].lost_gen << endl;
        
    }
    
}

void record_param(){
    
    // record
    // population size
    const int bufsize=1000000;
    char buf[bufsize];
    ofstream ofs;
    ofs.rdbuf()->pubsetbuf(buf,bufsize);

    stringstream FileName_pop;
    FileName_pop << output_folder_path << "/pop_size_" << gen << ".tsv";
    string filename_pop(FileName_pop.str());
    
    ofs.open(filename_pop, ios_base::trunc | ios_base::in);
    
    ofs << "gen\tpop\tpop_size" << endl;
    for (int k=0; k<pop_number; ++k){
        ofs << gen << "\t" << k << "\t" << pop_data[k] << endl;
    }
    
    // mutation rate
    ofstream ofs2;
    ofs2.rdbuf()->pubsetbuf(buf,bufsize);
    
    stringstream FileName_mut;
    FileName_mut << output_folder_path << "/mut_rate_" << gen << ".tsv";
    string filename_mut(FileName_mut.str());
    
    ofs2.open(filename_mut, ios_base::trunc | ios_base::in);
    ofs2 << "gen\tpop\tmut_rate" << endl;
    for (int k=0; k<pop_number; ++k){
        ofs2 << gen << "\t" << k << "\t" << u_data[k] << endl;
    }
    
    
}


vector <string> commonElement(vector<string> a, vector<string> b){
    vector <string> c;
    
    for (const auto& element:a){
        if (find(b.begin(), b.end(), element) != b.end()){
            c.push_back(element);
        }
    }
    
    return c;
}

vector <string> diffElement(vector<string> a, vector<string> b){ // only a
    vector <string> d;
    
    for (const auto& element:a){
        if (find(b.begin(), b.end(), element) == b.end()){
            d.push_back(element);
        }
    }
    
    return d;
}



