//
//  KmerCounter.hpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 2/2/23.
//  Copyright © 2023 USC_Mark. All rights reserved.
//

#ifndef KmerCounter_hpp
#define KmerCounter_hpp

#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstring>
#include <unordered_set>
#include <algorithm>
#include <thread>
#include <atomic>
#include <mutex>

#include "fasta.hpp"
#include "fastq.hpp"
#include "KmerStruct.hpp"

using namespace std;


using kmer64_dict = std::unordered_map<u128, uint, hash_128> ;
using kmer32_dict = std::unordered_map<ull, uint > ;

using kmer64_dict_nt = std::unordered_map<u128, uint8, hash_128> ;
using kmer32_dict_nt = std::unordered_map<ull, uint8>;

using kmer64_dict_mul = std::unordered_map<u128, uint*, hash_128> ;
using kmer32_dict_mul = std::unordered_map<ull, uint* > ;

template <int dictsize>
class kmer_counter
{
    using kmer_int = typename std::conditional<(dictsize>32), u128, ull>::type;
    using kmer_dict_type = typename std::conditional<(dictsize>32), kmer64_dict, kmer32_dict>::type;
    using kmer_dict_type_nt = typename std::conditional<(dictsize>32), kmer64_dict_nt, kmer32_dict_nt>::type;
    using kmer_dict_type_mul = typename std::conditional<(dictsize>32), kmer64_dict_mul, kmer32_dict_mul>::type;
    
    kmer_dict_type target_map;
    kmer_dict_type_nt target_map_nt;
    kmer_dict_type_mul target_map_mul;
    
    
    ull totalkmers = 1;
    kmer_int *kmer_records;
    
    int klen = 31 , knum = 0;
    bool iftarget = 0;
    vector<pair<ull,ull>> targetranges;
    uint targetindex = 0;
    std::atomic_uint restfileindex ;
    std::mutex Threads_lock;
    
    std::vector<std::thread*> threads;
    std::vector<std::string> inputfiles;
    std::vector<std::string> outputfiles;
    std::vector<std::string> prefixes;

public:
    
    kmer_counter (int kmersize): klen(31)
    {
        kmer_records = (kmer_int *) malloc(1000);
    };
    ~kmer_counter()
    {
        free(kmer_records);
    };
    template <class typefile>
    void read_counttarget(typefile &fastafile);
    
    void read_target(const char* infile);
        
    template <class typefile>
    void count_kmer(typefile &fastafile, uint8* samplevecs);
    
    void read_files(std::vector<std::string>& inputfiles, std::vector<std::string>& outputfiles, std::vector<std::string>& prefixes ,int numthread);
    
    void read_file();
            
    void write(const char* outfile, uint8*, pair<ull,ull> range);

};

//storage additional positions
template <typename T1, typename T2>
inline static void initiate_counter_mul(T2 &target_map, T1 &larger_kmer, uint index)
{
    
    typename T2::iterator map_find = target_map.find(larger_kmer);

    if (map_find == target_map.end() )
    {
        target_map[larger_kmer] = (uint*) malloc(sizeof(uint)*2);
        target_map[larger_kmer][0] = 1;
        target_map[larger_kmer][1] = index ;
    }
    
    else
    {
        uint* &data = map_find->second;
        if (data[0]%5 == 1)
        {
            data = (uint*) realloc(data, sizeof(uint)*(data[0]+1 + 5));
        }
            
        data[0]++;
        data[data[0]] = index ;
    }
}

template <typename T1, typename T2>
inline static void update_counter_mul(T2 &target_map_mul, T1 &larger_kmer, uint8* vec)
{
    
    typename T2::iterator map_find_mul = target_map_mul.find(larger_kmer);
    
    if (map_find_mul != target_map_mul.end())
    {
        uint* index_mul = map_find_mul->second;
        uint num_num = index_mul[0];
        
        for (uint i = 1 ; i < num_num + 1; ++i)
        {
            if ( vec[index_mul[i]] < 255)  vec[index_mul[i]] ++;
        }
    }
}


template <typename T1, typename T2, typename T3>
inline static void initiate_counter(T2 &target_map, T3 &target_map_mul, T1 &larger_kmer, ull &kindex, ull &start, T1 * &kmer_records)
{
    
    typename T2::iterator map_find = target_map.find(larger_kmer);

    if (map_find == target_map.end())
    {
        if ((kindex % 1000000) == 1)
        {
            kmer_records = (T1 *) realloc(kmer_records, sizeof(T1) * (kindex+1000000));
        }
        kmer_records[kindex] = larger_kmer;
        target_map[larger_kmer] = (int) kindex++;
    }
    
    else if (map_find->second < start)
    {
        initiate_counter_mul(target_map_mul, larger_kmer, map_find->second);
        
        if ((kindex % 1000000) == 1)
        {
            kmer_records = (T1 *) realloc(kmer_records, sizeof(T1) * (kindex+1000000));
        }
        kmer_records[kindex] = larger_kmer;
        target_map[larger_kmer] = (int) kindex++;
        
    }
    
}


template <typename T1, typename T2, typename T3>
inline static void update_counter(T2 &target_map,T3 &target_map_mul, T1 &larger_kmer, uint8* vec)
{
    
    typename T2::iterator map_find = target_map.find(larger_kmer);

    if (map_find != target_map.end())
    {
        uint index =map_find->second;
        
        if ( index == 0) return;
        
        if ( vec[index] < 255) vec[index] ++;
        
        update_counter_mul(target_map_mul, larger_kmer, vec);
    }
}


template <typename T>
static void kmer_read_c(char base, int klen, int &current_size, T &current_kmer, T &reverse_kmer)
{
    int converted = 0;
    T reverse_converted;
    
    if (base == '\n' || base == ' ') return;
    
    if (base_to_int(base, converted))
    {
        
        current_kmer <<= ( 8*sizeof(current_kmer) - 2*klen + 2 );
        current_kmer >>=  ( 8*sizeof(current_kmer) - 2*klen );
        current_kmer += converted;
        
        reverse_kmer >>= 2;
        reverse_converted = 0b11-converted;
        reverse_converted <<= (2*klen-2);
        reverse_kmer += reverse_converted;
        
        
    }
    
    else
    {
        current_size = -1;
        current_kmer = 0;
        reverse_kmer = 0;
    }
}

template <int dictsize>
template <class typefile>
void kmer_counter<dictsize>::read_counttarget(typefile &fastafile)
{
        
    int current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    iftarget = 1;
    targetranges.emplace_back();
    auto & targetrange = targetranges[targetranges.size()-1];
    targetrange.first = totalkmers;
    
    std::string StrLine;
    
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '@':  case '+': case '>':
                current_size = 0;
                continue;
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }
        
        for (auto base: StrLine)
        {
            if (base == '\0') break;
            
            if (base == '\n' || base == ' ') continue;

            kmer_read_c(base, klen, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen) continue;
                                
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer : reverse_kmer;
            
            initiate_counter(target_map, target_map_mul, larger_kmer, totalkmers, targetrange.first, kmer_records);
                
        }
    }
    
    fastafile.Close();
    
    targetrange.second = totalkmers ;
    
};
 
template <int dictsize>
void kmer_counter<dictsize>::read_target(const char* inputfile)
{
    
    fasta readsfile(inputfile);
    
    read_counttarget(readsfile);
    
}



template <int dictsize>
template <class typefile>
void kmer_counter<dictsize>::count_kmer(typefile &fastafile, uint8* samplevecs)
{
    
    int current_size = 0;
    
    kmer_int current_kmer = 0;
    kmer_int reverse_kmer = 0;
    
    //uint64_t ifmasked = 0;
    //int num_masked = 0 ;
    std::string StrLine;
    std::string Header;
    
    while (fastafile.nextLine(StrLine))
    {
        switch (StrLine[0])
        {
            case '@':  case '+': case '>':
                current_size = 0;
                continue;
            case ' ': case '\n': case '\t':
                continue;
            default:
                break;
        }
        
        for (auto base: StrLine)
        {
            if (base == '\0') break;
                        
            if (base == '\n' || base == ' ') continue;
 
            kmer_read_c(base, klen, current_size, current_kmer, reverse_kmer);
            
            if (++current_size < klen) continue;
            
            auto larger_kmer = (current_kmer >= reverse_kmer) ? current_kmer:reverse_kmer;
            
            if (iftarget)
            {
                update_counter(target_map, target_map_mul,larger_kmer, samplevecs);
            }
            else
            {
                totalkmers++;
                if (target_map_nt[larger_kmer]<255) target_map_nt[larger_kmer] ++;
                
            }
        }
    }
        
    fastafile.Close();
};

template <int dictsize>
void kmer_counter<dictsize>::read_files(std::vector<std::string>& inputs, std::vector<std::string>& outputs, std::vector<std::string>& prefs, int nthreads)
{
    
    inputfiles = inputs;
    outputfiles = outputs;
    prefixes = prefs;
    restfileindex = 0;
    
    std::vector<std::thread*> threads;
    
    for(int i=0; i< nthreads; ++i)
    {
        std::thread *newthread_ = new std::thread(&kmer_counter<dictsize>::read_file, this);
        threads.push_back(newthread_);
    }
    
    
    for(int i=0; i< nthreads; ++i)
    {
        threads[i]->join();
    }
    
}


template <int dictsize>
void kmer_counter<dictsize>::read_file()
{
    
    while (restfileindex < inputfiles.size())
    {
        
        Threads_lock.lock();
        
        int inputindex = restfileindex++ ;
        
        Threads_lock.unlock();
        
        if (inputindex >= inputfiles.size()) break;
        
        uint8* samplevecs = (uint8* )malloc(sizeof(uint8) * totalkmers);
        
        memset(samplevecs, 0, sizeof(uint8) * totalkmers);
        
        const char* inputfile = inputfiles[inputindex ].c_str();
        
        int pathlen = (int)strlen(inputfile);
        
        if ( pathlen > 2 && strcmp(inputfile+(pathlen-3),".gz") == 0 )
        {
            
            fastq readsfile(inputfile);
            
            count_kmer(readsfile, samplevecs);
            
        }
            
        else
        {
            
            fasta readsfile(inputfile);
            
            count_kmer(readsfile, samplevecs);
        }
        
        string prefix = "";
        
        if (prefixes.size() > inputindex) prefix = prefixes[inputindex];
        
        for (int j = 0; j <outputfiles.size(); ++j)
        {
            auto outputfile = outputfiles[j] + prefix;
            pair<ull,ull> range;
            if (iftarget)
            {
                range = targetranges[j];
            }
            else
            {
                range = make_pair(1,totalkmers);
            }
            
            write(outputfile.c_str(), samplevecs, range);
        }
        
        free(samplevecs);
        
    }
        
}




template <int dictsize>
void kmer_counter<dictsize>::write(const char * outputfile, uint8* counts, pair<ull,ull> range)
{
    
    FILE *fwrite=fopen(outputfile, "w");
    
    
    if (fwrite==NULL)
    {
        std::cerr << "ERROR: Cannot write file: " << outputfile << endl;
        
        std::_Exit(EXIT_FAILURE);
    }
    
    int code_bit = 1;
    
    int digit = (int)floor(1.0*klen/code_bit);
    
    char kmer_seq[digit+1];
    kmer_seq[digit] = '\0';
    
    if (iftarget)
    {
        fprintf(fwrite,"@targetsize\t%llu\n" , range.second -  range.first);
        
        for (auto i = range.first; i < range.second; ++i)
        {
            uint8 count = counts[i];
            kmer_int seq = kmer_records[i];
            
            if (count == 0) continue;
            
            for (int index = digit-1; index >= 0 ; --index)
            {
                kmer_seq[index] = "ACGT"[seq%4];
                seq /= 4;
            }
            
            fprintf(fwrite,"%s\t%d\n", kmer_seq, count);
        }
    }
    else
    {
        
        fprintf(fwrite,"@targetsize\t%llu\n", range.second -  range.first);
        
        for (const auto &kmer: target_map_nt)
        {
            
            auto seq = kmer.first;
            auto count = kmer.second;
            
            for (int index = digit-1; index >= 0 ; --index)
            {
                kmer_seq[index] = "ACGT"[seq%4];
                seq /= 4;
            }
            
            fprintf(fwrite,"%s\t%d\n", kmer_seq, count);
        }
    }
    
    fclose(fwrite);
    
    return ;
}



#endif /* KmerCounter_hpp */
