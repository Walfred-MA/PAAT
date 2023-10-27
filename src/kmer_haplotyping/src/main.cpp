//
//  main.cpp
//  kmer_haplotyping
//
//  Created by Wangfei MA on 10/13/21.
//  Copyright © 2021 USC_MarkLab. All rights reserved.
//

#include <iostream>
#include <string>
#include <filesystem>

#include "KmerSeacher.hpp"
#include "KmerCounter.hpp"
#include "KmerFilter.hpp"
#include "KmerCompare.hpp"


template <int dictsize>
void run(std::vector<std::string> &inputfiles, std::vector<std::string>& targetfiles, std::vector<std::string> &outputfiles, std::vector<std::string> &prefixes, const int kmer_size, const int nthreads, const int mode, const int cutoff)
{

    if (mode == 3)
    {
        kmer_filter<dictsize> filter(kmer_size);
      
    	std::vector<std::string> includes;
	std::vector<std::string> excludes;
	for (std::string &inputfile: inputfiles)
        {
            if (inputfile.size() && inputfile.find("_exclude.txt") == std::string::npos)
            {
                includes.push_back(inputfile.c_str());
            }
            else
            {
                excludes.push_back(inputfile.c_str());
            }
        }       
 
        filter.read_targets(includes, nthreads);
        filter.read_files(excludes, outputfiles, prefixes, nthreads);
        
        return;
    } 

    else if (mode == 4)
    {

        kmer_compare<dictsize> compare(kmer_size);
      
        for (std::string targetfile: targetfiles)
        {

            if (targetfile.size())
            {
                compare.read_target(targetfile.c_str());
            }
            
        }
        
        compare.read_files(inputfiles, outputfiles, prefixes, nthreads);       

        return;
    }

    if (mode == 1 or mode == 2)
    {
        kmer_map<dictsize> map(kmer_size, cutoff);
        
        for (std::string targetfile: targetfiles)
        {

            if (targetfile.size())
            {
                map.read_target(targetfile.c_str());
            }
            
        }
        
        map.read_files(inputfiles, outputfiles, prefixes, nthreads);
        
    }
    
    if (mode == 0 or mode == 2)
    {
        
        kmer_counter<dictsize> counter(kmer_size);
        
        for (std::string targetfile: targetfiles)
        {

            if (targetfile.size())
            {
                counter.read_target(targetfile.c_str());
            }
            
        }
        
        counter.read_files(inputfiles, outputfiles, prefixes, nthreads);
    }
    
}


int main(int argc, const char * argv[]) {
    
    
    std::vector<std::string> inputfiles;
    
    std::vector<std::string> prefixes;
    
    std::vector<std::string> targetfiles;
    
    std::vector<std::string> outputfiles;

    const char* Argument="";
        
    int mode = 0, kmer_size = 31,  nthreads = 1, cutoff =50;
    
    for (int i = 1; i < argc ; i++)
    {
        
        if (argv[i][0] == '-')
        {
            Argument = argv[i];
        }
        else if (strcmp(Argument, "-i")==0 or strcmp(Argument, "--input")==0)
        {

	    if (std::filesystem::is_directory(argv[i]))
            {
                for (const auto& entry : std::filesystem::directory_iterator(argv[i]))
                {
                    if (std::filesystem::is_regular_file(entry.path()))
                    {
                        inputfiles.push_back(entry.path().string());
                    }
                }
            }
            else
            {
                inputfiles.push_back(argv[i]);
            }

        }
        
        else if (strcmp(Argument, "-I")==0 or strcmp(Argument, "--Inputs")==0)
        {
            std::ifstream pathfile(argv[i]);
            std::string line;
            if(!pathfile)
            {
                std::cout<<"Error opening target file"<<std::endl;
                return -1;
            }
            while (std::getline(pathfile, line))
            {
                inputfiles.push_back(line);
            }
        }
        
        else if (strcmp(Argument, "-t")==0 or strcmp(Argument, "--target")==0)
        {
	    if (std::filesystem::is_directory(argv[i]))
            {
                for (const auto& entry : std::filesystem::directory_iterator(argv[i]))
                {
                    if (std::filesystem::is_regular_file(entry.path()))
                    {
                        targetfiles.push_back(entry.path().string());
                    }
                }
            }
            else
            {
                targetfiles.push_back(argv[i]);
            }
        }
        
        else if (strcmp(Argument, "-T")==0 or strcmp(Argument, "--Targets")==0)
        {
            
            std::ifstream pathfile(argv[i]);
            std::string line;
            if(!pathfile)
            {
                std::cout<<"Error opening target file"<<std::endl;
                return -1;
            }
            while (std::getline(pathfile, line))
            {
                targetfiles.push_back(line);
            }
        }
        
        else if (strcmp(Argument, "-p")==0 or strcmp(Argument, "--pref")==0)
        {
            prefixes.push_back(argv[i]);
        }
        
        else if (strcmp(Argument, "-P")==0 or strcmp(Argument, "--Prefs")==0)
        {
            std::ifstream pathfile(argv[i]);
            std::string line;
            if(!pathfile)
            {
                std::cout<<"Error opening target file"<<std::endl;
                return -1;
            }
            while (std::getline(pathfile, line))
            {
                prefixes.push_back(line);
            }
        }
        
        else if (strcmp(Argument, "-o")==0 or strcmp(Argument, "--output")==0)
        {
            outputfiles.push_back(argv[i]);
        }
        
        else if (strcmp(Argument, "-O")==0 or strcmp(Argument, "--Outputs")==0)
        {
            std::ifstream pathfile(argv[i]);

            if(!pathfile)
            {
                std::cout<<"Error opening output file"<<std::endl;
                return -1;
            }
            std::string line;
            while (std::getline(pathfile, line))
            {
                outputfiles.push_back(line);
            }

        }
        
        else if (strcmp(Argument, "-k")==0 or strcmp(Argument, "--kmer")==0)
        {
            kmer_size=(int)atoi(argv[i]);
        }
        
        else if (strcmp(Argument, "-s")==0 or strcmp(Argument, "--search")==0)
        {
            mode=(int)atoi(argv[i]);
        }
        
        else if (strcmp(Argument, "-n")==0 or strcmp(Argument, "--nthreads")==0)
        {
            nthreads=(int)atoi(argv[i]);
        }
        
        else if (strcmp(Argument, "-c")==0 or strcmp(Argument, "--cutoff")==0)
        {
            cutoff=(int)atoi(argv[i]);
        }
        
    }
    
    if (!inputfiles.size()) return 1;
    


 
    run<32>(inputfiles, targetfiles, outputfiles, prefixes,kmer_size, nthreads,mode,cutoff);
    
    /*
    else if (kmer_size<=64)
    {
        run<64>(inputfiles, targetfiles, outputfiles, prefixes, kmer_size, nthreads, mode);
    }
    
    
    else if (kmer_size<=128)
    {
        run<127>(inputfiles, targetfiles, outputfiles, kmer_size, mode);
    }
     */
    
    
    return 0;
}
