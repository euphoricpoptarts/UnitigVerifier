#include <iostream>
#include <fstream>
#include <functional>
#include <utility>
#include <filesystem>
#include <cstring>
#include <atomic>
#include <vector>
#include <thread>
#include <chrono>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

namespace fs = std::filesystem;
namespace chr = std::chrono;
using hc = chr::high_resolution_clock;

struct search_strings {
    uint64_t total_chars;
    std::string chars;
    uint64_t* offsets;
    uint32_t total_strings;
};

size_t getFileLength(const char *fname){
    // Obtain the size of the file.
    return fs::file_size(fname);
}

search_strings load_unitigs(char *fname) {

    std::ifstream infp(fname);
    if (!infp.is_open()) {
        printf("Error: Could not open input file. Exiting ...\n");
        exit(1);
    }

    auto sz = getFileLength(fname);

    // Create a buffer for file
    std::string s1(sz, '\0');
    // Read the whole file into the buffer.
    infp.read(s1.data(), sz);
    infp.close();

    uint32_t n = 0;
    uint64_t c_count = 0;
    const char* f = s1.c_str();
    for(long i = 0; i < sz; i++){
        if(f[i] != '\n') c_count++;
        else n++;
    }
    search_strings s;
    s.total_chars = c_count;
    s.total_strings = n;
    //n '\0' chars
    //i guess this would just be sz
    char* chars = new char[c_count + n];
    s.offsets = new uint64_t[n + 1];
    n = 0;
    s.offsets[0] = 0;
    n = 1;
    long write = 0;
    //basically strcpy but replace '\n' with '\0'
    for(long i = 0; i < sz; i++){
        if(f[i] == '\n') {
            chars[write++] = '\0';
            s.offsets[n++] = i + 1;
        } else {
            chars[write++] = f[i];
        }
    }
    std::string chars_s(chars, sz);
    s.chars = chars_s;
    delete[] chars;
    return s;
}

std::string load_fasta(char *fname){
    std::ifstream infp(fname);
    if (!infp.is_open()) {
        printf("Error: Could not open input file. Exiting ...\n");
        exit(1);
    }

    auto sz = getFileLength(fname);

    // Create a buffer for file
    std::string s1(sz, '\0');
    // Read the whole file into the buffer.
    infp.read(s1.data(), sz);
    infp.close();
    return s1;
}

void verify_unitigs(const search_strings& searches, seqan3::fm_index<char, seqan3::text_layout::single>& index, std::atomic<uint32_t>& found, long offset, long stride){
    for(long i = offset; i < searches.total_strings; i += stride){
        //const char* loc = strstr(s.c_str(), searches.chars + searches.offsets[i]);
        auto begin = searches.chars.begin();
        auto end = searches.chars.begin();
        std::advance(begin, searches.offsets[i]);
        std::advance(end, searches.offsets[i + 1] - 1);
        std::string x(begin, end);
        int hits = 0;
        for(auto && result : seqan3::search(x, index)){
            hits++;
        }
        if(hits > 0){
            uint32_t total_found = ++found;
            if(total_found % 100000 == 0) {
                printf("found %u\n", total_found);
            }
        } else {
            printf("Failed to find one!\n");
        }
    }
}

int main(int argc, char **argv) {

    if (argc < 3) {
        printf("You input %d args\n", argc);
        fprintf(stderr, "Usage: %s <unitig file> <fasta file>\n", argv[0]);
        return 1;
    }
    char *unitig_fname = argv[1];
    char *fasta_fname = argv[2];
    using tp = hc::time_point;
    tp t1 = hc::now();
    search_strings searches = load_unitigs(unitig_fname);
    std::string s = load_fasta(fasta_fname);
    tp t2 = hc::now();
    chr::duration<double> m1 = chr::duration_cast<chr::duration<double>>(t2 - t1);
    printf("loaded %u unitigs in %.3f seconds\n", searches.total_strings, m1.count());
    seqan3::fm_index index{s};
    tp t3 = hc::now();
    chr::duration<double> m2 = chr::duration_cast<chr::duration<double>>(t3 - t2);
    printf("built index in %.3f seconds\n", m2.count());
    std::atomic<uint32_t> found{0};
    std::vector<std::thread> threads;
    long total_threads = 32;
    for (long i = 0; i < total_threads; ++i){
        threads.push_back(std::thread(verify_unitigs, std::ref(searches), std::ref(index), std::ref(found), i, total_threads));
    }

    for (auto& t : threads){
        t.join();
    }
    tp t4 = hc::now();
    chr::duration<double> m3 = chr::duration_cast<chr::duration<double>>(t4 - t3);
    if(found == searches.total_strings){
        printf("Found every unitig\n");
    }
    printf("Done verifying unitigs in %.3f seconds!\n", m3.count());
    delete[] searches.offsets;
    return 0;
}
