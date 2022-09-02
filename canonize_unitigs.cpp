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
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/utility/type_list/type_list.hpp>
#include <seqan3/alphabet/views/all.hpp>
#include <seqan3/core/debug_stream.hpp>

namespace fs = std::filesystem;
namespace chr = std::chrono;
using hc = chr::high_resolution_clock;
using seqf_t = seqan3::sequence_file_input<seqan3::sequence_file_input_default_traits_dna>;
using seq_t = seqan3::dna5_vector;

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

void canonize_unitigs(const search_strings& searches){
    for(long i = 0; i < searches.total_strings; i ++){
        auto begin = searches.chars.begin();
        auto end = searches.chars.begin();
        std::advance(begin, searches.offsets[i]);
        std::advance(end, searches.offsets[i + 1] - 1);
        std::string y(begin, end);
        seqan3::dna5_vector x = y | seqan3::views::char_to<seqan3::dna5> | seqan3::views::to<std::vector>;
        seqan3::dna5_vector xr = x | std::views::reverse | seqan3::views::complement | seqan3::views::to<seqan3::dna5_vector>;
        if(x < xr){
            seqan3::debug_stream << x << '\n';
        } else {
            seqan3::debug_stream << xr << '\n';
        }
    }
}

int main(int argc, char **argv) {

    if (argc != 2) {
        printf("You input %d args\n", argc);
        fprintf(stderr, "Usage: %s <unitig file>\n", argv[0]);
        return 1;
    }
    char *unitig_fname = argv[1];
    using tp = hc::time_point;
    search_strings searches = load_unitigs(unitig_fname);
    canonize_unitigs(searches);
    delete[] searches.offsets;
    return 0;
}
