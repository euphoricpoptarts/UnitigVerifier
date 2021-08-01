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

std::vector<seq_t> load_fasta(char* fname){
    std::vector<seq_t> seqs;
    seqf_t fin{std::ifstream(fname), seqan3::format_fasta{}};
    for(auto & record : fin){
        seqs.push_back(record.sequence());
    }
    printf("%u total sequences\n", seqs.size());
    return seqs;
}

void verify_unitigs(const search_strings& searches, seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection>& index, std::atomic<uint32_t>& found, std::atomic<uint32_t>& not_found, std::atomic<uint32_t>& rcomp_found, int* missing, long offset, long stride){
    for(long i = offset; i < searches.total_strings; i += stride){
        //const char* loc = strstr(s.c_str(), searches.chars + searches.offsets[i]);
        auto begin = searches.chars.begin();
        auto end = searches.chars.begin();
        std::advance(begin, searches.offsets[i]);
        std::advance(end, searches.offsets[i + 1] - 1);
        std::string y(begin, end);
        seqan3::dna5_vector x = y | seqan3::views::char_to<seqan3::dna5> | seqan3::views::to<std::vector>;
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
            auto xr = x | std::views::reverse | seqan3::views::complement;
            for(auto && result : seqan3::search(xr, index)){
                hits++;
            }
            if(hits > 0){
                rcomp_found++;
                uint32_t total_found = ++found;
                if(total_found % 100000 == 0) {
                    printf("found %u\n", total_found);
                }
            }
            else {
                uint32_t t_n_f = not_found++;
                missing[t_n_f] = i;
            }
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
    std::vector<seq_t> seqs = load_fasta(fasta_fname);
    tp t2 = hc::now();
    chr::duration<double> m1 = chr::duration_cast<chr::duration<double>>(t2 - t1);
    printf("loaded %u unitigs in %.3f seconds\n", searches.total_strings, m1.count());
    seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> index{seqs};
    tp t3 = hc::now();
    chr::duration<double> m2 = chr::duration_cast<chr::duration<double>>(t3 - t2);
    printf("built index in %.3f seconds\n", m2.count());
    std::atomic<uint32_t> found{0};
    std::atomic<uint32_t> rcomp_found{0};
    std::atomic<uint32_t> not_found{0};
    std::vector<std::thread> threads;
    int* missing = new int[searches.total_strings];
    long total_threads = 32;
    for (long i = 0; i < total_threads; ++i){
        threads.push_back(std::thread(verify_unitigs, std::ref(searches), std::ref(index), std::ref(found), std::ref(not_found), std::ref(rcomp_found), missing, i, total_threads));
    }

    for (auto& t : threads){
        t.join();
    }
    tp t4 = hc::now();
    chr::duration<double> m3 = chr::duration_cast<chr::duration<double>>(t4 - t3);
    uint32_t total_found = found;
    uint32_t t_rcomp_found = rcomp_found;
    if(total_found == searches.total_strings){
        printf("Found every unitig; %u were found by reverse complement\n", t_rcomp_found);
    } else {
        printf("Found %u out of %u unitigs!\n", total_found, searches.total_strings);
        for(int u = 0; u < not_found; u++){
            int i = missing[u];
            auto begin = searches.chars.begin();
            auto end = searches.chars.begin();
            std::advance(begin, searches.offsets[i]);
            std::advance(end, searches.offsets[i + 1] - 1);
            std::string y(begin, end);
            seqan3::dna5_vector x = y | seqan3::views::char_to<seqan3::dna5> | seqan3::views::to<std::vector>;
            auto xr = x | std::views::reverse | seqan3::views::complement;
            seqan3::debug_stream << "Could not find: " << x << "\n or r-comp: " << xr << "\n";
        }
    }
    printf("Done verifying unitigs in %.3f seconds!\n", m3.count());
    delete[] searches.offsets;
    delete[] missing;
    return 0;
}
