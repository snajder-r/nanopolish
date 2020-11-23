//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_alphabet -- support for multiple alphabets
//
#ifndef NANOPOLISH_ALPHABET_H
#define NANOPOLISH_ALPHABET_H

#include <string>
#include <cstring>
#include <inttypes.h>
#include <assert.h>
#include <algorithm>
#include "nanopolish_iupac.h"
#include <iostream>

#define METHYLATED_C_SYMBOL 'Z'
#define METHYLATED_A_SYMBOL 'X'
#define METHYLATED_A_COMPLEMENT 'Q' // a T bound to a m6A

struct RecognitionMatch
{
    unsigned offset; // the matched position in the recognition site
    unsigned length; // the length of the match, 0 indicates no match
    bool covers_methylated_site; // does the match cover an M base?
};

class RecognitionMask
{
    public:
        bool* mask; //whether base is in a recognition site
        int* ri; //recognition site index of site starting here
        RecognitionMatch* matches; // offeset, length and met status of site

        RecognitionMask(size_t length){
            ri = new int[length];
            mask = new bool[length];
            matches = new RecognitionMatch[length];
        }
        ~RecognitionMask(){
            delete mask;
            delete ri;
            delete matches;
        }
};

// Check whether a recognition site starts at position i of str
inline RecognitionMatch match_to_site(const std::string& str, size_t i, const char* recognition, size_t rl)
{
    RecognitionMatch match;
    match.length = 0;
    match.offset = 0;
    match.covers_methylated_site = false;

    // Case 1: prefix of str is suffix of recognition
    if(i == 0) {
        for(int j = 0; j < rl; j++){
            size_t cl = std::min(rl-j, str.length());
            if(str.compare(0,cl,recognition,j,cl)==0){
                match.offset = j;
                match.length = cl;
                break;
            }
        }
    } else {
        // Case 2: the suffix str[i..n] is a prefix of recognition
        size_t cl = std::min(rl, str.length() - i);
        if(str.compare(i, cl, recognition, cl) == 0) {
            match.offset = 0;
            match.length = cl;
        }    
    }   

    //printf("Match site: %s %s %s %d %d\n", str.c_str(), str.substr(i).c_str(), recognition, match.offset, match.length);
    if(match.length > 0) {
        match.covers_methylated_site = 
            str.substr(i, match.length).find_first_of(METHYLATED_C_SYMBOL) != std::string::npos;
        match.covers_methylated_site = match.covers_methylated_site ||
            str.substr(i, match.length).find_first_of(METHYLATED_A_SYMBOL) != std::string::npos;
        match.covers_methylated_site = match.covers_methylated_site ||
            str.substr(i, match.length).find_first_of(METHYLATED_A_COMPLEMENT) != std::string::npos;

    }

    return match;
}


// Abstract base class for alphabets
class Alphabet
{
    public:
        // basic functions
        virtual uint8_t rank(char b) const = 0;
        virtual char base(uint8_t r) const = 0;
        virtual char complement(char b) const = 0;
        virtual uint32_t size() const = 0;
        virtual const char* get_name() const = 0;

        // support for methylated bases with recognition sequences
        virtual size_t num_recognition_sites() const = 0;
        virtual size_t recognition_length(size_t i) const = 0;
        virtual const char* get_recognition_site(size_t i) const = 0;
        virtual const char* get_recognition_site_methylated(size_t i) const = 0;
        virtual const char* get_recognition_site_methylated_complement(size_t i) const = 0;

        // return the lexicographic rank of the kmer amongst all strings of 
        // length k for this alphabet
        inline uint32_t kmer_rank(const char* str, uint32_t k) const
        {
            uint32_t p = 1;
            uint32_t r = 0;

            // from last base to first
            for(uint32_t i = 0; i < k; ++i) {
                r += rank(str[k - i - 1]) * p;
                p *= size();
            }
            return r;
        }
        
        // Increment the input string to be the next sequence in lexicographic order
        inline void lexicographic_next(std::string& str) const
        {
            int carry = 1;
            int i = str.size() - 1;
            do {
                uint32_t r = rank(str[i]) + carry;
                str[i] = base(r % size());
                carry = r / size();
                i -= 1;
            } while(carry > 0 && i >= 0);
        }

        // returns the number of unique strings of length l for this alphabet
        inline size_t get_num_strings(size_t l) const
        {
            size_t s = size();
            size_t n = 1;
            for(size_t i = 0; i < l; ++i) {
                n *= s;
            }
            return n;
        }

        inline RecognitionMask find_recognition_mask(const std::string& str, 
                const bool require_methylated, const bool require_full_length,
                const bool methylated_recognition_pattern) const
        {
            RecognitionMask ret(str.length());
            int remainder_of_site = 0;
            for(size_t i = 0; i < str.length(); i++) {
                int recognition_index = -1;
                RecognitionMatch match;

                // Does this location (partially) match a methylated recognition site?
                for(size_t j = 0; j < num_recognition_sites(); ++j) {
                    const char* recognition_pattern;
                    if(methylated_recognition_pattern){
                         recognition_pattern = get_recognition_site_methylated(j);
                    }else{
                         recognition_pattern = get_recognition_site(j); 
                    }

                    match = match_to_site(str, i, recognition_pattern, recognition_length(j));
                    if(match.length > 0 && (match.covers_methylated_site || !require_methylated) && 
                            ((match.length == recognition_length(j)) || !require_full_length)) {
                        recognition_index = j;
                        remainder_of_site = match.length;
                        break;
                    }
                }
                ret.matches[i] = match;
                ret.ri[i] = recognition_index;

                ret.mask[i] = (remainder_of_site--) > 0;
            }
            
            return ret;
        }

        // reverse-complement a string
        // when the string contains methylated bases, the methylation
        // symbol transfered to the output strand in the appropriate position
        virtual std::string reverse_complement(const std::string& str) const
        {
            std::string out(str.length(), 'A');
            RecognitionMask rmask = find_recognition_mask(str, true, false, true);

            size_t i = 0; // input
            int j = str.length() - 1; // output

            while(i < str.length()){
                // If this subsequence matched a methylated recognition site,
                // copy the complement of the site to the output
                if(rmask.ri[i] != -1) {
                    int im = i;
                    for(size_t k = rmask.matches[im].offset; k < rmask.matches[im].offset + rmask.matches[im].length; ++k) {
                        out[j--] = get_recognition_site_methylated_complement(rmask.ri[im])[k];
                        i += 1;
                        if(i < str.length()){
                           if(rmask.ri[i] != -1){
                               // If the next index is ALSO the starting point of a recognition
                               // site, use the new recognition site to fill in the rest of the
                               // methylated complement. This resolves overlapping recognition sites
                               // such as CGC (CG overlapping with GC)
                               break;
                           }
                        }
                    }
                } else {
                    // complement a single base
                    assert(str[i] != METHYLATED_C_SYMBOL);
                    assert(str[i] != METHYLATED_A_SYMBOL);
                    assert(str[i] != METHYLATED_A_COMPLEMENT);
                    out[j--] = complement(str[i++]);
                }
            }
            return out;
        }

        // return a new copy of the string with IUPAC ambiguity characters changed
        virtual std::string disambiguate(const std::string& str) const
        {
            RecognitionMask rmask = find_recognition_mask(str, false, false, true);
            // create output and convert lower case to upper case
            std::string out(str);
            std::transform(out.begin(), out.end(), out.begin(), ::toupper);

            for(size_t i = 0; i < out.length(); i++) {
                // disambiguate if not a recognition site
                if(!rmask.mask[i]){
                    assert(IUPAC::isValid(out[i]));
                    out[i] = IUPAC::getPossibleSymbols(out[i])[0];
                }
            }
            return out;
        }

        // If the alphabet supports methylated bases, convert str
        // to a methylated string using the recognition sites
        virtual std::string methylate(const std::string& str) const
        {
            std::string out(str);
            RecognitionMask rmask = find_recognition_mask(str, false, false, false);
            size_t i = 0;
            while(i < out.length()){
                // If this subsequence matched a recognition site,
                if(rmask.ri[i]!=-1) {
                    int im = i;
                    for(size_t k = rmask.matches[im].offset; k < rmask.matches[im].offset + rmask.matches[im].length; ++k) {
                        out[i] = get_recognition_site_methylated(rmask.ri[im])[k];
                        i += 1;
                        if(i < str.length()){
                           if(rmask.ri[i] != -1){
                               // If the next index is ALSO the starting point of a recognition
                               // site, use the new recognition site to fill in the rest
                               break;
                           }
                        }
                    }
                    assert(i!=im); // that would be a length 0 match
                }else{
                    i++;
                }
            }

            return out;
        }


        // Remove methylated bases according to the recognition site
        std::string unmethylate(const std::string& str) const
        {
            std::string out(str);
            RecognitionMask rmask = find_recognition_mask(str, false, false, true);

            size_t i = 0;
            while(i < out.length()){
                // If this subsequence matched a methylated recognition site,
                if(rmask.ri[i]!=-1) {
                    int im = i;
                    for(size_t k = rmask.matches[im].offset; k < rmask.matches[im].offset + rmask.matches[im].length; ++k) {
                        out[i] = get_recognition_site(rmask.ri[im])[k];
                        i += 1;
                        if(i < str.length()){
                           if(rmask.ri[i] != -1){
                               // If the next index is ALSO the starting point of a recognition
                               // site, use the new recognition site to fill in the rest
                               break;
                           }
                        }
                    }
                    assert(i!=im); // that would be a length 0 match
                }else{
                    i++;
                }
            }

            return out;
        }

        // does this alphabet contain all of the nucleotides in bases?
        virtual bool contains_all(const char *bases) const = 0;

        // check if the motif matches a recognition site
        bool is_motif_match(const std::string& str, const size_t i) const
        {
            RecognitionMatch match;
            for(size_t j = 0; j < num_recognition_sites(); ++j){
                match = match_to_site(str, i, get_recognition_site(j), recognition_length(j));
                if(match.length == recognition_length(j))
                    return true;
            }
            return false;
        }
};

#define BASIC_MEMBER_BOILERPLATE \
    static const uint8_t _rank[256]; \
    static const char* _name; \
    static const char* _base; \
    static const char* _complement; \
    static const uint32_t _size; 

#define BASIC_ACCESSOR_BOILERPLATE \
    virtual const char* get_name() const { return _name; } \
    virtual uint8_t rank(char b) const { return _rank[(int)b]; }        \
    virtual char base(uint8_t r) const { return _base[r]; } \
    virtual char complement(char b) const { return _complement[_rank[(int)b]]; } \
    virtual uint32_t size() const { return _size; } \

struct DNAAlphabet : public Alphabet
{
    // members
    BASIC_MEMBER_BOILERPLATE

    // functions
    BASIC_ACCESSOR_BOILERPLATE

    // no methylation in this alphabet
    virtual size_t num_recognition_sites() const { return 0; }
    virtual size_t recognition_length(size_t) const { return 0; }
    virtual const char* get_recognition_site(size_t) const { return NULL; }
    virtual const char* get_recognition_site_methylated(size_t) const { return NULL; }
    virtual const char* get_recognition_site_methylated_complement(size_t) const { 
        return NULL;
    }

    // does this alphabet contain all of the nucleotides in bases?
    virtual inline bool contains_all(const char *bases) const
    {
        return strspn(bases, _base) == strlen(bases);
    }
};

// An RNA alphabet where the Us have been converted to Ts
// This is used for direct RNA sequencing as most reference sequences
// contain Ts into of U
struct UtoTRNAAlphabet : public Alphabet
{
    // members
    BASIC_MEMBER_BOILERPLATE

    // functions
    BASIC_ACCESSOR_BOILERPLATE

    // no methylation in this alphabet
    virtual size_t num_recognition_sites() const { return 0; }
    virtual size_t recognition_length(size_t) const { return 0; }
    virtual const char* get_recognition_site(size_t) const { return NULL; }
    virtual const char* get_recognition_site_methylated(size_t) const { return NULL; }
    virtual const char* get_recognition_site_methylated_complement(size_t) const { 
        return NULL;
    }

    // does this alphabet contain all of the nucleotides in bases?
    virtual inline bool contains_all(const char *bases) const
    {
        return strspn(bases, _base) == strlen(bases);
    }
};


#define METHYLATION_MEMBER_BOILERPLATE \
    static const uint32_t _num_recognition_sites; \
    static const uint32_t _recognition_length[]; \
    static const char* _recognition_sites[]; \
    static const char* _recognition_sites_methylated[]; \
    static const char* _recognition_sites_methylated_complement[];

#define METHYLATION_ACCESSOR_BOILERPLATE \
    virtual size_t num_recognition_sites() const { return _num_recognition_sites; } \
    virtual size_t recognition_length(size_t i) const { return _recognition_length[i]; } \
    virtual const char* get_recognition_site(size_t i) const { return _recognition_sites[i]; } \
    virtual const char* get_recognition_site_methylated(size_t i) const { return _recognition_sites_methylated[i]; } \
    virtual const char* get_recognition_site_methylated_complement(size_t i) const { \
        return _recognition_sites_methylated_complement[i]; \
    } 

//
// methyl-cytosine in CG context
// 
struct MethylCpGAlphabet : public Alphabet
{
    // member variables, expanded by macrocs
    BASIC_MEMBER_BOILERPLATE
    METHYLATION_MEMBER_BOILERPLATE
    
    // member functions
    BASIC_ACCESSOR_BOILERPLATE
    METHYLATION_ACCESSOR_BOILERPLATE

    // does this alphabet contain all of the nucleotides in bases?
    virtual inline bool contains_all(const char *bases) const 
    {
        return strspn(bases, _base) == strlen(bases);
    }
};

struct M5CAlphabet : public Alphabet
{
    // member variables, expanded by macrocs
    BASIC_MEMBER_BOILERPLATE
    METHYLATION_MEMBER_BOILERPLATE

    // member functions
    BASIC_ACCESSOR_BOILERPLATE
    METHYLATION_ACCESSOR_BOILERPLATE

    // does this alphabet contain all of the nucleotides in bases?
    virtual inline bool contains_all(const char *bases) const
    {
        return strspn(bases, _base) == strlen(bases);
    }
};

struct NOMEAlphabet : public Alphabet
{
    // member variables, expanded by macrocs
    BASIC_MEMBER_BOILERPLATE
    METHYLATION_MEMBER_BOILERPLATE
    
    // member functions
    BASIC_ACCESSOR_BOILERPLATE
    METHYLATION_ACCESSOR_BOILERPLATE

    // does this alphabet contain all of the nucleotides in bases?
    virtual inline bool contains_all(const char *bases) const 
    {
        return strspn(bases, _base) == strlen(bases);
    }
};

//
// methyl-cytosine in GC context
//
struct MethylGpCAlphabet : public Alphabet
{
    // member variables, expanded by macrocs
    BASIC_MEMBER_BOILERPLATE
    METHYLATION_MEMBER_BOILERPLATE

    // member functions
    BASIC_ACCESSOR_BOILERPLATE
    METHYLATION_ACCESSOR_BOILERPLATE

    // does this alphabet contain all of the nucleotides in bases?
    virtual inline bool contains_all(const char *bases) const
    {
        return strspn(bases, _base) == strlen(bases);
    }
};

//
// Dam methylation: methyl-adenine in GATC context
// 
struct MethylDamAlphabet : public Alphabet
{
    // member variables, expanded by macrocs
    BASIC_MEMBER_BOILERPLATE
    METHYLATION_MEMBER_BOILERPLATE
    
    // member functions
    BASIC_ACCESSOR_BOILERPLATE
    METHYLATION_ACCESSOR_BOILERPLATE

    // does this alphabet contain all of the nucleotides in bases?
    virtual inline bool contains_all(const char *bases) const 
    {
        return strspn(bases, _base) == strlen(bases);
    }
};

//
// Dcm methylation: methyl-cytosine in CCAGG and CCTGG context
// 
struct MethylDcmAlphabet : public Alphabet
{
    // member variables, expanded by macrocs
    BASIC_MEMBER_BOILERPLATE
    METHYLATION_MEMBER_BOILERPLATE
    
    // member functions
    BASIC_ACCESSOR_BOILERPLATE
    METHYLATION_ACCESSOR_BOILERPLATE

    // does this alphabet contain all of the nucleotides in bases?
    virtual inline bool contains_all(const char *bases) const 
    {
        return strspn(bases, _base) == strlen(bases);
    }
};

// Global alphabet objects that can be re-used
extern DNAAlphabet gDNAAlphabet;
extern MethylCpGAlphabet gMCpGAlphabet;
extern MethylGpCAlphabet gMethylGpCAlphabet;
extern MethylDamAlphabet gMethylDamAlphabet;
extern MethylDcmAlphabet gMethylDcmAlphabet;
extern UtoTRNAAlphabet gUtoTRNAAlphabet;
extern M5CAlphabet gM5CAlphabet;
extern NOMEAlphabet gNOMEAlphabet;

std::vector<const Alphabet*> get_alphabet_list();
const Alphabet* best_alphabet(const char *bases);
const Alphabet* get_alphabet_by_name(const std::string& name);

#endif
