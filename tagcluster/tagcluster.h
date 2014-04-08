#ifndef __TAG_CLUSTER_H__
#define __TAG_CLUSTER_H__

#include <string>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>
#include "iodata.h"
#include "read.h"

struct matepair_t {
	int mate;
	std::string identifier;
	std::string sequence;
	
	matepair_t (int m, const char *pi, std::string s): 
		mate(m), identifier(pi), sequence(s) {};
	matepair_t (int m, std::stringstream &ssi, char *s): 
		mate(m), identifier(ssi.str()), sequence(s) {};

};


//typedef std::pair<std::string, std::string> matepair_t;
typedef std::unordered_map<std::string, std::vector<matepair_t>> tag_map_t;

typedef struct {
	tag_map_t &tag_cluster;
	char *sscs_path;
} data_t;

#endif
