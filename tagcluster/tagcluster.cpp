#include <cstring>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include "tagcluster.h"

static int max_cluster_size = 0;

void read_tagging(read_t *read, void *fdata)
{
	data_t *d = (data_t *) fdata;
	tag_map_t &tc = d->tag_cluster;

	char *ridenf = read->identifier;

	char *prid = strtok(ridenf, "/");
	char *ptag = strtok(NULL, "/");
	int mate_id = atoi(ptag + strlen(ptag) + 1);

	std::string rseq(read->sequence, read->length);
	matepair_t mp (mate_id, prid, rseq);
	
	std::string stag(ptag);
	tag_map_t::iterator exist = tc.find(stag);
	if (exist == tc.end()) {
		std::pair<tag_map_t::iterator, bool> pr = tc.insert(
				std::make_pair<std::string, std::vector<matepair_t>>(
					stag, std::vector<matepair_t>()));

		pr.first->second.push_back(mp);
	}
	else {
		exist->second.push_back(mp);

		int vsize = exist->second.size();
		if (vsize > max_cluster_size) {
			max_cluster_size = vsize;
		}
	}
}

void collapse_and_pair(data_t *d, int rlen)
{
	static char bases[5] = "ACTG";
	int sscs_id = 0;

	tag_map_t &tc = d->tag_cluster;
	std::ofstream ofs_sscs_mp1;
	std::ofstream ofs_sscs_mp2;

	std::unordered_map<std::string, matepair_t> r1map;
	std::unordered_map<std::string, matepair_t> r2map;

#if __cplusplus > 199711L
	r1map.reserve(tc.size() / 2);
	r2map.reserve(tc.size() / 2);
#endif

	if (d->sscs_path != NULL) {
		char out_fname[128];
		snprintf(out_fname, 127, "%s_1.fastq", d->sscs_path);
		ofs_sscs_mp1.open(out_fname, std::ofstream::out | std::ofstream::trunc);

		snprintf(out_fname, 127, "%s_2.fastq", d->sscs_path);
		ofs_sscs_mp2.open(out_fname, std::ofstream::out | std::ofstream::trunc);
	}

	char *pseudo_qual = new char[rlen + 1];
	memset(pseudo_qual, 'h', rlen);

	char *sscs = new char[rlen + 1];
	sscs[rlen] = '\0';
	int *nuc_tally = new int[rlen << 4];

	std::cerr << "Collaping reads into single strand consensus sequences (SSCSs) ..." << std::endl;

	// ----- collapse -----
	for (auto it = tc.begin(); it != tc.end(); ++it) {
		std::vector<matepair_t> *mpvec = &(it->second);

		if (mpvec->size() < 6) continue;

		memset(nuc_tally, 0, sizeof(*nuc_tally) * rlen * 8);

		//std::cerr << it->first << " # reads: " << mpvec->size() << std::endl;
		std::stringstream r1ident_ss;
		std::stringstream r2ident_ss;
		for (std::vector<matepair_t>::iterator mpit = mpvec->begin();
				mpit != mpvec->end(); ++mpit) {

			if (mpit->mate == 1) {
				r1ident_ss << mpit->identifier << "/";
			}
			else {
				r2ident_ss << mpit->identifier << "/";
			}

			const char *seq = mpit->sequence.c_str();
			int *t0 = nuc_tally + (mpit->mate - 1) * (rlen << 2);
			for (int i = 0; i < rlen; ++i) {
				int base = (seq[i] >> 1) & 3;
				++t0[(i << 2) + base];
			}
		}

		std::stringstream sr1;
		std::stringstream sr2;
		int n_reliable = 0;

		for (int r = 0; r < 2; r++) {

			int *t0 = nuc_tally + r * (rlen << 2);
			int n_ident = 0;
			for (int i = 0; i < rlen; ++i) {
				int n_nonzero = 0, n_max = 0, n_tie = 1;
				int argmax = 0;
				for (int j = 0; j < 4; ++j) {
					int tally_cnt = t0[(i << 2)+j];
					n_nonzero += (tally_cnt > 0);
					if (tally_cnt > n_max) {
						n_max = tally_cnt;
						n_tie = 1;
						argmax = j;
					}
					else if (tally_cnt == n_max) {
						++n_tie;
					}
				}
				//std::cout << std::endl;

				if (n_max >= 3 && n_tie < 2) {
					sscs[i] = bases[argmax];
				}
				else {
					sscs[i] = 'N';
				}

				n_ident += n_max; 
			}

			//std::cout << sscs << std::endl;

			double seq_ident = (double) n_ident / (rlen * mpvec->size() / 2);

			if (seq_ident > 0.9) {
				++n_reliable;

				matepair_t sscs_mp (r+1, r == 0 ? r1ident_ss : r2ident_ss, sscs);

				std::unordered_map<std::string, matepair_t> &rm = 
					(r == 0) ? r1map : r2map;
				std::string stag(it->first);
				rm.insert(std::make_pair<std::string, matepair_t>(
							stag, sscs_mp));

				if (d->sscs_path) {
					std::stringstream &sread = (r == 0 ? sr1 : sr2);
					sread << "@SSCS" << sscs_id << "/" << 
						it->first << "/" << (r+1) << std::endl 
						<< sscs << std::endl << "+" << std::endl
						<< pseudo_qual << std::endl;
				}
			}
		}

		/* both mate pair can generate reliable SSCS according to identity
		 * threshold (0.9) */
		if (n_reliable == 2) {
			++sscs_id;
			
			if (d->sscs_path) {
				ofs_sscs_mp1 << sr1.str();
				ofs_sscs_mp2 << sr2.str();
			}
		}
	}

	/* release all memory occupied by the original dictionary */
	tc.clear();

	std::cerr << "Pairing all SSCSs with identical tags..." << std::endl;

	// ----- pair -----
	std::cerr << r1map.size() << " / " << r2map.size() << std::endl;
	for (auto it = r1map.begin(); it != r1map.end(); ++it) {
		std::string revtag(it->first.substr(12, 12) + it->first.substr(0, 12));
		auto r2it = r2map.find(revtag);
		if (r2it != r2map.end()) {
			matepair_t *r1mp = &(it->second);
			matepair_t *r2mp = &(r2it->second);

			if (r1mp->sequence == r2mp->sequence) {
				//std::cout << "#" << it->first << "/" << revtag << std::endl;
				std::cout << "@" << r1mp->identifier << "1" << std::endl;
				std::cout << "@" << r2mp->identifier << "2" << std::endl;
				std::cout << r1mp->sequence << std::endl;
			}
		}
	}

#if 0
	for (auto it = tc.begin(); it != tc.end(); ++it) {
		std::vector<matepair_t> *mpvec = &(it->second);
		std::cerr << it->first << " # reads: " << mpvec->size() << std::endl;

		std::string revtag(it->first.substr(12, 12) + it->first.substr(0, 12));
		/*
		for (std::vector<matepair_t>::iterator mpit = mpvec->begin();
				mpit != mpvec->end(); ++mpit) {
			if (mpit->mate == 1) {
				std::cerr << mpit->sequence << std::endl;
			}
		}
		std::cerr << "------[ab]------" << std::endl;
		*/

		auto revit = tc.find(revtag);
		if (revit != tc.end()) {
			std::cout << mpvec->size() << " " << revit->second.size() << std::endl;
		}
		else {
			std::cout << mpvec->size() << " 0" << std::endl;
		}

		/*
		if (revit != tc.end()) {
			std::vector<matepair_t> *revmpvec = &(revit->second);
			std::cerr << revit->first << " # reads: " << revmpvec->size() << std::endl;

			for (std::vector<matepair_t>::iterator revmpit = revmpvec->begin();
					revmpit != revmpvec->end(); ++revmpit) {
				if (revmpit->mate == 2) {
					std::cerr << revmpit->sequence << std::endl;
				}
			}
		}
		std::cerr << "------[ba]------" << std::endl;
		getchar();
		*/
	}
#endif

	delete[] nuc_tally;
	delete[] sscs;
	delete[] pseudo_qual;

	if (d->sscs_path) {
		ofs_sscs_mp1.close();
		ofs_sscs_mp2.close();
	}
}

int main(int argc, char *argv[])
{
	tag_map_t tagmap;
	data_t d = {tagmap, NULL};

	char c;
	while ((c = getopt(argc, argv, "s:")) != -1) {
		switch (c) {
			case 's':
				d.sscs_path = (char *) calloc(strlen(optarg) + 1,
						sizeof(char));
				strncpy(d.sscs_path, optarg, strlen(optarg));
				break;

		}
	}

	io_iterate_reads_fastq(argv[optind], read_tagging, &d);

	std::cerr << "# Unique tags: " << tagmap.size() << std::endl;
	std::cerr << "Max cluster size: " << max_cluster_size << std::endl;

	collapse_and_pair(&d, 80);

	if (d.sscs_path) {
		free(d.sscs_path);
	}
}
