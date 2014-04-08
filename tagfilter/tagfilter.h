#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <deque>
#include <pcrecpp.h>

#define CACHE_LINESIZE 64
#define FASTQ_BUFFER_SIZE 1073741824

#define SKIP_TO_NEWLINE(p)      \
	for(; *(p) != '\n'; (p)++); \
	(p)++;                      \

typedef struct {
	int token;
	int position;
} token_t;

typedef struct {
	/* the fixed sequence after 12nt tag, used for validity. */
	char *fixed_seq;
	char *output_prefix;
	int fixed_seq_length;
	int tag_length;
	/* trim additional bases after 12nt tag and the fixed sequence. */
	int trim_5prime_length;
	/* threshold of consecutive nucleotides in the tag for
	 * discarding. */
	int tag_low_complexity_thresh; 
	FILE *out_fq1;
	FILE *out_fq2;
	/* detect homopolymer */
	pcrecpp::RE *regex;
} data;

void io_iterate_paired_reads(const char *fq1, const char *fq2,
		void (*callback)(char *seq1, char *seq2, char *qscore1, char *qscore2, 
			int rlen1,int rlen2, data *d),
		data *d);
size_t io_load_and_trim_buffer(char *buf, FILE *fp, size_t *buf_max, 
		size_t *off, int fully_load);
size_t io_file_size(const char *fname);
void paired_reads_parse(char *seq1, char *seq2, 
		char *qscore1, char *qscore2, int rlen1, int rlen2, data *d);
