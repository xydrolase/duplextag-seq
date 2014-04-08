#include "tagfilter.h"

size_t io_file_size(const char *fname) 
{
	int fd;
	struct stat fs;
	size_t file_size;
	if ((fd = open(fname, O_RDONLY)) == -1){
		return 0;
	}
	else {
		if (fstat(fd, &fs) != -1)
			file_size = fs.st_size;
		else {
			return 0;
		}
	}

	close(fd);

	return file_size;
}

size_t io_load_and_trim_buffer(char *buf, FILE *fp, size_t *buf_max, 
		size_t *off, int fully_load)
{
	size_t buf_size = (*buf_max + 1);
	fseek(fp, *off, SEEK_SET);
	fread(buf, buf_size, 1, fp);

	if (fully_load) {
		/* do not need to trim */
		return buf_size;
	}

	/* trimming: looking for @/+ or +/@ pairs */
	int trim_pos = 0;
	std::deque<token_t> global_tokens;
	for (int i = buf_size; i >= 0; i -= CACHE_LINESIZE) {
		std::deque<token_t> local_tokens;
		for (int j = i - CACHE_LINESIZE; j < i; j++) {
			if (buf[j] == '\n') {
				while (!global_tokens.empty()) {
					token_t &ftk = global_tokens.front();
					if (ftk.token != '\n' && 
							ftk.position - j > 1) {
						global_tokens.pop_front();
					}
					else {
						break;
					}
				}

				token_t _tk = {buf[j], j};
				local_tokens.push_back(_tk);
			}
			else if (buf[j] == '@' || buf[j] == '+') {
				if ((j > i - CACHE_LINESIZE) && buf[j-1] == '\n') {
					token_t _tk = {buf[j], j};
					local_tokens.push_back(_tk);
				}
				else if (j == i - CACHE_LINESIZE) {
					token_t _tk = {buf[j], j};
					local_tokens.push_back(_tk);
				}
			}
		}

		global_tokens.insert(global_tokens.begin(), local_tokens.begin(),
				local_tokens.end());

		/* process all tokens */
		char *token_series = new char[global_tokens.size() + 1];
		token_series[global_tokens.size()] = '\0';
		std::deque<token_t>::iterator it = global_tokens.begin();
		for (int i = 0; it != global_tokens.end(); it++, i++) {
			token_series[i] = it->token;	
		}

		const char *pos;
		if ((pos = strstr(token_series, "@\n\n+")) != NULL) {
			trim_pos = global_tokens[pos - token_series].position;
			break;
		}
		else if ((pos = strstr(token_series, "+\n\n@")) != NULL) {
			trim_pos = global_tokens[pos - token_series + 3].position;
			break;
		}

		delete token_series;
	}

	*off += trim_pos;
	*buf_max = trim_pos;

	return trim_pos;
}

void io_iterate_paired_reads(const char *fq1, const char *fq2,
		void (*callback)(char *seq1, char *seq2, char *qscore1, char *qscore2, 
			int rlen1,int rlen2, data *d),
		data *d)
{
	size_t fq1_size = io_file_size(fq1);
	size_t fq2_size = io_file_size(fq2);

	/* offset in FASTQ files for reading */
	size_t fq1_off = 0;
	size_t fq2_off = 0;

	/* allocate buffers */
	size_t buf_max_1 = (fq1_size < FASTQ_BUFFER_SIZE ? 
		fq1_size : FASTQ_BUFFER_SIZE) - 1;
	size_t buf_max_2 = (fq2_size < FASTQ_BUFFER_SIZE ? 
		fq2_size : FASTQ_BUFFER_SIZE) - 1;

	char *buf1, *buf2;
	posix_memalign((void **) &buf1, CACHE_LINESIZE, buf_max_1 + 1);
	posix_memalign((void **) &buf2, CACHE_LINESIZE, buf_max_2 + 1);

	char *p1 = buf1 + buf_max_1 + 1, *p2 = buf2 + buf_max_2 + 1;
	size_t off1 = 0, off2 = 0;
	FILE *fp1 = fopen(fq1, "r"), *fp2 = fopen(fq2, "r");

	while (fq1_size > 0 && fq2_size > 0) {
		/* reload buffer if necessary */
		if ((p1 - buf1) >= buf_max_1) {
			buf_max_1 = (fq1_size < FASTQ_BUFFER_SIZE ? 
					fq1_size : FASTQ_BUFFER_SIZE) - 1;
			fq1_size -= io_load_and_trim_buffer(buf1, fp1, &buf_max_1, &off1,
					fq1_size < FASTQ_BUFFER_SIZE);
			fprintf(stderr, "[%s] Remaining data to process: %.2f MB\n",
					fq1, (double) fq1_size / 1024 / 1024);
			p1 = buf1;
		}

		if ((p2 - buf2) >= buf_max_2) {
			buf_max_2 = (fq1_size < FASTQ_BUFFER_SIZE ? 
					fq1_size : FASTQ_BUFFER_SIZE) - 1;
			fq2_size -= io_load_and_trim_buffer(buf2, fp2, &buf_max_2, &off2,
					fq1_size < FASTQ_BUFFER_SIZE); 
			fprintf(stderr, "[%s] Remaining data to process: %.2f MB\n",
					fq2, (double) fq1_size / 1024 / 1024);
			p2 = buf2;
		}

		/* within-buffer data crunching */
		while ((p1 - buf1) < buf_max_1 &&
				(p2 - buf2) < buf_max_2) {
			/* fastq header */
			SKIP_TO_NEWLINE(p1);
			SKIP_TO_NEWLINE(p2);

			/* read sequence */
			char *pseq1 = p1, *pseq2 = p2;
			SKIP_TO_NEWLINE(p1);
			SKIP_TO_NEWLINE(p2);
			
			/* qscore header */
			SKIP_TO_NEWLINE(p1);
			SKIP_TO_NEWLINE(p2);
			
			/* don't care about qscore either */
			char *pqs1 = p1, *pqs2 = p2;
			SKIP_TO_NEWLINE(p1);
			SKIP_TO_NEWLINE(p2);

			callback(pseq1, pseq2, pqs1, pqs2, 
					p1 - pqs1 - 1, p2 - pqs2 - 1, d);
		}
	}

	free(buf1);
	free(buf2);
}

void paired_reads_parse(char *seq1, char *seq2, 
		char *qscore1, char *qscore2, int rlen1, int rlen2, data *d)
{
	static int read_id = 1;

	if (memchr(seq1, 'N', rlen1) != NULL || memchr(seq1, '.', rlen1) != NULL) return;
	if (memchr(seq2, 'N', rlen2) != NULL || memchr(seq2, '.', rlen2) != NULL) return;

	/* check for fixed sequence identity */
	if (strncmp(seq1 + d->tag_length, d->fixed_seq, 
				d->fixed_seq_length) != 0 ||
			strncmp(seq2 + d->tag_length, d->fixed_seq,
				d->fixed_seq_length) != 0) {
		return;
	}

	char *duplex_tag = new char[(d->tag_length << 1) + 1];
	duplex_tag[(d->tag_length << 1)] = '\0';
	memcpy(duplex_tag, seq1, d->tag_length);
	strncpy(duplex_tag + d->tag_length, seq2, d->tag_length);

	/* check for tag complexity */
	if (d->regex->PartialMatch(duplex_tag)) return;

	/* output trimmed sequences */
	int off = d->tag_length + d->fixed_seq_length + 
		d->trim_5prime_length;
	int trimmed_len = rlen1 - off;
	fprintf(d->out_fq1, "@%s.%d/%s\n%.*s\n+%s.%d\n%.*s\n",
			d->output_prefix, read_id, duplex_tag, 
			trimmed_len, seq1 + off, d->output_prefix, read_id, 
			trimmed_len, qscore1 + off);

	off = d->tag_length + d->fixed_seq_length + 
		d->trim_5prime_length;
	trimmed_len = rlen2 - off;
	fprintf(d->out_fq2, "@%s.%d/%s\n%.*s\n+%s.%d\n%.*s\n",
			d->output_prefix, read_id, duplex_tag,
			trimmed_len, seq2 + off, d->output_prefix, read_id, 
			trimmed_len, qscore2 + off);

	++read_id;

	delete duplex_tag;
}

void print_usage(char *synopsis)
{
	std::cout << "Usage: " << synopsis << " -p <prefix> [options] <p1.fastq> <p2.fastq>" 
		<< std::endl << std::endl;
	std::cout << "Options: " << std::endl;
	std::cout << "  -c INT  Low complexity threshold for duplex tags, " << std::endl 
		      << "          measured by a homopolymer of length <INT>." << std::endl;
	std::cout << "  -f STR  Fixed sequence (immediately appended after tag) " << std::endl 
		      << "          used for quality control." << std::endl;
	std::cout << "  -t INT  Length of tag sequence." << std::endl;
	std::cout << "  -T INT  Number of nucletides to trim after tag and" << std::endl
		      << "          the fixed sequence for quality control." << std::endl;
}

int main(int argc, char *argv[])
{
	data d = {"CAGTA", NULL, 5, 12, 4, 10, NULL, NULL, NULL};

	int arg_len;
	char *ptr_optarg;
	char c;
	while ((c = getopt (argc, argv, "hp:c:f:t:T:")) != -1) {
		switch (c) {
			case 'h':
				print_usage(argv[0]);
				exit(0);
				break;
			case 'p':
				arg_len = strlen(optarg);
				d.output_prefix = (char *) calloc(1, arg_len + 1);
				strncpy(d.output_prefix, optarg, arg_len);
				break;
			case 'c':
				d.tag_low_complexity_thresh = strtol(optarg, &ptr_optarg, 0);
				if (ptr_optarg == optarg) exit(1);
				break;
			// fix sequence
			case 'f':
				arg_len = strlen(optarg);
				d.fixed_seq_length = arg_len;
				d.fixed_seq = (char *) calloc(1, arg_len + 1);
				strncpy(d.fixed_seq, optarg, arg_len);
				break;
			// tag length
			case 't':
				d.tag_length = strtol(optarg, &ptr_optarg, 0);
				if (ptr_optarg == optarg) exit(1);
				break;
			case 'T':
				d.trim_5prime_length = strtol(optarg, &ptr_optarg, 0);
				if (ptr_optarg == optarg) exit(1);
				break;
		}
	}

	if (d.output_prefix == NULL) {
		print_usage(argv[0]);
		exit(1);
	}

	char out_fname[128];
	char homopol_pattern[32];
	snprintf(homopol_pattern, 31, "([ACGT])\\1{%d,}+", 
			d.tag_low_complexity_thresh-1);

	pcrecpp::RE regex_homo(homopol_pattern);
	d.regex = &regex_homo;

	snprintf(out_fname, 127, "%s_1.fastq", d.output_prefix);
	d.out_fq1 = fopen(out_fname, "w+");
	snprintf(out_fname, 127, "%s_2.fastq", d.output_prefix);
	d.out_fq2 = fopen(out_fname, "w+");

	fprintf(stderr, "Fixed sequence: %s\n", d.fixed_seq);
	fprintf(stderr, "Tag length: %d\n", d.tag_length);
	fprintf(stderr, "5' Trim length: %d\n", d.trim_5prime_length);
	fprintf(stderr, "Total trim length: %d\n", 
			d.fixed_seq_length + d.tag_length + d.trim_5prime_length);
	std::cerr << "Homopolymer RegEx pattern: " << d.regex->pattern() << std::endl;

	io_iterate_paired_reads(argv[optind], argv[optind+1], paired_reads_parse, &d);

	fclose(d.out_fq1);
	fclose(d.out_fq2);
}
