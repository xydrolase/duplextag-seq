#include <string.h>
#include <stdlib.h>
#include "iodata.h"
#include "read.h"

int io_iterate_reads_fastq(const char *fastq_file,
		void (*callback)(read_t *read, void *fdata), void *fdata)
{
	static char buf_ident[PMR_LINE_BUFFER_SIZE];

	int fd;
	struct stat fs;
	struct timeval tv;
	int64_t fastq_size;
	size_t size_ident;
	int read_idx = 0, n = 15;
	char *pnr = NULL;

	gettimeofday(&tv, NULL);

	if ((fd = open(fastq_file, O_RDONLY)) == -1){
		return -1;
	}
	else {
		if (fstat(fd, &fs) != -1)
			fastq_size = (int64_t) fs.st_size;
		else 
			return -1;
	}

	close(fd);
	
	/* create a temporary read */
	read_t *tmp_read = read_create();
	tmp_read->identifier = buf_ident;

	int chunk = 0;
	int chunks_to_read = (int) ceil((double) fastq_size / PMR_READ_BUFFER_SIZE);

	char *_buffer = (char *) malloc(
			fastq_size < PMR_READ_BUFFER_SIZE ? fastq_size :
			PMR_READ_BUFFER_SIZE);
	char *buf_loc = _buffer;

	/* process in blocks */
	FILE *fp = fopen(fastq_file, "r");
	while (fastq_size > 0) {
		size_t chunk_size = fastq_size < PMR_READ_BUFFER_SIZE ? fastq_size :
			PMR_READ_BUFFER_SIZE - (buf_loc - _buffer);

		size_t buf_size = fastq_size < PMR_READ_BUFFER_SIZE ? 
			fastq_size + (buf_loc - _buffer) :
			PMR_READ_BUFFER_SIZE;

		register size_t buf_max = buf_size - 1;

		char *p = _buffer;
		fread(buf_loc, chunk_size, 1, fp);
		fprintf(stderr, 
				"Reading in FASTQ data in chunk [%2d/%2d]...\n", ++chunk, 
				chunks_to_read);

		/* locate the last read in buffer */
		//size_t rem_buf_size = buf_size;
		while ((p - _buffer) < buf_size) {
			char *pread = p;
			buf_loc = _buffer;
			tmp_read->id = read_idx++;

			IO_LOCATE_NEWLINE;

			/* process header */
			char *p_ident = p + 1;
			for (; *p_ident != ' ' && p_ident <= pnr; p_ident++);
			size_ident = (size_t) (p_ident - p - 1);
			strncpy(buf_ident, p + 1, size_ident);
			buf_ident[size_ident] = '\0';

			/* read sequence */
			p = pnr + 1;
			tmp_read->sequence = p;

			IO_LOCATE_NEWLINE;

			int read_length = pnr - p;
			tmp_read->length = read_length;

			/* skip quality header line */
			p = pnr + 1;
			if (*p == '+') {
				IO_LOCATE_NEWLINE;
			}

			/* quality score */
			p = pnr + 1;
			tmp_read->qscore = p;
			IO_LOCATE_NEWLINE;

			/* invoke callback function */
			callback(tmp_read, fdata);

			if (read_idx >> n) {
				struct timeval _ntv;
				gettimeofday(&_ntv, NULL);

				fprintf(stderr, 
						"Processed %9d reads, time elapsed: %.3f\n", 
						read_idx,
						(0.0 + _ntv.tv_sec + (double) _ntv.tv_usec/1e6) - 
						(0.0 + tv.tv_sec + (double) tv.tv_usec/1e6));

				n++;
			}

			p = pnr + 1;
		}

		fastq_size -= buf_size;
	}

	struct timeval _ntv;
	gettimeofday(&_ntv, NULL);
	fprintf(stderr, 
			"All (%9d) reads processed in %.3f seconds.\n",
			read_idx, 
			(0.0 + _ntv.tv_sec + (double) _ntv.tv_usec/1e6) - 
			(0.0 + tv.tv_sec + (double) tv.tv_usec/1e6));

	read_destroy(tmp_read);

	fclose(fp);
	free(_buffer);
}
