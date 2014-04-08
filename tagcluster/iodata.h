#ifndef __PREMIER_IODATA_H__
#define __PREMIER_IODATA_H__

#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>

#include "read.h"

#define PMR_LINE_BUFFER_SIZE 1024
#define PMR_READ_BUFFER_SIZE 1073741824 //(1<<30) 

#define SKIP_TO_NEWLINE(p)      \
	for(; *(p) != '\n'; (p)++); \
	(p)++;                      \

#define IO_LOCATE_NEWLINE												\
	for(pnr = p; *pnr != '\n' && (pnr - _buffer) < buf_max; pnr++);		\
	if (*(pnr) != '\n') {												\
		size_t partial_read_size = (pnr - pread);						\
		memcpy(_buffer, pread, partial_read_size);						\
		buf_loc = _buffer + partial_read_size;							\
		break;															\
	}																	\

/*
int io_mmap_fastq_file(data *d, const char *filename);
void io_munmap_fastq_file(data *d);
*/
int io_iterate_reads_fastq(const char *fastq_file,
		void (*callback)(read_t *, void *),
		void *fdata);
		
#endif
