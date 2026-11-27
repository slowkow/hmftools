/*
 * The non-JNI part of gatk-bwamem-jni's jnibwa.c, for the Rust ports of hmftools (hmf-common-rs, module
 * bwa): opening a BWA index image (mapped from a file, or held in memory) and aligning a batch of
 * sequences into the packed buffer of alignments that the Java BwaMemAligner decodes.
 *
 * The option block is the mem_opt_t of mem_opt_init(), which the Rust side changes field by field, as
 * the Java aligners change it through a direct ByteBuffer over the same struct.
 *
 * Taken from esvee-rs's esvee_bwa.c (same code, hmf_bwa_ names).
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <errno.h>
#ifdef __linux__
#include <sys/vfs.h>
#endif
#include "bwamem.h"
#include "kstring.h"
#include "bntseq.h"
#include "bwt.h"
#include "utils.h"

static inline void kput32( int32_t val, kstring_t* str ) {
	kputsn((char*)&val, sizeof(int32_t), str);
}

static inline int cigarRefLen( int nCigar, uint32_t const* pCigar )
{
	uint32_t const* pEnd = pCigar + nCigar;
	int len = 0;
	while ( pCigar != pEnd ) {
		uint32_t lenOp = *pCigar++;
		int op = lenOp & 0xf;
		if ( !op || op == 2 )
			len += lenOp >> 4;
	}
	return len;
}

static void fmt_BAMish(mem_opt_t const* opt, bntseq_t const* bns, kstring_t *str, bseq1_t *s, int n, mem_aln_t const* list, int which, mem_aln_t const* p, mem_aln_t const* m) {
	if ( !which ) {
		size_t nInts = 12; // for space planning, assume mapped, unpaired reads with 3 cigar ops and an 8 character MD
		if ( p->flag & 0x1 ) nInts += 3; // if paired, add in enough space for mate info
		ks_resize(str, (n*nInts+1)*sizeof(int32_t));
		kput32(n, str);
	}
	int32_t flag_mapQ = p->flag;
	if ( p->flag & 0x10000 ) flag_mapQ |= 0x100;
	flag_mapQ = (flag_mapQ << 16) | (p->mapq & 0xff);
	kput32(flag_mapQ, str);
	if ( !(p->flag & 0x4) ) {
		kput32(p->rid, str);
		kput32(p->pos, str);
		kput32(p->NM, str);
		kput32(p->score, str);
		kput32(p->sub, str);
		int32_t nCig = p->n_cigar;
		kput32(nCig, str);
		uint32_t* pCig = p->cigar;
		while ( nCig-- ) {
			uint32_t lenOp = *pCig++;
			// op is encoded as MIDSH in a mem_aln_t, but as MIDNSH in a BAM
			if ( (lenOp & 0xf) > 2 ) ++lenOp;
			kput32(lenOp, str);
		}
		int32_t nMD = p->n_cigar ? strlen((char*)pCig) : 0;
		kput32(nMD, str);
		if ( nMD )
			kputsn((char*)pCig, (nMD+3)&~3, str);
		int32_t nXA = p->XA ? strlen(p->XA) : 0;
		kput32(nXA, str);
		if ( nXA ) {
			kputsn(p->XA, (nXA+3)&~3, str);
		}
	}
	if ( (p->flag & 0x9) == 1 ) {
		kput32(m->rid, str);
		kput32(m->pos, str);
		if ( (p->flag & 0x4) || p->rid != m->rid ) kput32(0, str);
		else {
			// get 5' end of + reads, and 3' end of - reads
			long p0 = p->pos;
			if ( p->is_rev ) p0 += cigarRefLen(p->n_cigar, p->cigar) - 1;
			long m0 = m->pos;
			if ( m->is_rev ) m0 += cigarRefLen(m->n_cigar, m->cigar) - 1;
			kput32(m0 - p0 + (p0 > m0 ? -1 : p0 < m0 ? 1 : 0), str);
		}
	}
}

static size_t bufLen( int32_t* pBuf ) {
	size_t totLen = 1;
	size_t nAligns = *pBuf++;
	while ( nAligns-- ) {
		int32_t flag = *pBuf++ >> 16;
		totLen += 1; // for flags
		if ( !(flag & 0x4) ) {
			totLen += 8; // refId, pos, NM, AS, XS, nCigOps, nMDchars, nXAchars
			pBuf += 5;
			int32_t nCig = *pBuf++;
			totLen += nCig;
			pBuf += nCig;
			int32_t nMD = (*pBuf++ + 3) >> 2;
			totLen += nMD;
			pBuf += nMD;
			int32_t nXA = (*pBuf++ + 3) >> 2;
			totLen += nXA;
			pBuf += nXA;
		}
		if ( (flag & 0x9) == 1 ) {
			totLen += 3; // mate rid, mate pos, tlen
			pBuf += 3;
		}
	}
	return totLen;
}

/* 1 if the file is on a FUSE (a cloud bucket mounted with gcsfuse, as on Google Batch) or NFS file
   system, where each page of a mapped file is fetched by its own small network read */
int hmf_bwa_on_network_fs( const char* path ) {
#ifdef __linux__
	struct statfs fs;
	if ( statfs(path, &fs) == 0 && (fs.f_type == 0x65735546 /* FUSE */ || fs.f_type == 0x6969 /* NFS */) ) return 1;
#endif
	(void)path;
	return 0;
}

static void set_formatter( void ) {
	mem_fmt_fnc = &fmt_BAMish;
	bwa_verbose = 0;
}

/* jnibwa_openIndex, taking a path: maps the image and sets the result formatter */
bwaidx_t* hmf_bwa_open_index( const char* path ) {
	int fd = open(path, O_RDONLY);
	if ( fd == -1 ) return 0;
	struct stat statBuf;
	if ( fstat(fd, &statBuf) == -1 ) { close(fd); return 0; }
	uint8_t* mem = mmap(0, statBuf.st_size, PROT_READ, MAP_SHARED, fd, 0);
	close(fd);
	if ( mem == MAP_FAILED ) return 0;
	bwaidx_t* pIdx = calloc(1, sizeof(bwaidx_t));
	bwa_mem2idx(statBuf.st_size, mem, pIdx);
	pIdx->is_shm = 1;
	set_formatter();
	return pIdx;
}

/* an index image held in memory by the caller (8-byte aligned), which must outlive the index */
bwaidx_t* hmf_bwa_open_index_mem( uint8_t* mem, int64_t size ) {
	bwaidx_t* pIdx = calloc(1, sizeof(bwaidx_t));
	bwa_mem2idx(size, mem, pIdx);
	pIdx->is_shm = 1;
	set_formatter();
	return pIdx;
}

/* jnibwa_destroyIndex; the memory of an image from hmf_bwa_open_index_mem stays with the caller */
int hmf_bwa_destroy_index( bwaidx_t* pIdx, int mapped ) {
	void* pMem = pIdx->mem;
	size_t memLen = pIdx->l_mem;
	bwa_idx_destroy(pIdx);
	return mapped ? munmap(pMem, memLen) : 0;
}

int hmf_bwa_contig_count( bwaidx_t* pIdx ) { return pIdx->bns->n_seqs; }
const char* hmf_bwa_contig_name( bwaidx_t* pIdx, int i ) { return pIdx->bns->anns[i].name; }
int hmf_bwa_contig_length( bwaidx_t* pIdx, int i ) { return pIdx->bns->anns[i].len; }

/* mem_opt_init(): the option block BwaMemAligner starts from */
mem_opt_t* hmf_bwa_default_opts( void ) { return mem_opt_init(); }
void hmf_bwa_destroy_opts( mem_opt_t* o ) { free(o); }
/* bwa_fill_scmat: the 5x5 matrix from a match score, a mismatch penalty and -1 for any N */
void hmf_bwa_fill_scmat( int a, int b, int8_t mat[25] ) { bwa_fill_scmat(a, b, mat); }

/* the layout check of the Rust mirror of mem_opt_t: its size and the offsets of a few fields */
size_t hmf_bwa_opts_layout( int which ) {
	switch ( which ) {
		case 0: return sizeof(mem_opt_t);
		case 1: return offsetof(mem_opt_t, max_mem_intv);
		case 2: return offsetof(mem_opt_t, n_threads);
		case 3: return offsetof(mem_opt_t, drop_ratio);
		case 4: return offsetof(mem_opt_t, mat);
		default: return 0;
	}
}

/*
 * jnibwa_createAlignments: pSeq is a count followed by that many NUL-terminated sequences; the result is a
 * malloc'ed buffer of 32-bit values (see fmt_BAMish) to be released with hmf_bwa_free.
 */
void* hmf_bwa_create_alignments( bwaidx_t* pIdx, mem_opt_t* pOpts, char* pSeq, size_t* pBufSize ) {
	char c = 0;
	char* emptyString = &c;
	uint32_t nSeqs = *(uint32_t*)pSeq;
	pSeq += sizeof(uint32_t);
	bseq1_t* pSeq1Beg = calloc(nSeqs, sizeof(bseq1_t));
	bseq1_t* pSeq1End = pSeq1Beg+nSeqs;
	bseq1_t* pSeq1;
	for ( pSeq1 = pSeq1Beg; pSeq1 != pSeq1End; ++pSeq1 ) {
		size_t seqLen = strlen(pSeq);
		pSeq1->l_seq = seqLen;
		pSeq1->seq = pSeq;
		pSeq1->name = emptyString;
		pSeq1->id = pSeq1-pSeq1Beg;
		pSeq += seqLen + 1;
	}
	mem_process_seqs(pOpts, pIdx->bwt, pIdx->bns, pIdx->pac, 0, nSeqs, pSeq1Beg, 0);
	size_t nInts = 0;
	for ( pSeq1 = pSeq1Beg; pSeq1 != pSeq1End; ++pSeq1 ) {
		if ( pSeq1->sam ) nInts += bufLen((int32_t*)pSeq1->sam);
	}
	int32_t* resultsBeg = malloc(nInts*sizeof(int32_t));
	int32_t* pResults = resultsBeg;
	for ( pSeq1 = pSeq1Beg; pSeq1 != pSeq1End; ++pSeq1 ) {
		int32_t* pBuf = (int32_t*)pSeq1->sam;
		if ( pBuf ) {
			size_t len = bufLen(pBuf);
			memcpy(pResults, pBuf, len*sizeof(int32_t));
			free(pBuf);
			pResults += len;
		}
	}
	free(pSeq1Beg);
	*pBufSize = nInts*sizeof(int32_t);
	return resultsBeg;
}

void hmf_bwa_free( void* p ) { free(p); }
