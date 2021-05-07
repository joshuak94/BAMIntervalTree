#include <htslib/sam.h>
#include <thread>

int sam_readrec(BGZF *ignored, void *fpv, void *bv, int *tid, hts_pos_t *beg, hts_pos_t *end)
{
    htsFile *fp = (htsFile *)fpv;
    bam1_t *b = (bam1_t *)bv;
    fp->line.l = 0;
    int ret = sam_read1(fp, fp->bam_header, b);
    if (ret >= 0) {
        *tid = b->core.tid;
        *beg = b->core.pos;
        *end = bam_endpos(b);
    }
    return ret;
}

int htslib_index(const char * input)
{
    return sam_index_build3(input, NULL, 0, std::thread::hardware_concurrency());
}

void htslib_overlap_file_offset(hts_idx_t * index, sam_hdr_t * header, char ** region, int count)
{
    hts_itr_t *iter;
    iter = sam_itr_regarray(index, header, region, count);
}

int htslib_overlap_records(hts_idx_t * index, sam_hdr_t * header, char ** region, int count, htsFile * in, htsFile * out)
{
    int result{};
    bam1_t * b = NULL;
    b = bam_init1();
    result = sam_hdr_write(out, header);
    hts_itr_t * iter = sam_itr_regarray(index, header, region, count);
    if (!iter)
    {
        bam_destroy1(b);
        return result;
    }
    while (sam_itr_next(in, iter, b) >= 0)
    {
        result = sam_write1(out, header, b);
    }

    hts_itr_destroy(iter);
    bam_destroy1(b);

    return result;
}
