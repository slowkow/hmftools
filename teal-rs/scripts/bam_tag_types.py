#!/usr/bin/env python3
"""usage: bam_tag_types.py <bam>: per record, its attributes as TAG:type (the binary type codes SAM text hides)"""
import gzip, struct, sys
data = gzip.open(sys.argv[1]).read()
p = 4
l_text = struct.unpack_from('<i', data, p)[0]; p += 4 + l_text
n_ref = struct.unpack_from('<i', data, p)[0]; p += 4
for _ in range(n_ref):
    l = struct.unpack_from('<i', data, p)[0]; p += 4 + l + 4
sizes = {'A': 1, 'c': 1, 'C': 1, 's': 2, 'S': 2, 'i': 4, 'I': 4, 'f': 4}
while p < len(data):
    bs = struct.unpack_from('<i', data, p)[0]
    r = data[p + 4:p + 4 + bs]; p += 4 + bs
    l_name, n_cig, l_seq = r[8], struct.unpack_from('<H', r, 12)[0], struct.unpack_from('<i', r, 16)[0]
    q = 32 + l_name + 4 * n_cig + (l_seq + 1) // 2 + l_seq
    out = []
    while q < len(r):
        tag, t = r[q:q + 2].decode(), chr(r[q + 2]); q += 3
        if t in sizes: q += sizes[t]
        elif t in 'ZH': q = r.index(b'\0', q) + 1
        elif t == 'B':
            st, n = chr(r[q]), struct.unpack_from('<i', r, q + 1)[0]; q += 5 + n * sizes[st]
        out.append(f"{tag}:{t}")
    print(r[32:32 + l_name - 1].decode(), " ".join(out))
