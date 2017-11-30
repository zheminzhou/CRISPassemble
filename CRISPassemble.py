# CRISPassemble

import sys, os, subprocess
from time import time
from json import dump, load

nucmer = '/home/zhemin/bin/nucmer --maxmatch -l 7 -c 11'
nucmer2 = '/home/zhemin/bin/nucmer --maxmatch -l 5 -c 7'
show_coords = '/home/zhemin/bin/show-coords -lqHT'
cap3 = '/home/zhemin/software/CAP3/cap3'
min_spacer = 20
max_spacer = 80

def extendEnd(part, maxdiff=3, edgeAware=False) :
    if part[2] < part[3] :
        end_diff = [min(part[0]-1, part[2]-1), min(part[7]-part[1], part[8]-part[3])]
        if sum(end_diff) <= maxdiff :
            part[0] -= end_diff[0]
            part[2] -= end_diff[0]
            part[1] += end_diff[1]
            part[3] += end_diff[1]
        else :
            end_diff = [-1, -1]
    else :
        end_diff = [min(part[0]-1, part[8]-part[2]), min(part[7]-part[1], part[3]-1)]
        if sum(end_diff) <= maxdiff :
            part[0] -= end_diff[0]
            part[2] += end_diff[1]
            part[1] += end_diff[1]
            part[3] -= end_diff[0]
        else :
            end_diff = [-1, -1]
    return end_diff

def iter_region(ends, i, j, conn = []) :
    res = []
    for i2 in i :
        for j2 in j :
            if ends[1][j2][1] - ends[0][i2][0] + 1 >= min_spacer and ends[1][j2][1] - ends[0][i2][0] + 1 <= max_spacer and ends[0][i2][1] - ends[1][j2][0] >= 7 :
                m1, m2 = len(ends[0][i2][2]), len(ends[1][j2][2])
                if float(m1*m2)/(m1+m2) < 1 :
                    continue
                n1, n2 = sum([x[5] for x in ends[0][i2][2]]), sum([x[5] for x in ends[1][j2][2]])
                c = [i2, j2, 2.*float(n1*n2)/(n1+n2)]
                res.extend(iter_region(ends, i - set([i2]), j - set([j2]), conn + [c]))
    if len(res) == 0 :
        res = [conn[:]]
    return res
def parse_spacers(save) :
    ends = [{}, {}]
    save.sort(key=lambda x:x[5], reverse=True)
    for part in save :
        if part[9].endswith('L') or part[9].endswith('W') :
            s = part[2]
            if s not in ends[0] :
                ends[0][s] = [part[2], part[3], [part]]
            else :
                if part[3] > ends[0][s][1] : ends[0][s][1] = part[3]
                ends[0][s][2].append(part)
        if part[9].endswith('R') or part[9].endswith('W') :
            e = part[3]
            if e not in ends[1] :
                ends[1][e] = [part[2], part[3], [part]]
            else :
                if part[2] < ends[1][e][0] : ends[1][e][0] = part[2]
                ends[1][e][2].append(part)
    ends = [ends[0].values(), ends[1].values()]
    regions = iter_region(ends, set(range(len(ends[0]))), set(range(len(ends[1]))), [])
    region = max([ [sum([r[2] for r in region]), region] for region in regions])[1]
    if len(regions) > 0 :
        ce = { e[0]:e[1:] for r in region for e in ((ends[0][r[0]][0], 0, r[0]), (ends[1][r[1]][1], 1, r[1])) }
        for id, end in enumerate(ends) :
            for e in end :
                if ce.get(e[id], [2])[0] != id :
                    for k in (-1, 1, -2, 2) :
                        if ce.get(e[id]+k, [2])[0] == id :
                            for ee in e[2] :
                                if ee[1] - ee[0] + 1 >= 11 :
                                    ee[id], ee[id+2], ee[5] = ee[id]+k, ee[id+2]+k, ee[5] - 4*abs(k)
                                    end[ce.get(e[id]+k)[1]][2].append(ee)
                            break
    regions = [ [save[0][10], ends[0][s][0], ends[1][e][1], depth, \
                 sorted(ends[0][s][2], key=lambda x:(x[9], x[0], x[1])), \
                 sorted(ends[1][e][2], key=lambda x:(x[9], x[0], x[1]))] for s, e, depth in region ]
    return regions


def parse_region(save) :
    regions = []
    save.sort(key=lambda x:x[5], reverse=True)
    for part in save :
        in_known_region = 0
        for reg in regions :
            in_s = max(part[2], reg[0])
            in_e = min(part[3], reg[1])
            if in_e - in_s+1 >= 0.4*(part[3]-part[2]+1) :
                in_known_region = 1
                if in_e - in_s + 1 >= 0.9 * (reg[1] - reg[0] + 1) and part[4] < reg[2][4] + 0.9 :
                    if part[2] < reg[0] : reg[0] = part[2]
                    if part[3] > reg[1] : reg[1] = part[3]
                    reg.append(part)
                break
        if not in_known_region :
            regions.append([part[2], part[3], part])
    for reg in regions :
        s, e = {}, {}
        for r in reg[2:] :
            s[r[2]], e[r[3]] = s.get(r[2], 0) + 1, e.get(r[3], 0) + 1
        reg[0] = int(max(s.iteritems(), key=lambda x:[x[1], -x[0]])[0])
        reg[1] = int(max(e.iteritems(), key=lambda x:[x[1], x[0]])[0])
    return {save[0][10]:sorted(regions)}

complement = {'A':'T', 'G':'C', 'C':'G', 'T':'A', 'N':'N'}
def rc(bases) :
    return ''.join([complement.get(b, 'N') for b in bases.upper()[::-1]])

def identify_positive_reads(prefix, DR, readfiles) :
    positive_region_file = '{0}.spacer.reads'.format(prefix)
    positive_link_file = '{0}.DR.links'.format(prefix)

    if debug :
        return positive_region_file, positive_link_file

    if os.path.isfile(positive_region_file) :
        os.unlink(positive_region_file)

    ori_t = time()
    if not debug :
        for id, read in enumerate(readfiles) :
            output_prefix = '{0}.{1}'.format(prefix, id)
            nucmer_cmd = '{nucmer} -p {2} {0} {1}'.format(DR, read, output_prefix, nucmer=nucmer)
            subprocess.Popen(nucmer_cmd.split(), stdout=subprocess.PIPE).wait()


    links = []
    for id, read in enumerate(readfiles) :
        reads = {}
        output_prefix = '{0}.{1}'.format(prefix, id)
        coords_cmd = '{show_coords} {0}.delta'.format(output_prefix, show_coords=show_coords)
        coords = subprocess.Popen(coords_cmd.split(), stdout=subprocess.PIPE)
        save = []
        for line in iter(coords.stdout.readline, r'') :
            part = line.strip().split()
            part[:9] = [float(x) for x in part[:9]]
            if part[4] < 14 :
                continue
            end_diff = sum(extendEnd(part))
            if end_diff >= 0 :
                m = max(part[4], part[5])
                part[6] = m*part[6]/(m+end_diff)
                if part[6] >= 85.0 :
                    part[4] = (m + end_diff) * (100. - part[6])/100.
                    part[5] = m * 2 - 3*part[4] + end_diff
                    if part[2] > part[3] :
                        part[:4] = [-part[1], -part[0], part[3], part[2]]
                    if len(save) > 0 and save[0][10] != part[10] :
                        reads.update(parse_region(save))
                        save = []
                    save.append(part)
        if len(save) > 0 :
            reads.update(parse_region(save))

        # get read sequences
        read_seq = {}
        with open(read) as fin:
            for line in fin :
                if line[0] == '>' :
                    name = line[1:].strip().split()[0]
                    if name in reads :
                        read_seq[name] = []
                elif name in read_seq :
                    read_seq[name].append(line.strip())
        for name in read_seq :
            read_seq[name] = ''.join(read_seq[name])

        # remove regions with only spacer but no DR
        with open(positive_region_file, 'a') as fout :
            for read, matches in reads.iteritems() :
                prev = 0
                for match in matches :
                    if match[2][0] > 0 :
                        direction = '+'
                        match[2:] = ['DR', '+', ''] # read_seq[read][(match[0]-1):match[1]]]
                    else :
                        direction = '-'
                        match[2:] = ['DR', '-', ''] # rc(read_seq[read][(match[0]-1):match[1]])]
                    if match[0] - prev - 1 >= 7 :
                        if prev > 0 :
                            if direction == '+' :
                                fout.write('>{0}:{1}:{2}:W\n{3}\n'.format(read, prev+1, match[0]-1, read_seq[read][prev:(match[0]-1)]))
                            else :
                                fout.write('>{0}:{2}:{1}:W\n{3}\n'.format(read, prev+1, match[0]-1, rc(read_seq[read][prev:(match[0]-1)])))
                        else :
                            if direction == '+' :
                                fout.write('>{0}:{1}:{2}:R\n{3}\n'.format(read, prev+1, match[0]-1, read_seq[read][prev:(match[0]-1)]))
                            else :
                                fout.write('>{0}:{2}:{1}:L\n{3}\n'.format(read, prev+1, match[0]-1, rc(read_seq[read][prev:(match[0]-1)])))
                    prev = match[1]
                if len(read_seq[read]) - prev >= 7 :
                    if direction == '+' :
                        fout.write('>{0}:{1}:{2}:L\n{3}\n'.format(read, prev+1, len(read_seq[read]), read_seq[read][prev:]))
                    else :
                        fout.write('>{0}:{2}:{1}:R\n{3}\n'.format(read, prev+1, len(read_seq[read]), rc(read_seq[read][prev:])))
                links.append([read, read_seq[read]] + [match for match in matches if match[2] == 'DR'])
    dump(links, open(positive_link_file, 'w'), indent=2)
    return positive_region_file, positive_link_file

def get_crispr(positive_region_file, positive_link_file) :
    # CAP3 assembly for spacers
    cap_run = subprocess.Popen('{0} {1} -a 11 -f 3 -i 21 -j 31 -o 16 -p 90 -s 251 -z1'.format(cap3, positive_region_file).split(), stdout=subprocess.PIPE)
    cap_run.communicate()

    asm = '{0}.cap.contigs'.format(positive_region_file)

    nucmer_cmd = '{nucmer} -p {0} {1} {0}'.format(asm, positive_region_file, nucmer=nucmer2)
    subprocess.Popen(nucmer_cmd.split(), stdout=subprocess.PIPE)
    coords_cmd = '{show_coords} {0}.delta'.format(asm, show_coords=show_coords)
    coords = subprocess.Popen(coords_cmd.split(), stdout=subprocess.PIPE)

    spacers, save = [], []
    for line in iter(coords.stdout.readline, r'') :
        part = line.strip().split('\t')
        part[:9] = [float(x) for x in part[:9]]
        if part[2] > part[3] or part[4] < 7 : continue
        if part[9].endswith('L') or part[9].endswith('W') :
            dist, d2 = part[0]-1, min(part[7]-part[1], part[8] - part[3])
            if (part[4] < 11 and max(dist, d2) > 0) or dist > 2 or part[2] - 1 < dist or d2 > 2 : continue
            part[0], part[2] = part[0] - dist, part[2] - dist
        if part[9].endswith('R') or part[9].endswith('W') :
            dist, d2 = part[7]-part[1], min(part[0]-1, part[2]-1)
            if (part[4] < 11 and max(dist, d2) > 0) or dist > 2 or part[8]-part[3] < dist or d2 > 2 : continue
            part[1], part[3] = part[1] + dist, part[3] + dist
        m1 ,m2 = max(part[4], part[5]), max(part[1]-part[0]+1, part[3]-part[2]+1)
        part[6] = m1*part[6]/m2
        if part[6] >= 80.0 :
            part[4] = m2 * (100. - part[6])/100.
            part[5] = (m1 * 2 - 3*part[4] + (m2 - m1))
            if len(save) > 0 and save[0][10] != part[10] :
                spacers.extend(parse_spacers(save))
                save = []
            save.append(part)
    if len(save) > 0 :
        spacers.extend(parse_spacers(save))
    spacer2 = []
    while len(spacers) > 0 :
        for spacer in spacers :
            m, n = sum([x[5] for x in spacer[4]]), sum([x[5] for x in spacer[5]])
            spacer[3] = m*n/float(m+n)

        s_id, s = max(enumerate(spacers), key=lambda x:x[1][3])
        spacer2.append(s)
        del spacers[s_id]
        regions = {x[9]:x[:2] for x in s[4] + s[5]}
        for s_id in xrange(len(spacers)-1, -1, -1) :
            spacer = spacers[s_id]
            overlap_cnt = [[0, 0], [0, 0]]
            for ss, sn in zip((spacer[4], spacer[5]), overlap_cnt) :
                for id in xrange(len(ss)-1, -1, -1) :
                    x = ss[id]
                    sn[1] += x[5]
                    if x[9] in regions :
                        in_s = max(regions[x[9]][0], x[0])
                        in_e = min(regions[x[9]][1], x[1])
                        if (in_e - in_s+1) >= 0.4 * (x[1]-x[0]+1) :
                            sn[0] += x[5]
                            del ss[id]
            if overlap_cnt[0][0]*2 >= overlap_cnt[0][1] or overlap_cnt[1][0]*2 >= overlap_cnt[1][1] :
                del spacers[s_id]
    spacers = spacer2
    # get contig sequences
    cont_seq = {}
    with open(asm) as fin:
        for line in fin :
            if line[0] == '>' :
                name = line[1:].strip().split()[0]
                cont_seq[name] = []
            else :
                cont_seq[name].append(line.strip())
    for name in cont_seq :
        cont_seq[name] = ''.join(cont_seq[name])
    blocks = {}
    for spacer in spacers :
        spacer[3] = cont_seq[spacer[0]][(int(spacer[1])-1):int(spacer[2])]
        used = {}
        for region in spacer[4] + spacer[5] :
            s, e = int(region[0] + spacer[1] - region[2]), int(region[1] + spacer[2] - region[3])
            tag = region[9].split(':')
            tag[1:3] = int(tag[1]), int(tag[2])
            if tag[1] < tag[2] :
                s, e, d = tag[1] + s - 1, tag[1] + e - 1, '+'
            else :
                s, e, d = tag[1] - e + 1, tag[1] - s + 1, '-'
            if (tag[0], s) not in used and (tag[0], e) not in used :
                used[(tag[0], s)] = used[(tag[0], e)] = 1
                if tag[0] not in blocks :
                    blocks[tag[0]] = [[s, e, tag[-1], d, spacer[3]]]
                else :
                    blocks[tag[0]].append([s, e, tag[-1], d, spacer[3]])
            else :
                print 'aa'
    edges = {e:{} for s in spacers for e in ((s[3], 'L'), (s[3], 'R')) }
    for r, b in blocks.iteritems() :
        if len(b) > 1 :
            b = sorted(b)
            for id in xrange(len(b)-1) :
                b1, b2 = b[id], b[id+1]
                if b1[3] == '+' :
                    if b1[2] in ('WR') and b2[2] in ('WL') :
                        edges[(b1[4], 'R')][(b2[4], 'L')] = edges[(b1[4], 'R')].get((b2[4], 'L'), 0) + 1
                        edges[(b2[4], 'L')][(b1[4], 'R')] = edges[(b1[4], 'R')][(b2[4], 'L')]
                else :
                    if b1[2] in ('WL') and b2[2] in ('WR') :
                        edges[(b1[4], 'L')][(b2[4], 'R')] = edges[(b1[4], 'L')].get((b2[4], 'R'), 0) + 1
                        edges[(b2[4], 'R')][(b1[4], 'L')] = edges[(b1[4], 'L')][(b2[4], 'R')]
    to_del = []
    for edge, conn in edges.iteritems() :
        if len(conn) > 1 :
            c = sorted(conn.items(), key=lambda x:x[1])
            for cc in c[1:] :
                if cc[1] * 4. <= c[0][1] :
                    to_del[edge, cc[0]]
    for m,n in to_del :
        edges[m].pop(n, None)
        edges[n].pop(m, None)
    #edges = { edge:conn for edge, conn in edges.iteritems() if len(conn) > 0 }

    linkage = load(open(positive_link_file))
    for link in linkage :
        if link[0] in blocks :
            b = blocks[link[0]]
            link.extend(b)
        link[2:] = sorted(link[2:])
        for lid in xrange(2, len(link)) :
            l = link[lid]
            if l[2] == 'DR' :
                if lid > 2 :
                    l[0] = link[lid-1][1] + 1
                if lid < len(link) - 1 :
                    l[1] = link[lid+1][0] - 1
                l[4] = link[1][(l[0]-1):l[1]] if l[3] == '+' else rc(link[1][(l[0]-1):l[1]])
    dump(linkage, open(positive_region_file+'.content', 'w'), indent=2)
    x = {x[0]:1 for x, y in edges.iteritems() if len(y) > 0}
    spacers = sorted([[s[3], 2.*len(s[4])*len(s[5])/float(len(s[4])+len(s[5]))] for s in spacers if s[3] in x], key=lambda x:x[1], reverse=True)
    edges = {e:c for e, c in edges.iteritems() if e[0] in x }
    return spacers, edges

def greedy_extension(spacer_conn, path_pool = [], current_path=None, maximum_occurence=3) :
    if len(path_pool) > 2048 :
        return
    if current_path == None :
        current_path = {'path':[], 'occurence':{'{0} {1}'.format(*sublist):0 for sublist in spacer_conn}}

    if len(current_path['path']) == 0 or current_path['path'][-1][:2] == 'L.' :
        if len(current_path['path']) > 0 :
            path_pool.append( deepcopy(current_path) )
        newRs = [conn for conn in spacer_conn if conn[0][:2] == 'R.' if current_path['occurence']['{0} {1}'.format(*conn)]==0]
        if len(newRs) > 0 :
            new_path = deepcopy(current_path)
            newR = newRs[0]
            newR_key = '{0} {1}'.format(*newR)
            if new_path['occurence'][newR_key] == 0 :
                new_path['path'].extend(newR[:2])
                new_path['occurence'][newR_key] += 1
                greedy_extension(spacer_conn, path_pool, new_path, maximum_occurence)
            return
    else :
        newRs = [conn for conn in spacer_conn if conn[0] == current_path['path'][-1]]
        if len(newRs) > 0 :
            for newR in newRs :
                new_path = deepcopy(current_path)
                newR_key = '{0} {1}'.format(*newR)
                if new_path['occurence'][newR_key] <= maximum_occurence :
                    new_path['path'].append(newR[1])
                    new_path['occurence'][newR_key] += 1
                    greedy_extension(spacer_conn, path_pool, new_path, maximum_occurence)
        else :
            new_path = deepcopy(current_path)
            newR = [new_path['path'][-1], 'L.new.?', 0]
            new_path['path'].append(newR[1])
            greedy_extension(spacer_conn, path_pool, new_path, maximum_occurence)

def evaluate_path(spacer_conn, path_pool) :
    read_count = float(sum([conn[2] for conn in spacer_conn]))
    conn_read = {'{0} {1}'.format(*sublist):sublist[2] for sublist in spacer_conn}
    for path in path_pool :
        occ = sum(path['occurence'].values())
        mean = read_count/occ

        var = sum([(conn_read[key] - mean*occ)**2 for key,occ in path['occurence'].iteritems()])/len(conn_read)
        #pdf = sum([math.log(normpdf(conn_read[key], mean*occ, var)) for key, occ in path['occurence'].iteritems()])
        pdf = sum([math.log(normpdf(conn_read[key], 0, var)) if occ < 1 else math.log(normpdf(conn_read[key]/occ, mean, var))*occ for key, occ in path['occurence'].iteritems()])
        #path['lk'] = 2*pdf - 2*occ
        path['lk'] = pdf
    max_pdf = max([x['lk'] for x in path_pool])
    return [x['path'] for x in path_pool if x['lk']==max_pdf]

def file_convert(prefix, reads) :
    import gzip
    outputs = []
    for rid, read in enumerate(reads) :
        fname = '{0}.{1}.fasta'.format(prefix, rid)
        outputs.append(fname)
        with gzip.open(read) as fin, open(fname, 'w') as fout :
            for id, line in enumerate(fin) :
                if id % 4 == 1 :
                    fout.write('>{0}.{2}\n{1}'.format(id/4, line, rid+1))
    return outputs

def extend_edges(edges) :
    out_edge = [e for e in edges.keys() if e[1] == 'R']
    while len(out_edge) :
        edge = out_edge.pop()
        if edge in edges and len(edges[edge]) == 1 :
            in_edge = edges[edge].keys()[0]
            if len(edges[in_edge]) == 1 :
                new_spacer = edge[0] + '-' + in_edge[0]
                edges.pop(edge)
                edges.pop(in_edge)
                edges[(new_spacer, 'R')] = edges.pop((in_edge[0], 'R'))
                out_edge.append((new_spacer, 'R'))
                for n in edges[(new_spacer, 'R')] :
                    edges[n][(new_spacer, 'R')] = edges[n].pop((in_edge[0], 'R'))
                edges[(new_spacer, 'L')] = edges.pop((edge[0], 'L'))
                for n in edges[(new_spacer, 'L')] :
                    edges[n][(new_spacer, 'L')] = edges[n].pop((edge[0], 'L'))
    return edges

def CRISPassemble(prefix, reads, DR) :
    if not debug :
        reads = file_convert(prefix, reads)
    else :
        reads = ['test.0.fasta', 'test.1.fasta']
    # initial screening of reads
    positive_regions, positive_links = identify_positive_reads(prefix, DR, reads)

    # define spacers in reads
    spacers, edges = get_crispr(positive_regions, positive_links)

    # find all possible links
    new_edges = extend_edges(edges)

    spacers = dict(spacers)
    for (edge, direction), conn in new_edges.iteritems() :
        if direction == 'R' :
            edges = edge.split('-')
            cov = sum([spacers[e] for e in edges])/len(edges)
            print 'CRISPR:\t{0}\t{1}'.format(edge, cov)
            if len(conn) > 0 :
                for e, c in conn.iteritems() :
                    print '## LINK:\t{2}\t{0}\t{1}'.format(edge, e[0], c)

debug = 1
if __name__ == '__main__' :
    prefix, DR = sys.argv[1:3]
    reads = sys.argv[3:]
    CRISPassemble(prefix, reads, DR)