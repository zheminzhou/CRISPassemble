# CRISPacer
# 1. combine a) confirmed DR; b) uncertain DR; c) CRT results
# 2. end randomizing for >50% instances
# 3. nucmer on all
# 4. identify reliable ends
# 5. remove suspect DRs
# 6. one-off clustering

import sys, os, subprocess, gzip, re
from json import dump, load

nucmer = '/home/zhemin/bin/nucmer --maxmatch -l 7 -c 9 -b 20 -g 10 --nosimplify '
show_coords = '/home/zhemin/bin/show-coords -lqHT '
min_DR = 25
max_DR = 55

cmpl = {
    'A':'T',
    'T':'A',
    'C':'G',
    'G':'C',
    'N':'N'
}
def rc(sequence) :
    return ''.join([cmpl.get(b, 'N') for b in sequence[::-1]])

def loadFasta(fname) :
    seq = []
    with open(fname) as fin :
        for line in fin :
            if line.startswith('>') :
                name = line[1:].strip()
                seq.append([])
            else :
                seq[-1].extend(line.strip().split())
    for n, s in enumerate(seq) :
        seq[n] = ''.join(s).upper()
    return seq

def loadPrediction(fname) :
    if fname.lower().endswith('.gz') :
        fin = gzip.open(fname)
    else :
        fin = open(fname)
    header = fin.readline()
    if header.startswith('Help on reading') :
        dr, crispr = loadPILERCR(fin)
    else :
        dr, crispr = loadCRT(fin)
    fin.close()
    return dr, crispr

def loadPILERCR(fin) :
    dr, crispr = {}, [[]]
    for line in fin :
        if line.startswith('DETAIL REP') :
            break
    for line in fin :
        if line.startswith('SUMMARY') :
            break
        else :
            part = line.strip().split()
            if len(part) > 0 and part[0][0] in ('0123456789') :
                if len(part) > 4 :
                    crispr[-1].extend([['DR', part[-2]], ['spacer', re.sub(r'[^ACGT]', 'N', part[-1])]])
                else :
                    cons = part[3]
                    for id, c in enumerate(crispr[-1]) :
                        if c[0] == 'DR' :
                            c[1] = re.sub(r'[^ACGT]', 'N', ''.join([ m if m != '.' else n for m,n in zip(list(c[1]), list(cons)) ]).replace('-', ''))
                            dr[c[1]] = dr.get(c[1], 0.) + 1
                    crispr.append([])
    return dr, [s for s in crispr if len(s)>0]

def loadCRT(fin) :
    dr, crispr = {}, []
    for line in fin :
        if line.startswith('--------') :
            crispr.append([])
        if line[0] in '0123456789' :
            part = line.split()
            s = re.sub(r'[^ACGT]', 'N', part[1].upper())
            dr[s] = dr.get(s, 0) + 1.
            crispr[-1].append(['DR', s])
            if len(part) >= 3 :
                s = re.sub(r'[^ACGT]', 'N', part[2].upper())
                crispr[-1].append(['spacer', s])
    return dr, [s for s in crispr if len(s)>0]

def loadCRISPRs(CRTfiles) :
    DRs, CRISPRs = {}, {}
    for CRTfile in CRTfiles :
        DR, CRISPR = loadPrediction(CRTfile)
        for s, c in DR.iteritems() :
            if len(s) >= min_DR and len(s) <= max_DR :
                DRs[s] = DRs.get(s, 0.) + c
        CRISPRs[CRTfile] = CRISPR
    return DRs, CRISPRs

def nucmerise(DRs) :
    confirmedDRs = {}
    with open('CRTparser.qry', 'w') as fout :
        for id, (DR, info) in enumerate(sorted(DRs.iteritems(), key=lambda x:x[1], reverse=True)) :
            DRname = '{1:05d}.{0}'.format(DR, id)
            confirmedDRs[DRname] = [0., 1]
            fout.write('>{1}\n{0}\n'.format(DR, DRname))
    nucmer_cmd = nucmer + '-p CRTparser CRTparser.qry CRTparser.qry'
    subprocess.Popen(nucmer_cmd.split()).communicate()
    show_coords_cmd = show_coords + 'CRTparser.delta'
    sc = subprocess.Popen(show_coords_cmd.split(), stdout=subprocess.PIPE)
    
    for line in iter(sc.stdout.readline, r'') :
        part = line.strip().split()
        part[:9] = [float(x) for x in part[:9]]
        idx1, idx2 = int(part[9].split('.')[0]), int(part[10].split('.')[0])
        if idx1 >= idx2 :
            continue
        score = max(part[4] * 2 - 5*part[4]*(100. - part[6])/100., 0.5)
        if part[10] not in confirmedDRs or confirmedDRs[part[10]][0] < score :
            d = -1 if part[2] > part[3] else 1
            d *= confirmedDRs[part[9]][1]
            confirmedDRs[part[10]] = [score, d]
    
    return {d.split('.', 1)[-1]:v[1] for d,v in confirmedDRs.iteritems() }

def CRTparser(*argv) :
    CRTs, confirmed = [], None
    for n in argv[1:] :
        if n.startswith('confirm') :
            confirmed = n.split('=', 1)[1]
        else :
            CRTs.append(n)
    DRs, CRISPRs = loadCRISPRs(CRTs)
    confirmedDRs = loadFasta(confirmed)
    max_occ = max(DRs.values())+1
    DRs.update({dr:max_occ+id for id, dr in enumerate(confirmedDRs[::-1]) })
    DRs = nucmerise(DRs)
    
    for fname, crispr in sorted(CRISPRs.iteritems()) :
        for c in crispr :
            positive = (sum([DRs.get(cc[1], 0) for cc in c if cc[0] == 'DR']) >= 0)
            if not positive :
                c = [ [cc[0], rc(cc[1])] for cc in c[::-1] ]
            c1 = '-'.join([ ('({0})'.format(cc[1]) if cc[0] == 'DR' else cc[1]) for cc in c ])
            print '{0}\t{1}'.format(fname, c1)

if __name__ == '__main__' :
    CRTparser(*sys.argv)