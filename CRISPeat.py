# CRISPacer
# 1. combine a) confirmed DR; b) uncertain DR; c) CRT results
# 2. end randomizing for >50% instances
# 3. nucmer on all
# 4. identify reliable ends
# 5. remove suspect DRs
# 6. one-off clustering

import sys, os, subprocess, gzip, re
from json import dump, load

nucmer = '/home/zhemin/bin/nucmer --maxmatch -l 5 -c 7 -b 10 -g 3 --nosimplify '
nucmer2 = '/home/zhemin/bin/nucmer --maxmatch -l 11 -c 11 -b 20 -g 10 --nosimplify -f'
show_coords = '/home/zhemin/bin/show-coords -lqHT '
min_DR = 25
max_DR = 55

cmpl = {
    'A':'T',
    'T':'A',
    'C':'G',
    'G':'C',
}
def rc(sequence) :
    return ''.join([cmpl[b] for b in sequence[::-1]])

def reverseSet(DRs) :
    for DR, info in DRs.items() :
        if info[1] < 1 :
            rDR = rc(DR)
            if rDR in DRs :
                DRs[rDR][2] += info[2]
            else :
                DRs[rDR] = [rDR, 0., info[2]]
    return DRs

def loadFasta(fname) :
    seq = {}
    with open(fname) as fin :
        for line in fin :
            if line.startswith('>') :
                name = line[1:].strip()
                seq[name] = []
            else :
                seq[name].extend(line.strip().split())
    for n, s in seq.iteritems() :
        seq[n] = ''.join(s).upper()
    return seq

def loadUncertain(fname) :
    seq_cnt = {}
    if fname.lower().endswith('.gz') :
        fin = gzip.open(fname)
    else :
        fin = open(fname)
    head = fin.readline()
    if head.startswith('>') :
        fin.close()
        seq = loadFasta(fname)
        for n, s in seq.iteritems() :
            nn = n.split()
            try :
                cnt = float(nn[1])
            except :
                cnt = 1.
            seq_cnt[s] = seq_cnt.get(s, 0) + cnt
    else :
        for line in fin :
            if line[0] in '0123456789' :
                s = line.split()[1].upper()
                if re.match(r'^[ACGT]+$', s) :
                    seq_cnt[s] = seq_cnt.get(s, 0) + 1.
        fin.close()
    return seq_cnt

def loadDRs(confirmed, uncertains) :
    seq = {}    
    if confirmed :
        confirmedSeq = loadFasta(confirmed)
        for n, s in confirmedSeq.iteritems() :
            nn = n.split()
            try :
                cnt = float(nn[1])
            except :
                cnt = 0.
            seq[s] = [nn[0], 2, cnt]

    uncertainSeq = [loadUncertain(uncertain) for uncertain in uncertains]
    
    for useq in uncertainSeq :
        for s, cnt in useq.iteritems() :
            if len(s) >= min_DR and len(s) <= max_DR :
                if s not in seq :
                    seq[s] = [s, 0, 0.]
                seq[s][2] += cnt
    return seq

def randomEnd(DRs) :
    DR_num = sum([ n[2] for n in DRs.itervalues() ])
    acc_cnt = 0
    new_drs = {}
    for DR, info in sorted(DRs.items(), key=lambda x:x[1][2], reverse=True) :
        info[1] += 1
        for u, v in ((0, 1), (-2, -1)) :
            for bu in set(cmpl) - set([DR[u]]) :
                for bv in set(cmpl) - set([DR[v]]) :
                    new_dr = ''.join([bu, bv, DR[(v+1):]]) if v > 0 else ''.join([DR[:u], bu, bv])
                    if new_dr in new_drs : continue
                    new_drs[new_dr] = 1
                    if new_dr not in DRs :
                        DRs[new_dr] = [new_dr, info[1], info[2]-0.5]
                    else :
                        DRs[new_dr][1] += 1
                        DRs[new_dr][2] = info[2]-0.5
        acc_cnt += info[2]
        if acc_cnt >= 0.75 * DR_num :
            break
    return DRs

def nucmerise(DRs) :
    confirmedDRs = {}
    with open('CRISPeat.ref', 'w') as fout1, open('CRISPeat.qry', 'w') as fout2 :
        for id, (DR, info) in enumerate(sorted(DRs.iteritems(), key=lambda x:x[1][1:], reverse=True)) :
            DRname = '{1:04d}.{0}'.format(DR, id)
            if id == 0 or info[1] >= 2 :
                confirmedDRs[DRname] = info + [ 0.5 ]
                fout1.write('>{1}\n{0}\n'.format(DR, DRname))
            else :
                confirmedDRs[DRname] = [ '', '', '', -0.5 ]
            fout2.write('>{1}\n{0}\n'.format(DR, DRname))
    nucmer_cmd = nucmer + '-p CRISPeat CRISPeat.ref CRISPeat.qry'
    subprocess.Popen(nucmer_cmd.split()).communicate()
    show_coords_cmd = show_coords + 'CRISPeat.delta'
    sc = subprocess.Popen(show_coords_cmd.split(), stdout=subprocess.PIPE)
    
    for line in iter(sc.stdout.readline, r'') :
        part = line.strip().split()
        part[:9] = [float(x) for x in part[:9]]
        idx1, idx2 = int(part[9][:4]), int(part[10][:4])
        if idx1 >= idx2 :
            continue
        score = max(part[4] * 2 - 5*part[4]*(100. - part[6])/100., 0.5)
        if part[2] > part[3] :
            if abs(confirmedDRs[part[10]][3]) < score :
                confirmedDRs[part[10]][3] = -score
            continue
        elif abs(confirmedDRs[part[10]][3]) < score :
            confirmedDRs[part[10]][3] = score
        
        if confirmedDRs[part[10]][0] == '' and part[4] >= min_DR-7 :
            editDist = [part[0]-1, part[7] - part[1], abs(part[2]-part[0]), abs( (part[8]-part[3]) - (part[7]-part[1]) )]
            if sum(editDist[:2]) <= 7 and max(editDist[:2]) <= 5 and sum(editDist[2:]) == 0 :
                confirmedDRs[part[10]][:3] = DRs[part[10].split('.', 1)[1]]
    
    confirmedDRs = {x:y for x,y in confirmedDRs.iteritems() if y[0] != '' and y[3] > 0}
    with open('CRISpeat.qry', 'w') as fout2 :
        for id, (DRname, info) in enumerate(sorted(confirmedDRs.iteritems(), key=lambda x:x[1][1:], reverse=True)) :
            DR = DRname.split('.', 1)[1]
            fout2.write('>{1}\n{0}\n'.format(DR, DRname))
    nucmer_cmd = nucmer + '-p CRISPeat CRISPeat.qry CRISPeat.qry'
    subprocess.Popen(nucmer_cmd.split()).communicate()
    show_coords_cmd = show_coords + 'CRISPeat.delta'
    sc = subprocess.Popen(show_coords_cmd.split(), stdout=subprocess.PIPE)
    
    for line in iter(sc.stdout.readline, r'') :
        part = line.strip().split()
        if part[9] not in confirmedDRs or part[10] not in confirmedDRs : continue
        part[:9] = [float(x) for x in part[:9]]
        idx1, idx2 = int(part[9][:4]), int(part[10][:4])
        if idx1 >= idx2 :
            continue
        editDist = [part[0]-1, part[7] - part[1], part[2]-1, part[8]-part[3]]
        if sum(editDist) == 0 and (1- part[6]/100.)*part[4] <= 1.9 :
            confirmedDRs.pop(part[10], None)
        
    return confirmedDRs

def CRISPace(*argv) :
    uncertains, confirmed = [], None
    for n in argv[1:] :
        if n.startswith('confirm') :
            confirmed = n.split('=', 1)[1]
        else :
            uncertains.append(n)
    DRs = loadDRs(confirmed, uncertains)
    reverseSet(randomEnd(DRs))
    DRs = nucmerise(DRs)
    for DR , info in sorted(DRs.iteritems(), key=lambda x:x[1][1:], reverse=True) :
        seq = DR.split('.', 1)[1]
        print '>{0}|{1}\n{2}'.format(info[0], info[2], seq)

debug = 1
if __name__ == '__main__' :
    CRISPace(*sys.argv)