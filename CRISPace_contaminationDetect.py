import sys, os

if __name__ == '__main__' :
    for fname in sys.argv[1:] :
        spacer = {}
        mean_depth = [0., 0.]
        with open(fname) as fin :
            for line in fin :
                part = line.strip().split()
                if part[0] == 'CRISPR:' :
                    depth = float(part[2])
                    for s in part[1].split('-') :
                        if not s.startswith('(') :
                            spacer[s] = depth
                            mean_depth[0] += depth
                            mean_depth[1] += 1
        if mean_depth[1] == 0 : continue
        md = mean_depth[0]/mean_depth[1]
        for s,d in spacer.iteritems() :
            if d * 2 < md :
                print '{0}\t{1}\t{2}\t{3}'.format(fname, s, md, d/md)
