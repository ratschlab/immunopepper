import numpy as np
import cPickle
import h5py
expr_path = '/Users/jiayu/PycharmProjects/CBM_RA/projects2018_immunopepper/tests/data/expr_ts.txt'
splicegraph_path = '/Users/jiayu/PycharmProjects/CBM_RA/projects2018_immunopepper/tests/data/posgraph/spladder/genes_graph_conf3.merge_graphs.pickle'
count_path = '/Users/jiayu/PycharmProjects/CBM_RA/projects2018_immunopepper/tests/data/posgraph/spladder/genes_graph_conf3.merge_graphs.count.hdf5'

count = h5py.File(count_path)
f = open(splicegraph_path,'r')
gene,metadata = cPickle.load(f)
gene0 = gene[0]

with open(expr_path,'r') as f:
    lines = f.readlines()
    expr_dict = {}
    for line in lines:
        item = line.split('\t')
        trans = item[0]
        expr = [int(iexpr) for iexpr in item[1:]]
        expr_dict[trans] = expr

expr1 = expr_dict['TRANS1.1']
expr3 = expr_dict['TRANS1.3']
expr4 = expr_dict['TRANS1.4']
expr5 = expr_dict['TRANS1.5']
expr2 = expr_dict['TRANS1.2']


# segment1 [4,26)
# trans1 (11,28,0),(38,49,0),(61,72,0)
# trans2 (60,74,0),(87,101,0)
# trans3 (11,28,1),(38,41,0) -> (12,28,0),(38,42)
# trans4 (14,25,2),(38,49,0),(66,74,0),(87,101,0) -> (16,25,0)
# trans5 (53,83,0)

dna = 'GTAATGTGTAAGATGACGCACGCATGGTGGTATTGGAGATGGGTTGCGGAGTAAGTTCGAGTTCAGGTACGTATAGTTGCTAAGTAGTCGACGTTCTGGTGGTAGACCTTTTGTTACCCAATTGGGTAAGGCAGCCATGTTGATATACAT'

ts_dict = {'TRANS1.1': 'GATGACGCACGCATGGTGATGGGTTGCGGATTCAGGTACGTA',
 'TRANS1.2': 'GTTCAGGTACGTATATCGACGTTCTGGTGG',
 'TRANS1.3': 'ATGACGCACGCATGGTGATGG',
 'TRANS1.4': 'CGCACGCATGATGGGTTGCGGAGTACGTATATCGACGTTCTGGTGG',
 'TRANS1.5': 'AGTTCGAGTTCAGGTACGTATAGTTGCTAAG'}

ts_to_cds={'TRANS1.1': [(11, 28, 0), (38, 49, 0), (61, 72, 0)],
 'TRANS1.2': [(60, 74, 0), (87, 101, 0)],
 'TRANS1.3': [(11, 28, 1), (38, 41, 0)],
 'TRANS1.4': [(14, 25, 2), (38, 49, 0), (66, 74, 0), (87, 101, 0)],
 'TRANS1.5': [(53, 83, 0)]}

# segment1 [4,26), 1287.2 = 1287.2
segment1 = np.zeros(26-4)
segment1[11-4:26-4] += expr1[:26-11]
segment1[12-4:26-4] += expr3[:25-11]
segment1[16-4:25-4] += expr4[:23-14]
assert np.mean(segment1) == count['segments'][0]

# segment2 [26,29),   # 2606 != 2528
segment2 = np.zeros(29-26)
assert(dna[26:29] == ts_dict['TRANS1.1'][15:18])
segment2 += expr1[26-11:29-11]
assert(dna[26:29] == ts_dict['TRANS1.3'][14:17])
segment2 += expr3[26-12:29-12]
print(np.mean(segment2))

# segment3 [38,49), 1495 != 1409.9
segment3 = np.zeros(49-38)
segment3 += expr1[18:29]
assert(dna[38:49] == ts_dict['TRANS1.1'][18:29])
segment3 += expr4[10:21]
assert(dna[38:49] == ts_dict['TRANS1.4'][10:21])
segment3[:4] += expr3[17:]
assert(dna[38:42] == ts_dict['TRANS1.3'][17:])
print(np.mean(segment3))

# segment4 [49,50) 1115 = 1115
segment4 = 0
assert(dna[49] == ts_dict['TRANS1.1'][29])
assert(dna[49] == ts_dict['TRANS1.4'][21])
segment4 += expr1[29]
segment4 += expr4[21]
print(np.mean(segment4))

# segment5 [50,53) 0 = 0

# segment6 [53,61) 297.15 = 297.15
segment = np.zeros(61-53)
assert(dna[53:61] == ts_dict['TRANS1.5'][:61-53])
assert(dna[60] == ts_dict['TRANS1.2'][0])
segment += expr5[:8]
segment[-1] += expr2[0]
print(np.mean(segment))

# segment7 [61,66) 1814 != 1763
segment = np.zeros(66-61)
assert(dna[61:66] == ts_dict['TRANS1.1'][-12:-7])
assert(dna[61:66] == ts_dict['TRANS1.2'][1:6])
assert(dna[61:66] == ts_dict['TRANS1.5'][8:13])
segment += expr1[-12:-7]
segment += expr2[1:6]
segment += expr5[8:13]
print(np.mean(segment))

# segment8 [66,75) 787 != 1717
segment = np.zeros(75-66)
#assert(dna[66:75] == ts_dict['TRANS1.4'][1:6])
assert(dna[66:75] == ts_dict['TRANS1.5'][13:22])
#segment += expr4[1:6]
segment += expr5[13:22]
print(np.mean(segment))


# segment9 [75,84) 244 = 244
assert(dna[75:84] == ts_dict['TRANS1.5'][22:])
print(np.mean(expr5[22:]))

# segment10 [84,87), 0 = 0

# segment11 [87,90) 874 != 807
segment = np.zeros(90-87)
assert(dna[89] == ts_dict['TRANS1.4'][-13])
assert(dna[87:90] == ts_dict['TRANS1.2'][15:18])
segment[-1] += expr4[-13]
segment += expr2[15:18]
print(np.mean(segment))

# segment12 [90,102) 363 = 363
segment = np.zeros(102-90)
assert(dna[90:101] == ts_dict['TRANS1.4'][-12:-1])
assert(dna[90:101] == ts_dict['TRANS1.2'][18:29])
segment[:-1] += expr4[-12:-1]
segment[:-1] += expr2[18:29]
print(np.mean(segment))



def get_relative_id(vstart,v_stop,v_pair_list):
    pass

def get_expected_expression(vstart, vstop, expr_dict, ts_dict, dna, ts_to_cds):
    for trans_name,trans_comp in enumerate(ts_to_cds.items()):
        read_frame = trans_comp[0][2]
        v_list = []
        for pair in trans_comp:
            pair_start = pair[0]+ read_frame
            pair_end = pair[1] + read_frame
            v_list.extend(list(set(range(vstart,vstop)).intersection(set(range(pair_start,pair_end)))))
            v_max = max(v_list)
            v_min = min(v_list)
            id = get_relative_id(v_min,v_max,trans_comp)



for idx in range(gene0.segmentgraph.segments.shape[1]):
    vstart = gene0.segmentgraph.segments[0,idx]
    vstop = gene0.segmentgraph.segments[1,idx]
    expe_expr = get_expected_expression(vstart,vstop,expr_dict,ts_dict,ts_to_cds)
    assert(int(count['segments'][idx]) == int(expe_expr))

