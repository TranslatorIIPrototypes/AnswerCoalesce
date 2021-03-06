import csv
from collections import defaultdict
import json

#A link is: link = (other_node, predicate, source)
# source is a boolean that is true if node is the source, and false if other_node is the source

def str2list(nlabels):
    x = nlabels[1:-1] #strip [ ]
    parts = x.split(',')
    labs = [ pi[1:-1] for pi in parts] #strip ""
    return labs

def add_labels(lfile,nid,nlabels,done):
    if nid in done:
        return
    labs = str2list(nlabels)
    lfile.write(f'{nid}\t{labs}\n')
    done.add(nid)

def go():
    wrote_labels = set()
    nodes_to_links = defaultdict(list)
    edgecounts = defaultdict(int)
    with open('everything.csv','r') as inf, open('nodelabels.txt','w') as labelfile:
        reader = csv.DictReader(inf)
        nl = 0
        for line in reader:
            nl += 1
            source_id = line['source_id']
            target_id = line['target_id']
            if source_id == '':
                continue
            pred = line['predicate']
            add_labels(labelfile,source_id,line['source_labels'],wrote_labels)
            add_labels(labelfile,target_id,line['target_labels'],wrote_labels)
            #Remove variant/anatomy edges from GTEX.  These are going away in the new load anyway
            if source_id.startswith('CAID'):
                if target_id.startswith('UBERON'):
                    if pred == 'affects_expression_of':
                        continue
            if source_id.startswith('UBERON'):
                if target_id.startswith('NCBIGene'):
                    if pred == 'expresses':
                        continue
            source_link = (target_id,pred,True)
            target_link = (source_id,pred,False)
            nodes_to_links[source_id].append(source_link)
            nodes_to_links[target_id].append(target_link)
            for tlabel in str2list(line['target_labels']):
                edgecounts[ (source_id, pred, True, tlabel) ] += 1
            for slabel in str2list(line['source_labels']):
                edgecounts[ (target_id, pred, False, slabel) ] += 1
            if nl % 1000000 == 0:
                print(nl)
    print('ate the whole thing')
    with open('links.txt','w') as outf:
        for node,links in nodes_to_links.items():
            outf.write(f'{node}\t{json.dumps(links)}\n')
    with open('backlinks.txt','w') as outf:
        for key,value in edgecounts.items():
            outf.write(f'{key}\t{value}\n')

if __name__ == '__main__':
    go()
