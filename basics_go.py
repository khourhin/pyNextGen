import sys

def ensemblGO_to_rtopGO(mart_export):
    """
    From a biomart export with the name of the transcript in one
    column and a GO term in the other, prepare for a R topGO input for
    the function readMappings.
    """
    GO_map = {}
    no_GOs = []
    
    with open(mart_export, 'r') as f:
        for line in f:
            line = line.strip()
            transcript_name = line.split()[0]
            
            if len(line.split()) == 1:
                no_GOs.append(transcript_name)
                continue
            
            go_term = line.split()[1]
            GO_map[transcript_name] = GO_map.get(transcript_name, []) + [go_term]
            
    for k in GO_map:
        sys.stdout.write('{0}\t{1}\n'.format(k, ', '.join(GO_map[k])))

    sys.stderr.write('{} transcripts did not have any GO term.\n'.format(len(no_GOs)))

    
ensemblGO_to_rtopGO('/home/ekornobis/Programming/pyNextGen/demo_data/ensembl_Tg_GOs.txt')
