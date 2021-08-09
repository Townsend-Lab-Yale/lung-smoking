import pandas as pd
from locations import gene_coordinates_file


genes = pd.read_csv(gene_coordinates_file, index_col='gene')

def number_mutations(data, mutations = None):
    numbers_per_mutation = {}

    numbers_per_mutation = {mutation: len(data[
            (data['Start_Position'] >= genes.loc[mutation, 'start']) &
                (data['Start_Position'] <= genes.loc[mutation, 'end']) &
                (data['Chromosome'] == genes.loc[mutation, 'chromosome'])]
            ['Sample ID'].unique())
                        for mutation in mutations}

    result = [(k,v) for k, v in sorted(numbers_per_mutation.items(), key=lambda item: item[1], reverse=True)]
    
    with open('ranked_mutations.txt', 'w') as fp:
        fp.write('\n'.join('%s %s' % x for x in result))
        fp.write('\n')
    
    return True

def produce_ranked_list(genes = 'ranked_mutations.txt', sep = "\n", header = None, set_length = None):
    dict = {}

    gene_list = pd.read_csv(genes, sep = sep, header = header)
    gene_list = gene_list.values.tolist()
    for i in gene_list:
        ind = i[0].find(' ')
        dict[i[0][:ind]] = int(i[0][ind+1:])
    if set_length == None:
        ranked_list = [k for k, v in sorted(dict.items(), key=lambda item: item[1], reverse=True)]
    else:
        ranked_list = [k for k, v in sorted(dict.items(), key=lambda item: item[1], reverse=True) 
                                    if v > (2**(set_length-1))-1]
    return ranked_list

def produce_prevalence_dict(genes_ = 'ranked_mutations.txt', sep_ = "\n", header_ = None):
    dict = {}

    gene_list = pd.read_csv(genes_, sep = sep_, header = header_)
    gene_list = gene_list.values.tolist()
    for i in gene_list:
        ind = i[0].find(' ')
        dict[i[0][:ind]] = int(i[0][ind+1:])
    return dict