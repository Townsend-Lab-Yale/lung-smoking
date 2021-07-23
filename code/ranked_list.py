import pandas as pd

def produce_ranked_list(genes_ = 'ranked_mutations.txt', sep_ = "\n", header_ = None, set_length_ = None):
    dict = {}

    gene_list = pd.read_csv(genes_, sep = sep_, header = header_)
    gene_list = gene_list.values.tolist()
    for i in gene_list:
        ind = i[0].find(' ')
        dict[i[0][:ind]] = int(i[0][ind+1:])
    if set_length_ == None:
        ranked_list_test5 = [k for k, v in sorted(dict.items(), key=lambda item: item[1])][::-1]
    else:
        ranked_list_test5 = [k for k, v in sorted(dict.items(), key=lambda item: item[1]) 
                                    if v > (2**(set_length_-1))-1][::-1]
    return ranked_list_test5
