import pandas as pd

def prepare_starting_molecules(filename):
    '''
    load starting molecules from csv file
    return: list of mol
    '''
    print('Loading starting molecules from %s' % filename)

    if filename[-3:] == 'csv':
        mols_list = list(pd.read_csv(filename)['mol'])
        starting_mols = set(mols_list)

    print('%d starting molecules loaded' % len(starting_mols))
    return starting_mols

def prepare_starting_molecules_subset(filename, saveto, nums=10):
    '''
    load starting molecules subset from csv file for test
    return: list of mol
    '''
    print('Loading starting molecules subset from %s' % filename)

    if filename[-3:] == 'csv':
        data = pd.read_csv(filename, nrows=nums+1)
        data.rename(columns={data.columns[0]: 'index'}, inplace=True)
        data.to_csv(saveto, index=False)

# prepare_starting_molecules_subset('../dataset/origin_dict.csv', '../dataset/origin_dict_subset.csv')

starting_mols = prepare_starting_molecules('../dataset/origin_dict.csv')
with open('../dataset/starting_mols.txt', 'w') as f:
    for mol in starting_mols:
        f.write('%s\n' % mol)
