"""
TODO: docstring
"""

from csv import reader
import pandas as pd


def process_data(fp):
    """TODO: doc"""
    with open(fp) as infile:
        read = reader(infile)

        table = {}
        for i, row in enumerate(read):
            # Full table, run through median machines
            if len(table) == 5:
                results = clean_chunk(table, i)
                table = {}
            elif row == '':
                continue
            else:
                table[i] = row

def clean_chunk(data_chunk, index) -> tuple:
    """TODO: doc"""
    jl_input = []
    mlab_input = []
    jl_dict = {[]}
    mlab_dict = {[]}

    ii = data_chunk[index - 5][0][9:]
    jj = data_chunk[index - 5][1][4:]

    # ALl jacked up
    for n in range(4, 1, -1):
        for _ in data_chunk:
            jl_dict[f'col{n}'] = append(data_chunk[index - n][0].replace('\t', ' '))
            mlab_dict[f'col{n}'].append(data_chunk[index - n][0].replace('\t', ' ').replace('im', 'i'))

    jl_input = [jl_dict]
    mlab_input = [mlab_dict]

    # print(mlab_input)
    mlab_df = pd.DataFrame.from_records(mlab_input)

    print(mlab_df)

    return 0,0
    


if __name__ == "__main__":
    FILE = "juliaOut/first_localfilt/tmp.csv"
    process_data(FILE)