import os
import pandas as pd
import numpy


import pandas as pd

def search_dataframe(target_string, df):
    # Obtém o número de linhas e colunas no DataFrame
    num_rows, num_cols = df.shape

    # Percorre cada linha do DataFrame
    for row_index in range(num_rows):
        row = df.iloc[row_index]  # Obtém a linha atual como uma Série

        # Percorre cada coluna na linha atual
        line = ""
        found = False
        for col_name in df.columns:
            cell_value = row[col_name]  # Obtém o valor da célula
            line += cell_value + "\t"
            # Verifica se o valor da célula é igual à string alvo
            if cell_value == target_string:
                # Se for igual, imprima a linha inteira
                #print(row)
                found = True
        if found:
            print(line)

def process_files(input_directory, string_list):
    new_directory = input_directory + "_net_act"
    os.makedirs(new_directory, exist_ok=True)

    input_files = [f for f in os.listdir(input_directory) if f.endswith(".tsv")]

    for input_file in input_files:
        print(f"Table:::::::::::::-->  {input_file}")
        input_path = os.path.join(input_directory, input_file)
        output_path = os.path.join(new_directory, input_file)

        with open(input_path, 'r') as file:
            header = file.readline()

        relevant_lines = []

        df = pd.read_csv(input_path, sep='\t', dtype=str)
        df = df.fillna('')
        for string in string_list:
            search_dataframe(string, df)
        #continue

        for string in string_list:
            mask = df.apply(lambda row: string in ' '.join(map(str, row)), axis=1)
            relevant_lines.extend(df[mask].values.tolist())

            #relevant_lines.extend(df[df.apply(lambda row: string in ' '.join(map(str, row)), axis=1)].values.tolist())
            #df = df[df.apply(lambda row: string in ' '.join(map(str, row)), axis=1)]
#        relevant_lines = df[df.apply(lambda row: any(string in ' '.join(map(str, row)) for string in string_list), axis=1)]
        #df = pd.DataFrame(relevant_lines, columns=df.columns)
        #print(df)
        with open(output_path, 'w') as output_file:
            output_file.write(header)
            for line in relevant_lines:
                output_file.write('\t'.join(map(str, line)) + '\n')


# Example usage:
string_list = ["Su(var)205", "Top3beta", "Mef2", "Clk", "Dref", "TfIIB", "Myc", "AGO2", "Nipped-B", "Cp190", "TfIIA-L",
               "Trl", "ash1", "Raf", "Abd-B", "Orc2", "Rbf", "mof", "msl-1", "Hmr"]
input_directory = "/home/saulo/snet/hyperon/das/das/flybase2metta/fb_data/2023_05/precomputed"
process_files(input_directory, string_list)
