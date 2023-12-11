import pandas as pd

def separate_and_write_datasets(input_file, column_name="SQL Table"):
    try:
        df = pd.read_csv(input_file, sep='\t')
        df = df.fillna("NONE")
        # Filtrar linhas com valores vazios na coluna especificada
        not_sql_columns = df[df[column_name] == 'NONE']
        sql_columns = df[df[column_name] != 'NONE']

        # Obter as colunas originais
        columns = df.columns.tolist()

        # Escrever os dois datasets em arquivos distintos, mantendo as mesmas colunas originais
        not_sql_columns.to_csv("mapped_not_sql_columns.tsv", sep='\t', index=False, columns=columns)
        sql_columns.to_csv("mapped_sql_columns.tsv", sep='\t', index=False, columns=columns)

        print("Datasets separados e escritos com sucesso.")
    except Exception as e:
        print("Erro ao separar e escrever os datasets:", e)

def main():
    input_file = "fb_precomputed_all_columns.tsv"
    column_name = "SQL Table"  # column to split the dataset...
    separate_and_write_datasets(input_file, column_name)

if __name__ == "__main__":
    main()
