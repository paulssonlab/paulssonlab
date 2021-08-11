import pandas as pd


def add_list_to_column(df, list_to_add, column_name):
    df = df.repartition(partition_size="25MB").persist()
    df["index"] = 1
    idx = df["index"].cumsum()
    df["index"] = idx

    list_to_add = pd.DataFrame(list_to_add)
    list_to_add["index"] = idx
    df = df.join(list_to_add.set_index("index"), how="left", on="index")

    df = df.drop(["index"], axis=1)

    df.columns = df.columns.tolist()[:-1] + [column_name]

    return df


def set_new_aligned_index(df, index_column):
    ### Sets a column to be the new index

    ### Fast, but assumes index is both sorted and division aligned to the primary index

    first_indices = df.loc[list(df.divisions)].compute()
    new_index_divisions = first_indices[index_column].to_list()
    output_df = df.reset_index(drop=False).set_index(
        index_column,
        drop=True,
        sorted=True,
        npartitions=df.npartitions,
        divisions=new_index_divisions,
    )
    return output_df
