import pandas as pd


def add_list_to_column(df, list_to_add, column_name, repartition=False):
    if repartition:
        df = df.repartition(partition_size="25MB").persist()

    index_name = df.index.name
    divisions = df.divisions

    df["index"] = 1
    idx = df["index"].cumsum()
    df["index"] = idx
    df = df.reset_index().set_index("index", sorted=True)

    list_to_add = pd.DataFrame(list_to_add)
    list_to_add["index"] = idx
    list_to_add = list_to_add.set_index("index")

    df = df.join(list_to_add, how="left", on="index")

    if divisions[0] != None:
        df = df.set_index(index_name, sorted=True, divisions=divisions)

    else:
        df = df.set_index(index_name)

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
