# fmt: off
from .utils import pandas_hdf5_handler
from .trcluster import dask_controller

import h5py
import os

import skimage as sk
import pandas as pd
import numpy as np
import dask.dataframe as dd
import dask.delayed as delayed

from distributed.client import futures_of
from time import sleep

from matplotlib import pyplot as plt
## HERE

class regionprops_extractor:
    def __init__(self,headpath,segmentationdir,intensity_channel_list=None,include_background=False,props=['centroid','area','mean_intensity'],unpack_dict={'centroid':["centroid_y","centroid_x"]}):
        self.headpath = headpath
        self.intensity_channel_list = intensity_channel_list
        self.intensity_channel_dict = {channel:i for i,channel in enumerate(intensity_channel_list)}
        self.include_background = include_background
        self.kymographpath = headpath + "/kymograph"
        self.segmentationpath = headpath + "/" + segmentationdir
        self.metapath = self.kymographpath + "/metadata"
        self.analysispath = headpath + "/analysis"
        self.props = props
        self.unpack_dict = unpack_dict

    def get_file_regionprops(self,file_idx):
        segmentation_file = self.segmentationpath + "/segmentation_" + str(file_idx) + ".hdf5"
        kymograph_file = self.kymographpath + "/kymograph_" + str(file_idx) + ".hdf5"

        with h5py.File(segmentation_file,"r") as segfile:
            seg_arr = segfile["data"][:]
        if self.intensity_channel_list is not None:
            kymo_arr_list = []
            with h5py.File(kymograph_file,"r") as kymofile:
                for intensity_channel in self.intensity_channel_list:
                    kymo_arr_list.append(kymofile[intensity_channel][:])
        all_props_list = []
        for k in range(seg_arr.shape[0]):
            for t in range(seg_arr.shape[1]):
                labels = sk.measure.label(seg_arr[k,t])
                ## Measure regionprops of background pixels; will always be marked as the first object
                if self.include_background:
                    labels += 1
                if self.intensity_channel_list is not None:
                    for i,intensity_channel in enumerate(self.intensity_channel_list):
                        rps = sk.measure.regionprops(labels, kymo_arr_list[i][k,t])
                        props_list = []
                        for obj,rp in enumerate(rps):
                            props_entry = [file_idx, k, t, obj, intensity_channel]
                            for prop_key in self.props:
                                if prop_key in self.unpack_dict.keys():
                                    prop_split = self.unpack_dict[prop_key]
                                    prop_output = rp[prop_key]
                                    props_entry += [prop_output[i] for i in range(len(prop_split))]
                                else:
                                    props_entry += [rp[prop_key]]
                            props_list.append(props_entry)
                        all_props_list+=props_list

#                         for prop_key in self.props:
#                             props_list.append([file_idx, k, t, obj, intensity_channel])
#                         props_list = [[file_idx, k, t, obj, intensity_channel]+[getattr(rp, prop_key) for prop_key in self.props] for obj,rp in enumerate(rps)]
#                         all_props_list+=props_list
                else:
                    rps = sk.measure.regionprops(labels)
                    props_list = []
                    for obj,rp in enumerate(rps):
                        props_entry = [file_idx, k, t, obj]
                        for prop_key in self.props:
                            if prop_key in self.unpack_dict.keys():
                                prop_split = self.unpack_dict[prop_key]
                                prop_output = rp[prop_key]
                                props_entry += [prop_output[i] for i in range(len(prop_split))]
                            else:
                                props_entry += [rp[prop_key]]
                        props_list.append(props_entry)
                    all_props_list+=props_list




#                     props_list = [[file_idx, k, t, obj]+[getattr(rp, prop_key) for prop_key in self.props] for obj,rp in enumerate(rps)]
#                     all_props_list+=props_list

        property_columns = [self.unpack_dict[prop] if prop in self.unpack_dict.keys() else [prop] for prop in self.props]
        property_columns = [item for sublist in property_columns for item in sublist]

        if self.intensity_channel_list is not None:
            column_list = ['File Index','File Trench Index','timepoints','Objectid','Intensity Channel'] + property_columns
            df_out = pd.DataFrame(all_props_list, columns=column_list).reset_index()
        else:
            column_list = ['File Index','File Trench Index','timepoints','Objectid'] + property_columns
            df_out = pd.DataFrame(all_props_list, columns=column_list).reset_index()

        file_idx = df_out.apply(lambda x: int(f"{x['File Index']:04}{x['File Trench Index']:04}{x['timepoints']:04}{x['Objectid']:02}{self.intensity_channel_dict[x['Intensity Channel']]:02}"), axis=1)

        df_out["File Parquet Index"] = [item for item in file_idx]
        df_out = df_out.set_index("File Parquet Index").sort_index()
        del df_out["index"]

        return df_out

    def analyze_all_files(self,dask_cont):
        df = dd.read_parquet(self.metapath)
        file_list = df["File Index"].unique().compute().tolist()
#         kymo_meta = dd.read_parquet(self.metapath)
#         file_list = kymo_meta["File Index"].unique().tolist()

        delayed_list = []
        for file_idx in file_list:
            df_delayed = delayed(self.get_file_regionprops)(file_idx)
            delayed_list.append(df_delayed.persist())

        ## filtering out non-failed dataframes ##
        all_delayed_futures = []
        for item in delayed_list:
            all_delayed_futures += futures_of(item)
        while any(future.status == "pending" for future in all_delayed_futures):
            sleep(0.1)

        good_delayed = []
        for item in delayed_list:
            if all([future.status == "finished" for future in futures_of(item)]):
                good_delayed.append(item)

        ## compiling output dataframe ##
        df_out = dd.from_delayed(good_delayed).persist()
        df_out["File Parquet Index"] = df_out.index
        df_out = df_out.set_index("File Parquet Index", drop=True, sorted=False)
        df_out = df_out.repartition(partition_size="25MB").persist()

        kymo_df = dd.read_parquet(self.metapath)
        kymo_df["File Merge Index"] = kymo_df["File Parquet Index"]
        kymo_df = kymo_df.set_index("File Merge Index", sorted=True)
        kymo_df = kymo_df.drop(["File Index","File Trench Index","timepoints","File Parquet Index"], axis=1)

        df_out["File Merge Index"] = df_out.apply(lambda x: int(f'{x["File Index"]:04}{x["File Trench Index"]:04}{x["timepoints"]:04}'), axis=1)
        df_out = df_out.reset_index(drop=False)
        df_out = df_out.set_index("File Merge Index", sorted=True)

        df_out = df_out.join(kymo_df)
        df_out = df_out.set_index("File Parquet Index",sorted=True)

        dd.to_parquet(
            df_out,
            self.analysispath,
            engine="fastparquet",
            compression="gzip",
            write_metadata_file=True,
        )

    def export_all_data(self,n_workers=20,memory='8GB'):

        dask_cont = dask_controller(walltime='01:00:00',local=False,n_workers=n_workers,memory=memory,working_directory=self.headpath+"/dask")
        dask_cont.startdask()
        dask_cont.displaydashboard()
        dask_cont.futures = {}

        try:
            self.analyze_all_files(dask_cont)
            dask_cont.shutdown(delete_files=False)
        except:
            dask_cont.shutdown(delete_files=False)
            raise

class kymograph_viewer:
    def __init__(self,headpath,channel,segdir):
        self.headpath = headpath
        self.kymopath = headpath + "/kymograph"
        self.segpath = headpath + "/" + segdir
        self.channel = channel

    def get_kymograph_data(self,file_idx,trench_idx):
        with h5py.File(self.kymopath + "/kymograph_"+str(file_idx)+".hdf5", "r") as infile:
            kymodat = infile[self.channel][trench_idx]
        with h5py.File(self.segpath+"/segmentation_"+str(file_idx)+".hdf5", "r") as infile:
            segdat = infile["data"][trench_idx]
        segdat = np.array([sk.morphology.label(segdat[t],connectivity=1) for t in range(segdat.shape[0])])
        return kymodat,segdat

    def plot_kymograph(self,kymograph):
        """Helper function for plotting kymographs. Takes a kymograph array of
        shape (y_dim,x_dim,t_dim).

        Args:
            kymograph (array): kymograph array of shape (y_dim,x_dim,t_dim).
        """
        list_in_t = [kymograph[t,:,:] for t in range(kymograph.shape[0])]
        img_arr = np.concatenate(list_in_t,axis=1)
        plt.imshow(img_arr)

    def plot_kymograph_data(self,kymodat,segdat,x_size=20,y_size=6):
        fig=plt.figure(figsize=(x_size, y_size))
        fig.add_subplot(2, 1, 1)
        self.plot_kymograph(kymodat)
        fig.add_subplot(2, 1, 2)
        self.plot_kymograph(segdat)
        plt.show()

    def inspect_trench(self,file_idx,trench_idx,x_size=20,y_size=6):
        kymodat,segdat = self.get_kymograph_data(file_idx,trench_idx)
        self.plot_kymograph_data(kymodat,segdat,x_size=x_size,y_size=y_size)

def get_image_measurements(
    kymographpath, channels, file_idx, output_name, img_fn, *args, **kwargs
):

    df = dd.read_parquet(kymographpath + "/metadata")
    df = df.set_index("File Parquet Index",sorted=True)

    start_idx = int(str(file_idx) + "00000000")
    end_idx = int(str(file_idx) + "99999999")

    working_dfs = []

    proc_file_path = kymographpath + "/kymograph_" + str(file_idx) + ".hdf5"
    with h5py.File(proc_file_path, "r") as infile:
        working_filedf = df.loc[start_idx:end_idx].compute()
        trench_idx_list = working_filedf["File Trench Index"].unique().tolist()
        for trench_idx in trench_idx_list:
            trench_df = working_filedf[working_filedf["File Trench Index"] == trench_idx]
            for channel in channels:
                kymo_arr = infile[channel][trench_idx]
                fn_out = [
                img_fn(kymo_arr[i], *args, **kwargs)
                    for i in range(kymo_arr.shape[0])
                ]
                trench_df[channel + " " + output_name] = fn_out
            working_dfs.append(trench_df)


    out_df = pd.concat(working_dfs)
    return out_df


def get_all_image_measurements(
    headpath, output_path, channels, output_name, img_fn, *args, **kwargs
):
    kymographpath = headpath + "/kymograph"
    df = dd.read_parquet(kymographpath + "/metadata")

    file_list = df["File Index"].unique().compute().tolist()

    delayed_list = []
    for file_idx in file_list:
        df_delayed = delayed(get_image_measurements)(
            kymographpath, channels, file_idx, output_name, img_fn, *args, **kwargs
        )
        delayed_list.append(df_delayed.persist())

    ## filtering out non-failed dataframes ##
    all_delayed_futures = []
    for item in delayed_list:
        all_delayed_futures += futures_of(item)
    while any(future.status == "pending" for future in all_delayed_futures):
        sleep(0.1)

    good_delayed = []
    for item in delayed_list:
        if all([future.status == "finished" for future in futures_of(item)]):
            good_delayed.append(item)

    ## compiling output dataframe ##
    df_out = dd.from_delayed(good_delayed).persist()
    df_out["FOV Parquet Index"] = df_out.index
    df_out = df_out.set_index("FOV Parquet Index", drop=True, sorted=False)
    df_out = df_out.repartition(partition_size="25MB").persist()

    dd.to_parquet(
        df_out,
        output_path,
        engine="fastparquet",
        compression="gzip",
        write_metadata_file=True,
    )
