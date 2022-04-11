import numpy as np
import pandas as pd
import matplotlib
import MMCorePy
from IPython.display import clear_output
import matplotlib.pyplot as plt
import time

import h5py


class scopeCore:
    def __init__(
        self,
        configpath,
        logpath,
        camera_name="BSI Prime",
        shutter_name="SpectraIII",
        xystage_name="XYStage",
        focus_name="ZDrive",
    ):
        self.mmc = MMCorePy.CMMCore()
        self.mmc.loadSystemConfiguration(configpath)
        self.mmc.setPrimaryLogFile(logpath)
        self.mmc.setCameraDevice(camera_name)

        self.camera_name = camera_name
        self.shutter_name = shutter_name
        self.xystage_name = xystage_name
        self.focus_name = focus_name

    def snap_image(self, img_size=(12, 12)):
        self.mmc.snapImage()
        im1 = self.mmc.getImage()
        return im1

    def auto_contrast(self, img, low_percentile=0, high_percentile=100):
        low = np.percentile(img, low_percentile)
        high = np.percentile(img, high_percentile)
        return low, high

    def plot_img(self, img, low, high, img_size=(12, 12)):
        clear_output(wait=True)
        plt.figure(figsize=img_size)
        plt.imshow(im1, interpolation="None", vmin=low, vmax=high)
        plt.show()

    def liveview(self, img_size=(12, 12), low=None, high=None):  # W,interval=0.5):
        while True:
            try:
                while self.mmc.deviceBusy(self.camera_name):
                    time.sleep(0.005)

                im1 = self.snap_image()
                clear_output(wait=True)
                plt.figure(figsize=img_size)
                if low == None or high == None:
                    plt.imshow(im1, interpolation="None", cmap="gray")
                else:
                    plt.imshow(
                        im1, interpolation="None", cmap="gray", vmin=low, vmax=high
                    )
                plt.show()
            except KeyboardInterrupt:
                break
        while self.mmc.deviceBusy(self.camera_name):
            time.sleep(0.01)

    def set_grid(self, num_col, num_row, col_step=333.0, row_step=686.0):
        grid_coords = []

        x_ori, y_ori = self.mmc.getXYPosition()
        start_left = True

        for row in range(num_row):
            y_disp = row * row_step

            if start_left:
                for col in range(num_col):
                    x_disp = col * (-col_step)
                    current_coord = (x_ori + x_disp, y_ori + y_disp)
                    grid_coords.append(current_coord)
                start_left = False
            else:
                for col in range(num_col):
                    x_disp = (num_col - col - 1) * (-col_step)
                    current_coord = (x_ori + x_disp, y_ori + y_disp)
                    grid_coords.append(current_coord)
                start_left = True

        return grid_coords

    def multipoint_aq(
        self,
        grid_coords,
        config_list,
        timepoint,
        output_folder="./",
        group_name="FISH_channels",
    ):

        ### Make sure configs are valid ###
        undefined_configs = []
        for config in config_list:
            config_defined = self.mmc.isConfigDefined(group_name, config)
            if not config_defined:
                undefined_configs.append(config)

        if len(undefined_configs) > 0:
            raise ValueError(
                "The following configs are undefined: " + ", ".join(undefined_configs)
            )

        ### Gather basic metadata ###

        t_start = time.time()
        x_dim = self.mmc.getProperty(self.camera_name, "X-dimension")
        y_dim = self.mmc.getProperty(self.camera_name, "Y-dimension")

        ## Note change the write to disk later ##

        imgs = []
        imgs_metadata = []

        x_coord, y_coord = grid_coords[0]
        self.mmc.setXYPosition(x_coord, y_coord)

        for fov_num, (x_coord, y_coord) in enumerate(grid_coords):
            while self.mmc.deviceBusy(self.xystage_name):
                time.sleep(0.1)
                pass
            self.mmc.setXYPosition(x_coord, y_coord)

            for config in config_list:
                while self.mmc.systemBusy():
                    time.sleep(0.1)
                    pass

                self.mmc.setConfig(group_name, config)

                #                 ### put write here because it is likely the slow step ###

                #                 for img_num in range(len(imgs)):
                #                     img = imgs[img_num]
                #                     metadata_entry = imgs_metadata[img_num]

                #                     with h5py.File(output_folder + "fov=" + str(metadata_entry["fov"]) + "_config=" + str(metadata_entry["config"]) + "_t=" + str(timepoint),"w") as h5pyfile:
                #                         hdf5_dataset = h5pyfile.create_dataset("data", data=img, chunks=(128,128), dtype='uint16')
                #                         all_metadata.append(metadata_entry)

                while self.mmc.systemBusy():
                    time.sleep(0.1)
                    pass

                if "noPFS" not in config:
                    self.mmc.setProperty("PFS", "FocusMaintenance", "On")
                    time.sleep(0.25)

                self.mmc.setShutterOpen(self.shutter_name, True)
                self.mmc.snapImage()
                self.mmc.setShutterOpen(self.shutter_name, False)

                read_x_coord, read_y_coord = self.mmc.getXYPosition(self.xystage_name)
                read_z_coord = self.mmc.getPosition(self.focus_name)
                current_time = time.time() - t_start

                metadata_entry = {
                    "fov": fov_num,
                    "config": config,
                    "x": read_x_coord,
                    "y": read_y_coord,
                    "z": read_z_coord,
                    "t": current_time,
                }
                img = self.mmc.getImage()

                imgs.append(img)
                imgs_metadata.append(metadata_entry)
            self.mmc.setConfig(group_name, config_list[0])
            if "noPFS" not in config_list[0]:
                self.mmc.setProperty("PFS", "FocusMaintenance", "On")
                time.sleep(0.25)
        self.mmc.setConfig(group_name, config_list[0])
        if "noPFS" not in config_list[0]:
            self.mmc.setProperty("PFS", "FocusMaintenance", "On")
            time.sleep(0.25)
        x_coord, y_coord = grid_coords[0]
        self.mmc.setXYPosition(x_coord, y_coord)

        for img_num in range(len(imgs)):
            img = imgs[img_num]
            metadata_entry = imgs_metadata[img_num]

            with h5py.File(
                output_folder
                + "fov="
                + str(metadata_entry["fov"])
                + "_config="
                + str(metadata_entry["config"])
                + "_t="
                + str(timepoint),
                "w",
            ) as h5pyfile:
                hdf5_dataset = h5pyfile.create_dataset(
                    "data", data=img, chunks=(128, 128), dtype="uint16"
                )

        metadata = pd.DataFrame.from_dict(imgs_metadata)
        metadata.to_hdf(
            output_folder + "metadata_" + str(timepoint) + ".hdf5", key="data", mode="w"
        )
