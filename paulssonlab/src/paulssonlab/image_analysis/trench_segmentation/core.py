# import numpy as np
# import pandas as pd
# import zarr
# import skimage
# from functools import partial
# from geometry import get_image_limits
# from util import tqdm_auto, map_collections, iterate_get_collection_value

# def trenchwise_apply(img_stack, trench_set_points, trench_idx, func, channel_slice=slice(None), time_slice=slice(None)):
#     x_lim, y_lim = get_img_limits(img_stack.shape[-2:])
#     ul, lr = get_trench_bbox(trench_set_points, trench_idx, x_lim, y_lim)
#     if len(img_stack.shape) == 4:
#         thumb_stack = (img_stack.oindex if isinstance(img_stack, zarr.Array) else img_stack)[channel_slice,time_slice,ul[1]:lr[1],ul[0]:lr[0]]
#     else:
#         thumb_stack = img_stack[...,ul[1]:lr[1],ul[0]:lr[0]]
#     res = func(thumb_stack)
#     return res

# def trenchwise_map(img_stack, trench_points, func, progress_bar=tqdm_auto, preload=True, **kwargs):
#     if preload:
#         img_stack = (img_stack.oindex if isinstance(img_stack, zarr.Array) else img_stack)[kwargs.get('channel_slice', slice(None)),kwargs.get('time_slice', slice(None)),:,:]
#         del kwargs['channel_slice']
#         del kwargs['time_slice']
#     obj = {trench_set_idx: {trench_idx: trenchwise_apply(img_stack, trench_set_points, trench_idx, func, **kwargs) for trench_idx in progress_bar(range(len(trench_set_points[0])))}
#                for trench_set_idx, trench_set_points in progress_bar(trench_points.items())}
#     representative_obj = iterate_get_collection_value(obj, 2)
#     if isinstance(representative_obj, (pd.Series, pd.DataFrame)):
#         df = map_collections(partial(pd.concat, axis=1), obj, max_level=2)
#         df.columns.set_names('trench_set_idx', level=0, inplace=True)
#         df.columns.set_names('trench_idx', level=1, inplace=True)
#         return df
#     else:
#         return obj

# def positionwise_trenchwise_map(img_group, trench_points_pos, func, positions=None, progress_bar=tqdm_auto, **kwargs):
#     if positions is None:
#         positions = trench_points_pos.keys()
#     obj = {pos: trenchwise_map(img_group[pos], trench_points_pos[pos], func, progress_bar=lambda x: x, **kwargs)
#                for pos in progress_bar(positions)}
#     if isinstance(iterate_get_collection_value(obj, 1), (pd.Series, pd.DataFrame)):
#         df = pd.concat(obj, axis=1)
#         df.columns.set_names('position', level=0, inplace=True)
#         return df
#     else:
#         return obj
