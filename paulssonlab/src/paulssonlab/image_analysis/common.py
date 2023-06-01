from functools import wraps

import holoviews as hv

# make sure holoviews backend is selected
# (TOOD: allow user to set using hv.extension in the notebook, transmit to dask workers),
hv.extension("bokeh")
# always pickle options with holoviews elements
hv.Store.save_option_state = True

orig_setstate = hv.core.dimension.LabelledData.__setstate__


# monkey-patch hv.Dimension to always load custom options if they are pickled using hv.Store.save_option_state = True
@wraps(orig_setstate)
def new_setstate(self, d):
    hv.Store.load_counter_offset = hv.StoreOptions.id_offset()
    orig_setstate(self, d)
    hv.Store.load_counter_offset = None


hv.core.dimension.LabelledData.__setstate__ = new_setstate
