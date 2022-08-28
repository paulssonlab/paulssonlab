from tqdm.auto import tqdm

# FROM: https://gist.github.com/drs251/c0c348e88a1cc95e3e70205bff0ed91b
def holomap_to_video(hmap, out=None, fps=10, dpi=200, size=None, progress_bar=tqdm):
    if out is None:
        out = io.BytesIO()
    renderer = hv.renderer("matplotlib")
    if dpi:
        old_dpi = renderer.dpi
        renderer.dpi = dpi
    if size:
        old_size = renderer.size
        renderer.size = size
    writer = imageio.get_writer(out, fps=fps, format="mp4")
    if isinstance(hmap, hv.Layout) or isinstance(hmap, hv.NdLayout):
        kdim = hmap[hmap.keys()[0]].kdims[0].name  # TODO: make more elegant/robust
        keys = hmap[hmap.keys()[0]].keys()
        key_items = [(k, hmap.select(**{kdim: k})) for k in keys]
    else:
        key_items = hmap.items()
    if progress_bar:
        key_items = progress_bar(key_items)
    for key, item in key_items:
        canvas = renderer.get_plot(item).state.canvas
        canvas.draw()
        ary = np.array(canvas.get_renderer()._renderer)
        writer.append_data(ary)
    writer.close()
    if dpi:
        renderer.dpi = old_dpi
    if size:
        renderer.size = old_size
    return out


def trench_movie(trench_bboxes, key, channel, ts):
    get_trench_frame = lambda t: ui._trench_img(
        workflow.get_trench_image(
            trench_bboxes, key[0], key[1], channel, t, key[2], key[3]
        )
    )
    movie = hv.HoloMap({t: get_trench_frame(t) for t in tqdm_notebook(ts)})
    return movie
