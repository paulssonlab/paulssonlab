import numpy as np


def point_linspace(anchor0, anchor1, num_points):
    for s in np.linspace(0, 1, num_points)[1:-1]:
        anchor = (1 - s) * anchor0 + s * anchor1
        yield anchor


def coords_along(x0, x1):
    # TODO: is floor or trunc correct?
    # do we care, since coords should be positive?
    x0 = np.floor(x0)
    x1 = np.floor(x1)
    # TODO: are ceil, length=0 -> 1 correct??
    length = int(np.ceil(np.sqrt(np.sum((x1 - x0) ** 2))))
    if length == 0:
        length = 1
    xs = np.linspace(x0[0], x1[0], length).astype(np.int_)  # [1:-1]
    ys = np.linspace(x0[1], x1[1], length).astype(np.int_)  # [1:-1]
    return xs, ys


def edge_point(x0, theta, x_lim, y_lim):
    x_min, x_max = x_lim
    y_min, y_max = y_lim
    # TODO: hack to fix getting edge points where x0 is on the border
    if x0[0] in x_lim and x0[1] in y_lim:
        return x0
    theta = theta % (2 * np.pi)
    if 0 <= theta < np.pi / 2:
        print("a")
        corner_x, corner_y = x_min, y_max
    elif np.pi / 2 <= theta < np.pi:
        print("b")
        corner_x, corner_y = x_max, y_max
    elif np.pi <= theta < 3 / 2 * np.pi:
        print("c")
        corner_x, corner_y = x_max, y_min
    elif 3 / 2 * np.pi <= theta:
        print("d")
        corner_x, corner_y = x_min, y_min
    angle_to_corner = np.arctan2(corner_y - x0[1], x0[0] - corner_x) % (2 * np.pi)
    print("theta", theta, "angle_to_corner", angle_to_corner)
    for expr in [
        "theta >= angle_to_corner",
        "0 < theta <= np.pi / 2",
        "theta < angle_to_corner",
        "np.pi / 2 <= theta < np.pi",
        "theta >= angle_to_corner",
        "np.pi < theta <= 3 / 2 * np.pi",
        "theta < angle_to_corner",
        "3 / 2 * np.pi <= theta < 2 * np.pi",
    ]:
        print("***", expr, eval(expr))
    if (
        (theta >= angle_to_corner and 0 < theta <= np.pi / 2)
        or (theta < angle_to_corner and np.pi / 2 <= theta < np.pi)
        or (theta >= angle_to_corner and np.pi < theta <= 3 / 2 * np.pi)
        or (theta < angle_to_corner and 3 / 2 * np.pi <= theta < 2 * np.pi)
    ):
        print("X1a")
        # top/bottom
        x1 = np.array([x0[0] - (corner_y - x0[1]) / np.tan(theta), corner_y])
    else:
        print("X1b")
        # left/right
        x1 = np.array([corner_x, x0[1] - (corner_x - x0[0]) * np.tan(theta)])
    print(
        np.array([x0[0] - (corner_y - x0[1]) / np.tan(theta), corner_y]),
        "|",
        np.array([corner_x, x0[1] - (corner_x - x0[0]) * np.tan(theta)]),
    )
    return x1
