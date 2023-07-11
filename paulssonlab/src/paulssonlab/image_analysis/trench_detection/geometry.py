import numpy as np
from skspatial.objects import Line


def angled_line(point, angle):
    return Line(point, (np.sin(angle), -np.cos(angle)))


def intersect_line_with_segment(line, line_segment):
    line2 = Line.from_points(line_segment.point_a, line_segment.point_b)
    try:
        intersection = line.intersect_line(line2)
    except:
        return None
    if line_segment.contains_point(intersection):
        return intersection
    else:
        return None


def trench_anchors(angle, anchor_rho, rho_min, rho_max, x_lim, y_lim):
    x_min, x_max = x_lim
    y_min, y_max = y_lim
    if angle < 0:
        anchor_rho = rho_max - anchor_rho
    anchors = (
        anchor_rho[:, np.newaxis]
        * np.array((np.cos(angle), np.sin(angle)))[np.newaxis, :]
    )
    if angle < 0:
        upper_right = np.array((x_max, 0))
        anchors = upper_right - anchors
    return anchors


def edge_point(x0, theta, x_lim, y_lim):
    x_min, x_max = x_lim
    y_min, y_max = y_lim
    theta = theta % (2 * np.pi)
    if 0 <= theta < np.pi / 2:
        corner_x, corner_y = x_min, y_max
    elif np.pi / 2 <= theta < np.pi:
        corner_x, corner_y = x_max, y_max
    elif np.pi <= theta < 3 / 2 * np.pi:
        corner_x, corner_y = x_max, y_min
    elif 3 / 2 * np.pi <= theta:
        corner_x, corner_y = x_min, y_min
    angle_to_corner = np.arctan2(corner_y - x0[1], x0[0] - corner_x) % (2 * np.pi)
    if (
        (theta >= angle_to_corner and 0 < theta <= np.pi / 2)
        or (theta < angle_to_corner and np.pi / 2 <= theta < np.pi)
        or (theta >= angle_to_corner and np.pi < theta <= 3 / 2 * np.pi)
        # need to shift so that when angle_to_corner=0, theta < angle_to_corner is evaluated correctly
        or (
            (
                (theta - 3 / 2 * np.pi) % (2 * np.pi)
                < (angle_to_corner - 3 / 2 * np.pi) % (2 * np.pi)
            )
            and 3 / 2 * np.pi <= theta < 2 * np.pi
        )
    ):
        # top/bottom
        x1 = np.array([x0[0] - (corner_y - x0[1]) / np.tan(theta), corner_y])
    else:
        # left/right
        x1 = np.array([corner_x, x0[1] - (corner_x - x0[0]) * np.tan(theta)])
    return x1
