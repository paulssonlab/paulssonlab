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
