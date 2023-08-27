import numpy as np


def angled_vector(angle):
    return np.array([np.sin(angle), -np.cos(angle)])


def _vdot(a, b):
    return np.einsum("ij,ij->i", a, b)


def _cosine_similarity(a, b):
    return _vdot(a, b) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))


def intersect_lines_with_segment(points, directions, line_segment):
    line_segment_direction = (line_segment[1] - line_segment[0])[np.newaxis, :]
    perp = np.cross(directions, line_segment_direction)
    vector_ab = np.asarray(line_segment[0]) - points
    num = np.cross(vector_ab, line_segment_direction) * perp
    denom = perp**2
    intersections = points + (num / denom)[:, np.newaxis] * directions
    vector_ia = line_segment[0] - intersections
    vector_ib = line_segment[1] - intersections
    similarity = _cosine_similarity(vector_ia, vector_ib)
    mask = np.logical_or(
        np.isclose(similarity, -1),
        np.logical_or(
            np.all(np.isclose(vector_ia, 0), axis=1),
            np.all(np.isclose(vector_ib, 0), axis=1),
        ),
    )
    return np.where(
        mask[:, np.newaxis],
        intersections,
        np.nan,
    )
