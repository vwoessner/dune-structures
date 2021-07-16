import numpy as np
import scipy.integrate as sint


def angle_between(line1, line2):
    # NOTE: `direction` is normed
    return np.arccos(np.dot(line1.direction, line2.direction))


def rotate(line, angle):
    """Counter-clockwise rotation of a line.

    Returns:
        Line: Rotated line
    """
    R = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    start_new = np.dot(R, line.start)
    end_new = np.dot(R, line.end)

    return Line(start_new, end_new)


def three_line_integral(line_ref, line1, line2):
    # Shift lines to origin
    ref_start = line_ref.start
    line1 = Line(line1.start - ref_start, line1.end - ref_start)
    line2 = Line(line2.start - ref_start, line2.end - ref_start)
    line_ref = Line(line_ref.start - ref_start, line_ref.end - ref_start)

    # Find rotation angle to match reference line with x-axis
    def rotate_cc(line):
        if line.start[1] > line.end[1]:
            return True
        else:
            return False

    x_axis = Line([0.0, 0.0], [1.0, 0.0])
    rot_angle = angle_between(line_ref, x_axis)
    rot_angle = rot_angle if rotate_cc(line_ref) is True else -rot_angle

    # Rotate lines
    line_ref, line1, line2 = (
        rotate(line_ref, rot_angle),
        rotate(line1, rot_angle),
        rotate(line2, rot_angle),
    )

    # Integrate absolute areas
    line_ref, line1, line2 = (
        LineFunc.from_line(line_ref, absolute=True),
        LineFunc.from_line(line1, absolute=True),
        LineFunc.from_line(line2, absolute=True),
    )

    def integrate_min_range(line_func1, line_func2):
        low = max(line_func1.line.start[0], line_func2.line.start[0])
        high = min(line_func1.line.end[0], line_func2.line.end[0])

        return sint.quad(lambda x: np.abs(line_func1(x) - line_func2(x)), low, high)[0]

    return integrate_min_range(line_ref, line1) + integrate_min_range(line_ref, line2)


class Line:
    def __init__(self, start, end):
        start = np.array(start)
        end = np.array(end)
        if start.shape != end.shape:
            raise RuntimeError("Input shapes do not match!")

        # Assert orientation
        if start[0] > end[0]:
            cache = np.copy(end)
            end = np.copy(start)
            start = cache

        self.start = start
        self.end = end

    @property
    def length(self):
        return np.sqrt(np.sum((self.end - self.start) ** 2))

    @property
    def direction(self):
        return (self.end - self.start) / self.length

    def func(self):
        """Return a functional representation of the line in the fashion ``y = f(x)``.
        Returns ``None`` if ``x`` is out of bounds of the line.
        """
        # x_0, x = np.zeros_like(self.start), np.zeros_like(self.start)
        # x[0] = 1.0
        # ang = angle_between(self, Line(start=x_0, end=x))
        # return (
        #     lambda x, ang=ang: self.start[1] + np.arctan(ang) * x
        #     if self.start[0] < x < self.end[0]
        #     else None
        # )
        return LineFunc.from_line(self)


class LineFunc:
    def __init__(self, start, end, absolute=False):
        self.line = Line(start, end)
        self._abs = absolute

        # Determine slope
        x_0, x = np.zeros_like(self.line.start), np.zeros_like(self.line.start)
        x[0] = 1.0

        slope_angle = angle_between(self.line, Line(start=x_0, end=x))
        slope_angle = (
            slope_angle if self.line.start[1] < self.line.end[1] else -slope_angle
        )
        self._slope = np.tan(slope_angle)

    @classmethod
    def from_line(cls, line, absolute=False):
        return LineFunc(line.start, line.end, absolute)

    def __call__(self, x):
        def maybe_abs(val):
            if self._abs:
                return np.abs(val)
            else:
                return val

        x = np.asarray(x)
        return np.where(
            np.logical_and(self.line.start[0] <= x, x <= self.line.end[0]),
            maybe_abs(self.line.start[1] + self._slope * (x - self.line.start[0])),
            np.nan,
        )
