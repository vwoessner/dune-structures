import pytest

import math
import numpy as np

from dune.structures.fibergrowth.line import (
    Line,
    LineFunc,
    angle_between,
    rotate,
    three_line_integral,
)

# --- Fixtures ---


@pytest.fixture
def x_axis_line():
    return Line([0.0, 0.0], [1.0, 0.0])


@pytest.fixture
def y_axis_line():
    return Line([0.0, 0.0], [0.0, 1.0])


# --- Tests ---


def test_line_init():
    line = Line([0.0], [1.0])
    assert np.allclose(line.start, [0.0])
    assert np.allclose(line.end, [1.0])

    # Flipping
    line = Line([1.0], [0.0])
    assert np.allclose(line.start, [0.0])
    assert np.allclose(line.end, [1.0])

    with pytest.raises(RuntimeError):
        Line([0.0, 1.0], [0.0])


def test_line_length(x_axis_line):
    assert math.isclose(x_axis_line.length, 1.0)
    assert math.isclose(Line([0.0, 0.0], [1.0, 1.0]).length, np.sqrt(2))


def test_line_direction(x_axis_line, y_axis_line):
    assert np.allclose(x_axis_line.direction, x_axis_line.end)
    assert np.allclose(y_axis_line.direction, y_axis_line.end)
    assert np.allclose(Line([1.0, 0.0], [3.0, 0.0]).direction, [1.0, 0.0])


def test_angle(x_axis_line, y_axis_line):
    assert math.isclose(angle_between(x_axis_line, x_axis_line), 0.0)
    assert math.isclose(angle_between(x_axis_line, y_axis_line), math.pi / 2)
    assert math.isclose(angle_between(y_axis_line, x_axis_line), math.pi / 2)


def test_line_func_init(x_axis_line):
    linef = LineFunc(x_axis_line.start, x_axis_line.end)
    assert np.allclose(linef.line.start, x_axis_line.start)
    assert np.allclose(linef.line.end, x_axis_line.end)

    linef = LineFunc.from_line(x_axis_line)
    assert np.allclose(linef.line.start, x_axis_line.start)
    assert np.allclose(linef.line.end, x_axis_line.end)


def test_line_func_eval(x_axis_line):
    linef = LineFunc.from_line(x_axis_line)
    assert np.allclose(linef([0.0, 1.0, 2.0]), [0.0, 0.0, np.nan], equal_nan=True)

    linef = LineFunc([0.0, 0.0], [2.0, 1.0])
    assert np.allclose(
        linef([-1.0, 0.0, 1.0, 2.0, 3.0]),
        [np.nan, 0.0, 0.5, 1.0, np.nan],
        equal_nan=True,
    )

    # Check line NOT starting at origin
    linef = LineFunc([1.0, 0.0], [3.0, 1.0])
    assert np.allclose(
        linef([0.0, 1.0, 2.0, 3.0, 4.0]),
        [np.nan, 0.0, 0.5, 1.0, np.nan],
        equal_nan=True,
    )

    linef = LineFunc([2.0, 1.0], [3.0, -1.0], absolute=True)
    assert np.allclose(linef([2.0, 2.5, 3.0]), [1.0, 0.0, 1.0])


def test_line_func_eval_abs():
    linef = LineFunc([0.0, -1.0], [1.0, 1.0])
    assert np.allclose(linef([0.0, 0.5, 1.0]), [-1.0, 0.0, 1.0])

    linef = LineFunc([0.0, -1.0], [1.0, 1.0], absolute=True)
    assert np.allclose(linef([0.0, 0.5, 1.0]), [1.0, 0.0, 1.0])


def test_line_rotate(x_axis_line, y_axis_line):
    new = rotate(x_axis_line, math.pi / 2)
    assert np.allclose([new.start, new.end], [y_axis_line.start, y_axis_line.end])

    new = rotate(y_axis_line, -math.pi / 2)
    assert np.allclose([new.start, new.end], [x_axis_line.start, x_axis_line.end])


def test_three_line_integral():
    line_ref = Line([0.0, 0.0], [3.0, 0.0])
    line1 = Line([0.0, 0.0], [1.0, 1.0])
    line2 = Line([2.0, 1.0], [3.0, -1.0])

    assert math.isclose(three_line_integral(line_ref, line1, line2), 1.0)

    # Rotate by 45 deg
    line_ref, line1, line2 = (
        rotate(line_ref, math.pi / 4),
        rotate(line1, math.pi / 4),
        rotate(line2, math.pi / 4),
    )
    assert math.isclose(three_line_integral(line_ref, line1, line2), 1.0)

    # Rotate by -45 deg
    line_ref, line1, line2 = (
        rotate(line_ref, -math.pi / 4),
        rotate(line1, -math.pi / 4),
        rotate(line2, -math.pi / 4),
    )
    assert math.isclose(three_line_integral(line_ref, line1, line2), 1.0)
