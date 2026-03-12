import importlib.util
import sys
import types
from pathlib import Path


MODULE_PATH = (
    Path(__file__).resolve().parents[1] / "docs" / "md_interface_analysis.py"
)


class _FakeArray(list):
    def mean(self):
        return sum(self) / len(self) if self else 0.0

    def std(self):
        if not self:
            return 0.0
        mu = self.mean()
        return (sum((value - mu) ** 2 for value in self) / len(self)) ** 0.5


def _load_md_interface_analysis_module():
    fake_numpy = types.ModuleType("numpy")
    fake_numpy.ndarray = _FakeArray
    fake_numpy.array = lambda values: _FakeArray(values)
    fake_numpy.argwhere = lambda values: []
    fake_numpy.column_stack = lambda values: values
    fake_numpy.savetxt = lambda *args, **kwargs: None

    fake_mda = types.ModuleType("MDAnalysis")
    fake_mda.Universe = object
    fake_analysis = types.ModuleType("MDAnalysis.analysis")
    fake_analysis.distances = types.SimpleNamespace(distance_array=lambda *args, **kwargs: [])

    previous = {
        "numpy": sys.modules.get("numpy"),
        "MDAnalysis": sys.modules.get("MDAnalysis"),
        "MDAnalysis.analysis": sys.modules.get("MDAnalysis.analysis"),
    }
    sys.modules["numpy"] = fake_numpy
    sys.modules["MDAnalysis"] = fake_mda
    sys.modules["MDAnalysis.analysis"] = fake_analysis

    try:
        spec = importlib.util.spec_from_file_location(
            "test_md_interface_analysis_module",
            MODULE_PATH,
        )
        module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(module)
        return module
    finally:
        for name, original in previous.items():
            if original is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = original


class _FakeTs:
    def __init__(self, frame, time):
        self.frame = frame
        self.time = time


class _FakeTrajectory:
    def __init__(self, times):
        self._frames = [_FakeTs(index, time) for index, time in enumerate(times)]

    def __getitem__(self, item):
        return self._frames[item]


class _FakeAtomGroup:
    def __init__(self, resids):
        self.resids = resids


class _FakeUniverse:
    def __init__(self, times, partner_resids):
        self.trajectory = _FakeTrajectory(times)
        self._partner_resids = partner_resids

    def select_atoms(self, selection):
        return _FakeAtomGroup(self._partner_resids)


def test_compute_residue_occupancy_respects_start_time_ps():
    module = _load_md_interface_analysis_module()
    universe = _FakeUniverse(
        times=[0.0, 50_000.0, 100_000.0, 150_000.0],
        partner_resids=[101, 102],
    )
    contacts_by_frame = {
        0: {(1, 101)},
        1: set(),
        2: {(1, 102)},
        3: {(1, 102)},
    }

    def fake_get_residue_contact_pairs(universe_obj, receptor_sel, partner_sel, cutoff, frame_idx=0):
        return contacts_by_frame[frame_idx]

    module.get_residue_contact_pairs = fake_get_residue_contact_pairs

    full_occupancy = module.compute_residue_occupancy(
        universe,
        "segid A",
        "segid B",
        stride=1,
    )
    windowed_occupancy = module.compute_residue_occupancy(
        universe,
        "segid A",
        "segid B",
        stride=1,
        start_time_ps=100_000.0,
    )

    assert abs(full_occupancy[101] - 0.25) < 1e-9
    assert abs(full_occupancy[102] - 0.5) < 1e-9
    assert abs(windowed_occupancy[101] - 0.0) < 1e-9
    assert abs(windowed_occupancy[102] - 1.0) < 1e-9
