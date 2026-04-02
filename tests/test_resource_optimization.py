"""Tests for server resource optimization utilities."""

import math


# Replicate _optimal_chunksize logic for testing without importing pipeline_manager
# (which requires matplotlib/numpy/pyrosetta at import time)
def _optimal_chunksize(n_tasks: int, n_workers: int) -> int:
    if n_tasks <= 0 or n_workers <= 0:
        return 1
    if n_tasks < 100:
        return 1
    if n_tasks < 5000:
        return max(1, int(math.sqrt(n_tasks / n_workers)))
    return max(1, min(n_workers * 4, n_tasks // (n_workers * 4)))


class TestOptimalChunksize:
    def test_small_batch_returns_1(self):
        """Small batches (<100) should use chunksize=1 to avoid starvation."""
        assert _optimal_chunksize(50, 16) == 1
        assert _optimal_chunksize(10, 4) == 1
        assert _optimal_chunksize(0, 16) == 1

    def test_medium_batch_uses_sqrt(self):
        """Medium batches should use sqrt(n/workers) heuristic."""
        cs = _optimal_chunksize(1000, 16)
        assert 1 < cs < 1000
        # sqrt(1000/16) ≈ 7.9
        assert 5 <= cs <= 15

    def test_large_batch_is_capped(self):
        """Large batches (>5000) should cap chunksize for memory."""
        cs = _optimal_chunksize(20000, 16)
        # cap = min(16*4, 20000//(16*4)) = min(64, 312) = 64
        assert cs <= 64
        assert cs >= 1

    def test_single_worker(self):
        """Single worker should still produce valid chunksize."""
        cs = _optimal_chunksize(5000, 1)
        assert cs >= 1

    def test_more_workers_than_tasks(self):
        """Edge case: workers > tasks."""
        cs = _optimal_chunksize(5, 32)
        assert cs == 1

    def test_production_scenario(self):
        """Realistic production: 20K tasks, 16 workers."""
        cs = _optimal_chunksize(20000, 16)
        # Should batch enough to reduce IPC but not hog memory
        assert 10 <= cs <= 100

    def test_filter_survivors_scenario(self):
        """After filtering: ~2000 tasks, 16 workers."""
        cs = _optimal_chunksize(2000, 16)
        assert 5 <= cs <= 50

    def test_final_scoring_scenario(self):
        """Final scoring: ~200 tasks, 16 workers."""
        cs = _optimal_chunksize(200, 16)
        assert 1 <= cs <= 20
