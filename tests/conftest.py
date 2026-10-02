"""Shared test configuration."""
import pytest


@pytest.fixture(autouse=True)
def _no_ncbi_network(monkeypatch):
    """
    Tests must never reach NCBI: an online BLASTn test that forgets to mock the
    client would otherwise send real searches (counting against the user's
    daily NCBI quota) and sleep through NCBI's 20–60 s polling intervals.
    """
    def _blocked(*args, **kwargs):
        raise RuntimeError("network access to NCBI is blocked in tests — mock NcbiBlastClient")

    # Blocking the transport is enough: the client fails on its first request
    # (the submission), before any polling sleep. time.sleep is NOT patched —
    # that would be global and break tests that measure real delays.
    monkeypatch.setattr("viralquest.blastn.urlopen", _blocked)
