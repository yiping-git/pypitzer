import hashlib
import sqlite3

import pytest

from pypitzer._database import DATABASE_SHA256, ParameterRepository


def test_database_snapshot_checksum_and_counts():
    repository = ParameterRepository()
    assert hashlib.sha256(repository.path.read_bytes()).hexdigest() == DATABASE_SHA256
    with repository._connect() as connection:
        counts = {
            table: connection.execute(f"SELECT count(*) FROM {table}").fetchone()[0]
            for table in ("binary_params", "ternary_params", "reaction")
        }
    assert counts == {"binary_params": 96, "ternary_params": 123, "reaction": 70}


def test_database_connection_is_read_only():
    repository = ParameterRepository()
    with repository._connect() as connection:
        with pytest.raises(sqlite3.OperationalError):
            connection.execute("DELETE FROM binary_params")


def test_phase_lookup_resolves_internal_reaction():
    reaction = ParameterRepository().phase_reaction("sylvite", 298.15)
    assert reaction.equation == "KCl(s) = Cl-(aq) + K+(aq)"
    assert reaction.stoichiometry["KCl(s)"] == -1.0
    assert reaction.stoichiometry["K+(aq)"] == 1.0
