"""Tests for the GDAL bootstrapper.

Everything here is offline and installs nothing.  The parts that would touch
the network or a package manager are exercised through :func:`make_plan`,
which decides what to do without doing any of it -- that separation is the
main reason the module is shaped the way it is.
"""

from __future__ import annotations

import io
import json
import os
import subprocess
import sys
import tarfile

import pytest

from TopoAnalysis import _gdal_source, gdal_setup


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _isolate_gdal_environment():
    """Undo every environment change these tests can cause.

    ``activate_environment`` assigns to ``os.environ`` directly, so
    ``monkeypatch`` only restores the variables it was itself asked to touch.
    ``PROJ_DATA`` in particular was being left behind pointing at a fixture's
    deliberately-invalid ``proj.db``, which then broke unrelated tests --
    invisibly on GDAL 3.4, which ignores ``PROJ_DATA``, and visibly on 3.13,
    which does not.
    """
    names = ("PROJ_LIB", "PROJ_DATA", "GDAL_DATA", "TOPOANALYSIS_NO_GDAL_ENV")
    saved = {name: os.environ.get(name) for name in names}
    yield
    for name, value in saved.items():
        if value is None:
            os.environ.pop(name, None)
        else:
            os.environ[name] = value


def make_probe(**fields):
    """A Probe with the given raw fields and its problems worked out."""
    config = fields.pop("_config", None)
    found = gdal_setup.Probe(**fields)
    found.problems = tuple(gdal_setup._diagnose(found, config))
    return found


WORKING = {
    "bindings_version": "3.9.3",
    "has_gdal_array": True,
    "projection_works": True,
    "gdal_config": "/opt/env/bin/gdal-config",
    "_config": {"version": "3.9.3", "prefix": "/nonexistent-prefix", "path": "/opt/env/bin/gdal-config"},
}


# ---------------------------------------------------------------------------
# Version handling
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "text, expected",
    [
        ("3.9.3", (3, 9, 3)),
        ("3.4.1", (3, 4, 1)),
        ("3.10", (3, 10, 0)),
        ("3.13.3dev", (3, 13, 3)),
        ("", None),
        (None, None),
        ("not a version", None),
    ],
)
def test_parse_version(text, expected):
    assert gdal_setup._parse_version(text) == expected


def test_bindings_requirements_prefer_an_exact_match():
    """The bindings must match the C library; exact first, then the series."""
    assert gdal_setup._bindings_requirements("3.4.1") == ["gdal==3.4.1", "gdal==3.4.*"]


def test_bindings_requirements_without_a_version():
    assert gdal_setup._bindings_requirements(None) == ["gdal"]


# ---------------------------------------------------------------------------
# Diagnosis
# ---------------------------------------------------------------------------


def test_a_working_install_has_no_problems():
    found = make_probe(**WORKING)
    assert found.ok
    assert found.problems == ()
    assert "GDAL is installed and working." in found.describe()


def test_missing_bindings_and_library_are_both_reported():
    found = make_probe(bindings_error="ModuleNotFoundError: osgeo", _config=None)
    assert set(found.codes) == {"no-bindings", "no-library"}


def test_missing_numpy_bridge_is_caught():
    """The trap this module exists for: osgeo imports but ReadAsArray fails."""
    fields = dict(WORKING, has_gdal_array=False, gdal_array_error="ImportError: numpy")
    found = make_probe(**fields)
    assert found.has("no-gdal-array")
    assert any("ReadAsArray" in message for _, message in found.problems)


def test_version_mismatch_between_bindings_and_library():
    fields = dict(WORKING, bindings_version="3.4.1")
    found = make_probe(**fields)
    assert found.has("version-mismatch")


def test_patch_level_differences_are_not_a_mismatch():
    """GDAL does not break its C API within a minor version."""
    fields = dict(WORKING, bindings_version="3.9.1")
    found = make_probe(**fields)
    assert not found.has("version-mismatch")


def test_broken_projection_lookup_is_reported():
    fields = dict(WORKING, projection_works=False, projection_error="proj.db missing")
    found = make_probe(**fields)
    assert found.has("no-proj-data")


def test_library_older_than_supported():
    fields = dict(
        WORKING,
        bindings_version="2.4.4",
        _config={"version": "2.4.4", "prefix": "/nonexistent-prefix", "path": "/x/gdal-config"},
    )
    found = make_probe(**fields)
    assert found.has("old-library")


def test_describe_survives_a_completely_empty_probe():
    found = make_probe()
    assert any("not installed" in line for line in found.describe())


# ---------------------------------------------------------------------------
# Choosing a strategy
# ---------------------------------------------------------------------------


def test_a_working_install_plans_nothing():
    plan = gdal_setup.make_plan(make_probe(**WORKING))
    assert plan.strategy == "none"
    assert plan.steps == ()


def test_an_existing_library_is_reused_rather_than_reinstalled(monkeypatch):
    """Only the bindings are missing, so only the bindings get built."""
    monkeypatch.setattr(
        gdal_setup,
        "_gdal_config_info",
        lambda *a, **k: {
            "version": "3.9.3",
            "prefix": "/opt/env",
            "path": "/opt/env/bin/gdal-config",
        },
    )
    found = make_probe(
        bindings_error="ModuleNotFoundError: osgeo",
        gdal_config="/opt/env/bin/gdal-config",
        _config={"version": "3.9.3", "prefix": "/nonexistent-prefix", "path": "/opt/env/bin/gdal-config"},
    )
    plan = gdal_setup.make_plan(found)
    assert plan.strategy == "bindings"
    joined = " ".join(" ".join(step.command or []) for step in plan.steps)
    assert "gdal==3.9.3" in joined
    # The whole point: build against the runtime numpy, not an isolated one.
    assert "--no-build-isolation" in joined


def test_an_architecture_mismatch_does_not_go_back_to_the_package_manager(monkeypatch):
    """Homebrew supplied the wrong architecture; asking it again would loop."""
    # Windows short-circuits to the prebuilt-wheel route before any of this,
    # so pin the platform rather than let the assertion depend on the runner.
    monkeypatch.setattr(sys, "platform", "darwin")
    monkeypatch.setattr(gdal_setup, "_conda_executable", lambda: None)
    monkeypatch.setattr(
        gdal_setup, "_system_manager", lambda: ("/opt/homebrew/bin/brew", ["install", "gdal"], False)
    )
    found = make_probe(**WORKING)
    found.problems = (("arch-mismatch", "library is arm64 but Python is x86_64"),)
    strategy, notes = gdal_setup._choose_strategy(found)
    assert strategy == "source"
    assert any("same architecture" in note for note in notes)


def test_a_conda_environment_is_preferred_over_a_source_build(tmp_path, monkeypatch):
    monkeypatch.setattr(gdal_setup, "_conda_executable", lambda: "/opt/conda/bin/mamba")
    # conda-meta is what marks a prefix as a conda environment.
    (tmp_path / "conda-meta").mkdir()
    found = make_probe(bindings_error="no osgeo", sys_prefix=str(tmp_path))
    assert found.is_conda
    strategy, _ = gdal_setup._choose_strategy(found)
    assert strategy == "conda"


def test_root_is_not_assumed_for_system_packages(monkeypatch):
    monkeypatch.setattr(sys, "platform", "linux")
    monkeypatch.setattr(gdal_setup, "_conda_executable", lambda: None)
    monkeypatch.setattr(
        gdal_setup, "_system_manager", lambda: ("/usr/bin/apt-get", ["install", "-y", "libgdal-dev"], True)
    )
    monkeypatch.setattr(gdal_setup, "_is_root", lambda: False)
    found = make_probe(bindings_error="no osgeo")
    strategy, notes = gdal_setup._choose_strategy(found)
    assert strategy == "source"
    assert any("administrator rights" in note for note in notes)


def test_allow_root_makes_the_package_manager_usable(monkeypatch):
    """Otherwise --allow-root would still fall through to a half-hour build."""
    monkeypatch.setattr(sys, "platform", "linux")
    monkeypatch.setattr(gdal_setup, "_conda_executable", lambda: None)
    monkeypatch.setattr(
        gdal_setup,
        "_system_manager",
        lambda: ("/usr/bin/apt-get", ["install", "-y", "libgdal-dev"], True),
    )
    monkeypatch.setattr(gdal_setup, "_is_root", lambda: False)
    found = make_probe(bindings_error="no osgeo")
    strategy, _ = gdal_setup._choose_strategy(found, allow_root=True)
    assert strategy == "system"


def test_windows_uses_a_prebuilt_wheel(monkeypatch):
    monkeypatch.setattr(sys, "platform", "win32")
    monkeypatch.setattr(gdal_setup, "_conda_executable", lambda: None)
    found = make_probe(bindings_error="no osgeo")
    strategy, _ = gdal_setup._choose_strategy(found)
    assert strategy == "wheel"


def test_apt_refreshes_its_package_lists_first(monkeypatch):
    """A stock Debian image has no package lists; install alone would fail."""
    monkeypatch.setattr(
        gdal_setup,
        "_system_manager",
        lambda: ("/usr/bin/apt-get", ["install", "-y", "libgdal-dev"], True),
    )
    monkeypatch.setattr(gdal_setup, "_is_root", lambda: True)
    found = make_probe(bindings_error="no osgeo")
    plan = gdal_setup.make_plan(found, strategy="system")
    assert plan.steps[0].command == ["/usr/bin/apt-get", "update"]
    assert "install" in plan.steps[1].command


def test_repairing_broken_bindings_forces_a_rebuild(monkeypatch):
    """pip would otherwise say 'already satisfied' and change nothing.

    The commonest repair is a missing gdal_array, which leaves osgeo at
    exactly the right version -- and pip's wheel cache holds the bad build.
    """
    monkeypatch.setattr(
        gdal_setup,
        "_gdal_config_info",
        lambda *a, **k: {
            "version": "3.9.3",
            "prefix": "/opt/env",
            "path": "/opt/env/bin/gdal-config",
        },
    )
    fields = dict(WORKING, has_gdal_array=False, gdal_array_error="ImportError: numpy")
    found = make_probe(**fields)
    plan = gdal_setup.make_plan(found, strategy="bindings")
    command = plan.steps[-1].command
    assert "--force-reinstall" in command
    assert "--no-cache-dir" in command


def test_a_first_install_of_the_bindings_does_not_force(monkeypatch):
    """Nothing to displace, so leave pip's cache alone."""
    monkeypatch.setattr(
        gdal_setup,
        "_gdal_config_info",
        lambda *a, **k: {
            "version": "3.9.3",
            "prefix": "/opt/env",
            "path": "/opt/env/bin/gdal-config",
        },
    )
    found = make_probe(
        bindings_error="ModuleNotFoundError: osgeo",
        gdal_config="/opt/env/bin/gdal-config",
        _config={
            "version": "3.9.3",
            "prefix": "/nonexistent-prefix",
            "path": "/opt/env/bin/gdal-config",
        },
    )
    plan = gdal_setup.make_plan(found, strategy="bindings")
    assert "--force-reinstall" not in plan.steps[-1].command


# ---------------------------------------------------------------------------
# Plans
# ---------------------------------------------------------------------------


def test_a_source_plan_builds_the_library_then_the_bindings():
    found = make_probe(bindings_error="no osgeo")
    plan = gdal_setup.make_plan(found, strategy="source", prefix="/tmp/example-prefix")
    summaries = [step.summary for step in plan.steps]
    assert any("proj" in text for text in summaries)
    assert any("gdal" in text for text in summaries)
    # The bindings are always last: they need the library to exist first.
    assert "bindings" in summaries[-1]
    assert "/tmp/example-prefix" in " ".join(plan.notes)


def test_a_source_plan_can_be_described_without_running_anything():
    found = make_probe(bindings_error="no osgeo")
    plan = gdal_setup.make_plan(found, strategy="source", prefix="/tmp/example-prefix")
    text = "\n".join(plan.describe())
    assert "strategy: source" in text
    assert "estimated time" in text


def test_a_source_plan_refuses_to_start_without_a_compiler(monkeypatch, capsys):
    """Better a sentence naming the missing tool than cmake failing later."""
    monkeypatch.setattr(
        _gdal_source, "missing_build_tools", lambda: ["cmake", "a C/C++ compiler"]
    )
    found = make_probe(bindings_error="no osgeo")
    plan = gdal_setup.make_plan(found, strategy="source", prefix="/tmp/example-prefix")
    assert plan.blockers
    assert "cmake" in plan.blockers[0]
    # The steps are still listed, so `plan` can show what would happen.
    assert plan.steps
    assert gdal_setup.apply_plan(plan) is False
    assert "cannot run on this machine" in capsys.readouterr().out


def test_a_source_plan_has_no_blockers_when_the_tools_are_present(monkeypatch):
    monkeypatch.setattr(_gdal_source, "missing_build_tools", lambda: [])
    found = make_probe(bindings_error="no osgeo")
    plan = gdal_setup.make_plan(found, strategy="source", prefix="/tmp/example-prefix")
    assert plan.blockers == ()


def test_forcing_the_bindings_strategy_without_a_library_is_an_error(monkeypatch):
    monkeypatch.setattr(gdal_setup, "_gdal_config_info", lambda *a, **k: None)
    found = make_probe(bindings_error="no osgeo")
    with pytest.raises(ValueError, match="gdal-config"):
        gdal_setup.make_plan(found, strategy="bindings")


def test_a_system_plan_marks_the_step_that_needs_root(monkeypatch):
    monkeypatch.setattr(
        gdal_setup, "_system_manager", lambda: ("/usr/bin/apt-get", ["install", "-y", "libgdal-dev"], True)
    )
    monkeypatch.setattr(gdal_setup, "_is_root", lambda: False)
    found = make_probe(bindings_error="no osgeo")
    plan = gdal_setup.make_plan(found, strategy="system")
    assert plan.needs_root
    assert plan.steps[0].command[0] == "sudo"


def test_a_plan_needing_root_refuses_to_run_without_permission(capsys):
    plan = gdal_setup.Plan(
        "system",
        [gdal_setup.Step("install things", command=["sudo", "apt-get", "install"], needs_root=True)],
    )
    assert gdal_setup.apply_plan(plan, allow_root=False) is False
    printed = capsys.readouterr().out
    assert "sudo apt-get install" in printed
    assert "--strategy source" in printed


def test_dry_run_executes_nothing(capsys):
    def explode(log):
        raise AssertionError("a dry run must not run anything")

    plan = gdal_setup.Plan("source", [gdal_setup.Step("build the world", action=explode)])
    assert gdal_setup.apply_plan(plan, dry_run=True) is True
    assert "build the world" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# The prefix manifest
# ---------------------------------------------------------------------------


@pytest.fixture
def built_gdal(tmp_path, monkeypatch):
    """A recorded source-built GDAL belonging to a fake environment."""
    prefix = tmp_path / "gdal"
    (prefix / "share" / "proj").mkdir(parents=True)
    (prefix / "share" / "proj" / "proj.db").write_bytes(b"not really a database")
    (prefix / "share" / "gdal").mkdir(parents=True)

    environment = tmp_path / "env"
    environment.mkdir()
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.delenv("TOPOANALYSIS_NO_GDAL_ENV", raising=False)
    # activate_environment only acts for the interpreter the bindings were
    # built for, so the test has to be that interpreter.
    monkeypatch.setattr(sys, "prefix", str(environment))
    gdal_setup.write_manifest(str(prefix), str(environment))
    return prefix, environment


def test_manifest_round_trip(built_gdal):
    prefix, environment = built_gdal
    manifest = gdal_setup.read_manifest()
    assert manifest["prefix"] == str(prefix)
    assert manifest["python_prefix"] == str(environment)
    assert manifest["proj_lib"].endswith(os.path.join("share", "proj"))


def test_the_manifest_is_written_inside_the_environment(built_gdal):
    """So two virtualenvs can each have their own build."""
    _, environment = built_gdal
    assert (environment / "share" / "topoanalysis" / "gdal-prefix.json").exists()


def test_activate_environment_sets_data_paths(built_gdal, monkeypatch):
    prefix, _ = built_gdal
    monkeypatch.delenv("PROJ_LIB", raising=False)
    monkeypatch.delenv("GDAL_DATA", raising=False)
    gdal_setup.activate_environment()
    assert os.environ["PROJ_LIB"] == str(prefix / "share" / "proj")
    assert os.environ["PROJ_DATA"] == str(prefix / "share" / "proj")
    assert os.environ["GDAL_DATA"] == str(prefix / "share" / "gdal")


def test_activate_environment_displaces_a_stale_setting(built_gdal, monkeypatch):
    """The case this exists for: conda exported another GDAL's data paths.

    A PROJ 9 library handed PROJ 8's proj.db refuses every lookup, so leaving
    the inherited value in place is not the conservative choice -- it is the
    broken one.
    """
    prefix, _ = built_gdal
    monkeypatch.setenv("PROJ_LIB", "/opt/conda/envs/other/share/proj")
    monkeypatch.setenv("GDAL_DATA", "/opt/conda/envs/other/share/gdal")
    gdal_setup.activate_environment()
    assert os.environ["PROJ_LIB"] == str(prefix / "share" / "proj")
    assert os.environ["GDAL_DATA"] == str(prefix / "share" / "gdal")


def test_activate_environment_can_be_switched_off(built_gdal, monkeypatch):
    monkeypatch.setenv("TOPOANALYSIS_NO_GDAL_ENV", "1")
    monkeypatch.setenv("PROJ_LIB", "/somewhere/the/user/chose")
    gdal_setup.activate_environment()
    assert os.environ["PROJ_LIB"] == "/somewhere/the/user/chose"


def test_a_manifest_from_another_environment_is_ignored(built_gdal, monkeypatch):
    """Those data files belong to a library this interpreter is not running."""
    monkeypatch.setattr(sys, "prefix", "/some/entirely/different/env")
    monkeypatch.setenv("PROJ_LIB", "/opt/conda/envs/other/share/proj")
    gdal_setup.activate_environment()
    assert os.environ["PROJ_LIB"] == "/opt/conda/envs/other/share/proj"


def test_activate_environment_ignores_a_prefix_without_a_database(
    tmp_path, monkeypatch
):
    """An empty share/proj is not a usable replacement, so leave things be."""
    prefix = tmp_path / "gdal"
    (prefix / "share" / "proj").mkdir(parents=True)  # no proj.db
    environment = tmp_path / "env"
    environment.mkdir()
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setattr(sys, "prefix", str(environment))
    gdal_setup.write_manifest(str(prefix), str(environment))
    monkeypatch.setenv("PROJ_LIB", "/opt/conda/envs/other/share/proj")
    gdal_setup.activate_environment()
    assert os.environ["PROJ_LIB"] == "/opt/conda/envs/other/share/proj"


def test_activate_environment_is_silent_without_a_manifest(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path))
    monkeypatch.setattr(sys, "prefix", str(tmp_path / "nothing-here"))
    gdal_setup.activate_environment()  # must not raise


# ---------------------------------------------------------------------------
# Source pins
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("component", sorted(_gdal_source.SOURCES))
def test_every_source_is_pinned_by_checksum_over_https(component):
    spec = _gdal_source.SOURCES[component]
    assert spec["url"].startswith("https://")
    assert len(spec["sha256"]) == 64
    assert set(spec["sha256"]) <= set("0123456789abcdef")


@pytest.mark.parametrize("component", sorted(_gdal_source.SOURCES))
def test_the_pinned_url_matches_the_pinned_version(component):
    """A version bumped without its URL would download the old tarball."""
    spec = _gdal_source.SOURCES[component]
    major, minor, patch = gdal_setup._parse_version(spec["version"])
    plain = "%d.%d.%d" % (major, minor, patch)
    # SQLite encodes 3.53.4 in its filenames as 3530400.
    encoded = "%d%02d%02d00" % (major, minor, patch)
    assert plain in spec["url"] or encoded in spec["url"], spec["url"]


def test_library_directory_prefers_what_gdal_config_reports(tmp_path):
    """On Fedora and RHEL a 64-bit libgdal lives in lib64, not lib."""
    lib64 = tmp_path / "lib64"
    lib64.mkdir()
    (tmp_path / "lib").mkdir()
    suffix = ".dylib" if sys.platform == "darwin" else ".so"
    (lib64 / ("libgdal" + suffix)).write_bytes(b"")

    assert gdal_setup._library_directory(
        {"libs": "-L%s -lgdal" % lib64, "prefix": str(tmp_path)}
    ) == str(lib64)


def test_library_directory_finds_the_library_when_there_is_no_minus_l(tmp_path):
    """Guessing 'lib' would point the RPATH at an empty directory."""
    lib64 = tmp_path / "lib64"
    lib64.mkdir()
    (tmp_path / "lib").mkdir()  # exists but holds nothing
    suffix = ".dylib" if sys.platform == "darwin" else ".so"
    (lib64 / ("libgdal" + suffix)).write_bytes(b"")

    assert gdal_setup._library_directory(
        {"libs": "-lgdal", "prefix": str(tmp_path)}
    ) == str(lib64)


def test_a_too_old_cmake_is_caught_before_the_build_starts(monkeypatch):
    """PROJ 9 and GDAL 3.13 both declare a 3.16 floor."""
    monkeypatch.setattr(_gdal_source, "_cmake_version", lambda cmake: (3, 10))
    monkeypatch.setattr(_gdal_source, "_c_compiler", lambda: "/usr/bin/cc")
    monkeypatch.setattr(_gdal_source, "_cxx_compiler", lambda: "/usr/bin/c++")
    monkeypatch.setattr(_gdal_source.shutil, "which", lambda name: "/usr/bin/" + name)
    missing = _gdal_source.missing_build_tools()
    assert any("3.16" in item and "3.10" in item for item in missing)


def test_a_recent_enough_cmake_is_accepted(monkeypatch):
    monkeypatch.setattr(_gdal_source, "_cmake_version", lambda cmake: (3, 16))
    monkeypatch.setattr(_gdal_source, "_c_compiler", lambda: "/usr/bin/cc")
    monkeypatch.setattr(_gdal_source, "_cxx_compiler", lambda: "/usr/bin/c++")
    monkeypatch.setattr(_gdal_source.shutil, "which", lambda name: "/usr/bin/" + name)
    assert _gdal_source.missing_build_tools() == []


def test_an_unreadable_cmake_version_does_not_block_the_build(monkeypatch):
    """An unknown answer must not veto a build that would have worked."""
    monkeypatch.setattr(_gdal_source, "_cmake_version", lambda cmake: None)
    monkeypatch.setattr(_gdal_source, "_c_compiler", lambda: "/usr/bin/cc")
    monkeypatch.setattr(_gdal_source, "_cxx_compiler", lambda: "/usr/bin/c++")
    monkeypatch.setattr(_gdal_source.shutil, "which", lambda name: "/usr/bin/" + name)
    assert _gdal_source.missing_build_tools() == []


def test_a_source_build_needs_a_cxx_compiler_not_just_a_c_one(monkeypatch):
    """PROJ and GDAL are C++17; a bare gcc without g++ is a real machine."""
    monkeypatch.setattr(_gdal_source, "_cxx_compiler", lambda: None)
    monkeypatch.setattr(_gdal_source, "_c_compiler", lambda: "/usr/bin/cc")
    monkeypatch.setattr(_gdal_source.shutil, "which", lambda name: "/usr/bin/" + name)
    assert "a C++ compiler" in _gdal_source.missing_build_tools()


def test_a_stale_proj_variable_is_diagnosed_not_reinstalled():
    """Reinstalling GDAL would leave the offending variable exactly where it is."""
    found = make_probe(
        bindings_version="3.9.3",
        has_gdal_array=True,
        projection_works=False,
        projection_error="RuntimeError: proj.db lacks DATABASE.LAYOUT.VERSION",
        gdal_config="/opt/env/bin/gdal-config",
        proj_env={"PROJ_LIB": "/opt/conda/envs/other/share/proj"},
        _config={
            "version": "3.9.3",
            "prefix": "/nonexistent-prefix",
            "path": "/opt/env/bin/gdal-config",
        },
    )
    assert found.codes == ("no-proj-data",)
    plan = gdal_setup.make_plan(found)
    assert plan.strategy == "proj-data"
    assert plan.steps == ()
    assert any("unset PROJ_LIB" in blocker for blocker in plan.blockers)
    # The offending value is named, because that is the whole diagnosis.
    assert any("/opt/conda/envs/other/share/proj" in m for _, m in found.problems)


def test_build_order_covers_every_source():
    assert set(_gdal_source.BUILD_ORDER) == set(_gdal_source.SOURCES)


def test_overriding_a_pin_requires_its_checksum(monkeypatch):
    monkeypatch.setenv("TOPOANALYSIS_GDAL_SOURCE_GDAL_URL", "https://example.org/g.tar.gz")
    monkeypatch.delenv("TOPOANALYSIS_GDAL_SOURCE_GDAL_SHA256", raising=False)
    with pytest.raises(_gdal_source.BuildError, match="must be set together"):
        _gdal_source.apply_environment_overrides()


def test_components_needed_is_a_subset_of_the_build_order():
    needed = _gdal_source.components_needed()
    assert set(needed) <= set(_gdal_source.BUILD_ORDER)
    # PROJ and GDAL are never optional.
    assert "proj" in needed and "gdal" in needed


def test_extracting_a_hostile_archive_is_refused(tmp_path):
    """A tarball that writes outside its directory must not be unpacked."""
    archive = tmp_path / "evil.tar.gz"
    payload = io.BytesIO(b"pwned")
    with tarfile.open(archive, "w:gz") as tar:
        info = tarfile.TarInfo("../escaped.txt")
        info.size = len(payload.getvalue())
        tar.addfile(info, payload)

    with pytest.raises((_gdal_source.BuildError, tarfile.TarError)):
        _gdal_source._extract(str(archive), str(tmp_path / "into"))
    assert not (tmp_path / "escaped.txt").exists()


def test_a_corrupt_download_is_rejected(tmp_path):
    """A cached file with the wrong checksum is discarded, not built."""
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "thing.tar.gz").write_bytes(b"not what was pinned")
    with pytest.raises(_gdal_source.BuildError):
        # The URL is unreachable, so the only way this could succeed is by
        # trusting the bad cached copy.
        _gdal_source._download(
            "https://localhost:1/thing.tar.gz", "0" * 64, str(cache), lambda line: None
        )


# ---------------------------------------------------------------------------
# The command line
# ---------------------------------------------------------------------------


def test_doctor_reports_on_the_running_interpreter(capsys):
    code = gdal_setup.main(["doctor"])
    printed = capsys.readouterr().out
    assert "GDAL bindings" in printed
    assert code in (0, 1)


def test_plan_prints_a_plan_without_installing(capsys):
    code = gdal_setup.main(["plan", "--strategy", "source", "--prefix", "/tmp/x"])
    printed = capsys.readouterr().out
    assert code == 0
    assert "strategy: source" in printed


def test_installer_dry_run_is_harmless(tmp_path):
    """``install.py --dry-run`` must not create or modify anything."""
    import TopoAnalysis

    root = os.path.dirname(os.path.abspath(TopoAnalysis.__file__))
    script = os.path.join(root, "install.py")
    if not os.path.exists(script):  # an installed copy may not ship it
        pytest.skip("install.py is not part of this installation")

    before = set(os.listdir(root))
    result = subprocess.run(
        [sys.executable, script, "--dry-run"],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        cwd=root,
    )
    assert result.returncode == 0, result.stdout.decode("utf8", "replace")
    assert "stopping before anything is installed" in result.stdout.decode("utf8", "replace")
    assert set(os.listdir(root)) == before
