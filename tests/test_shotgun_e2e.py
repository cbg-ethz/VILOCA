import subprocess
import filecmp
import os
import pytest
import shutil

cwd = "./data_1"
f_path = "work/raw_reads/w-HXB2-2804-3004.extended-ref.fas"
base_files = []

@pytest.fixture(scope="session", autouse=True)
def base_files():
    return os.listdir(cwd)

@pytest.fixture(autouse=True)
def cleanup(base_files):
    def run(base_files):
        for f in set(os.listdir(cwd)) - set(base_files):
            dirpath = os.path.join(cwd, f)
            if os.path.exists(dirpath) and os.path.isdir(dirpath):
                shutil.rmtree(dirpath)
            if os.path.isfile(dirpath) and not os.path.isdir(dirpath):
                os.remove(dirpath)

    run(base_files)
    yield
    run(base_files)

def test_e2e_shorah():
    cwd_test_e2e = "./data_1/test_e2e_shorah"
    # Check if the directory exists, if not, create it
    if not os.path.exists(cwd_test_e2e):
        os.makedirs(cwd_test_e2e)

    original = subprocess.run(
        "../shotgun_test.sh", shell=True, check=True, cwd=cwd_test_e2e
    )
    assert original.returncode == 0
    assert os.path.isfile(os.path.join(cwd, f_path)) == False

    assert filecmp.cmp(
        "./data_1/test.csv",
        "./data_1/test_e2e_shorah/work/snv/SNVs_0.010000_final.csv",
        shallow=False
    )

def test_e2e_shorah_with_extended_window_mode():
    cwd_test_e2e_extended_window_mode = "./data_1/test_e2e_extended_window_mode"
    # Check if the directory exists, if not, create it
    if not os.path.exists(cwd_test_e2e_extended_window_mode):
        os.makedirs(cwd_test_e2e_extended_window_mode)
    p = subprocess.run(
        "../shotgun_test.sh --extended_window_mode", shell=True, check=True, cwd=cwd_test_e2e_extended_window_mode
    )

    assert p.returncode == 0
    assert os.path.isfile(os.path.join(cwd_test_e2e_extended_window_mode, f_path)) == True
