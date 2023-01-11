import subprocess
import filecmp
import os.path

cwd = "./data_1"
f_path = "raw_reads/w-HXB2-2804-3004.extended-ref.fas.gz"

def test_e2e_shorah():
    original = subprocess.run(
        "./shotgun_test.sh", shell=True, check=True, cwd=cwd
    )
    assert original.returncode == 0
    assert os.path.isfile(os.path.join(cwd, f_path)) == False

    assert filecmp.cmp(
        "./data_1/test.csv",
        "./data_1/snv/SNVs_0.010000_final.csv",
        shallow=False
    )

# TODO do all full cleanup in between

def test_e2e_shorah_with_extended_window_mode():
    p = subprocess.run(
        "./shotgun_test.sh --extended_window_mode", shell=True, check=True, cwd=cwd
    )

    assert p.returncode == 0
    assert os.path.isfile(os.path.join(cwd, f_path)) == True
