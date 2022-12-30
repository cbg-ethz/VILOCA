import subprocess
import filecmp
import os.path

f_path = "raw_reads/w-HXB2-2804-3004.extended-ref.fas.gz"

def test_e2e_shorah():
    original = subprocess.run(
        "./shotgun_test.sh", shell=True, check=True, cwd="./data_1"
    )
    assert original.returncode == 0
    assert os.path.isfile(f_path) == False

    assert filecmp.cmp(
        "./data_1/test.csv",
        "./data_1/snv/SNVs_0.010000_final.csv",
        shallow=False
    )

def test_e2e_shorah_with_extended_window_mode():
    p = subprocess.run(
        "echo ./shotgun_test.sh --extended_window_mode", shell=True, check=True, cwd="./data_1"
    )

    assert p.returncode == 0
    assert os.path.isfile(f_path) == True
