from glob import glob
from pymol import cmd
import os
import sys

def align(args):
    # default
    num_args = len(args) - 1
    if num_args == 2:
        protein_name = args[1]
        mobile_path = args[2]
        tgtpdb_path = glob(f"./dat/base_receptor/{protein_name}/*.pdb")[0]
        ligand_path = glob(f"./dat/base_ligand/{protein_name}/*.pdb")[0]
        print(tgtpdb_path, mobile_path)
        try:
            cmd.reinitialize("everything")
            cmd.load(mobile_path, object="target")
            cmd.load(tgtpdb_path, object="reference")
            cmd.load(ligand_path, object="ligand")
            res=cmd.align(mobile="target", target="reference", cycles=0)
            # ligandから8Å以内の残基をそれぞれ選択
            cmd.select("seleT", 'byres target within 8.0 of ligand')
            cmd.select("seleR", 'byres reference within 8.0 of ligand')
            # align
            res=cmd.align(mobile="seleT", target="seleR", cycles=0)

            output_path = f"{os.path.splitext(mobile_path)[0]}_align.pdb"
            cmd.save(output_path, "target")
            print("saved", output_path, res[0])
            return output_path
        except:
            output_path = f"{os.path.dirname(mobile_path)}/{os.path.splitext(os.path.basename(mobile_path))[0]}_align.pdb"
            print(f"[ERROR] {output_path}")
    elif num_args == 3:
        protein_name = args[1]
        pdb1_path = args[2]
        pdb2_path = args[3]

        tgtpdb_path = glob(f"./dat/base_receptor/{protein_name}/*.pdb")[0]
        ligand_path = glob(f"./dat/base_ligand/{protein_name}/*.pdb")[0]

        try:
            cmd.reinitialize("everything")
            cmd.load(pdb1_path, object="pdb1")
            cmd.load(tgtpdb_path, object="reference")
            cmd.load(ligand_path, object="ligand")
            res=cmd.align(mobile="pdb1", target="reference", cycles=0)
            cmd.load(pdb2_path, object="pdb2")
            res=cmd.align(mobile="pdb2", target="reference", cycles=0)
            # ligandから8Å以内の残基をそれぞれ選択
            cmd.select("seleT", 'byres pdb1 within 8.0 of ligand')
            cmd.select("seleR", 'byres pdb2 within 8.0 of ligand')
            res=cmd.align(mobile="seleT", target="seleR", cycles=0)
            cmd.reinitialize("everything")
            return res[0]
        except:
            return None


if __name__ == '__main__':
    args = sys.argv
    align(args)