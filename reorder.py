import os
import shutil
import glob

j = os.path.join
xs = [
    'AcNH2AcNH2', 'BenzeneMeNH2NHpi', 'EthynePentane', 'NeopentanePentane',
    'PyridinePyridineTS', 'WaterPyridine', 'AcNH2Uracil', 'BenzeneMeOHOHpi',
    'EthyneWaterCHO', 'PentaneAcNH2', 'PyridinePyridinepipi', 'WaterWater',
    'AcOHAcOH', 'BenzeneNeopentane', 'MeNH2MeNH2', 'PentaneAcOH',
    'PyridineUracilpipi', 'UracilCyclopentane', 'UracilNeopentane',
    'AcOHUracil', 'BenzenePeptideNHpi', 'MeNH2MeOH', 'PentanePentane',
    'BenzenePyridineTS', 'MeNH2Peptide', 'PeptideEthene', 'UracilEthene',
    'BenzeneAcNH2NHpi', 'BenzenePyridinepipi', 'MeNH2Pyridine',
    'PeptideMeNH2', 'UracilEthyne', 'PeptideWater', 'UracilUracilpipi',
    'BenzeneAcOH', 'BenzeneUracilpipi', 'MeNH2Water', 'PeptideMeOH',
    'BenzeneAcOHOHpi', 'BenzeneWaterOHpi', 'MeOHMeNH2',
    'PeptidePentane', 'UracilPentane', 'WaterMeOH',
    'BenzeneBenzeneTS', 'CyclopentaneCyclopentane', 'MeOHMeOH',
    'PeptidePeptide', 'UracilUracilBP', 'WaterMeNH2',
    'BenzeneBenzenepipi', 'CyclopentaneNeopentane', 'MeOHPeptide',
    'BenzeneCyclopentane  EthenePentane', 'MeOHPyridine', 'PyridineEthene',
    'BenzeneEthene', 'EthyneAcOHOHpi', 'MeOHWater', 'PyridineEthyne',
    'BenzeneEthyneCHpi', 'EthyneEthyneTS', 'NeopentaneNeopentane',
    'PyridinePyridineCHN', 'WaterPeptide']


ds = ['0.9', '0.95', '1.0', '1.05', '1.1', '1.25', '1.5', '2.0']
files = ["INCAR", "POSCAR", "KPOINTS", "POTCAR"]
root = '/nfs/slac/g/suncatfs/ksb/s66x8'
sub = '/nfs/slac/g/suncatfs/ksb/beefjobs/bulks/7500/Al/latopt/subVASP.sh'
beef = '/nfs/slac/g/suncatfs/yasheng/msurf_new/candidate_1/s66x8/BEEFCAR'


def submit(pth: str) -> None:
    os.system(';'.join([
        "cd %s" % pth,
        "chmod 755 %s/subVASP.sh" % pth,
        "bsub -n 16 -q suncat3 -W 40:00 %s/subVASP.sh" % pth]))


def copy2files(p: str) -> None:
    shutil.copyfile(sub, j(p, 'subVASP.sh'))
    shutil.copyfile(beef, j(p, 'BEEFCAR'))
    print('\tcopied BEEFCAR to ', p)
    print('\tcopied subVASP.sh to ', p)


def main() -> None:
    flag = True
    for x in xs:
        print(x)
        for dist in ds:
            p = j(root, x, dist)
            if not os.path.exists(p):
                os.mkdir(p)
                print('created ', p)
            copy2files(p)
            for f in files:
                name = f + '_' + x + "_complex-" + dist
                fi = j(root, x, name)
                shutil.copyfile(fi, j(p, f))
                print('\tmoved ', name, 'to ', p)
            submit(p)

        for d in '01':
            suf = '_mol' + d
            p = j(root, x, suf[1:])
            if not os.path.exists(p):
                os.mkdir(p)
                print('created ', p)
            copy2files(p)
            for f in files:
                name = f + '_' + x + suf
                fi = j(root, x, name)
                shutil.copyfile(fi, j(p, f))
                print('\tmoved ', name, 'to ', p)
            submit(p)
        if flag:
            import pdb; pdb.set_trace()
            flag = False


def dbh() -> None:
    root = '/nfs/slac/g/suncatfs/ksb/dbh24_scan/dbh24_r%i/'
    for i in range(1, 13):
        r = root % i
        mols = [x[x.rfind('24_') + 3:] for x in glob.glob(r + 'OUTCAR*')]
        for mol in mols:
            os.mkdir(r + mol)
            for f in files:
                shutil.copyfile(j(r, f + '_dbh24_' + mol), j(r + mol, f))

            shutil.copyfile('/nfs/slac/g/suncatfs/ksb/re42_scan/re42_gas/run_gamma.bash', j(r + mol, 'run_gamma.bash'))

            os.chdir(r + mol)
            os.system("sed -i 's/BF/SCAN/g' INCAR")
            os.system("chmod 755 run_gamma.bash")

            os.system('/usr/local/bin/bsub -n 16 -q suncat3 -W 40:00 -o out.log -e err.log `pwd`/run_gamma.bash')





if __name__ == '__main__':
    dbh()
