import os,sys,subprocess,glob
import shutil
import re
import matplotlib.pyplot as plt

def pdos_ecalj(cif_name):
    #to make pdos from cif using 'ecalj'
    try:
        cif2cell_c='cif2cell '+cif_name+' -p vasp --vasp-cartesian --vasp-format=5'
        subprocess.run(cif2cell_c,shell=True)
        subprocess.run('vasp2ctrl POSCAR',shell=True)
        subprocess.run('cp ctrls.POSCAR.vasp2ctrl ctrls.tmp',shell=True)
        subprocess.run('ctrlgenPDOStest.py tmp',shell=True)
        subprocess.run('cp ctrlgenPDOStest.ctrl.tmp ctrl.tmp',shell=True)
        subprocess.run('lmfa tmp',shell=True)
        subprocess.run('job_pdos tmp -np 4 NoGnuplot',shell=True)
        #subprocess.run('cif2cell '+cif_name+' >>ciftext.txt',shell=True)
    except:
        print('error in pdos_ecalj')

def listup_cif(directory):
    """
    with open(cifdir+'/cif_list.txt',mode='r') as f:
        ciflist=[re.search(r"([^/]*?)$",s.strip()).group() for s in f.readlines()]
    ciflist[i].replace('.cif',"")
    """
    if not os.path.isdir(directory):
        print('No such directory')
        sys.exit()
    files=glob.glob(directory+"/*.cif")
    for i,file in enumerate(files):
        files[i]=file.replace("./"+directory+"/",'')
    
    with open(directory+'/cif_list.txt',mode='w') as f:
        f.write('\n'.join(files))
    return

def make_pdos_in_di(cif_directory):
    '''sample program
    make_pdos_in_di('~/Directory where cif files are grouped together')
    '''
    get_now_directory=os.getcwd()
    listup_cif(cif_directory)
    with open(cif_directory+'/cif_list.txt',mode='r') as f:
        cif_list=[s.strip() for s in f.readlines()]
    result_di=cif_directory+'/result'
    try:
        os.mkdir(result_di)
    except:
        shutil.rmtree(result_di)
        os.mkdir(result_di)

    for cif_file in cif_list:
        try:
            directory_name=cif_file.replace('.cif','')
            directory_name=result_di+'/'+re.search(r"([^/]*?)$",directory_name).group()
            #subprocess.run('cp '+cif_file+' '+directory_name+'/'+cif_number+'.cif',shell=True)
            os.mkdir(directory_name)
            os.chdir(directory_name)
            pdos_ecalj(cif_file)
            os.chdir('..')
        except:
            print('error in for cif')
    os.chdir(get_now_directory)
