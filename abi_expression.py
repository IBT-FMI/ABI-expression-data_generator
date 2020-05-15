import argparse
import copy
import json
import os
import sys
import shutil
import urllib
import urllib.request
import zipfile
import numpy
import nibabel
import re
import socket
import tarfile
from socket import timeout
from collections import defaultdict
from mhd_utils_3d import *
from nipype.interfaces.ants import ApplyTransforms
#from nipype.interfaces.fsl import fslorient

API_SERVER = "http://api.brain-map.org/"
API_DATA_PATH = API_SERVER + "api/v2/data/"

#set default timeout for urllib.request
socket.setdefaulttimeout(100)

def GetGeneNames(startRow=0,numRows=2000,totalRows=-1):
    """
    Queries the Allen Mouse Brain Institute website for all gene expression data available for download.

    Parameters
    -----------
    startRow: int
      Starting row
    numRows: int
      Number of rows shown per query. Max is 2000.
    totalRows: int
      Number of total rows to query. If set to -1, all available rows will be displayed.

    Returns
    --------
    info: defaultdict
        key: genename, value: list of corresponding SectionDataSetID
        (SectionDataSet: see "http://help.brain-map.org/display/api/Data+Model")
        ID needed to specify download target.

    """
    startRow = startRow
    numRows = numRows
    totalRows = totalRows
    rows = []
    GeneNames = []
    SectionDataSetID = []
    info = defaultdict(list)

    done = False

    while not done:
        pagedUrl = API_DATA_PATH +"query.json?criteria=model::SectionDataSet,rma::criteria,[failed$eqfalse],products[abbreviation$eq'Mouse'],treatments[name$eq'ISH'],rma::include,genes,specimen(donor(age)),plane_of_section" + '&startRow={0}&numRows={1}'.format(startRow,numRows)

        source = urllib.request.urlopen(pagedUrl).read()
        response = json.loads(source)
        rows += response['msg']
        for x in response['msg']:

            if x['failed'] == False:
                info[x['genes'][0]['acronym']].append(x['id'])

        if totalRows < 0:
            totalRows = int(response['total_rows'])

        startRow += len(response['msg'])

        if startRow >= totalRows:
            done = True

    return info

def sub_name(s):
   """replaces brackets with '-' and removes ',* """

   #replace brackets with '_' and remove ',*
   #s = re.sub('\W', '',s_n)  dont do that, will replace '-' as well
   s_n = re.sub('[()]',"-",s)
   s_n = re.sub("'","",s_n)
   s_n = re.sub('\*',"",s_n)
   return s_n

def download_all_ISH(info,folder_name="ABI-expression-data-9999"):
    """
    Downloads all datasets corresponding to SectionDataSetID given, converts data format from mhd/raw to nii and transforms data to dsurqec template.

    Parameters
   -----------
    info: defaultdict
        key: genename, value: list of corresponding SectionDataSetID (SectionDataSet: see "http://help.brain-map.org/display/api/Data+Model")
        ID needed to specify download target.

    """
    failed_downloads = list()
    data_path = folder_name
    #TODO: test timeout, and if timeout occurs, repeat attempt download
    if not os.path.exists(data_path): os.mkdir(data_path)
    download_url = "http://api.brain-map.org/grid_data/download/"
    for gene in info.items():
        gene_name = gene[0]
        gene_ids = gene[1]
        gene_r = sub_name(gene_name)
        path_to_gene = os.path.join(data_path,gene_r)
        #TODO: check if already downloaded,right file exists
        if not os.path.exists(path_to_gene) : os.mkdir(path_to_gene)
        for gene_id in gene_ids:
            url = download_url + str(gene_id)
            try:
                fh = urllib.request.urlretrieve(url)
            except timeout:
                print("timeout with " + str(gene_id)+ gene_r)
                failed_downloads.append(gene_id)
                shutil.rmtree(path_to_gene)
                continue

            zf = zipfile.ZipFile(fh[0])
            filename = str.split((fh[1]._headers[6][1]),'filename=')[1]
            filename = str.split(filename,'.zip')[0]
            filename = sub_name(filename)
            path_to_folder = os.path.join(path_to_gene,filename)
            zf.extractall(os.path.join(path_to_gene,filename))
            zf.close()

            #some datasets without energy.mhd file (not available on API) Skip and delete folder
            if not os.path.isfile(os.path.join(path_to_folder,"energy.mhd")):
                    print("removing" + str(gene_id)+gene_r)
                    shutil.rmtree(path_to_folder)
                    continue

            path_to_mhd = os.path.join(path_to_folder,"energy.mhd")
            path_to_raw = os.path.join(path_to_folder,"energy.raw")
            path_to_nifti = convert_raw_to_nii(path_to_mhd,filename)
            apply_composite(path_to_nifti)
            os.remove(path_to_nifti)
            os.remove(path_to_mhd)
            os.remove(path_to_raw)

    if len(failed_downloads) > 0:
        print("failed: ")
        for item in failed_downloads:
            print(str(item))

def struc_unionize(my_id):
   """Queries the ABI API for an expression summary for structure-id = 997 (the entire brain) and returns
   expression density, expression energy.

   Parameters
   ----------

   id: int
      Unique identifier for ABI SectionDataSetID

   Returns
   -------

   density: float
      Expression density (fraction of voxels with signal detected).
   energy: float
      Expression energy, expression density modulated by signal intensity.


   """
   url = "http://api.brain-map.org/api/v2/data/SectionDataSet/{}.json?include=structure_unionizes[structure_id$eq997]".format(str(my_id))
   source = urllib.request.urlopen(url).read()
   response = json.loads(source)
   i = 0
   density = 0
   energy = 0
   for x in response['msg']:
      density = x['structure_unionizes'][0]['expression_density']
      energy = x['structure_unionizes'][0]['expression_energy']

   return density, energy

def convert_raw_to_nii(input_file,output_file):
    """
    Converts mhd/raw format to NIfTI and orients data matrix in RAS-space.

    Parameters
    -----------
    input_file : str
      path to .mhd file
    output_file : str
      filename prefix

    Returns
    ---------
    output_path : str
      path to generated NIfTI - file

    """
    path = os.path.abspath('.')
    image_array, meta_header = load_raw_data_with_mhd(input_file)

    #Read header infomormation and create affine matrix
    dims = numpy.array(meta_header["ElementSpacing"].split(" "),dtype=numpy.float)
    affine_matrix = numpy.zeros((4,4),dtype=numpy.float)
    affine_matrix[0,0] = dims[0]
    affine_matrix[1,1] = dims[1]
    affine_matrix[2,2] = dims[2]

    #Orient in RAS
    image_array = numpy.swapaxes(image_array,1,2)

    image_array = image_array[::-1,:,:]
    image_array = image_array[:,:,::-1]
    image_array = image_array[:,::-1,:]

    #Bring to the right units
    affine_matrix = affine_matrix*0.001

    affine_matrix[3,3] = 1

    img = nibabel.Nifti1Image(image_array,affine_matrix)
    name = output_file + '.nii.gz'
    output_path = os.path.join(os.path.dirname(input_file),name)
    nibabel.save(img,output_path)

    return output_path

def apply_composite(file):
    """
    Uses ANTS ApplyTransforms to transform image to target space.

    Parameters :
    ------------

    file : str
        path to image

    """
    at = ApplyTransforms()
    at.inputs.dimension = 3
    at.inputs.input_image = file
    at.inputs.reference_image = '/usr/share/mouse-brain-atlases/dsurqec_200micron_masked.nii'
    name = str.split(os.path.basename(file),'.nii')[0] + '_2dsurqec.nii.gz'
    at.inputs.interpolation = 'NearestNeighbor' #TODO: Sure?? Yes, avoiding values between -1 and 0 (ckeck)
    at.inputs.output_image = os.path.join(os.path.dirname(file),name)
    at.inputs.transforms = '/usr/share/mouse-brain-atlases/abi2dsurqec_Composite.h5'
    at.run()

    #TODO sform to qform

#TODO: possibly info where no dataset is available

def save_dens_energy(info, folder_name):
   """saves the informatin obtained from the Allen Mouse Brain structure unionize module"""
   if not os.path.isdir(folder_name):os.mkdir(folder_name)
   file_path = os.path.join(folder_name,"density_energy.csv")
   f = open(file_path,"w+")
   f.write("acronym,id,density,energy")
   for gene in info:
      for gene_id in info[gene]:
         d,e = struc_unionize(gene_id)
         f.write('\n')
         f.write(gene + "," + str(gene_id) + "," + str(d) + "," + str(e))

def save_info(info,folder_name):
   """saves the information about genename and correspoding SectionDataSetID as csv """
   file_path=os.path.join(folder_name,"ABI-genes-datasetid.csv")
   f = open(file_path,"w+")
   for gene in info:
        f.write('\n')
        f.write(gene)
        for gene_id in info[gene]:
            f.write("," + str(gene_id))

def create_archive(folder_name):
   """creates .tar.xz archive """
   tar_name = folder_name + ".tar.xz"
   with tarfile.open(tar_name, "w:xz") as tar_handle:
      for root,dirs,files in os.walk(folder_name):
         for file in files:
            print(file)
            tar_handle.add(os.path.join(root,file))

def main():
   parser = argparse.ArgumentParser(description="ABI-expression",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('--package_name','-n',type=str,default="ABI-expression-data")
   parser.add_argument('--package_version','-v',type=str,default="9999")
   parser.add_argument('--startRow','-s',type=int,default=0)
   parser.add_argument('--numRows','-r',type=int,default=2000)
   parser.add_argument('--totalRows','-t',type=int,default=-1)
   args=parser.parse_args()

   folder_name = args.package_name + "-" + args.package_version
   info=GetGeneNames(startRow=args.startRow,numRows=args.numRows,totalRows=args.totalRows)
   save_dens_energy(info,folder_name)
   download_all_ISH(info,folder_name=folder_name)
   save_info(info,folder_name)
   create_archive(folder_name)

if __name__ == "__main__":
    main()
